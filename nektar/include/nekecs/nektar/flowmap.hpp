// Nektar++ backend for nekecs: a vector space over velocity fields and a
// deterministic flow map u -> Phi_T(u) built on the incompressible
// velocity-correction scheme (VCS).
//
// State vector
//   FieldVector::u[i]  physical-space values of velocity component i at the
//                      quadrature points (pressure is NOT part of the state:
//                      it is a slave variable recomputed by the solver).
//   FieldVector::extra scalar unknowns appended by the problem (period T,
//                      shifts, ...). Empty for equilibria.
//
// Multi-level states (StartMode::History): u holds [u_n, u_{n-1}, ..]
// components, k = IMEX order levels; FieldVector::u.size() = k * nvel.
//
// Inner product: <a, b> = sum_i \int_Omega a_i b_i dx  (quadrature-weighted,
// MPI-reduced by ExpList::Integral, over all levels)  +  sum_j a.extra[j] b.extra[j].
//
// Determinism (see DESIGN.md "Hidden history"): before every flow-map call the
// VCS high-order pressure-BC history is regenerated, the pressure field is
// zeroed and time/step counters are reset, so Phi_T depends on u only.
//
// Time stepping uses our own loop (FlowMapVCS::Integrate) so the multistep
// start-up can be chosen per call (StartMode; DESIGN.md "Start-up artefact").
#pragma once

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <SpatialDomains/MeshGraphIO.h>

#include <filesystem>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace nekecs
{
namespace nektar
{

using Nektar::Array;
using Nektar::NekDouble;
using Nektar::OneD;

struct FieldVector
{
    std::vector<Array<OneD, NekDouble>> u;
    std::vector<double> extra;
};

// Vector space over FieldVector for the nekecs core algorithms.
class NektarSpace
{
public:
    using Vector = FieldVector;

    NektarSpace(Array<OneD, Nektar::MultiRegions::ExpListSharedPtr> fields,
                int nvel)
        : m_fields(fields), m_nvel(nvel)
    {
    }

    Vector zeros_like(const Vector &x) const
    {
        Vector z;
        for (const auto &a : x.u)
        {
            z.u.emplace_back(a.size(), 0.0); // new storage
        }
        z.extra.assign(x.extra.size(), 0.0);
        return z;
    }

    void copy(const Vector &src, Vector &dst) const
    {
        for (std::size_t i = 0; i < src.u.size(); ++i)
        {
            Vmath::Vcopy(src.u[i].size(), src.u[i], 1, dst.u[i], 1);
        }
        dst.extra = src.extra;
    }

    double dot(const Vector &a, const Vector &b) const
    {
        // All arrays: velocity components of every stored time level
        // (layout [u_n components, u_{n-1} components, ...]).
        double s = 0.0;
        for (std::size_t i = 0; i < a.u.size(); ++i)
        {
            const int n = a.u[i].size();
            Array<OneD, NekDouble> ab(n);
            Vmath::Vmul(n, a.u[i], 1, b.u[i], 1, ab, 1);
            s += m_fields[i % m_nvel]->Integral(ab);
        }
        for (std::size_t j = 0; j < a.extra.size(); ++j)
        {
            s += a.extra[j] * b.extra[j];
        }
        return s;
    }

    void axpy(double a, const Vector &x, Vector &y) const
    {
        for (std::size_t i = 0; i < x.u.size(); ++i)
        {
            Vmath::Svtvp(x.u[i].size(), a, x.u[i], 1, y.u[i], 1,
                                 y.u[i], 1);
        }
        for (std::size_t j = 0; j < x.extra.size(); ++j)
        {
            y.extra[j] += a * x.extra[j];
        }
    }

    void scale(double a, Vector &x) const
    {
        for (auto &arr : x.u)
        {
            Vmath::Smul(arr.size(), a, arr, 1, arr, 1);
        }
        for (auto &e : x.extra)
        {
            e *= a;
        }
    }

private:
    Array<OneD, Nektar::MultiRegions::ExpListSharedPtr> m_fields;
    int m_nvel;
};

// How a flow-map call starts the IMEX multistep integrator.
enum class StartMode
{
    // Nektar++'s own start-up: first step(s) with lower-order IMEX phases and
    // lower-order pressure-BC extrapolation. Reproduces DoSolve().
    Standard,
    // All k history levels set to u (values y_{n-j} = u, explicit terms
    // dt N(y_{n-j}) = dt N(u)) and the pressure-BC history filled from u;
    // every step uses the order-k scheme. A steady state of the regular
    // scheme is then an exact fixed point for any number of steps.
    ConstantHistory,
    // The k-1 older levels u_{n-1} .. u_{n-k+1} are given explicitly (part of
    // the state). Integrator and pressure-BC histories are rebuilt from them,
    // so a call is exactly nsteps regular order-k steps: the map on
    // (u_n, .., u_{n-k+1}) is a semigroup (exact periodic orbits, exact
    // eigenvalues of the time-stepper).
    History,
};

// VelocityCorrectionScheme with access to the pressure-BC history and its
// own time loop (instead of UnsteadySystem::DoSolve, which always runs the
// start-up phases and also drives filters, IO and CFL control).
class FlowMapVCS : public Nektar::VelocityCorrectionScheme
{
public:
    FlowMapVCS(const Nektar::LibUtilities::SessionReaderSharedPtr &s,
               const Nektar::SpatialDomains::MeshGraphSharedPtr &g)
        : UnsteadySystem(s, g), VelocityCorrectionScheme(s, g)
    {
    }

    void ResetHistory()
    {
        m_extrapolation->GenerateHOPBCMap(m_session);
    }

    using Levels = std::vector<Nektar::Array<OneD, Nektar::Array<OneD, NekDouble>>>;

    int Order() const
    {
        return static_cast<int>(m_intScheme->GetOrder());
    }

    // Change the time step (periodic orbits: dt = T / N). Besides m_timestep
    // this updates the pressure-BC extrapolation, which keeps its own copy of
    // TimeStep from the session (used in the du/dt boundary term), and drops
    // the cached viscous Helmholtz matrices (lambda = 1/(a_ii dt nu)); the
    // pressure Poisson matrix does not depend on dt and is kept.
    void SetDt(NekDouble dt)
    {
        if (dt == m_timestep)
        {
            return;
        }
        m_session->SetParameter("TimeStep", dt);
        m_timestep = dt;
        SetUpExtrapolation();
        for (int i = 0; i < static_cast<int>(m_velocity.size()); ++i)
        {
            m_fields[m_velocity[i]]->ClearGlobalLinSysManager();
        }
    }

    // Change the kinematic viscosity (continuation in Re). nu enters the
    // viscous Helmholtz factor lambda = 1/(a_ii dt nu) and the forcing scale
    // 1/nu (m_diffCoeff), and the pressure BCs (m_kinvis, passed to the
    // extrapolation at every step; it caches no copy). Cached velocity
    // matrices are dropped. Not supported: per-field DiffusionCoefficient
    // functions and SVV (its kernels are scaled by 1/nu at set-up).
    void SetKinvis(NekDouble nu)
    {
        if (nu == m_kinvis)
        {
            return;
        }
        if (m_useSpecVanVisc || m_useHomo1DSpecVanVisc)
        {
            throw std::runtime_error(
                "nekecs: changing Kinvis with SVV is not supported");
        }
        for (int n = 0; n < m_nConvectiveFields; ++n)
        {
            if (m_session->DefinesFunction("DiffusionCoefficient",
                                           m_session->GetVariable(n)))
            {
                throw std::runtime_error(
                    "nekecs: changing Kinvis with a DiffusionCoefficient "
                    "function is not supported");
            }
        }
        m_kinvis = nu;
        for (int n = 0; n < m_nConvectiveFields; ++n)
        {
            m_diffCoeff[n] = nu;
        }
        m_session->SetParameter("Kinvis", nu);
        std::ostringstream os;
        os << std::setprecision(17) << nu;
        m_fieldMetaDataMap["Kinvis"] = os.str();
        for (int i = 0; i < static_cast<int>(m_velocity.size()); ++i)
        {
            m_fields[m_velocity[i]]->ClearGlobalLinSysManager();
        }
    }

    NekDouble Kinvis() const
    {
        return m_kinvis;
    }

    // Re-evaluate every boundary condition from its expression (session
    // parameters are bound at evaluation time, so a changed parameter takes
    // effect). Time 0: flow-map calls always start at t = 0.
    void ReevaluateBCs()
    {
        for (int i = 0; i < static_cast<int>(m_fields.size()); ++i)
        {
            m_fields[i]->EvaluateBoundaryConditions(0.0,
                                                    m_session->GetVariable(i));
        }
    }

    // Advance m_fields by nsteps fixed steps of m_timestep from m_time.
    // Mirrors the core of UnsteadySystem::v_DoSolve (pre-integrate, step,
    // copy back, forward transform, post-integrate) without filters, IO,
    // CFL or abort handling.
    //
    // mode == History: `prev` holds u_{n-1} .. u_{n-k+1} (k = IMEX order) on
    // input and the corresponding levels at the end on output.
    void Integrate(int nsteps, StartMode mode, Levels *prev = nullptr)
    {
        using namespace Nektar;
        const int nvar = static_cast<int>(m_intVariables.size());
        Array<OneD, Array<OneD, NekDouble>> fields(nvar);
        for (int i = 0; i < nvar; ++i)
        {
            fields[i] = m_fields[m_intVariables[i]]->UpdatePhys();
            m_fields[m_intVariables[i]]->SetPhysState(false);
        }

        // Builds the final-phase solution vector at u and evaluates the ODE
        // RHS once at u (which also performs the first pressure-BC update).
        m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

        auto deep = [](const Array<OneD, Array<OneD, NekDouble>> &a) {
            Array<OneD, Array<OneD, NekDouble>> b(a.size());
            for (std::size_t j = 0; j < a.size(); ++j)
            {
                b[j] = Array<OneD, NekDouble>(a[j].size());
                Vmath::Vcopy(a[j].size(), a[j], 1, b[j], 1);
            }
            return b;
        };
        const int k = Order();

        std::size_t phase_offset = 0;
        if (mode != StartMode::Standard)
        {
            if (m_intScheme->GetName() != "IMEX")
            {
                throw std::runtime_error(
                    "nekecs: history starts need an IMEX scheme");
            }
            const auto &sv = m_intScheme->GetSolutionVector();
            if (static_cast<int>(sv.size()) != 2 * k)
            {
                throw std::runtime_error(
                    "nekecs: unexpected IMEX solution-vector layout");
            }
            // History levels u_{n-1} .. u_{n-k+1}.
            Levels hist;
            if (mode == StartMode::History)
            {
                if (!prev || static_cast<int>(prev->size()) != k - 1)
                {
                    throw std::runtime_error(
                        "nekecs: History start needs k-1 previous levels");
                }
                hist = *prev;
            }
            else
            {
                for (int j = 1; j < k; ++j)
                {
                    hist.push_back(deep(fields));
                }
            }

            // Solution vector layout: [y_n .. y_{n-k+1}, dtN_n .. dtN_{n-k+1}].
            // Explicit terms dt N(u_{n-j}) come from the ODE RHS; evaluating
            // it oldest -> newest, then at u_n, also leaves the VCS pressure-BC
            // history (m_pressureHBCs, m_iprodnormvel) holding exactly the
            // levels u_{n-k+1} .. u_n in order, as after regular stepping.
            Array<OneD, Array<OneD, NekDouble>> rhs(nvar);
            for (int i = 0; i < nvar; ++i)
            {
                rhs[i] = Array<OneD, NekDouble>(fields[i].size());
            }
            for (int j = k - 1; j >= 1; --j)
            {
                const auto &y = hist[j - 1]; // u_{n-j}
                m_ode.DoOdeRhs(y, rhs, m_time - j * m_timestep);
                auto dtN = deep(rhs);
                for (int i = 0; i < nvar; ++i)
                {
                    Vmath::Smul(dtN[i].size(), m_timestep, dtN[i], 1, dtN[i], 1);
                }
                m_intScheme->SetSolutionVector(j, deep(y));
                m_intScheme->SetSolutionVector(k + j, dtN);
            }
            m_ode.DoOdeRhs(fields, rhs, m_time); // pressure-BC level for u_n
            // Use the final (order-k) phase from the first step.
            phase_offset = 1000;
        }

        for (int step = 0; step < nsteps; ++step)
        {
            v_PreIntegrate(step);
            fields = m_intScheme->TimeIntegrate(phase_offset + step, m_timestep);
            m_time += m_timestep;
            for (int i = 0; i < nvar; ++i)
            {
                auto &f = m_fields[m_intVariables[i]];
                f->SetPhys(fields[i]);
                fields[i] = f->UpdatePhys();
                if (v_RequireFwdTrans())
                {
                    f->FwdTransLocalElmt(f->GetPhys(), f->UpdateCoeffs());
                }
                f->SetPhysState(false);
            }
            v_PostIntegrate(step);
        }

        if (mode == StartMode::History)
        {
            const auto &sv = m_intScheme->GetSolutionVector();
            for (int j = 1; j < k; ++j)
            {
                for (int i = 0; i < nvar; ++i)
                {
                    Vmath::Vcopy(sv[j][i].size(), sv[j][i], 1, (*prev)[j - 1][i], 1);
                }
            }
        }
    }
};

// Flow map Phi_{N dt}(u) for a VCS session.
class VCSFlowMap
{
public:
    VCSFlowMap(const Nektar::LibUtilities::SessionReaderSharedPtr &session,
               const Nektar::SpatialDomains::MeshGraphSharedPtr &graph)
        : m_session(session)
    {
        if (session->GetSolverInfo("SolverType") != "VelocityCorrectionScheme")
        {
            throw std::runtime_error(
                "nekecs: only SolverType=VelocityCorrectionScheme is "
                "supported (other schemes have not been checked for hidden "
                "time-stepping history)");
        }
        m_eq = Nektar::MemoryManager<FlowMapVCS>::AllocateSharedPtr(session,
                                                                     graph);
        std::static_pointer_cast<Nektar::SolverUtils::EquationSystem>(m_eq)
            ->InitObject();
        m_eq->DoInitialise();
        m_nvel = static_cast<int>(session->GetVariables().size()) - 1;
        if (session->GetVariables().back() != "p")
        {
            throw std::runtime_error("nekecs: expected variables u,v[,w],p");
        }
    }

    int nvel() const
    {
        return m_nvel;
    }
    NekDouble timestep() const
    {
        return m_eq->GetTimeStep();
    }
    int calls() const
    {
        return m_calls;
    }
    const Nektar::LibUtilities::SessionReaderSharedPtr &session() const
    {
        return m_session;
    }

    // Diagnostics only: disabling either makes Phi depend on solver state
    // left over from the previous call.
    void set_zero_pressure(bool on)
    {
        m_zero_pressure = on;
    }
    void set_reset_history(bool on)
    {
        m_reset_history = on;
    }

    // Integrator start-up for every flow-map call (see StartMode).
    void set_start_mode(StartMode m)
    {
        m_start = m;
    }
    StartMode start_mode() const
    {
        return m_start;
    }
    // Regression/diagnostics: use Nektar++'s DoSolve instead of our loop
    // (equivalent to StartMode::Standard, plus filters and IO).
    void set_use_dosolve(bool on)
    {
        m_use_dosolve = on;
    }

    NektarSpace space() const
    {
        return NektarSpace(m_eq->UpdateFields(), m_nvel);
    }

    // Velocity currently held by the solver (after DoInitialise: the
    // session's initial condition).
    FieldVector state() const
    {
        FieldVector x;
        auto &f = m_eq->UpdateFields();
        for (int i = 0; i < m_nvel; ++i)
        {
            Array<OneD, NekDouble> a(f[i]->GetNpoints());
            Vmath::Vcopy(a.size(), f[i]->GetPhys(), 1, a, 1);
            x.u.push_back(a);
        }
        return x;
    }

    int order() const
    {
        return m_eq->Order();
    }

    // Coordinates of the quadrature points (2D).
    void coords(Array<OneD, NekDouble> &x, Array<OneD, NekDouble> &y) const
    {
        m_eq->UpdateFields()[0]->GetCoords(x, y);
    }

    // Time step for subsequent calls (periodic orbits: dt = T / N).
    void set_timestep(NekDouble dt)
    {
        m_eq->SetDt(dt);
    }

    // Continuation parameters.
    //  "Re"      sets Re and Kinvis = 1/Re (viscosity; see FlowMapVCS::SetKinvis)
    //  "Kinvis"  sets the viscosity directly
    //  other     a session parameter used in boundary-condition expressions
    //            (e.g. an inflow angle: u = cos(AoA), v = sin(AoA)); all BCs
    //            are re-evaluated. The parameter must appear in the
    //            expressions itself: parameters DERIVED from it in the
    //            session (P2 = 2*P1) were evaluated once at load and do not
    //            follow. Forcing and other functions are not re-evaluated.
    // Geometric parameters (gap ratio, ...) change the mesh and are not
    // handled here.
    void set_parameter(const std::string &name, double value)
    {
        if (name == "Re")
        {
            if (value != m_session->GetParameter("Re"))
            {
                m_session->SetParameter("Re", value);
            }
            m_eq->SetKinvis(1.0 / value);
            return;
        }
        if (name == "Kinvis")
        {
            m_eq->SetKinvis(value);
            return;
        }
        if (!m_session->DefinesParameter(name))
        {
            throw std::runtime_error("nekecs: unknown session parameter " + name);
        }
        if (value == m_session->GetParameter(name))
        {
            return;
        }
        m_session->SetParameter(name, value);
        m_eq->ReevaluateBCs();
    }

    double parameter(const std::string &name) const
    {
        if (name == "Kinvis")
        {
            return m_eq->Kinvis();
        }
        return m_session->GetParameter(name);
    }

    // Velocity component `comp` of a state at the point (px, py), by
    // spectral/hp interpolation (level 0 of a multi-level state). NaN if the
    // point is outside the mesh.
    double probe(const FieldVector &x, int comp, double px, double py) const
    {
        auto &f = m_eq->UpdateFields()[comp];
        Array<OneD, NekDouble> c(f->GetCoordim(0), 0.0); // must match coord dim
        c[0] = px;
        if (c.size() > 1)
        {
            c[1] = py;
        }
        try
        {
            return f->PhysEvaluate(c, x.u[comp]);
        }
        catch (const std::exception &)
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Multi-level state (u, u, ..) with all k levels equal to u.
    FieldVector with_history(const FieldVector &u) const
    {
        FieldVector x;
        const int k = order();
        for (int j = 0; j < k; ++j)
        {
            for (int i = 0; i < m_nvel; ++i)
            {
                Array<OneD, NekDouble> a(u.u[i].size());
                Vmath::Vcopy(a.size(), u.u[i], 1, a, 1);
                x.u.push_back(a);
            }
        }
        x.extra = u.extra;
        return x;
    }

    // Level 0 (u_n) of a state, without extras.
    FieldVector level0(const FieldVector &x) const
    {
        FieldVector u;
        for (int i = 0; i < m_nvel; ++i)
        {
            Array<OneD, NekDouble> a(x.u[i].size());
            Vmath::Vcopy(a.size(), x.u[i], 1, a, 1);
            u.u.push_back(a);
        }
        return u;
    }

    // out <- Phi_{nsteps dt}(u). `out` must have the shape of u.
    // Single-level u: start according to set_start_mode(). Multi-level u
    // (k * nvel arrays): StartMode::History, all levels advanced.
    void advance(const FieldVector &u, int nsteps, FieldVector &out)
    {
        const int nlev = static_cast<int>(u.u.size()) / m_nvel;
        set_state(u);
        if (m_reset_history)
        {
            m_eq->ResetHistory();
        }
        m_eq->SetTime(0.0);
        m_eq->SetInitialStep(0);

        FlowMapVCS::Levels prev;
        if (nlev > 1)
        {
            if (nlev != order())
            {
                throw std::runtime_error(
                    "nekecs: multi-level state must have IMEX-order levels");
            }
            auto &f = m_eq->UpdateFields();
            for (int j = 1; j < nlev; ++j)
            {
                Array<OneD, Array<OneD, NekDouble>> lev(m_nvel);
                for (int i = 0; i < m_nvel; ++i)
                {
                    // Same projection as level 0 (Dirichlet values imposed).
                    // ContField::FwdTrans takes the Dirichlet values from the
                    // output array, so seed it with the field's coefficients
                    // (set by set_state, which carry the BC data); a zero
                    // array would impose homogeneous Dirichlet data.
                    Array<OneD, NekDouble> c(f[i]->GetNcoeffs());
                    Vmath::Vcopy(c.size(), f[i]->GetCoeffs(), 1, c, 1);
                    lev[i] = Array<OneD, NekDouble>(f[i]->GetNpoints());
                    f[i]->FwdTrans(u.u[j * m_nvel + i], c);
                    f[i]->BwdTrans(c, lev[i]);
                }
                prev.push_back(lev);
            }
            m_eq->Integrate(nsteps, StartMode::History, &prev);
        }
        else if (m_use_dosolve)
        {
            m_eq->SetSteps(nsteps);
            m_eq->DoSolve();
        }
        else
        {
            m_eq->Integrate(nsteps, m_start);
        }
        ++m_calls;

        auto &f = m_eq->UpdateFields();
        for (int i = 0; i < m_nvel; ++i)
        {
            Vmath::Vcopy(out.u[i].size(), f[i]->GetPhys(), 1, out.u[i], 1);
        }
        for (int j = 1; j < nlev; ++j)
        {
            for (int i = 0; i < m_nvel; ++i)
            {
                Vmath::Vcopy(out.u[j * m_nvel + i].size(), prev[j - 1][i], 1,
                             out.u[j * m_nvel + i], 1);
            }
        }
    }

    // Write u to a .fld file. The pressure is taken from one time step started
    // at u (pressure is not part of the state, so this is the solver's
    // pressure for velocity ~u); the velocity written is u itself.
    void write(const FieldVector &u, const std::string &file)
    {
        FieldVector tmp = space().zeros_like(u);
        advance(u, 1, tmp);
        auto &f = m_eq->UpdateFields();
        for (int i = 0; i < m_nvel; ++i)
        {
            Vmath::Vcopy(f[i]->GetNpoints(), u.u[i], 1,
                                 f[i]->UpdatePhys(), 1);
            f[i]->FwdTrans(f[i]->GetPhys(), f[i]->UpdateCoeffs());
        }
        write_fld(file);
    }

    // Write a perturbation field (e.g. an eigenvector) to a .fld file. Unlike
    // write(), no boundary values are imposed (FwdTrans on a ContField would
    // put the base-flow Dirichlet data into it): elemental projection only,
    // pressure zero.
    void write_perturbation(const FieldVector &v, const std::string &file)
    {
        auto &f = m_eq->UpdateFields();
        for (int i = 0; i < m_nvel; ++i)
        {
            Vmath::Vcopy(f[i]->GetNpoints(), v.u[i], 1, f[i]->UpdatePhys(), 1);
            f[i]->FwdTransLocalElmt(f[i]->GetPhys(), f[i]->UpdateCoeffs());
        }
        auto &p = f[m_nvel];
        Vmath::Zero(p->GetNpoints(), p->UpdatePhys(), 1);
        Vmath::Zero(p->GetNcoeffs(), p->UpdateCoeffs(), 1);
        write_fld(file);
    }

    // Pseudo-random field of the shape of the state (deterministic seed).
    // Not admissible itself (ignores BCs); map it through the flow-map
    // Jacobian once to get an admissible start vector.
    FieldVector random_field(unsigned seed) const
    {
        FieldVector x = state();
        unsigned s    = seed * 2654435761u + 12345u;
        for (auto &a : x.u)
        {
            for (std::size_t j = 0; j < a.size(); ++j)
            {
                s    = s * 1664525u + 1013904223u; // LCG, portable
                a[j] = (s >> 8) * (1.0 / 16777216.0) - 0.5;
            }
        }
        return x;
    }

private:
    // Write the solver fields to `file`, replacing it. (Nektar++ would rename
    // an existing file to <name>.bak<N>.fld first, which leaves dozens of
    // backups of files rewritten at every Newton step.)
    void write_fld(const std::string &file)
    {
        auto comm = m_session->GetComm();
        if (comm->GetRank() == 0)
        {
            std::error_code ec;
            std::filesystem::remove_all(file, ec); // .fld is a directory in parallel
        }
        comm->Block();
        m_eq->WriteFld(file);
    }

    // Load velocity into the solver (projected onto the discrete space with
    // the Dirichlet values imposed) and zero the pressure.
    void set_state(const FieldVector &u)
    {
        auto &f = m_eq->UpdateFields();
        for (int i = 0; i < m_nvel; ++i)
        {
            Vmath::Vcopy(f[i]->GetNpoints(), u.u[i], 1,
                                 f[i]->UpdatePhys(), 1);
            f[i]->FwdTrans(f[i]->GetPhys(), f[i]->UpdateCoeffs());
            f[i]->BwdTrans(f[i]->GetCoeffs(), f[i]->UpdatePhys());
        }
        if (m_zero_pressure)
        {
            auto &p = f[m_nvel];
            Vmath::Zero(p->GetNpoints(), p->UpdatePhys(), 1);
            Vmath::Zero(p->GetNcoeffs(), p->UpdateCoeffs(), 1);
        }
    }

    Nektar::LibUtilities::SessionReaderSharedPtr m_session;
    std::shared_ptr<FlowMapVCS> m_eq;
    int m_nvel  = 0;
    int m_calls = 0;
    bool m_zero_pressure = true;
    bool m_reset_history = true;
    bool m_use_dosolve   = false;
    StartMode m_start    = StartMode::Standard;
};

} // namespace nektar
} // namespace nekecs
