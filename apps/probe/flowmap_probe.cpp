// Out-of-tree probe: drive Nektar++'s incompressible NS EquationSystem as a
// flow map u -> Phi_T(u) and check the two properties nekecs relies on:
//   1. it links and runs outside the Nektar++ source tree;
//   2. Phi_T is a deterministic function of u (same input -> same output,
//      i.e. no hidden multistep history leaks between calls).
//
// Usage: flowmap_probe [--nsteps N] mesh.xml session.xml

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SolverUtils/EquationSystem.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>

using namespace Nektar;

// The velocity-correction scheme keeps a multistep history for the
// high-order pressure BCs (Extrapolate::m_pressureHBCs, m_iprodnormvel,
// m_pressureCalls) that UnsteadySystem::DoSolve does not reset. Re-running
// GenerateHOPBCMap zeroes it, so every flow-map call starts identically.
class FlowMapVCS : public VelocityCorrectionScheme
{
public:
    FlowMapVCS(const LibUtilities::SessionReaderSharedPtr &s,
               const SpatialDomains::MeshGraphSharedPtr &g)
        : UnsteadySystem(s, g), VelocityCorrectionScheme(s, g)
    {
    }

    void ResetHistory()
    {
        m_extrapolation->GenerateHOPBCMap(m_session);
    }
};

using State = std::vector<Array<OneD, NekDouble>>;

static State GetState(SolverUtils::EquationSystemSharedPtr eq)
{
    State s;
    for (auto &f : eq->UpdateFields())
    {
        Array<OneD, NekDouble> copy(f->GetNpoints());
        Vmath::Vcopy(f->GetNpoints(), f->GetPhys(), 1, copy, 1);
        s.push_back(copy);
    }
    return s;
}

static void SetState(SolverUtils::EquationSystemSharedPtr eq, const State &s)
{
    auto &fields = eq->UpdateFields();
    for (size_t i = 0; i < fields.size(); ++i)
    {
        Vmath::Vcopy(fields[i]->GetNpoints(), s[i], 1, fields[i]->UpdatePhys(),
                     1);
        fields[i]->FwdTrans(fields[i]->GetPhys(), fields[i]->UpdateCoeffs());
        fields[i]->BwdTrans(fields[i]->GetCoeffs(), fields[i]->UpdatePhys());
    }
}

// L2 norm over the domain (quadrature-weighted), summed over the velocity
// components only (pressure is a slave variable, not part of the state).
static NekDouble L2(SolverUtils::EquationSystemSharedPtr eq, const State &s,
                    int nvel)
{
    auto &fields = eq->UpdateFields();
    NekDouble sum = 0.0;
    for (int i = 0; i < nvel; ++i)
    {
        Array<OneD, NekDouble> sq(fields[i]->GetNpoints());
        Vmath::Vmul(sq.size(), s[i], 1, s[i], 1, sq, 1);
        sum += fields[i]->Integral(sq);
    }
    return std::sqrt(sum);
}

static State Diff(const State &a, const State &b)
{
    State d;
    for (size_t i = 0; i < a.size(); ++i)
    {
        Array<OneD, NekDouble> x(a[i].size());
        Vmath::Vsub(x.size(), a[i], 1, b[i], 1, x, 1);
        d.push_back(x);
    }
    return d;
}

static State FlowMap(std::shared_ptr<FlowMapVCS> eq, const State &u0,
                     int nsteps, bool reset)
{
    if (reset)
    {
        eq->ResetHistory();
    }
    SetState(eq, u0);
    eq->SetTime(0.0);
    eq->SetInitialStep(0);
    eq->SetSteps(nsteps);
    eq->DoSolve();
    return GetState(eq);
}

int main(int argc, char *argv[])
{
    // Strip our own "--nsteps N" option; everything else goes to Nektar++.
    int nsteps = 50;
    for (int i = 1; i + 1 < argc; ++i)
    {
        if (std::string(argv[i]) == "--nsteps")
        {
            nsteps = std::atoi(argv[i + 1]);
            for (int j = i; j + 2 < argc; ++j)
            {
                argv[j] = argv[j + 2];
            }
            argc -= 2;
            break;
        }
    }

    auto session = LibUtilities::SessionReader::CreateInstance(argc, argv);
    auto graph   = SpatialDomains::MeshGraphIO::Read(session);
    // Construct our VelocityCorrectionScheme subclass directly (the factory
    // would build the stock class, which has no history reset).
    ASSERTL0(session->GetSolverInfo("SolverType") == "VelocityCorrectionScheme",
             "probe supports SolverType=VelocityCorrectionScheme only");
    auto eq = MemoryManager<FlowMapVCS>::AllocateSharedPtr(session, graph);
    std::static_pointer_cast<SolverUtils::EquationSystem>(eq)->InitObject();
    eq->DoInitialise();

    const int nvel = session->GetVariables().size() - 1; // u,v[,w],p
    const NekDouble dt = eq->GetTimeStep();
    const State u0 = GetState(eq);

    std::cout << "[probe] variables=" << session->GetVariables().size()
              << " nvel=" << nvel << " npts=" << u0[0].size()
              << " dt=" << dt << " nsteps=" << nsteps
              << " T=" << dt * nsteps << std::endl;

    // Without reset: second call inherits the first call's pressure-BC history.
    const State a1 = FlowMap(eq, u0, nsteps, false);
    const State a2 = FlowMap(eq, u0, nsteps, false);
    // With reset before every call.
    const State b1 = FlowMap(eq, u0, nsteps, true);
    const State b2 = FlowMap(eq, u0, nsteps, true);

    const NekDouble n0 = L2(eq, u0, nvel);
    const NekDouble res     = L2(eq, Diff(a1, u0), nvel) / n0;
    const NekDouble noreset = L2(eq, Diff(a2, a1), nvel) / n0;
    const NekDouble withres = L2(eq, Diff(b2, b1), nvel) / n0;
    const NekDouble fresh   = L2(eq, Diff(b1, a1), nvel) / n0;

    std::cout.precision(6);
    std::cout << std::scientific;
    std::cout << "[probe] ||u0||                        = " << n0 << "\n"
              << "[probe] ||Phi_T(u0) - u0|| / ||u0||   = " << res << "\n"
              << "[probe] repeat diff, no reset         = " << noreset << "\n"
              << "[probe] repeat diff, with reset       = " << withres << "\n"
              << "[probe] reset call vs first call      = " << fresh
              << std::endl;

    const bool ok = withres <= 1e-14 && fresh <= 1e-14;
    std::cout << (ok ? "[probe] PASS: flow map is deterministic with reset"
                     : "[probe] FAIL: flow map still depends on hidden state")
              << std::endl;

    session->Finalise();
    return ok ? 0 : 2;
}
