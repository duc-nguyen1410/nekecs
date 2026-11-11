// shed: DNS from a perturbed equilibrium to a saturated limit cycle, to make
// a periodic-orbit initial guess.
//
// Usage: shed mesh.xml session.xml   (IC: the equilibrium)
//
// Session parameters (defaults in brackets):
//   SHED_Time     total simulated time                     [300]
//   SHED_Amp      perturbation amplitude / |u|             [1e-3]
//   SHED_ProbeX, SHED_ProbeY   probe location for v         [2, 0]
//   SHED_Chunk    steps per flow-map call                  [10]
//
// Time stepping uses the two-level (History) state, so chunked calls are
// exactly continuous IMEX-2 stepping. Writes <session>_probe.dat (t, v) and
// <session>_shed.fld (final state), and prints the periods between upward
// zero crossings of v at the probe.
#include <nekecs/nektar/flowmap.hpp>
#include <nekecs/space.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using namespace nekecs;
using namespace nekecs::nektar;

int main(int argc, char *argv[])
{
    auto session = Nektar::LibUtilities::SessionReader::CreateInstance(argc, argv);
    auto graph   = Nektar::SpatialDomains::MeshGraphIO::Read(session);
    NekDouble t_end, amp, px, py;
    int chunk;
    session->LoadParameter("SHED_Time", t_end, 300.0);
    session->LoadParameter("SHED_Amp", amp, 1e-3);
    session->LoadParameter("SHED_ProbeX", px, 2.0);
    session->LoadParameter("SHED_ProbeY", py, 0.0);
    session->LoadParameter("SHED_Chunk", chunk, 10);
    std::string name = session->GetFilenames().back();
    name = name.substr(name.find_last_of('/') + 1);
    name = name.substr(0, name.find('.'));

    VCSFlowMap fm(session, graph);
    const NektarSpace S = fm.space();
    const double dt     = fm.timestep();

    // Probe: quadrature point nearest (px, py).
    std::size_t probe = 0;
    {
        const int np = fm.state().u[0].size();
        Array<OneD, NekDouble> X(np), Y(np);
        fm.coords(X, Y);
        double best = std::numeric_limits<double>::infinity();
        for (int i = 0; i < np; ++i)
        {
            const double d = (X[i] - px) * (X[i] - px) + (Y[i] - py) * (Y[i] - py);
            if (d < best)
            {
                best  = d;
                probe = i;
            }
        }
        std::cout << "[shed] probe at (" << X[probe] << ", " << Y[probe] << ")"
                  << std::endl;
    }

    FieldVector u = fm.state();
    FieldVector w = fm.random_field(7);
    S.scale(amp * norm(S, u) / norm(S, w), w);
    S.axpy(1.0, u, w);
    FieldVector x = fm.with_history(w);
    FieldVector y = fm.space().zeros_like(x);

    std::ofstream out(name + "_probe.dat");
    out << std::setprecision(10);
    std::vector<double> crossings;
    double t = 0.0, vprev = x.u[1][probe];
    const int nchunks = static_cast<int>(std::ceil(t_end / (chunk * dt)));
    for (int c = 0; c < nchunks; ++c)
    {
        fm.advance(x, chunk, y);
        std::swap(x, y);
        t += chunk * dt;
        const double v = x.u[1][probe];
        out << t << " " << v << "\n";
        if (vprev < 0.0 && v >= 0.0)
        {
            // linear interpolation of the crossing time
            crossings.push_back(t - chunk * dt * v / (v - vprev));
        }
        vprev = v;
    }
    out.close();

    std::cout << std::setprecision(8) << "[shed] t = " << t << ", "
              << crossings.size() << " upward crossings of v at the probe\n";
    const std::size_t n = crossings.size();
    for (std::size_t i = (n > 12 ? n - 12 : 1); i < n; ++i)
    {
        std::cout << "[shed]   period " << i << ": "
                  << crossings[i] - crossings[i - 1] << "\n";
    }
    if (n >= 6)
    {
        double s = 0.0;
        for (std::size_t i = n - 5; i < n; ++i)
        {
            s += crossings[i] - crossings[i - 1];
        }
        std::cout << "[shed] mean of last 5 periods: T = " << s / 5.0
                  << ", St = " << 5.0 / s << "\n";
    }
    fm.write(fm.level0(x), name + "_shed.fld");
    std::cout << "[shed] wrote " << name << "_shed.fld" << std::endl;
    session->Finalise();
    return 0;
}
