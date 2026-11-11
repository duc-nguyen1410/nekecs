// Check that changing a parameter in place (VCSFlowMap::set_parameter) gives
// the same flow map as a solver built from scratch at that value.
//
// Usage: diag_param mesh.xml session.xml NAME VALUE
//   NAME = Re      in-place SetKinvis vs a fresh session with -P Kinvis=1/VALUE
//   NAME = other   in-place BC re-evaluation vs a fresh session with -P NAME=VALUE
//
// For each: Phi (two-level state, ECS_Steps steps) from the session IC after
// the change, and after changing back to the original value (must reproduce
// the original map exactly).
#include <nekecs/nektar/flowmap.hpp>
#include <nekecs/space.hpp>

#include <iomanip>
#include <iostream>
#include <sstream>

using namespace nekecs;
using namespace nekecs::nektar;

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        std::cerr << "usage: diag_param mesh.xml session.xml NAME VALUE\n";
        return 2;
    }
    const std::string pname = argv[3];
    const double value      = std::stod(argv[4]);

    // Session 1: as in the file.
    char *args1[] = {argv[0], argv[1], argv[2], nullptr};
    auto s1 = Nektar::LibUtilities::SessionReader::CreateInstance(3, args1);
    auto g1 = Nektar::SpatialDomains::MeshGraphIO::Read(s1);
    VCSFlowMap fm1(s1, g1);
    int nsteps;
    s1->LoadParameter("ECS_Steps", nsteps, 100);
    const NektarSpace S1 = fm1.space();
    const FieldVector x  = fm1.with_history(fm1.state());
    const double orig    = fm1.parameter(pname);

    FieldVector a0 = S1.zeros_like(x), a1 = S1.zeros_like(x),
                a2 = S1.zeros_like(x);
    fm1.advance(x, nsteps, a0);           // original value
    fm1.set_parameter(pname, value);
    fm1.advance(x, nsteps, a1);           // changed in place
    fm1.set_parameter(pname, orig);
    fm1.advance(x, nsteps, a2);           // changed back

    // Session 2: fresh solver at the new value.
    std::ostringstream ov;
    ov << std::setprecision(17)
       << (pname == "Re" ? "Kinvis" : pname) << "="
       << (pname == "Re" ? 1.0 / value : value);
    const std::string ovs = ov.str();
    std::string p2        = "-P";
    char *args2[] = {argv[0], &p2[0], const_cast<char *>(ovs.c_str()), argv[1],
                     argv[2], nullptr};
    auto s2 = Nektar::LibUtilities::SessionReader::CreateInstance(5, args2);
    auto g2 = Nektar::SpatialDomains::MeshGraphIO::Read(s2);
    VCSFlowMap fm2(s2, g2);
    const NektarSpace S2 = fm2.space();
    FieldVector b1       = S2.zeros_like(x);
    fm2.advance(x, nsteps, b1); // same input state (same mesh, same layout)

    auto rel = [&](const FieldVector &p, const FieldVector &q) {
        FieldVector d = clone(S1, p);
        S1.axpy(-1.0, q, d);
        return norm(S1, d) / norm(S1, q);
    };
    std::cout << std::scientific << std::setprecision(3)
              << "[diag_param] " << pname << ": " << orig << " -> " << value
              << " (fresh session: -P " << ovs << "), " << nsteps
              << " steps, two-level state\n"
              << "[diag_param] |Phi_new - Phi_orig| / |Phi_orig|        = "
              << rel(a1, a0) << "   (size of the parameter effect)\n"
              << "[diag_param] |Phi_inplace - Phi_fresh| / |Phi_fresh| = "
              << rel(a1, b1) << "\n"
              << "[diag_param] |Phi_back - Phi_orig| / |Phi_orig|       = "
              << rel(a2, a0) << std::endl;
    s1->Finalise();
    return 0;
}
