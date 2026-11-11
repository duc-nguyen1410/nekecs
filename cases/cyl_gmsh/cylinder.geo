// Cylinder of diameter 1 in a rectangular box, quadrilateral mesh.
// Make the mesh with:   gmsh -2 cylinder.geo -o cylinder.msh
// then convert it:      NekMesh cylinder.msh mesh.xml
//
// Physical groups become Nektar++ composites with the same number, used in
// session.xml:
//   C[1] cylinder wall    C[2] inflow + sides    C[3] outflow    C[10] fluid
SetFactory("OpenCASCADE");

// Domain size (in cylinder diameters) and mesh sizes.
xin = -10; xout = 25; ylim = 10;
h_cyl = 0.12;   // element size on the cylinder
h_wake = 0.6;   // in the wake, at the outflow
h_far = 2.0;    // far field

Rectangle(1) = {xin, -ylim, 0, xout - xin, 2 * ylim};
Disk(2) = {0, 0, 0, 0.5};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };

// Sizes: fine on the cylinder, coarser downstream and far away.
Field[1] = Distance;      Field[1].CurvesList = {5};  // the circle (checked below)
Field[2] = Threshold;     Field[2].InField = 1;
Field[2].SizeMin = h_cyl; Field[2].SizeMax = h_far;
Field[2].DistMin = 0.2;   Field[2].DistMax = 12;
Field[3] = Box;           Field[3].VIn = h_wake; Field[3].VOut = h_far;
Field[3].XMin = 0; Field[3].XMax = xout; Field[3].YMin = -3; Field[3].YMax = 3;
Field[3].Thickness = 4;
Field[4] = Min;           Field[4].FieldsList = {2, 3};
Background Field = 4;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Quadrilaterals, second order (curved cylinder edges for Nektar++).
Mesh.Algorithm = 6;          // Frontal-Delaunay, then recombined into quads
Mesh.RecombineAll = 1;
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;  // put mid-edge nodes on the true circle
Mesh.MshFileVersion = 2.2;

// Boundaries, picked by position.
eps = 1e-3;
wall[]   = Curve In BoundingBox{-0.5 - eps, -0.5 - eps, -eps, 0.5 + eps, 0.5 + eps, eps};
inflow[] = Curve In BoundingBox{xin - eps, -ylim - eps, -eps, xin + eps, ylim + eps, eps};
bottom[] = Curve In BoundingBox{xin - eps, -ylim - eps, -eps, xout + eps, -ylim + eps, eps};
top[]    = Curve In BoundingBox{xin - eps, ylim - eps, -eps, xout + eps, ylim + eps, eps};
outflow[] = Curve In BoundingBox{xout - eps, -ylim - eps, -eps, xout + eps, ylim + eps, eps};
Field[1].CurvesList = {wall[]};

Physical Curve(1) = {wall[]};
Physical Curve(2) = {inflow[], bottom[], top[]};
Physical Curve(3) = {outflow[]};
Physical Surface(10) = {3};
