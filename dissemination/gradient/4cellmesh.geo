// Gmsh project created on Tue Jul 14 14:56:25 2026
//+
Point(1) = {0, 0, 0, 10.0};
//+
Point(2) = {7, 0, 0, 10.0};
//+
Point(3) = {3, 7, 0, 10.0};
//+
Point(4) = {-2, 9, 0, 10.0};
//+
Point(5) = {9, 2, 0, 10.0};
//+
Point(6) = {5, -6, 0, 10.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {1, 3};
//+
Line(4) = {2, 5};
//+
Line(5) = {5, 3};
//+
Line(6) = {3, 4};
//+
Line(7) = {4, 1};
//+
Line(8) = {6, 1};
//+
Line(9) = {6, 2};
//+
Curve Loop(1) = {7, 3, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -2, -1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, -2, 4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {1, -9, 8};
//+
Plane Surface(4) = {4};
