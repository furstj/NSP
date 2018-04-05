//+
SetFactory("Built-in");

h  = 0.8;
hf = h/8;
a  = 2;
b  = 0.5;

//+
Point(1) = {-5, -5, 0, h};
//+
Point(2) = {15, -5, 0, h};
//+
Point(3) = {15, 5, 0, h};
//+
Point(4) = {-5, 5, 0, h};
//+
Point(5) = {a, 0, 0, hf/2};
//+
Point(6) = {0, b, 0, hf};
//+
Point(7) = {0, -b, 0, hf};
//+
Point(8) = {-a, 0, 0, hf/2};
//+
Point(9) = {0, 0, 0, hf};
//+
Ellipse(1) = {8, 9, 8, 6};
//+
Ellipse(2) = {5, 9, 5, 6};
//+
Ellipse(3) = {8, 9, 8, 7};
//+
Ellipse(4) = {5, 9, 5, 7};
//+
Line(5) = {1, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 4};
//+
Line(8) = {4, 1};
//+
Line Loop(1) = {8, 5, 6, 7};
//+
Line Loop(2) = {3, -4, 2, -1};
//+
Plane Surface(1) = {1, 2};
//+
Physical Line("INLET") = {8};
//+
Physical Line("OUTLET") = {6};
//+
Physical Line("BOTTOM") = {5};
//+
Physical Line("TOP") = {7};
//+
Physical Line("PROFILE") = {3, 4, 2, 1};

Physical Surface("DOMAIN") = {1};
//+
//+
Field[1] = BoundaryLayer;
//+
Field[1].EdgesList = {1, 2, 3, 4};
//+
Field[1].hwall_n = 0.01;
//+
Field[1].thickness = 0.25;
//+
BoundaryLayer Field = 1;
