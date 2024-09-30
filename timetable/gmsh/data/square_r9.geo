/******************************************************************************
* Copyright (C) Siarhei Uzunbajakau, 2023.
*
* This program is free software. You can use, modify, and redistribute it under
* the terms of the GNU Lesser General Public License as published by the Free
* Software Foundation, either version 3 or (at your option) any later version.
* This program is distributed without any warranty.
*
* Refer to COPYING.LESSER for more details.
******************************************************************************/

r = 9; //@1

a = 1.0;

Point(1) = {-a,-a, 0};
Point(2) = {-a, a, 0};
Point(3) = { a, a, 0};
Point(4) = { a,-a, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface(100) = {1};

Physical Line(1) = {2, 3};
Physical Line(2) = {1, 4};

// Physical Line(1) = {1, 2, 3, 4};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = 2*r - 1;

