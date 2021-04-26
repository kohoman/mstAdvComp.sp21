// ... Inputs ...
boxdim = 1;
gridsize = boxdim/10;

// ... original Vertices ...
//Point(1) = {0,0,0,gridsize};
//Point(2) = {boxdim,0,0,gridsize};
//Point(3) = {boxdim,boxdim,0,gridsize};
//Point(4) = {0,boxdim,0,gridsize};

// ... revised Vertices ...
Point(1) = {0,0,0};
Point(2) = {boxdim,0,0};
Point(3) = {boxdim,boxdim,0};
Point(4) = {0,boxdim,0};


// ... Edges ...
Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

// ... Block ...
Line Loop(9) = {5,6,7,8};
Plane Surface(10) = 9;
Recombine Surface{10};
