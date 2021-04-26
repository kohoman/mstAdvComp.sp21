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
Transfinite Line{5,6,7,8} = boxdim/gridsize;
Transfinite Surface{10};
Recombine Surface{10};
// ... generate 3D by extrusion ...
depth = 0.1;
newEntities[]=
Extrude { 0,0,depth }
{
   Surface{10};
   Layers{1};
   Recombine;
};
// ... face names ...
Physical Surface("movingWall") = {newEntities[4]};
Physical Surface("fixedWalls") = {newEntities[{2,3,5}]};
Physical Surface("frontAndBack") = {10,newEntities[0]};

// ... volume identification ...
Physical Volume(100) = {newEntities[1]};
