
//Inputs
boxdim = 1;
gridsize = boxdim/10;

Point(1) = {0,0,0,gridsize};
Point(2) = {boxdim,0,0,gridsize};
Point(3) = {boxdim,boxdim,0,gridsize};
Point(4) = {0,boxdim,0,gridsize};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(9) = {5,6,7,8};
Plane Surface(10) = 9;

Transfinite Line{5,6,7,8} = boxdim/gridsize;
Transfinite Surface{10};
Recombine Surface{10};

//Transfinite Line{5,6,7,8} = boxdim/gridsize Using Bump 0.25;
//Transfinite Line{5,6} = boxdim/gridsize Using Progression 1.2;
//Transfinite Line{-7,-8} = boxdim/gridsize Using Progression 1.2;

//Generate 3D by extrusion
newEntities[]=
Extrude { 0,0,1 }
{
   Surface{10};
   Layers{boxdim/gridsize};
   Recombine;
};

Physical Surface("back") = {10};
Physical Surface("front")  = {newEntities[0]};
Physical Surface("bottom") = {newEntities[2]};
Physical Surface("right")  = {newEntities[3]};
Physical Surface("top")    = {newEntities[4]};
Physical Surface("left")   = {newEntities[5]};
Physical Volume(100) = {newEntities[1]};
