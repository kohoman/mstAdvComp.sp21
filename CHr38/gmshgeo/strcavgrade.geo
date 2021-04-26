<str2dboxuniform>>
Transfinite Line{5,6,7,8} = boxdim/gridsize Using Bump 0.25;
Transfinite Line{5,6} = boxdim/gridsize Using Progression 1.2;
Transfinite Line{-7,-8} = boxdim/gridsize Using Progression 1.2;
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
