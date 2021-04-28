Class Hour 39
-----

## geo files in OpenFOAM compatible format

located in gmshgeo.alt/ directory

Specific files (all 3D):
unscavtri.geo    - unstructured, tri mesh
unscavquad.geo   - unstructured, quad mesh
strcavuniform.geo
strcavgrade.geo

## Example cases based on blockMesh

cpm000 - uniform grid via OpenFOAM blockMesh
------
OpenFOAM solution for 2D driven cavity (laminar)

0/U - velocity vector bcs
  p - pressure field bcs

constant/transportProperties - fluid props

system/blockMeshDict - block mesh parameters
       controlDict   - stepping controls
       fvSchemes     - discretization schemes
       fvSolution   - solution algorithms/params

Execution steps:
$ make obm
$ make run


cpm010 - non-uniform structured grid
------   via OpenFOAM blockMesh

Execution steps: as for cpm000.

Notice output in the time-level subdirs and
in the postProcessing subdir.

cpm307 - writes point data to file
------

cpm308 - writes data along line to file
------


## Example cases based on gmsh

cpm320 - gmsh based unstructured tri mesh
------

Execution steps:
$ make gmo
$ make run


cpm321 - gmsh based unstructured quad mesh
------

cpm322 - structured uniform mesh
------

cpm323 - structured non-uniform mesh
------

