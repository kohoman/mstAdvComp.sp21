# ... GNU Make for Bourne Shell ...
SHELL = /bin/bash
# ... variable definitions ...
ENV  = /opt/openfoam6/etc/bashrc
SLVR = icoFoam
GEO  = case.geo
MSH  = case.msh
CDIR = $(notdir $(shell pwd))
PYC  =
# ... OpenFoam initial condition dir ...
startDir = 0
# ... target rules and commands ...
.PHONY : clean help
help :
	@echo "make       -> help"
	@echo "make obm   -> execute OpenFOAM blockMesh"
	@echo "make gmo   -> convert gmsh to OpenFOAM"
	@echo "make run   -> execute $(SLVR) solver"
	@echo "make pps   -> execute post-processing ops"
	@echo "make all   -> execute msh and run"
	@echo "make clean -> delete msh, run and all"

clean :
	-rm -f obm gmo run pps all

all :
	make obm
	make run
	make pps
	touch all

obm : system/blockMeshDict
	source $(ENV) && blockMesh
	touch obm

gmo :
	gmsh $(GEO) -3 -o $(MSH)
	source $(ENV) && gmshToFoam $(MSH) -case ../$(CDIR)	
	touch gmo

run : obm $(startDir)/* constant/* system/*
	source $(ENV) && $(SLVR)
	touch run

pps : run
	source $(ENV) && postProcess -func CourantNo
	source $(ENV) && pisoFoam -postProcess -func wallShearStress
	touch pps
	touch run
