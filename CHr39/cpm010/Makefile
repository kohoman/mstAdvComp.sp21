# ... GNU Make for Bourne Shell ...
SHELL = /bin/bash
# ... variable definitions ...
ENV  = /opt/openfoam5/etc/bashrc
SLVR = icoFoam
GEO  = case.geo
MSH  = case.msh
CDIR = $(notdir $(shell pwd))
PYC  =
# ... target rules and commands ...
.PHONY : clean help
help :
	@echo "make       -> help"
	@echo "make obm   -> execute OpenFOAM blockMesh"
	@echo "make gmo   -> convert gmsh to OpenFOAM"
	@echo "make run   -> execute $(SLVR) solver"
	@echo "make all   -> execute msh and run"
	@echo "make clean -> delete msh, run and all"
obm :
	source $(ENV) && blockMesh
	touch obm

gmo :
	gmsh $(GEO) -3 -o $(MSH)
	source $(ENV) && gmshToFoam $(MSH) -case ../$(CDIR)	
	touch gmo

run :
	source $(ENV) && $(SLVR)
	touch run
all :
	make msh
	make run
	touch all
clean :
	-rm -f msh run all
