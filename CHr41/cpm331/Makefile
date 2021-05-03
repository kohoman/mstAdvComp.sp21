# ... GNU Make for Bourne Shell ...
SHELL = /bin/bash
# ... variable definitions ...
ENV  = /opt/openfoam8/etc/bashrc
SLVR = pimpleFoam
PYC  =
# ... target rules and commands ...
.PHONY : clean help
help :
	@echo "make       -> help"
	@echo "make msh   -> execute blockMesh"
	@echo "make run   -> execute $(SLVR)"
	@echo "make all   -> execute msh and run"
	@echo "make clean -> delete msh, run and all"
msh :
	source $(ENV) && blockMesh
	touch msh
run :
	source $(ENV) && $(SLVR)
	touch run
all :
	make msh
	make run
	touch all
clean :
	-rm -f msh run all
