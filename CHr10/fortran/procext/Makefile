SHELL=/bin/sh

SOURCE = main.f90
EXBL   = main.exe
F90    = /usr/bin/gfortran

.PHONY : clean help

exe : $(EXBL)

$(EXBL) : $(SOURCE)
	$(F90) -o $(EXBL) $(SOURCE)

run : $(EXBL)
	./$(EXBL)
	touch run

all : exe run
	make exe
	make run

clean :
	-rm -f exe run $(EXBL)

help :
	@echo "make       -> compile source"
	@echo "make run   -> run executable"
	@echo "make clean -> remove executable"
