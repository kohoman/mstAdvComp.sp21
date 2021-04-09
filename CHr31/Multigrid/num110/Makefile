# ... GNU Make for Bourne Shell ...
SHELL = /bin/sh
# ... variable declarations ...
SOURCE = main.f
EXBL   = main.exe
INCS   = incs.com
INPS   = main.in
LFLAGS = 
LIBS   =
F77    = /usr/bin/gfortran
# ... compile rules ...
.PHONY : exe clean help
exe : $(EXBL)
$(EXBL): $(SOURCE) $(INCS) Makefile
	$(F77) $(LFLAGS) -o $(EXBL) $(SOURCE) $(LIBS)
run : $(EXBL) $(INPS)
	-rm -f *.out
	./$(EXBL)
	touch run
clean :
	-rm -f exe run *.exe *.txt *.out
help :
	@echo "make     -> compile source"
	@echo "make exe -> compile source"
	@echo "make run -> compile and execute"
	@echo "make clean -> rm generated files"
