# ... my first Makefile ...

SHELL=/bin/bash

foo=simple

$(foo).pdf : $(foo).tex
	pdflatex $(foo).tex

.PHONY : clean help

clean :
	-rm -f $(foo).{aux,log,pdf}

help :
	@echo "make       -> compile $(foo).tex"
	@echo "make clean -> remove all generated files"
