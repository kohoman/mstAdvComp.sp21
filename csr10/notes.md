# Fortran, C++ and CodeBlocks (IDE)
---

## Fortran

Formula Translation -> Fortran

Fortran is the original numerically-intensive computing language

Fortran I - 1953
also II, II, IV, 66

Fortran 77 - ANSI Standard
           - fixed source form - defined column format
	   - standard file extension - .f
	   - assumed variable declaration by first letter: I-N integer
	   - accepted practice: implicit none
	   - comments - C or c in first column

Fortran 90 - introduced internal procedures
	   - free source form
	   - standard file extension - .f90
	   - comments - anything following the bang (!)

* Fortran is NOT case-sensitive

Compile command:
$ gfortran main.f90


## C++

Also a compiled language, like Fortran

Typical file extension - .cpp

Basic building blocks are functions

Comments - short comments: anything after //
         - C-style comments: anything between /* and */,
	   including multi-line content

Compile command:
$ g++ main.cpp
