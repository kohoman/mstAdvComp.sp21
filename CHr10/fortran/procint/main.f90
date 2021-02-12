! fortran 90 - free source format
program main
  call sub()
contains
  subroutine sub()
    implicit none
    write(*,*) "My first fortran subroutine"
    return
  end subroutine sub      
end program main
