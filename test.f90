program test

  implicit none

  real(8), dimension(50,3,4,0:1) :: prova

  prova(1,1,1,0) = 3.d0

CONTAINS
  subroutine mira(vec)
    implicit none
    real(8), dimension(50), intent(in) :: vec

    print *, vec


  end subroutine mira


end program test
