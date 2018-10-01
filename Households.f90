!===== SINGLES ====================================================================================
subroutine Singles(mysex)
  implicit none

  integer, intent(in) :: mysex

  call iniSingles(mysex)
  call VFSingles()
  call SimSingles(mysex,.false.)
end subroutine Singles

!===== MARRIED ====================================================================================
subroutine MarriedHH()
  implicit none

  call iniMarried()
  call VFMarried()
  call SimMarried(.false.)

end subroutine MarriedHH
