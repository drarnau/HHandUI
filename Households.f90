!===== SINGLES ====================================================================================
subroutine Singles(mysex)
  implicit none

  integer, intent(in) :: mysex

  call iniSingles(mysex)
  call VFSingles()
  call SimSingles(mysex,.FALSE.)
end subroutine Singles

!===== MARRIED ====================================================================================
subroutine MarriedHH()
  implicit none

  call iniMarried()
  call VFMarried()
  call SimMarried(.FALSE.)

end subroutine MarriedHH
