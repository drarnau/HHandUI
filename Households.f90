!===== SINGLES ====================================================================================
subroutine Singles(mysex)
  implicit none

  integer, intent(in) :: mysex
  character(len=2) :: identity

  if (mysex.eq.1) then
    identity = "sm"
  elseif (mysex.eq.2) then
    identity = "sf"
  end if

  call iniSingles(identity)
  call VFSingles()
  call SimSingles(mysex)
end subroutine Singles

!===== MARRIED ====================================================================================
subroutine MarriedHH()
  implicit none

  call iniMarried()
  call VFMarried()
  call SimMarried()

end subroutine MarriedHH
