!===== SINGLES ====================================================================================
subroutine Singles(mysex)
  use Globals ! TO BE DELATED WHEN RESULTS NOT PRINTED
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

  ! Print results
  print *, "====================================="
  print *, "Sex", mysex, "Transitions:"
  print '(3f7.4)', transitions(0,mysex,1,:)
  print '(3f7.4)', transitions(0,mysex,2,:)
  print '(3f7.4)', transitions(0,mysex,3,:)
  print *,""
  print '(a,f7.4)', " Employment rate:", Erate(0,mysex)
  print '(a,f7.4)', " Unemployment rate:", Urate(0,mysex)
  print '(a,f7.4)', " Share OLF:", Nrate(0,mysex)
  print *, "====================================="
  print *, ""
end subroutine Singles

!===== MARRIED ====================================================================================
subroutine MarriedHH()
  implicit none

  call iniMarried()
  call VFMarried()
  call SimMarried()

end subroutine MarriedHH
