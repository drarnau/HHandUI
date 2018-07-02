program calibration
  use Globals
  use Utils
  implicit none
  
  real(8) :: t_start, t_finish

! INITIALISATION: Read parameters from outside world, compute grids, tranistion matrices, and initialise equilibrium variables
  call cpu_time(t_start)

  call initialisation(shock_z, shock_q, shock_g)

  call cpu_time(t_finish)
  print *, "Initialisation time:"
  call mytime(t_finish-t_start)

! COMPUTE DECISION RULES
  call cpu_time(t_start)

  call VFiteration()

  call cpu_time(t_finish)
  print *, "Time to compute decision rules:"
  call mytime(t_finish-t_start)

! SIMULATE THE ECONOMY
  call cpu_time(t_start)




  call cpu_time(t_finish)
  print *, "Time to simulate:"
  call mytime(t_finish-t_start)

end program calibration
