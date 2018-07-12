program calibration
  use Globals
  use Utils
  implicit none

  integer, parameter :: maxIter = 1
  integer :: iter
  real(8), parameter :: adj_KL = 0.5d0, adj_T = 0.5d0, adj_avgz = 0.5d0, &
                        tol_KL = 0.01, tol_T = 1.0d-4, tol_avgz = 0.001
  real(8) :: start, finish, t_start, t_finish, error_KL, error_T, error_avgz

  call cpu_time(start)
! INITIALISATION: Read global parameters from outside world and initialise equilibrium variables
  call cpu_time(t_start)

  call initialisation()

  call cpu_time(t_finish)
  print *, "Initialisation time:"
  call mytime(t_finish-t_start)

  do iter = 1, maxIter
  ! SOLVE THE MODEL FOR SINGLE MALES
    call cpu_time(t_start)

    call Singles(1)

    call cpu_time(t_finish)
    print *, "Time to solve single males:"
    call mytime(t_finish-t_start)

  ! SOLVE THE MODEL FOR SINGLE FEMALES
    call cpu_time(t_start)

    ! call Singles(2)

    call cpu_time(t_finish)
    print *, "Time to solve single females:"
    call mytime(t_finish-t_start)

  ! SOLVE THE MODEL FOR MARRIED HOUSEHOLDS
    call cpu_time(t_start)

    call MarriedHH()

    call cpu_time(t_finish)
    print *, "Time to solve married households:"
    call mytime(t_finish-t_start)

  ! ! COMPUTE DECISION RULES
  !   call cpu_time(t_start)
  !
  !   call VFiteration()
  !
  !   call cpu_time(t_finish)
  !   print *, "Time to compute decision rules:"
  !   call mytime(t_finish-t_start)
  !
  ! ! SIMULATE THE ECONOMY
  !   call cpu_time(t_start)
  !
  !   call simulation()
  !
  !   call cpu_time(t_finish)
  !   print *, "Time to simulate:"
  !   call mytime(t_finish-t_start)
  !
  ! ! CHECK EQUILIBRIUM IS REACHED
  !   ! Compute errors
  !   error_KL = abs(KLratio-new_KLratio)
  !   error_T = abs(T-new_T)
  !   error_avgz = abs(average_z-new_average_z)
  !
  !   ! Print current situation
  !   print *, "Current equilibirum errors, at iteration:", iter
  !   print '(a,3f9.4)', "  KL ratio (old new error): ", KLratio, new_KLratio, error_KL
  !   print '(a,3f9.4)', "  T (old new error):        ", T, new_T, error_T
  !   print '(a,3f9.4)', "  Average z (old new error):", average_z, new_average_z, error_avgz
  !   print *, ""
  !
  !   ! Update equilibrium values
  !   KLratio = adj_KL*new_KLratio + (1.d0-adj_KL)*KLratio
  !   T = adj_T*new_T + (1.d0-adj_T)*T
  !   average_z = adj_avgz*new_average_z + (1.d0-adj_avgz)*average_z
  !
  !   ! Check errors and tolerance
  !   if ((error_KL.lt.tol_KL).and.(error_T.lt.tol_T).and.(error_avgz.lt.tol_avgz)) then
  !     print *, "Equilibirum reached at iteration:", iter
  !     print *, ""
  !     exit
  !   end if
  end do

! PRINT RESULTS


  call cpu_time(finish)
  print *, "Total time to solve the model:"
  call mytime(finish-start)
end program calibration
