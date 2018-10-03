program calibration
  use Globals
  use Utils
  implicit none

  integer, parameter :: maxIter = 100
  integer :: iter, myms, mysex
  real(8), parameter :: adj_KL = 0.5d0, adj_T = 0.5d0, adj_avgz = 0.5d0, &
                        tol_KL = 0.1 , tol_T = 0.1, tol_avgz = 0.1
                        ! tol_KL = 0.01 , tol_T = 1.0d-4, tol_avgz = 0.001
  real(8) :: start, finish, t_start, t_finish, error_KL, error_T, error_avgz
  character(len=1024) :: filename

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

    call Singles(2)

    call cpu_time(t_finish)
    print *, "Time to solve single females:"
    call mytime(t_finish-t_start)

  ! SOLVE THE MODEL FOR MARRIED HOUSEHOLDS
    call cpu_time(t_start)

    call MarriedHH()

    call cpu_time(t_finish)
    print *, "Time to solve married households:"
    call mytime(t_finish-t_start)

  ! CHECK EQUILIBRIUM IS REACHED
    ! Compute aggregates
    new_KLratio = aggregate3(aux_KLratio)
    new_T = aggregate3(aux_T)
    new_average_z = aggregate3(aux_average_z)

    ! Compute errors
    error_KL = abs(KLratio-new_KLratio)
    error_T = abs(T-new_T)
    error_avgz = abs(average_z-new_average_z)

    ! Print current situation
    print *, "Current equilibirum errors, at iteration:", iter
    print '(a,3f9.4)', "  KL ratio (old new error): ", KLratio, new_KLratio, error_KL
    print '(a,3f9.4)', "  T (old new error):        ", T, new_T, error_T
    print '(a,3f9.4)', "  Average z (old new error):", average_z, new_average_z, error_avgz
    print *, ""

    ! Update equilibrium values
    KLratio = adj_KL*new_KLratio + (1.d0-adj_KL)*KLratio
    T = adj_T*new_T + (1.d0-adj_T)*T
    average_z = adj_avgz*new_average_z + (1.d0-adj_avgz)*average_z
    call Prices(KLratio, int_rate, wage)

    ! Check errors and tolerance
    if ((error_KL.lt.tol_KL).and.(error_T.lt.tol_T).and.(error_avgz.lt.tol_avgz)) then
      print *, "Equilibirum reached at iteration:", iter
      print *, ""
      exit
    end if
  end do

! PRINT RESULTS
  do myms = 0, 1
  do mysex = 1, 2
    ! Screen
    print *, "====================================="
    print *, "Marital Status:", myms
    print *, "Sex", mysex
    print *, "Transitions:"
    print '(3f7.4)', transitions(myms,mysex,1,:)
    print '(3f7.4)', transitions(myms,mysex,2,:)
    print '(3f7.4)', transitions(myms,mysex,3,:)
    print *,""
    print '(a,f7.4)', " Employment rate:", Erate(myms,mysex)
    print '(a,f7.4)', " Unemployment rate:", Urate(myms,mysex)
    print '(a,f7.4)', " Share OLF:", Nrate(myms,mysex)
    print *, "====================================="
    print *, ""

    ! Txt file
    write (filename, "(A10,I1,A4,I1,A4)") "results_ms", myms, "_sex", mysex, ".txt"
    open(unit=12, file=trim(filename))
    write(12,'(f10.7,a)') Erate(myms,mysex), "      ! Employment rate"
    write(12,'(f10.7,a)') Urate(myms,mysex), "      ! Unmployment rate"
    write(12,'(f10.7,a)') Nrate(myms,mysex), "      ! Share OLF"
    write(12,'(f10.7,a)') transitions(myms,mysex,1,1), "      ! EE Transition"
    write(12,'(f10.7,a)') transitions(myms,mysex,1,2), "      ! EU Transition"
    write(12,'(f10.7,a)') transitions(myms,mysex,1,3), "      ! EN Transition"
    write(12,'(f10.7,a)') transitions(myms,mysex,2,1), "      ! UE Transition"
    write(12,'(f10.7,a)') transitions(myms,mysex,2,2), "      ! UU Transition"
    write(12,'(f10.7,a)') transitions(myms,mysex,2,3), "      ! UN Transition"
    write(12,'(f10.7,a)') transitions(myms,mysex,3,1), "      ! NE Transition"
    write(12,'(f10.7,a)') transitions(myms,mysex,3,2), "      ! NU Transition"
    write(12,'(f10.7,a)') transitions(myms,mysex,3,3), "      ! NN Transition"
    close(12)
  end do
  end do


  call cpu_time(finish)
  print *, "Total time to solve the model:"
  call mytime(finish-start)
end program calibration
