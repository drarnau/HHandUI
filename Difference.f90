subroutine Difference()
  use Globals
  implicit none

  integer :: myms, mysex, td, tw
  real(8) :: dint_rate, dmarsing_ratio, marsing_ratio
  real(8), dimension(0:1,1:2) :: dErate, dUrate
  real(8), dimension(0:1,1:2,1:3,1:3) :: dtrans

  ! Aggregate data
  dint_rate = 0.00327                       ! Interest rate
  dmarsing_ratio = 2.8481675393d0           ! Married-single wealth ratio

  ! Load data from txt file
  open(unit = 11, file = "data.txt")
  do myms = 0, 1
  do mysex = 1, 2
    do td = 1, 3
    do tw = 1, 3
      read (11, *) dtrans(myms,mysex,td,tw)
    end do
    end do
    read(11,*) dErate(myms,mysex)
    read(11,*) dUrate(myms,mysex)

    ! Adjust transitions to sum up to one
    ! EU Transition
    dtrans(myms,mysex,1,2) = 1.d0 -(dtrans(myms,mysex,1,1)+dtrans(myms,mysex,1,3))
    ! UN Transition
    dtrans(myms,mysex,2,3) = 1.d0 -(dtrans(myms,mysex,2,1)+dtrans(myms,mysex,2,2))
    ! NU Transition
    dtrans(myms,mysex,3,2) = 1.d0 -(dtrans(myms,mysex,3,1)+dtrans(myms,mysex,3,3))
  end do
  end do
  close(11)

  ! Model Married-single wealth ratio
  marsing_ratio = (weights(married,male)+weights(married,female))*wealth(male+female) &
                /((wealth(male)*weights(single,male))+(wealth(female)*weights(single,female)))

  ! Print results to screen
    print *, "*************************************************************************"
    print '(a,f10.7,a,f10.7)', " Interest rate (Data vs. Model):", dint_rate, "    ", int_rate
    print '(a,f10.7,a,f10.7)', " Married-single wealth ratio (Data vs. Model):", &
                                dmarsing_ratio, "   ", marsing_ratio
    print *, "*************************************************************************"
    do myms = 0, 1
    do mysex = 1, 2
      ! Screen
      print *, "========================================================================="
      print *, "Marital Status:", myms
      print *, "Sex", mysex
      print *, "Transitions (Data vs. Model):"
      print '(3f7.4,a,3f7.4)',dtrans(myms,mysex,1,:), "   ", transitions(myms,mysex,1,:)
      print '(3f7.4,a,3f7.4)',dtrans(myms,mysex,2,:), "   ", transitions(myms,mysex,2,:)
      print '(3f7.4,a,3f7.4)',dtrans(myms,mysex,3,:), "   ", transitions(myms,mysex,3,:)
      print *,""
      print '(a,f7.4,a,f7.4)', " Employment rate:", dErate(myms,mysex), "   ", Erate(myms,mysex)
      print '(a,f7.4,a,f7.4)', " Unemployment rate:", dUrate(myms,mysex), "   ", Urate(myms,mysex)
      print *, "========================================================================="
      print *, ""
    end do
    end do

  ! Print results to txt file
  open(unit = 12, file = "model.txt")
  do myms = 0, 1
  do mysex = 1, 2
    do td = 1, 3
    do tw = 1, 3
      write (12, *) transitions(myms,mysex,td,tw)
    end do
    end do
    write(12,*) Erate(myms,mysex)
    write(12,*) Urate(myms,mysex)
  end do
  end do
  close(12)

  ! Print aggregate variables
  open(unit = 13, file = "aggvars.txt")
  write(13,*) int_rate
  write(13,*) tau
  write(13,*) KLratio
  write(13,*) average_z
  write(13,*) wage
  close(13)

end subroutine Difference
