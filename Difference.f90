subroutine Difference()
  use Globals
  implicit none

  integer :: myms, mysex
  real(8) :: dint_rate, dmarsing_ratio, mydiff, marsing_ratio
  real(8), dimension(0:1,1:2) :: dErate, dUrate
  real(8), dimension(0:1,1:2,1:3,1:3) :: dtrans

  ! Aggregate data
  dint_rate = 0.00327                       ! Interest rate
  dmarsing_ratio = 2.8481675393d0           ! Married-single wealth ratio

  ! Data single females
  dtrans(single,female,1,1) = 0.9681d0      ! EE Transition
  dtrans(single,female,1,3) = 0.0209d0      ! EN Transition

  dtrans(single,female,2,1) = 0.2472d0      ! UE Transition
  dtrans(single,female,2,2) = 0.5361d0      ! UU Transition

  dtrans(single,female,3,1) = 0.0310d0      ! NE Transition
  dtrans(single,female,3,3) = 0.9516d0      ! NN Transition


  dtrans(single,female,1,2) = &             ! EU Transition
      1.d0 -(dtrans(single,female,1,1)+dtrans(single,female,1,3))
  dtrans(single,female,2,3) = &             ! UN Transition
      1.d0 -(dtrans(single,female,2,1)+dtrans(single,female,2,2))
  dtrans(single,female,3,2) = &             ! NU Transition
      1.d0 -(dtrans(single,female,3,1)+dtrans(single,female,3,3))

  dErate(single,female) = 0.6076d0
  dUrate(single,female) = 0.0487d0

  ! Data single males
  dtrans(single,male,1,1) = 0.9629d0        ! EE Transition
  dtrans(single,male,1,3) = 0.0216d0        ! EN Transition

  dtrans(single,male,2,1) = 0.2586d0        ! UE Transition
  dtrans(single,male,2,2) = 0.5733d0        ! UU Transition

  dtrans(single,male,3,1) = 0.0358d0        ! NE Transition
  dtrans(single,male,3,3) = 0.9450d0        ! NN Transition


  dtrans(single,male,1,2) = &               ! EU Transition
      1.d0 -(dtrans(single,male,1,1)+dtrans(single,male,1,3))
  dtrans(single,male,2,3) = &               ! UN Transition
      1.d0 -(dtrans(single,male,2,1)+dtrans(single,male,2,2))
  dtrans(single,male,3,2) = &               ! NU Transition
      1.d0 -(dtrans(single,male,3,1)+dtrans(single,male,3,3))

  dErate(single,male) = 0.6397d0
  dUrate(single,male) = 0.0588d0

  ! Data married females
  dtrans(married,female,1,1) = 0.9618d0      ! EE Transition
  dtrans(married,female,1,3) = 0.0305d0      ! EN Transition

  dtrans(married,female,2,1) = 0.2678d0      ! UE Transition
  dtrans(married,female,2,2) = 0.4702d0      ! UU Transition

  dtrans(married,female,3,1) = 0.0330d0      ! NE Transition
  dtrans(married,female,3,3) = 0.9557d0      ! NN Transition


  dtrans(married,female,1,2) = &             ! EU Transition
      1.d0 -(dtrans(married,female,1,1)+dtrans(married,female,1,3))
  dtrans(married,female,2,3) = &             ! UN Transition
      1.d0 -(dtrans(married,female,2,1)+dtrans(married,female,2,2))
  dtrans(married,female,3,2) = &             ! NU Transition
      1.d0 -(dtrans(married,female,3,1)+dtrans(married,female,3,3))

  dErate(married,female) = 0.5532d0
  dUrate(married,female) = 0.0324d0

  ! Data married males
  dtrans(married,male,1,1) = 0.9760d0       ! EE Transition
  dtrans(married,male,1,3) = 0.0148d0       ! EN Transition

  dtrans(married,male,2,1) = 0.2988d0       ! UE Transition
  dtrans(married,male,2,2) = 0.5658d0       ! UU Transition

  dtrans(married,male,3,1) = 0.0337d0       ! NE Transition
  dtrans(married,male,3,3) = 0.9547d0       ! NN Transition


  dtrans(married,male,1,2) = &              ! EU Transition
      1.d0 -(dtrans(married,male,1,1)+dtrans(married,male,1,3))
  dtrans(married,male,2,3) = &              ! UN Transition
      1.d0 -(dtrans(married,male,2,1)+dtrans(married,male,2,2))
  dtrans(married,male,3,2) = &              ! NU Transition
      1.d0 -(dtrans(married,male,3,1)+dtrans(married,male,3,3))

  dErate(married,male) = 0.7286d0
  dUrate(married,male) = 0.0296d0

  ! Model Married-single wealth ratio
  marsing_ratio = wealth(male+female)/((wealth(male)*weights(single,male))&
                                      +(wealth(female)*weights(single,female)))

  ! Compute distance between data and model
  mydiff = 0.d0
  mydiff = mydiff + fdiff(dint_rate,int_rate)             ! Interest rate
  mydiff = mydiff + fdiff(dmarsing_ratio,marsing_ratio)   ! Married-single wealth ratio

  do myms = 0, 1
  do mysex = 1, 2
    ! EE Transition
    mydiff = mydiff + weights(myms,mysex)*fdiff(dtrans(myms,mysex,1,1),transitions(myms,mysex,1,1))
    ! EN Transition
    mydiff = mydiff + weights(myms,mysex)*fdiff(dtrans(myms,mysex,1,3),transitions(myms,mysex,1,3))
    ! UE Transition
    mydiff = mydiff + weights(myms,mysex)*fdiff(dtrans(myms,mysex,2,1),transitions(myms,mysex,2,1))
    ! UU Transition
    mydiff = mydiff + weights(myms,mysex)*fdiff(dtrans(myms,mysex,2,2),transitions(myms,mysex,2,2))
    ! NE Transition
    mydiff = mydiff + weights(myms,mysex)*fdiff(dtrans(myms,mysex,3,1),transitions(myms,mysex,3,1))
    ! NN Transition
    mydiff = mydiff + weights(myms,mysex)*fdiff(dtrans(myms,mysex,3,3),transitions(myms,mysex,3,3))
    ! Employment rate
    mydiff = mydiff + weights(myms,mysex)*fdiff(dErate(myms,mysex),Erate(myms,mysex))
    ! Unemployment rate
    mydiff = mydiff + weights(myms,mysex)*fdiff(dErate(myms,mysex),Erate(myms,mysex))
  end do
  end do

  ! PRINT RESULTS
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

    open(unit = 11, file = "distance.txt")
      write(11,*) mydiff
    close(11)

CONTAINS
  real(8) function fdiff(data,model)
    implicit none
    real(8) :: data, model

    fdiff = (log(model/data))**2.d0
    ! fdiff = ((1.d0-model)/data)**2.d0
  end function

end subroutine Difference
