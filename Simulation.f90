subroutine realise_shocks(agents, periods, trans_matrix, gp_tm, shocks)
  ! Computes the realised shock implied by trans_matrix, and returns it in shocks

  implicit none

  integer, intent(in) :: agents, periods, gp_tm
  integer :: ind_ag, ind_p, ind_sh, rep
  integer, dimension(agents, periods), intent(out) :: shocks
  real(8) :: aux_sum
  real(8), dimension(gp_tm,gp_tm), intent(in) :: trans_matrix
  real(8), dimension(agents, periods) :: aux_shocks

  call random_number(aux_shocks)

  ! Initial guess
  shocks(:,1) = 1

  ! Repeat it twice to have output start from stationary
  do rep = 1, 2
    do ind_ag = 1, agents
    do ind_p = 2, periods
      ! Find out shock
      aux_sum = 0.d0
      do ind_sh = 1, gp_tm
        aux_sum  = aux_sum + trans_matrix(shocks(ind_ag,ind_p-1),ind_sh)
        if (aux_sum.ge.aux_shocks(ind_ag,ind_p)) exit
      end do
      shocks(ind_ag,ind_p) = ind_sh
    end do
    end do
    if (rep.eq.1) then
      ! Shocks of initial period take the value of final shocks of first repetition
      shocks(:,1) = shocks(:,periods)
    end if
  end do
end subroutine realise_shocks

subroutine random_integers(rows, columns, lb, ub, int_matrix)
  ! Returns a int_matrix(rows, columns) of random integer number numbers in [lb,ub]
  integer, intent(in) :: rows, columns, lb, ub
  integer, dimension(rows, columns), intent(out) :: int_matrix
  real(8), dimension(rows, columns) :: aux_mat

  call random_number(aux_mat)

  int_matrix = lb + floor(real(ub+1-lb)*aux_mat)
end subroutine random_integers

! subroutine simulation
!   use Globals
!   use Utils
!
!   implicit none
!
! end subroutine simulation
