subroutine initialisation()
  use Globals
  use Utils

  implicit none

  real(8), dimension(gp_q, gp_q) :: aux_q

  allocate(shock_z(agents,periods))
  allocate(shock_q(agents,periods))
  allocate(shock_g(agents,periods))
  allocate(shock_mu(agents,periods))
  allocate(shock_lm(agents,periods))

  ! Set seed for fair price iteration
  call setseed(12345)

  ! Assigned
  open(unit=1, file='assigned.txt')
  read(1,*) mu    ! Average duration UI
  read(1,*) b_0   ! Default replacement ratio
  read(1,*) b_bar ! Benefits cap
  read(1,*) theta ! Capital share of output in aggregate production function
  read(1,*) delta ! Capital depreciation
  read(1,*) tau   ! Proportional tax on labor income
  close(1)

  ! Calibrated
  open(unit=2, file='calibrated.txt')
  read(2,*) alpha           ! Utility cost of working
  read(2,*) beta            ! Discount factor
  read(2,*) gamma_bar       ! Average search cost
  read(2,*) epsilon_gamma   ! Standard deviation search cost
  read(2,*) rho_z           ! Persistence productivity process
  read(2,*) sigma_epsilon   ! Standard deviation productivity process
  read(2,*) sigma_q         ! Standard deviation match quality process
  read(2,*) lambda_e        ! Probability of finding another job for employed agents
  read(2,*) lambda_u        ! Probability of finding a job for unemployed agents
  read(2,*) lambda_n        ! Probability of finding a job for OLD agents
  read(2,*) sigma           ! Probability of losing a job for employed agents
  close(2)

  ! Grid for assets
  a_values = loggrid(min_a, max_a, gp_a)

  ! Productivity process
  call rouwenhorst(rho_z, 0.d0, sigma_epsilon, gp_z, cover_z, z_values, z_trans)
  ! call tauchen(rho_z, sigma_epsilon, cover_z, gp_z, z_values, z_trans)
  z_values = exp(z_values)
  z_ssdist = my_ss(z_trans,gp_z)
  call realise_shocks(agents, periods, z_trans, gp_z, shock_z)

  ! Match quality process
  call rouwenhorst(0.d0, 0.d0, sigma_q, gp_q, cover_q, q_values, aux_q)
  q_values = exp(q_values)
  q_trans = aux_q(1,:)
  call realise_shocks(agents, periods, aux_q, gp_q, shock_q)

  ! Search cost process
  gamma_values(1) = gamma_bar - epsilon_gamma
  gamma_values(2) = gamma_bar
  gamma_values(3) = gamma_bar + epsilon_gamma
  gamma_trans = 1.d0/real(gp_gamma)

  ! call random_integers(agents,periods,1,3,shock_g)

  ! UI process, 1: Entilted, 0: Not entailted
  IB_trans(1,1) = 1.d0 - mu
  IB_trans(1,0) = mu
  IB_trans(0,1) = 0.d0
  IB_trans(0,0) = 1.d0

  ! Guess inital prices
  KLratio = 129.314057
  call Prices(KLratio, int_rate, wage)

  ! Guess average z and T
  average_z = 2.5674
  T = 1.40182

  ! Generate shocks
  call random_number(shock_lm)
  call realise_iid(agents, periods, mu, shock_mu)
  call random_integers(agents, periods, 1, gp_gamma, shock_g)
end subroutine initialisation

!===== REALISE SHOCKS =============================================================================
subroutine realise_shocks(agents, periods, trans_matrix, gp_tm, shocks)
  ! Computes the realised shock implied by trans_matrix, and returns it in shocks
  use Utils
  implicit none

  integer, intent(in) :: agents, periods, gp_tm
  integer :: ind_ag, ind_p, ind_sh
  integer, dimension(agents, periods), intent(out) :: shocks
  real(8) :: aux_sum
  real(8), dimension(gp_tm) :: ss_dist, aux_vec
  real(8), dimension(gp_tm,gp_tm), intent(in) :: trans_matrix
  real(8), dimension(agents, periods) :: aux_shocks

  call random_number(aux_shocks)

  ! Find steady state distribution
  ss_dist = my_ss(trans_matrix, gp_tm)

  ! Simulate for all periods
  do ind_ag = 1, agents
  do ind_p = 1, periods
    ! Chose auxiliary vector for first period
    if (ind_p.eq.1) then
      aux_vec = ss_dist
    else
      aux_vec = trans_matrix(shocks(ind_ag,ind_p-1),:)
    end if
    ! Find out shock
    aux_sum = 0.d0
    do ind_sh = 1, gp_tm
      aux_sum  = aux_sum + aux_vec(ind_sh)
      if (aux_sum.ge.aux_shocks(ind_ag,ind_p)) exit
    end do
    shocks(ind_ag,ind_p) = ind_sh
  end do
  end do
end subroutine realise_shocks

!===== REALISE SIMPLE IID SHOCKS ==================================================================
subroutine realise_iid(agents, periods, param_value, shocks)
  ! Computes the realised shock of an iid process

  implicit none

  integer, intent(in) :: agents, periods
  integer :: ind_ag, ind_p
  integer, dimension(agents, periods), intent(out) :: shocks
  real(8), intent(in) :: param_value
  real(8), dimension(agents, periods) :: aux_shocks

  shocks = 0

  call random_number(aux_shocks)

  do ind_ag = 1, agents
  do ind_p = 1, periods
    if (aux_shocks(ind_ag,ind_p).le.param_value) then
      shocks(ind_ag,ind_p) = 1
    end if
  end do
  end do
end subroutine realise_iid

!===== RANDOM INTEGERS ============================================================================
subroutine random_integers(rows, columns, lb, ub, int_matrix)
  ! Returns a int_matrix(rows, columns) of random integer number numbers in [lb,ub]
  integer, intent(in) :: rows, columns, lb, ub
  integer, dimension(rows, columns), intent(out) :: int_matrix
  real(8), dimension(rows, columns) :: aux_mat

  call random_number(aux_mat)

  int_matrix = lb + floor(real(ub+1-lb)*aux_mat)
end subroutine random_integers
