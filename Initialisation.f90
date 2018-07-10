subroutine initialisation()
  use Globals
  use Utils

  implicit none

  allocate(shock_z(agents,periods))
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
  read(1,*) weights(single,male) ! Share of single males in the economy
  read(1,*) weights(single,female) ! Share of single females in the economy
  read(1,*) weights(married,male) ! Share of married males in the economy
  read(1,*) weights(married,female) ! Share of married females in the economy
  close(1)

  ! Check weights add up to 1
  if (abs(sum(weights)-1.d0).gt.tiny) then
    print *, "ERROR: Weights in assigned.txt o NOT add up to 1"
  end if

  ! UI process, 1: Entilted, 0: Not entailted
  IB_trans(1,1) = 1.d0 - mu
  IB_trans(1,0) = mu
  IB_trans(0,1) = 0.d0
  IB_trans(0,0) = 1.d0

  ! Guess inital prices
  KLratio = 119.60
  call Prices(KLratio, int_rate, wage)

  ! Guess average z and T
  average_z = 2.33
  T = 1.28

  ! Generate shocks
  call random_number(shock_z)
  call random_number(shock_g)
  call realise_iid(agents, periods, mu, shock_mu)
  call random_number(shock_lm)
end subroutine initialisation

!===== SINGLES INITIALISATION =====================================================================
subroutine iniSingles(myidentity)
  use Globals
  use GlobalsSingles
  use Utils

  implicit none

  character(len=2), intent(in) :: myidentity
  character(len=1024) :: format_string, aux_name

  ! Calibrated
  ! Name file
  format_string = "(A11,A2,A4)"
  write (aux_name,format_string), "calibrated_", myidentity, ".txt"
  open(unit=2, file = aux_name)
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
  ! call rouwenhorst(rho_z, 0.d0, sigma_epsilon, gp_z, cover_z, z_values, z_trans)
  call tauchen(rho_z, sigma_epsilon, cover_z, gp_z, z_values, z_trans)
  z_values = exp(z_values)
  z_ssdist = my_ss(z_trans,gp_z)

  ! Search cost process
  gamma_values(1) = gamma_bar - epsilon_gamma
  gamma_values(2) = gamma_bar
  gamma_values(3) = gamma_bar + epsilon_gamma
  gamma_trans = 1.d0/real(gp_gamma)
end subroutine iniSingles

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
