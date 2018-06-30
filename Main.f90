program calibration
  use Globals
  use Utils
  implicit none

  real(8) :: t_start, t_finish, myx, myf
  real(8), dimension(gp_q, gp_q) :: aux_q

! INITIALISATION: Read parameters from outside world, compute grids, tranistion matrices, and initialise equilibrium variables
  call cpu_time(t_start)
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
  call loggrid(min_a, max_a, gp_a, a_values)

  ! Productivity process
  call tauchen(rho_z, sigma_epsilon, cover_z, gp_z, z_values, z_trans)
  z_values = exp(z_values)

  ! Match quality process
  call tauchen(0.d0, sigma_q, cover_q, gp_q, q_values, aux_q)
  q_values = exp(q_values)
  q_trans = aux_q(1,:)

  ! Search cost process
  gamma_values(1) = gamma_bar - epsilon_gamma
  gamma_values(2) = gamma_bar
  gamma_values(3) = gamma_bar + epsilon_gamma
  gamma_trans = 1.d0/real(gp_gamma)

  ! UI process, 1: Entilted, 2: Not entailted
  IB_values(1) = 1.d0
  IB_values(2) = 0.d0
  IB_trans(1,1) = mu
  IB_trans(1,2) = 1.d0 - mu
  IB_trans(2,1) = 0.d0
  IB_trans(2,2) = 1.d0

  ! Guess inital prices
  KLratio = 129.314057
  call Prices(KLratio, int_rate, wage)

  ! Guess average z and T
  average_z = 2.5674
  T = 1.40182
  call cpu_time(t_finish)
  print *, "Initalisation time:"
  call mytime(t_finish-t_start)

! COMPUTE DECISION RULES
  call cpu_time(t_start)
  call golden_method(myfun, -100.d0, 100.d0, myx, myf)
  print *, myx, myf, myfun(3.d0)
  call cpu_time(t_finish)
  print *, "Time to compute decision rules:"
  call mytime(t_finish-t_start)

CONTAINS
  real(8) function myfun(x)
    real(8), intent(in) :: x
    myfun = -(x-5)**2.d0
  end function myfun

end program calibration
