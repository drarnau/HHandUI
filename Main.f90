program calibration
  use Globals
  use Utils
  implicit none

  integer, parameter :: maxIter = 10000, showError = 300
  real(8), parameter :: tolerance = 1.0d-4
  integer :: iter, ind_a, ind_z, ind_q, ind_g, ind_b

  real(8) :: t_start, t_finish, error_N_vf, error_U_vf, error_W_vf, error_J_vf, error_V_vf, &
            error_N_pf, error_U_pf, error_W_pf
  real(8), dimension(gp_a) :: aux_exp
  real(8), dimension(gp_q, gp_q) :: aux_q
  real(8), dimension(gp_z,gp_a) :: new_N_vf, new_N_pf, exp_N
  real(8), dimension(1:2,gp_z,gp_a) :: exp_U
  real(8), dimension(gp_q,gp_z,gp_a) :: new_W_vf, new_W_pf, exp_W
  real(8), dimension(1:2,gp_gamma,gp_z,gp_a) :: new_U_vf, new_J_vf, new_U_pf
  real(8), dimension(1:2,gp_gamma,gp_q,gp_z,gp_a) :: new_V_vf

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
  print *, "Initialisation time:"
  call mytime(t_finish-t_start)

! COMPUTE DECISION RULES
  call cpu_time(t_start)

  ! Guess a value for global value functions
  call random_number(N_vf) ! Assign random numbers to N_vf
  call random_number(U_vf)
  call random_number(W_vf)
  call ValueFunctions(N_vf,U_vf,W_vf,J_vf,V_vf)

  ! Initialise new value functions
  new_N_vf = 0.d0
  new_U_vf = 0.d0
  new_W_vf = 0.d0
  new_J_vf = 0.d0
  new_V_vf = 0.d0

  ! Initialise global policy functions
  N_pf = 0.d0
  U_pf = 0.d0
  W_pf = 0.d0

  ! Initialise new policy functions
  new_N_pf = 1.d0
  new_U_pf = 1.d0
  new_W_pf = 1.d0

  ! Value function iteration
  do iter = 1, maxIter
    ! Compute expected values
    call ExpectedValues(exp_N, exp_U, exp_W)

    ! Iterate over all states TODAY
    do ind_a = 1, gp_a
    do ind_z = 1, gp_z
      aux_exp = exp_N(ind_z, :)
      call value_N(ind_a, aux_exp, new_N_pf(ind_z,ind_a), new_N_vf(ind_z,ind_a))

      do ind_q = 1, gp_q
        aux_exp = exp_W(ind_q, ind_z, :)
        call value_W(ind_a, ind_z, ind_q, aux_exp, &
                    new_W_pf(ind_q,ind_z,ind_a), new_W_vf(ind_q,ind_z,ind_a))
      end do

      do ind_b = 1,2
        aux_exp = exp_U(ind_b, ind_z, :)
        do ind_g = 1, gp_gamma
          call value_U(ind_a, ind_z, ind_g, ind_b, aux_exp, &
                      new_U_pf(ind_b,ind_g,ind_z,ind_a), new_U_vf(ind_b,ind_g,ind_z,ind_a))
        end do
      end do
    end do
    end do

    ! Compute new J and new V
    call ValueFunctions(new_N_vf,new_U_vf,new_W_vf,new_J_vf,new_V_vf)

    ! Compute errors
    error_N_vf = maxval(abs(new_N_vf-N_vf))
    error_U_vf = maxval(abs(new_U_vf-U_vf))
    error_W_vf = maxval(abs(new_W_vf-W_vf))
    error_J_vf = maxval(abs(new_J_vf-J_vf))
    error_V_vf = maxval(abs(new_V_vf-V_vf))

    error_N_pf = maxval(abs(new_N_pf-N_pf))
    error_U_pf = maxval(abs(new_U_pf-U_pf))
    error_W_pf = maxval(abs(new_W_pf-W_pf))

    ! Update value functions and decision rules
    N_vf = new_N_vf
    U_vf = new_U_vf
    W_vf = new_W_vf
    J_vf = new_J_vf
    V_vf = new_V_vf

    N_pf = new_N_pf
    U_pf = new_U_pf
    W_pf = new_W_pf

    if ((error_N_vf+error_U_vf+error_W_vf+error_J_vf+error_V_vf.lt.tolerance).and.&
        (error_N_pf+error_U_pf+error_W_pf.lt.tolerance)) then
      print *, "Value function iteration finished at iteration:", iter
      exit
    elseif (mod(iter,showError).eq.0) then
      print *, "Current value function iteration errors, at iteration:", iter
      print *, "Error for N is:", error_N_vf
      print *, "Error for U is:", error_U_vf
      print *, "Error for W is:", error_W_vf
      print *, "Error for J is:", error_J_vf
      print *, "Error for V is:", error_V_vf
      print *, "Error for N's policy function is:", error_N_pf
      print *, "Error for U's policy function is:", error_U_pf
      print *, "Error for W's policy function is:", error_W_pf
    end if
  end do
  call cpu_time(t_finish)
  print *, "Time to compute decision rules:"
  call mytime(t_finish-t_start)

end program calibration
