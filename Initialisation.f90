subroutine initialisation()
  use Globals
  use Utils

  implicit none

  allocate(shock_z(agents,periods))
  allocate(shock_g(1:2,agents,periods))
  allocate(shock_mu(1:2,agents,periods))
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
  call random_number(shock_mu)
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
  read(2,*) c_min           ! Consumption floor
  read(2,*) gamma_bar       ! Average search cost
  read(2,*) epsilon_gamma   ! Standard deviation search cost
  read(2,*) rho_z           ! Persistence productivity process
  read(2,*) sigma_epsilon   ! Standard deviation productivity process
  read(2,*) lambda_u        ! Probability of finding a job for unemployed agents
  read(2,*) lambda_n        ! Probability of finding a job for OLF agents
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

!===== MARRIED INITIALISATION =====================================================================
subroutine iniMarried()
  use Globals
  use GlobalsMarried
  use Utils

  implicit none

  integer :: sex, ind_m, ind_f, ind_mp, ind_fp, fila, columna
  real(8) :: auxm_z_values(gp_z), auxm_z_trans(gp_z,gp_z), auxf_z_values(gp_z), &
  auxf_z_trans(gp_z,gp_z)

  open(unit = 3, file = "calibrated_married.txt")
  read(3,*) alpha(male)             ! Utility cost of working male
  read(3,*) alpha(female)           ! Utility cost of working female
  read(3,*) alpha(3)                ! Utility cost of joint work
  read(3,*) beta                    ! Discount factor
  read(3,*) c_min                   ! Consumption floor
  read(3,*) chi                     ! Adult-equivalent scale
  read(3,*) gamma_bar(male)         ! Average search cost male
  read(3,*) gamma_bar(female)       ! Average search cost female
  read(3,*) epsilon_gamma(male)     ! Standard deviation search cost male
  read(3,*) epsilon_gamma(female)   ! Standard deviation search cost female
  read(3,*) rho_z(male)             ! Persistence productivity male
  read(3,*) rho_z(female)           ! Persistence productivity female
  read(3,*) sigma_epsilon(male)     ! Standard deviation productivity male
  read(3,*) sigma_epsilon(female)   ! Standard deviation productivity female
  read(3,*) lambda_u(male)          ! Probability of finding a job for unemployed male
  read(3,*) lambda_u(female)        ! Probability of finding a job for unemployed female
  read(3,*) lambda_n(male)          ! Probability of finding a job for OLF male
  read(3,*) lambda_n(female)        ! Probability of finding a job for OLF female
  read(3,*) sigma(male)             ! Probability of losing a job for employed male
  read(3,*) sigma(female)           ! Probability of losing a job for employed female
  close(3)

  ! Grid for assets
  a_values = loggrid(min_a, max_a, gp_a)

  ! TEMPORARY Z PROCESS
  call tauchen(rho_z(male), sigma_epsilon(male), cover_z, gp_z, auxm_z_values, auxm_z_trans)
  call tauchen(rho_z(female), sigma_epsilon(female), cover_z, gp_z, auxf_z_values, auxf_z_trans)

  fila = 0

  do ind_m = 1, gp_z
  do ind_f = 1, gp_z
    fila = fila + 1
    columna = 0
    z_values(male,fila) = exp(auxm_z_values(ind_m))
    z_values(female,fila) = exp(auxf_z_values(ind_f))
  do ind_mp = 1, gp_z
  do ind_fp = 1, gp_z
    columna = columna + 1
    z_trans(fila,columna) = auxm_z_trans(ind_m,ind_mp)*auxf_z_trans(ind_f,ind_fp)
  end do
  end do
  end do
  end do

  z_ssdist = my_ss(z_trans,gp_z2)

  ! Search cost process
  do sex = 1, 2
    gamma_values(sex,1) = gamma_bar(sex) - epsilon_gamma(sex)
    gamma_values(sex,2) = gamma_bar(sex)
    gamma_values(sex,3) = gamma_bar(sex) + epsilon_gamma(sex)
    gamma_trans(sex,:) = 1.d0/real(gp_gamma)
  end do
end subroutine iniMarried
