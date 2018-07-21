subroutine initialisation()
  use Globals
  use Utils

  implicit none

  integer :: ind_param

  allocate(shock_z(agents,periods))
  allocate(shock_g(1:2,agents,periods))
  allocate(shock_mu(1:2,agents,periods))
  allocate(shock_lm(agents,periods))

  ! Set seed for fair price iteration
  call setseed(12345)

  ! Assigned parameters
  weights(single,male) = 0.25     ! Share of single males in the economy
  weights(single,female) = 0.25   ! Share of single females in the economy
  weights(married,male) = 0.25    ! Share of married males in the economy
  weights(married,female) = 0.25  ! Share of married females in the economy

  ! Calibrated parameters
  open(unit = 1, file = "calibrated.txt")
  read(1,*) beta      ! Discount factor
  do ind_param = 2, 26
    read(1,*) aux_param(ind_param)
  end do
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
  T = 1.3

  ! Generate shocks
  call random_number(shock_z)
  call random_number(shock_g)
  call random_number(shock_mu)
  call random_number(shock_lm)
end subroutine initialisation

!===== SINGLES INITIALISATION =====================================================================
subroutine iniSingles(mysex)
  use Globals
  use GlobalsSingles
  use Utils

  implicit none

  integer, intent(in) :: mysex

  ! Calibrated
  ! Name file
  if (mysex.eq.male) then
    alpha = aux_param(2)          ! Utility cost of working
    c_min = aux_param(3)          ! Consumption floor
    epsilon_gamma = aux_param(4)  ! Standard deviation search cost
    lambda_u = aux_param(5)       ! Probability of finding a job for unemployed agents
    lambda_n = aux_param(6)       ! Probability of finding a job for OLF agents
    sigma = aux_param(7)          ! Probability of losing a job for employed agents
  elseif (mysex.eq.female) then
    alpha = aux_param(8)          ! Utility cost of working
    c_min = aux_param(9)         ! Consumption floor
    epsilon_gamma = aux_param(10) ! Standard deviation search cost
    lambda_u = aux_param(11)      ! Probability of finding a job for unemployed agents
    lambda_n = aux_param(12)      ! Probability of finding a job for OLF agents
    sigma = aux_param(13)         ! Probability of losing a job for employed agents
  else
    print *, "Error reading sex in iniSingles"
  end if

  ! Grid for assets
  a_values = loggrid(min_a, max_a, gp_a)

  ! Productivity process
  call rouwenhurst(rho_z, 0.d0, sigma_epsilon, gp_z, z_values, z_trans)
  ! call tauchen(rho_z, sigma_epsilon, cover_z, gp_z, z_values, z_trans)
  z_values = exp(z_values)
  z_ssdist = my_ss(z_trans,gp_z)

  ! Search cost process
  gamma_bar = (3.5d0/40.d0)*alpha
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

  integer :: sex !, ind_m, ind_f, ind_mp, ind_fp, fila, columna
  ! real(8) :: auxm_z_values(gp_z), auxm_z_trans(gp_z,gp_z), auxf_z_values(gp_z), &
  ! auxf_z_trans(gp_z,gp_z)

  real(8), dimension(1:2,1:2) :: aux_rho, aux_sigma

  alpha(male) = aux_param(14)           ! Utility cost of working male
  alpha(female) = aux_param(15)         ! Utility cost of working female
  alpha(male+female) = aux_param(16)    ! Utility cost of joint work
  c_min = aux_param(17)                 ! Consumption floor
  chi = aux_param(18)                   ! Adult-equivalent scale
  epsilon_gamma(male) = aux_param(19)   ! Standard deviation search cost male
  epsilon_gamma(female) = aux_param(20) ! Standard deviation search cost female
  lambda_u(male) = aux_param(21)        ! Probability of finding a job for unemployed male
  lambda_u(female) = aux_param(22)      ! Probability of finding a job for unemployed female
  lambda_n(male) = aux_param(23)        ! Probability of finding a job for OLF male
  lambda_n(female) = aux_param(24)      ! Probability of finding a job for OLF female
  sigma(male) = aux_param(25)           ! Probability of losing a job for employed male
  sigma(female) = aux_param(26)         ! Probability of losing a job for employed female

  ! Grid for assets
  a_values = loggrid(min_a, max_a, gp_a)

  ! Tauchen for Z process - Independent
  ! call tauchen(rho_z, sigma_epsilon, cover_z, gp_z, auxm_z_values, auxm_z_trans)
  ! call tauchen(rho_z, sigma_epsilon, cover_z, gp_z, auxf_z_values, auxf_z_trans)
  !
  ! fila = 0
  ! do ind_m = 1, gp_z
  ! do ind_f = 1, gp_z
  !   fila = fila + 1
  !   columna = 0
  !   z_values(male,fila) = exp(auxm_z_values(ind_m))
  !   z_values(female,fila) = exp(auxf_z_values(ind_f))
  !   do ind_mp = 1, gp_z
  !   do ind_fp = 1, gp_z
  !     columna = columna + 1
  !     z_trans(fila,columna) = auxm_z_trans(ind_m,ind_mp)*auxf_z_trans(ind_f,ind_fp)
  !   end do
  !   end do
  ! end do
  ! end do

  ! VAR(1) for Z process
  aux_rho = 0.d0
  aux_rho(male,male) = rho_z
  aux_rho(female,female) = rho_z
  aux_sigma = 0.d0
  aux_sigma(male,male) = sigma_epsilon
  aux_sigma(female,female) = sigma_epsilon
  aux_sigma = aux_sigma**2.d0
  call discretize2vars(aux_rho,aux_sigma,gp_z,0,z_trans,z_values)
  z_values = exp(z_values)

  ! Stationary distribution Z process
  z_ssdist = my_ss(z_trans,gp_z2)

  ! Search cost process
  do sex = 1, 2
    gamma_bar(sex) = (3.5d0/40.d0)*alpha(sex)
    gamma_values(sex,1) = gamma_bar(sex) - epsilon_gamma(sex)
    gamma_values(sex,2) = gamma_bar(sex)
    gamma_values(sex,3) = gamma_bar(sex) + epsilon_gamma(sex)
    gamma_trans(sex,:) = 1.d0/real(gp_gamma)
  end do
end subroutine iniMarried
