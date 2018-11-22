module Globals
  implicit none

  integer, parameter :: single = 0, married = 1, male = 1, female = 2
  integer, parameter :: agents = 10000, periods = 5001, printed = 500

  ! Real parameters
  real(8), parameter :: cover_z = 2.d0 , cover_q = 2.d0, tiny = 1.0d-10

  ! Assigned parameters
  real(8), parameter :: mu = 1.d0/6.d0          ! Average duration UI
  real(8), parameter :: b_0 = 0.0               ! Default replacement ratio 0.23
  real(8), parameter :: b_bar = 0.0             ! Benefits cap 0.465
  real(8), parameter :: theta = 0.3             ! Capital share of output in production function
  real(8), parameter :: delta = 0.0067          ! Capital depreciation
  real(8), parameter :: T = 1.75                ! Household lump-sum transfer
  real(8), parameter :: rho_z = 0.99617         ! Persistence productivity process
  real(8), parameter :: sigma_epsilon = 0.09565 ! Standard deviation productivity process

  ! Real variables
  real(8) :: beta, int_rate, wage, KLratio, average_z, tau, new_KLratio, new_average_z, new_tau

  ! Real vectors
  real(8), dimension(1:3) :: aux_KLratio, aux_average_z, aux_labincome, aux_bpaid, wealth
  real(8), dimension(1:3) :: aux_tot_z
  real(8), dimension(1:26) :: aux_param

  ! Real matrices
  real(8), dimension(:,:), allocatable :: shock_lm, shock_z
  real(8), dimension(:,:,:), allocatable :: shock_g, shock_mu
  real(8), dimension(0:1, 0:1) :: IB_trans
  real(8), dimension(0:1,1:2) :: Erate, Urate, Nrate, weights
  real(8), dimension(0:1,1:2,1:3,1:3) :: transitions

CONTAINS

  ! Computes equilibrium prices
  subroutine Prices(my_KL,my_r,my_w)
    implicit none
    real(8), intent(in) :: my_KL
    real(8), intent(out) ::  my_r, my_w

    my_w = (1.d0-theta)*(my_KL**theta)
    my_r = (theta*(my_KL**(theta-1.d0))) - delta
  end subroutine Prices

  ! Computes unemployment benefits
  real(8) function benefits(productivity)
    implicit none
    real(8), intent(in) :: productivity

    if ((productivity*b_0).lt.(average_z*b_bar)) then
      benefits = productivity*b_0*wage
    else
      benefits = average_z*b_bar*wage
    end if
  end function benefits

  ! Aggregates variables according to module weights
  real(8) function aggregate4(matrix)
    implicit none
    integer :: sex, mstatus
    real(8) :: aux_agg
    real(8), dimension(0:1,1:2) :: matrix

    aux_agg = 0.d0
    do mstatus = 0,1
    do sex = 1, 2
      aux_agg = aux_agg + weights(mstatus,sex)*matrix(mstatus,sex)
    end do
    end do
    aggregate4 = aux_agg
  end function aggregate4

  real(8) function aggregate3(vector)
    implicit none
    real(8), dimension(1:3) :: vector

    aggregate3 = (vector(1)*weights(0,1)) + (vector(2)*weights(0,2)) + &
                 (vector(3)*(weights(1,1)+weights(1,2)))
  end function aggregate3
end module Globals

!===== SINGLES GLOBALS ============================================================================
module GlobalsSingles
  implicit none

  ! Integer parameters
  integer, parameter :: gp_a = 55, gp_z = 17, gp_gamma = 3

  ! Real parameters
  real(8), parameter :: min_a = 0, max_a = 1440.0

  ! Real variables
  real(8) :: alpha, c_min, gamma_bar, epsilon_gamma, lambda_u, lambda_n, sigma

  ! Real vectors
  real(8), dimension(gp_a) :: a_values
  real(8), dimension(gp_z) :: z_values, z_ssdist
  real(8), dimension(gp_gamma) :: gamma_values, gamma_trans

  ! Real matrices
  real(8), dimension(gp_z,gp_z) :: z_trans
  real(8), dimension(gp_a,gp_z) :: N_vf, N_pf, W_vf, W_pf
  real(8), dimension(gp_a,gp_z,gp_gamma,0:1) :: U_vf, J_vf, U_pf, V_vf

CONTAINS
  ! Computes singles utility
  real(8) function u(consumption)
    implicit none
    real(8), intent(in) :: consumption

    u = log(max(1.0d-11,consumption-c_min))
  end function u

  ! Creates the value functions J and V given N, U, and W
  subroutine ValueFunctions(N,U,W,J,V)
    implicit none
    integer :: ind_a, ind_z, ind_g, ind_b
    real(8), dimension(gp_a,gp_z), intent(in) :: N, W
    real(8), dimension(gp_a,gp_z,gp_gamma,0:1), intent(in) :: U
    real(8), dimension(gp_a,gp_z,gp_gamma,0:1), intent(out) :: J, V

    do ind_a = 1, gp_a
    do ind_z = 1, gp_z
    do ind_g = 1, gp_gamma
    do ind_b = 0, 1
      J(ind_a,ind_z,ind_g,ind_b) = max(U(ind_a,ind_z,ind_g,ind_b), N(ind_a,ind_z))
      V(ind_a,ind_z,ind_g,ind_b) = max(W(ind_a,ind_z),J(ind_a,ind_z,ind_g,ind_b))
    end do
    end do
    end do
    end do
  end subroutine ValueFunctions
end module GlobalsSingles

!===== MARRIED GLOBALS ============================================================================
module GlobalsMarried
  implicit none

  ! Integer parameters
  integer, parameter :: gp_a = 55, gp_z = 9, gp_gamma = 3, gp_z2 = gp_z**2

  ! Real parameters
  real(8), parameter :: min_a = 0, max_a = 2880.0

  ! Real variables
  real(8) :: c_min, chi

  ! Real vectors
  real(8), dimension(1:3) :: alpha
  real(8), dimension(1:2) :: gamma_bar, epsilon_gamma, lambda_u, lambda_n, sigma
  real(8), dimension(gp_a) :: a_values
  real(8), dimension(gp_z2) :: z_ssdist

  ! Real matrices
  real(8), dimension(1:2,gp_z2) :: z_values
  real(8), dimension(1:2,gp_gamma) :: gamma_values, gamma_trans
  real(8), dimension(gp_z2,gp_z2) :: z_trans
  real(8), dimension(gp_a,gp_z2) :: NN_vf, NN_pf, WW_vf, WW_pf, WN_vf, WN_pf, NW_vf, NW_pf
  real(8), dimension(gp_a,gp_z2,gp_gamma,0:1) :: WU_vf, WU_pf, UW_vf, UW_pf, NU_vf, NU_pf, &
                                                      UN_vf, UN_pf
  real(8), dimension(gp_a,gp_z2,gp_gamma,gp_gamma,0:1,0:1) :: UU_vf, UU_pf, JJ_vf, VJ_vf, &
                                                                    JV_vf, VV_vf

CONTAINS

  ! Computes married utility
  real(8) function u(consumption)
    implicit none
    real(8), intent(in) :: consumption

    u = log(max(1.0d-11,((consumption-c_min)/chi)))
  end function u

  ! Create value functions JJ, VJ, JV, and VV, given NN, NU, NW, UN, UU, UW, WW, WU, WN
  subroutine ValueFunctions(WW, WU, WN, UW, UU, UN, NW, NU, NN, JJ, VJ, JV, VV)
    implicit none
    integer :: ind_a, ind_z, ind_g_m, ind_g_f, ind_b_m, ind_b_f
    real(8), dimension(gp_a,gp_z2), intent(in) :: NN, WW, WN, NW
    real(8), dimension(gp_a,gp_z2,gp_gamma,0:1), intent(in) :: WU, UW, NU, UN
    real(8), dimension(gp_a,gp_z2,gp_gamma,gp_gamma,0:1,0:1), intent(in) :: UU
    real(8), dimension(gp_a,gp_z2,gp_gamma,gp_gamma,0:1,0:1), intent(out) :: JJ, VJ, JV, VV

    do ind_a = 1, gp_a
    do ind_z = 1, gp_z2
    do ind_g_m = 1, gp_gamma
    do ind_g_f = 1, gp_gamma
    do ind_b_m = 0, 1
    do ind_b_f = 0, 1
      ! None has a job offer
      JJ(ind_a,ind_z,ind_g_m,ind_g_f,ind_b_m,ind_b_f) = &
        max(NN(ind_a,ind_z), NU(ind_a,ind_z,ind_g_f,ind_b_f), UN(ind_a,ind_z,ind_g_m,ind_b_m), &
            UU(ind_a,ind_z,ind_g_m,ind_g_f,ind_b_m,ind_b_f))

      ! Male has a job offer
      VJ(ind_a,ind_z,ind_g_m,ind_g_f,ind_b_m,ind_b_f) = &
        max(WU(ind_a,ind_z,ind_g_f,ind_b_f), WN(ind_a,ind_z), &
            JJ(ind_a,ind_z,ind_g_m,ind_g_f,ind_b_m,ind_b_f))

      ! Female has a job offer
      JV(ind_a,ind_z,ind_g_m,ind_g_f,ind_b_m,ind_b_f) = &
        max(UW(ind_a,ind_z,ind_g_m,ind_b_m), NW(ind_a,ind_z), &
            JJ(ind_a,ind_z,ind_g_m,ind_g_f,ind_b_m,ind_b_f))

      ! Both have a job offer
      VV(ind_a,ind_z,ind_g_m,ind_g_f,ind_b_m,ind_b_f) = max(WW(ind_a,ind_z), &
            WU(ind_a,ind_z,ind_g_f,ind_b_f), WN(ind_a,ind_z), &
            UW(ind_a,ind_z,ind_g_m,ind_b_m), NW(ind_a,ind_z), &
            JJ(ind_a,ind_z,ind_g_m,ind_g_f,ind_b_m,ind_b_f))
    end do
    end do
    end do
    end do
    end do
    end do
  end subroutine ValueFunctions
end module GlobalsMarried
