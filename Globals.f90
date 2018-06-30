module globals
  implicit none

  ! Integer parameters
  integer, parameter :: gp_a = 48, gp_z = 20, gp_q = 7, gp_gamma = 3

  ! Real parameters
  real(8), parameter :: cover_z = 2.d0 , cover_q = 2.d0, min_a = 0, max_a = 1440.0, &
                        tiny = 1.0d-10 

  ! Real variables
  real(8) :: alpha, beta, gamma_bar, epsilon_gamma, rho_z, sigma_epsilon, sigma_q, lambda_e, &
  lambda_u, lambda_n, sigma, mu, b_0, b_bar, theta, delta, tau, int_rate, wage, KLratio, average_z, T

  ! Real vectors
  real(8), dimension(gp_a) :: a_values
  real(8), dimension(gp_z) :: z_values
  real(8), dimension(gp_q) :: q_values, q_trans
  real(8), dimension(gp_gamma) :: gamma_values, gamma_trans
  real(8), dimension(1:2) :: IB_values

  ! Real matrices
  real(8), dimension(gp_z,gp_z) :: z_trans
  real(8), dimension(1:2, 1:2) :: IB_trans
  real(8), dimension(gp_z,gp_a) :: N_vf, N_pf
  real(8), dimension(gp_q,gp_z,gp_a) :: W_vf, W_pf
  real(8), dimension(1:2,gp_gamma,gp_z,gp_a) :: U_vf, J_vf, U_pf
  real(8), dimension(1:2,gp_gamma,gp_q,gp_z,gp_a) :: V_vf

CONTAINS
  real(8) function u(consumption)
    implicit none
    real(8), intent(in) :: consumption

    u = log(max(1.0d-320,consumption))
  end function u

  ! Creates the value functions J and V given N, U, and W
  subroutine ValueFunctions(N,U,W,J,V)
    implicit none
    integer :: ind_a, ind_z, ind_g, ind_b, ind_q
    real(8), dimension(gp_z,gp_a), intent(in) :: N
    real(8), dimension(gp_q,gp_z,gp_a), intent(in) :: W
    real(8), dimension(1:2,gp_gamma,gp_z,gp_a), intent(in) :: U
    real(8), dimension(1:2,gp_gamma,gp_z,gp_a), intent(out) :: J
    real(8), dimension(1:2,gp_gamma,gp_q,gp_z,gp_a), intent(out) :: V

    do ind_a = 1, gp_a
    do ind_z = 1, gp_z
    do ind_g = 1, gp_gamma
    do ind_b = 1, 2
      J(ind_b,ind_g,ind_z,ind_a) = max(U(ind_b,ind_g,ind_z,ind_a), N(ind_z,ind_a))
      do ind_q = 1, gp_q
        V(ind_b,ind_g,ind_q,ind_z,ind_a) = max(W(ind_q,ind_z,ind_a),J(ind_b,ind_g,ind_z,ind_a))
      end do
    end do
    end do
    end do
    end do
  end subroutine ValueFunctions

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
end module globals
