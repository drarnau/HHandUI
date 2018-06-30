subroutine ExpectedValues(exp_N, exp_U, exp_W)
  ! Computes auxiliary vector to store the expected value for each level of assets tomorrow
  use Globals
  implicit none
  integer :: ind_z, ind_q, ind_b, ind_ap, ind_zp, ind_qp, ind_gp, ind_bp
  real(8) :: prob_N, prob_U, prob_W
  real(8), dimension(gp_z,gp_a), intent(out) :: exp_N
  real(8), dimension(1:2,gp_z,gp_a), intent(out) :: exp_U
  real(8), dimension(gp_q,gp_z,gp_a), intent(out) :: exp_W

  ! Initialise arrays to 0
  exp_N = 0.d0
  exp_U = 0.d0
  exp_W = 0.d0

  ! Iterate over values of z TODAY
  do ind_z = 1, gp_z
    ! Iterate over all values of assets, z, q, and gammas tomorrow
    do ind_ap = 1, gp_a
    do ind_zp = 1, gp_z
    do ind_qp = 1, gp_q
    do ind_gp = 1, gp_gamma
      prob_N = z_trans(ind_z,ind_zp)*q_trans(ind_qp)*gamma_trans(ind_gp)
      exp_N(ind_z,ind_ap) = exp_N(ind_z,ind_ap) + prob_N*(&
                      (lambda_n*V_vf(2,ind_gp,ind_qp,ind_zp,ind_ap)) + &
                      ((1.d0-lambda_n)*J_vf(2,ind_gp,ind_zp,ind_ap)))

      ! Iterate over values of match quality TODAY
      prob_W = prob_N
      do ind_q = 1, gp_q
        exp_W(ind_q,ind_z,ind_ap) = exp_W(ind_q,ind_z,ind_ap) + prob_W*(&
                        ((1.d0-sigma-lambda_e)*V_vf(2,ind_gp,ind_qp,ind_zp,ind_ap)) + &
                        (lambda_e*V_vf(2,ind_gp,max(ind_q,ind_qp),ind_zp,ind_ap)) + &
                        (sigma*(1.d0-lambda_u)*J_vf(1,ind_gp,ind_zp,ind_ap)) + &
                        (sigma*lambda_u*V_vf(1,ind_gp,ind_qp,ind_zp,ind_ap)))
      end do

      ! Iterate over values of UI TODAY and TOMORROW
      do ind_b = 1, 2
      do ind_bp = 1, 2
        prob_U = prob_W*IB_trans(ind_b,ind_bp)
        exp_U(ind_b,ind_z,ind_ap) = exp_U(ind_b,ind_z,ind_ap) + prob_U*(&
                        (lambda_u*V_vf(ind_bp,ind_gp,ind_qp,ind_zp,ind_ap)) + &
                        ((1.d0-lambda_u)*J_vf(ind_bp,ind_gp,ind_zp,ind_ap)))
      end do
      end do
    end do
    end do
    end do
    end do
  end do
end subroutine ExpectedValues
