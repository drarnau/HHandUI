!===== VF SINGLES =================================================================================
subroutine VFSingles()
  use Globals
  use GlobalsSingles
  use Utils

  implicit none

  integer, parameter :: maxIter = 20000, showError = 500
  real(8), parameter :: tolerance = 1.0d-4
  integer :: iter, ind_a, ind_z, ind_g, ind_b

  real(8) :: error_N_vf, error_U_vf, error_W_vf, error_J_vf, error_V_vf, error_N_pf, error_U_pf, error_W_pf

  real(8), dimension(gp_a,gp_z) :: new_N_vf, new_N_pf, exp_N, new_W_vf, new_W_pf, exp_W
  real(8), dimension(gp_a,gp_z,0:1) :: exp_U
  real(8), dimension(gp_a,gp_z,gp_gamma,0:1) :: new_U_vf, new_J_vf, new_U_pf, new_V_vf

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
    call ExpectedValuesS(exp_N, exp_U, exp_W)

    ! Iterate over all states TODAY
    do ind_a = 1, gp_a
    do ind_z = 1, gp_z
      call value_N(ind_a, exp_N(:,ind_z), new_N_pf(ind_a,ind_z), new_N_vf(ind_a,ind_z))
      call value_W(ind_a, ind_z, exp_W(:,ind_z), &
                    new_W_pf(ind_a,ind_z), new_W_vf(ind_a,ind_z))

      do ind_b = 0, 1
        do ind_g = 1, gp_gamma
          call value_U(ind_a, ind_z, ind_g, ind_b, exp_U(:,ind_z,ind_b), &
                      new_U_pf(ind_a,ind_z,ind_g,ind_b), new_U_vf(ind_a,ind_z,ind_g,ind_b))
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
      print *, ""
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
      print *, ""
    end if
  end do

end subroutine VFSingles

!===== EXPECTED VALUES SINGLES ====================================================================
subroutine ExpectedValuesS(exp_N, exp_U, exp_W)
  ! Computes auxiliary vector to store the expected value for each level of assets tomorrow
  use Globals
  use GlobalsSingles
  implicit none

  integer :: ind_z, ind_b, ind_ap, ind_zp, ind_gp, ind_bp
  real(8) :: prob_N, prob_U, prob_W
  real(8), dimension(gp_a,gp_z), intent(out) :: exp_N, exp_W
  real(8), dimension(gp_a,gp_z,0:1), intent(out) :: exp_U

  ! Initialise arrays to 0
  exp_N = 0.d0
  exp_U = 0.d0
  exp_W = 0.d0

  ! Iterate over values of z TODAY
  do ind_z = 1, gp_z
    ! Iterate over all values of assets, z, and gammas tomorrow
    do ind_ap = 1, gp_a
    do ind_zp = 1, gp_z
    do ind_gp = 1, gp_gamma
      prob_N = z_trans(ind_z,ind_zp)*gamma_trans(ind_gp)
      exp_N(ind_ap,ind_z) = exp_N(ind_ap,ind_z) + prob_N*(&
                      (lambda_n*V_vf(ind_ap,ind_zp,ind_gp,0)) + &
                      ((1.d0-lambda_n)*J_vf(ind_ap,ind_zp,ind_gp,0)))

      prob_W = prob_N
      exp_W(ind_ap,ind_z) = exp_W(ind_ap,ind_z) + prob_W*(&
                        ((1.d0-sigma)*V_vf(ind_ap,ind_zp,ind_gp,0)) + &
                        (sigma*(1.d0-lambda_u)*J_vf(ind_ap,ind_zp,ind_gp,1)) + &
                        (sigma*lambda_u*V_vf(ind_ap,ind_zp,ind_gp,1)))

      ! Iterate over values of UI TODAY and TOMORROW
      do ind_b = 0, 1
      do ind_bp = 0, 1
        prob_U = prob_W*IB_trans(ind_b,ind_bp)
        exp_U(ind_ap,ind_z,ind_b) = exp_U(ind_ap,ind_z,ind_b) + prob_U*(&
                        (lambda_u*V_vf(ind_ap,ind_zp,ind_gp,ind_bp)) + &
                        ((1.d0-lambda_u)*J_vf(ind_ap,ind_zp,ind_gp,ind_bp)))
      end do
      end do
    end do
    end do
    end do
  end do
end subroutine ExpectedValuesS

!===== OLF ========================================================================================
subroutine value_N(ind_a, aux_exp, pf, vf)
  use Globals
  use GlobalsSingles
  use Utils
  implicit none

  integer, intent(in) :: ind_a
  real(8) :: aux_max_a, income
  real(8), dimension(gp_a), intent(in) :: aux_exp
  real(8), intent(out) :: pf, vf

  ! Compute income
  income = (1.d0+int_rate)*a_values(ind_a) + T

  ! Compute upper bound for assets
  aux_max_a = min(max_a-tiny, max(min_a, income))

  ! Find level of assets that maximised Bellman Equation
  call golden_method(valor_N, min_a, aux_max_a, pf, vf)

CONTAINS
  real(8) function valor_N(ap)
    implicit none
    real(8), intent(in) :: ap
    real(8) :: aux_inter, cons

    ! Compute interpolated value for tomorrow
    aux_inter = my_inter(a_values, aux_exp, gp_a, ap)

    ! Compute consumption
    cons = income - ap

    valor_N = u(cons) + beta*aux_inter
  end function valor_N
end subroutine value_N

!===== EMPLOYED ===================================================================================
subroutine value_W(ind_a, ind_z, aux_exp, pf, vf)
  use Globals
  use GlobalsSingles
  use Utils
  implicit none

  integer, intent(in) :: ind_a, ind_z
  real(8) :: aux_max_a, income
  real(8), dimension(gp_a), intent(in) :: aux_exp
  real(8), intent(out) :: pf, vf

  ! Compute income
  income = (1.d0+int_rate)*a_values(ind_a) + T + (1.d0-tau)*wage*z_values(ind_z)

  ! Compute upper bound for assets
  aux_max_a = min(max_a-tiny, max(min_a, income))

  ! Find level of assets that maximised Bellman Equation
  call golden_method(valor_W, min_a, aux_max_a, pf, vf)

CONTAINS
  real(8) function valor_W(ap)
    implicit none
    real(8), intent(in) :: ap
    real(8) :: aux_inter, cons

    ! Compute interpolated value for tomorrow
    aux_inter = my_inter(a_values, aux_exp, gp_a, ap)

    ! Compute consumption
    cons = income - ap

    valor_W = u(cons) - alpha + beta*aux_inter
  end function valor_W
end subroutine value_W

!===== EMPLOYED ===================================================================================
subroutine value_U(ind_a, ind_z, ind_g, ind_b, aux_exp, pf, vf)
  use Globals
  use GlobalsSingles
  use Utils
  implicit none

  integer, intent(in) :: ind_a, ind_z, ind_g, ind_b
  real(8) :: aux_max_a, income
  real(8), dimension(gp_a), intent(in) :: aux_exp
  real(8), intent(out) :: pf, vf

  ! Compute income
  income = (1.d0+int_rate)*a_values(ind_a) + T + (1.d0-tau)*real(ind_b)*benefits(z_values(ind_z))

  ! Compute upper bound for assets
  aux_max_a = min(max_a-tiny, max(min_a, income))

  ! Find level of assets that maximised Bellman Equation
  call golden_method(valor_U, min_a, aux_max_a, pf, vf)

CONTAINS
  real(8) function valor_U(ap)
    implicit none
    real(8), intent(in) :: ap
    real(8) :: aux_inter, cons

    ! Compute interpolated value for tomorrow
    aux_inter = my_inter(a_values, aux_exp, gp_a, ap)

    ! Compute consumption
    cons = income - ap

    valor_U = u(cons) - gamma_values(ind_g) + beta*aux_inter
  end function valor_U
end subroutine value_U

!===== VF MARRIED =================================================================================
subroutine VFMarried()
  use Globals
  use GlobalsMarried
  use Utils

  implicit none

  real(8), dimension(gp_a,gp_z**2) :: new_NN_vf, new_NN_pf, new_WW_vf, new_WW_pf, new_WN_vf, &
                                      new_WN_pf, new_NW_vf, new_NW_pf, &
                                      exp_WW, exp_WN, exp_NN, exp_NW
  real(8), dimension(gp_a,gp_z**2,0:1) :: exp_WU, exp_UW, exp_NU, exp_UN
  real(8), dimension(gp_a,gp_z**2,0:1,0:1) :: exp_UU
  real(8), dimension(gp_a,gp_z**2,gp_gamma,0:1) :: new_WU_vf, new_WU_pf, new_UW_vf, new_UW_pf, &
                                                      new_NU_vf, new_NU_pf, new_UN_vf, new_UN_pf
  real(8), dimension(gp_a,gp_z**2,gp_gamma,gp_gamma,0:1,0:1) :: new_UU_vf, new_UU_pf

  ! Guess a value for global value functions
  call random_number(WW_vf)
  call random_number(WU_vf)
  call random_number(WN_vf)
  call random_number(UW_vf)
  call random_number(UU_vf)
  call random_number(UN_vf)
  call random_number(NW_vf)
  call random_number(NU_vf)
  call random_number(NN_vf)
  call ValueFunctions(WW_vf,WU_vf,WN_vf,UW_vf,UU_vf,UN_vf,NW_vf,NU_vf,NN_vf,JJ_vf,VJ_vf,JV_vf,VV_vf)

  ! Initialise new value functions
  new_WW_vf = 0.d0
  new_WU_vf = 0.d0
  new_WN_vf = 0.d0
  new_UW_vf = 0.d0
  new_UU_vf = 0.d0
  new_UN_vf = 0.d0
  new_NW_vf = 0.d0
  new_NU_vf = 0.d0
  new_NN_vf = 0.d0

  ! Initialise global policy functions
  WW_pf = 0.d0
  WU_pf = 0.d0
  WN_pf = 0.d0
  UW_pf = 0.d0
  UU_pf = 0.d0
  UN_pf = 0.d0
  NW_pf = 0.d0
  NU_pf = 0.d0
  NN_pf = 0.d0

  ! Initialise new policy functions
  new_WW_pf = 1.d0
  new_WU_pf = 1.d0
  new_WN_pf = 1.d0
  new_UW_pf = 1.d0
  new_UU_pf = 1.d0
  new_UN_pf = 1.d0
  new_NW_pf = 1.d0
  new_NU_pf = 1.d0
  new_NN_pf = 1.d0

  ! Compute expected values
  call ExpectedValuesM(exp_WW,exp_WU,exp_WN,exp_UW,exp_UU,exp_UN,exp_NW,exp_NU,exp_NN)



end subroutine VFMarried

!===== EXPECTED VALUES MARRIED ====================================================================
subroutine ExpectedValuesM(exp_WW,exp_WU,exp_WN,exp_UW,exp_UU,exp_UN,exp_NW,exp_NU,exp_NN)
  use Globals
  use GlobalsMarried
  implicit none

  integer :: ind_z, ind_b_f, ind_ap, ind_zp, ind_gp_m, ind_gp_f, ind_b, ind_bp, ind_bp_f
  real(8) :: prob_WW, prob_WU, prob_WN, prob_UW, prob_UU, prob_UN, prob_NW, prob_NU, prob_NN
  real(8), dimension(gp_a,gp_z**2), intent(out) :: exp_WW, exp_WN, exp_NN, exp_NW
  real(8), dimension(gp_a,gp_z**2,0:1), intent(out) :: exp_WU, exp_UW, exp_NU, exp_UN
  real(8), dimension(gp_a,gp_z**2,0:1,0:1), intent(out) :: exp_UU

  ! Initialise arrays to 0
  exp_WW = 0.d0
  exp_WU = 0.d0
  exp_WN = 0.d0
  exp_UW = 0.d0
  exp_UU = 0.d0
  exp_UN = 0.d0
  exp_NW = 0.d0
  exp_NU = 0.d0
  exp_NN = 0.d0

  ! Iterate over values of z TODAY
  do ind_z = 1, gp_z**2
    ! Iterate over values of assets, z, and gammas TOMORROW
    do ind_ap = 1, gp_a
    do ind_zp = 1, gp_z**2
    do ind_gp_m = 1, gp_gamma
    do ind_gp_f = 1, gp_gamma
      prob_NN = z_trans(ind_z,ind_zp)*gamma_trans(male,ind_gp_m)*gamma_trans(female,ind_gp_f)
      exp_NN(ind_ap,ind_z) = exp_NN(ind_ap,ind_z) + prob_NN*( &
        (lambda_n(male)*lambda_n(female)*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,0))+&
        (lambda_n(male)*(1.d0-lambda_n(female))*VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,0))+&
        ((1.d0-lambda_n(male))*lambda_n(female)*JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,0))+&
        ((1.d0-lambda_n(male))*(1.d0-lambda_n(female))*JJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,0)))

      prob_WN = prob_NN
      exp_WN(ind_ap,ind_z) = exp_WN(ind_ap,ind_z) + prob_WN*( &
        ((1.d0-sigma(male))*lambda_n(female)*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,0))+&
        ((1.d0-sigma(male))*(1.d0-lambda_n(female))*VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,0))+&
        (sigma(male)*lambda_u(male)*lambda_n(female)*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,0))+&
        (sigma(male)*lambda_u(male)*(1.d0-lambda_n(female))*&
                                                    VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,0))+&
        (sigma(male)*(1.d0-lambda_u(male))*lambda_n(female)*&
                                                    JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,0))+&
        (sigma(male)*(1.d0-lambda_u(male))*(1.d0-lambda_n(female))*&
                                                    JJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,0)))

      prob_NW = prob_WN
      exp_NW(ind_ap,ind_z) = exp_NW(ind_ap,ind_z) + prob_NW*( &
        ((1.d0-sigma(female))*lambda_n(male)*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,0))+&
        ((1.d0-sigma(female))*(1.d0-lambda_n(male))*VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,0))+&
        (sigma(female)*lambda_u(female)*lambda_n(male)*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,1))+&
        (sigma(female)*lambda_u(female)*(1.d0-lambda_n(male))*&
                                                    VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,1))+&
        (sigma(female)*(1.d0-lambda_u(female))*lambda_n(male)*&
                                                    JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,1))+&
        (sigma(female)*(1.d0-lambda_u(female))*(1.d0-lambda_n(male))*&
                                                    JJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,1)))

      prob_WW = prob_NW
      exp_WW(ind_ap,ind_z) = exp_WW(ind_ap,ind_z) + prob_WW*( &
        ((1.d0-sigma(male))*(1.d0-sigma(female))*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,0))+&
        ((1.d0-sigma(male))*sigma(female)*lambda_u(female)*&
                                                    VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,1))+&
        ((1.d0-sigma(male))*sigma(female)*(1.d0-lambda_u(female))*&
                                                    VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,1))+&
        (sigma(male)*lambda_u(male)*(1.d0-sigma(female))*&
                                                    VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,0))+&
        (sigma(male)*lambda_u(male)*sigma(female)*lambda_u(female)*&
                                                    VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,1))+&
        (sigma(male)*lambda_u(male)*sigma(female)*(1.d0-lambda_u(female))*&
                                                    VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,1))+&
        (sigma(male)*(1.d0-lambda_u(male))*(1.d0-sigma(female))*&
                                                    JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,0))+&
        (sigma(male)*(1.d0-lambda_u(male))*sigma(female)*lambda_u(female)*&
                                                    JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,1))+&
        (sigma(male)*(1.d0-lambda_u(male))*sigma(female)*(1.d0-sigma(female))*&
                                                    JJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,1)))

      ! Iterate over values of UI TODAY and TOMORROW
      do ind_b = 0, 1
      do ind_bp = 0, 1
        prob_NU = prob_WW*IB_trans(ind_b,ind_bp)
        exp_NU(ind_ap,ind_z,ind_b) = exp_NU(ind_ap,ind_z,ind_b) + prob_NU*( &
          (lambda_n(male)*lambda_u(female)*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,ind_bp))+&
          (lambda_n(male)*(1.d0-lambda_u(female))*&
                                          VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,ind_bp))+&
          ((1.d0-lambda_n(male))*lambda_u(female)*&
                                          JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,ind_bp))+&
          ((1.d0-lambda_n(male))*(1.d0-lambda_u(female))*&
                                          JJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,ind_bp)))

        prob_UN = prob_NU
        exp_UN(ind_ap,ind_z,ind_b) = exp_UN(ind_ap,ind_z,ind_b) + prob_UN*( &
          (lambda_n(female)*lambda_u(male)*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,0))+&
          (lambda_n(female)*(1.d0-lambda_u(male))*&
                                          VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,0))+&
          ((1.d0-lambda_n(female))*lambda_u(male)*&
                                          JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,0))+&
          ((1.d0-lambda_n(female))*(1.d0-lambda_u(male))*&
                                          JJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,0)))

        prob_WU = prob_NU
        exp_WU(ind_ap,ind_z,ind_b) = exp_WU(ind_ap,ind_z,ind_b) + prob_WU*( &
          ((1.d0-sigma(male))*lambda_u(female)*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,ind_bp))+&
          ((1.d0-sigma(male))*(1.d0-lambda_u(female))*&
                                            VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,0,ind_bp))+&
          (sigma(male)*lambda_u(male)*lambda_u(female)*&
                                            VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,ind_bp))+&
          (sigma(male)*lambda_u(male)*(1.d0-lambda_u(female))*&
                                            VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,ind_bp))+&
          (sigma(male)*(1.d0-lambda_u(male))*lambda_u(female)*&
                                            JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,ind_bp))+&
          (sigma(male)*(1.d0-lambda_u(male))*(1.d0-lambda_u(female))*&
                                            JJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,1,ind_bp)))

        prob_UW = prob_UN
        exp_UW(ind_ap,ind_z,ind_b) = exp_UW(ind_ap,ind_z,ind_b) + prob_UW*( &
          ((1.d0-sigma(female))*lambda_u(male)*VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,0))+&
          ((1.d0-sigma(female))*(1.d0-lambda_u(male))*&
                                            VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,0))+&
          (sigma(female)*lambda_u(female)*lambda_u(male)*&
                                            VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,1))+&
          (sigma(female)*lambda_u(female)*(1.d0-lambda_u(male))*&
                                            VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,1))+&
          (sigma(female)*(1.d0-lambda_u(female))*lambda_u(male)*&
                                            JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,1))+&
          (sigma(female)*(1.d0-lambda_u(female))*(1.d0-lambda_u(male))*&
                                            JJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,1)))

        ! Iterate over value of UI TODAY and TOMORROW for both
        do ind_b_f = 0, 1
        do ind_bp_f = 0, 1
          prob_UU = prob_NU*IB_trans(ind_b_f,ind_bp_f)
          exp_UU(ind_ap,ind_z,ind_b,ind_b_f) = exp_UU(ind_ap,ind_z,ind_b,ind_b_f) + prob_UU*( &
            (lambda_u(male)*lambda_u(female)*&
                                          VV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,ind_bp_f))+&
            (lambda_u(male)*(1.d0-lambda_u(female))*&
                                          VJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,ind_bp_f))+&
            ((1.d0-lambda_u(male))*lambda_u(female)*&
                                          JV_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,ind_bp_f))+&
            ((1.d0-lambda_u(male))*(1.d0-lambda_u(female))*&
                                          JJ_vf(ind_ap,ind_zp,ind_gp_m,ind_gp_f,ind_bp,ind_bp_f)))
        end do
        end do
      end do
      end do
    end do
    end do
    end do
    end do
  end do



end subroutine ExpectedValuesM
