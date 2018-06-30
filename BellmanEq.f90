!===== OLF ========================================================================================
subroutine value_N(ind_a, aux_exp, pf, vf)
  use Globals
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
    call my_inter(a_values, aux_exp, gp_a, ap, aux_inter)

    ! Compute consumption
    cons = income - ap

    valor_N = u(cons) + beta*aux_inter
  end function valor_N
end subroutine value_N

!===== EMPLOYED ===================================================================================
subroutine value_W(ind_a, ind_z, ind_q, aux_exp, pf, vf)
  use Globals
  use Utils
  implicit none

  integer, intent(in) :: ind_a, ind_z, ind_q
  real(8) :: aux_max_a, income
  real(8), dimension(gp_a), intent(in) :: aux_exp
  real(8), intent(out) :: pf, vf

  ! Compute income
  income = (1.d0+int_rate)*a_values(ind_a) + T + (1.d0-tau)*wage*z_values(ind_z)*q_values(ind_q)

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
    call my_inter(a_values, aux_exp, gp_a, ap, aux_inter)

    ! Compute consumption
    cons = income - ap

    valor_W = u(cons) - alpha + beta*aux_inter
  end function valor_W
end subroutine value_W

!===== EMPLOYED ===================================================================================
subroutine value_U(ind_a, ind_z, ind_g, ind_b, aux_exp, pf, vf)
  use Globals
  use Utils
  implicit none

  integer, intent(in) :: ind_a, ind_z, ind_g, ind_b
  real(8) :: aux_max_a, income
  real(8), dimension(gp_a), intent(in) :: aux_exp
  real(8), intent(out) :: pf, vf

  ! Compute income
  income = (1.d0+int_rate)*a_values(ind_a) + T + (1.d0-tau)*IB_values(ind_b)*benefits(z_values(ind_z))

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
    call my_inter(a_values, aux_exp, gp_a, ap, aux_inter)

    ! Compute consumption
    cons = income - ap

    valor_U = u(cons) - gamma_values(ind_g) + beta*aux_inter
  end function valor_U
end subroutine value_U
