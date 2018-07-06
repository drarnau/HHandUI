subroutine simulation()
  use Globals
  use Utils

  implicit none

  integer, parameter :: sim_gp_a = 1000
  integer :: ind_ag, ind_p, ind_z, ind_q, ind_g, ind_b, reps
  integer, dimension(agents) :: LMstatus, new_LMstatus, entitled, new_entitled, assets, new_assets
  integer, dimension(sim_gp_a,gp_z) :: sim_N_pf
  integer, dimension(sim_gp_a,gp_z,gp_q) :: sim_W_pf
  integer, dimension(sim_gp_a,gp_z,gp_gamma,0:1) :: sim_U_pf
  real(8) :: employed, unemployed, OLF, tot_labincome, tot_income, tot_taxrev, tot_bpaid, &
             tot_assets, tot_z, tot_q, income, a, a_income, labincome, taxrev, bpaid, z, q, ap
  real(8), dimension(sim_gp_a) :: sim_a_values
  real(8), dimension(3,3) :: trans

  ! Create grid of assets for the simulation (easier to determine if distribution is statinary)
  sim_a_values = loggrid(min_a, max_a, sim_gp_a)

  ! Reshape decision rules with new grid
  do ind_z = 1, gp_z
    sim_N_pf(:,ind_z) = regrid(gp_a,a_values,N_pf(:,ind_z),sim_gp_a,sim_a_values)
    do ind_q = 1, gp_q
      sim_W_pf(:,ind_z,ind_q) = regrid(gp_a,a_values,W_pf(:,ind_z,ind_q),sim_gp_a,sim_a_values)
    end do
    do ind_g = 1, gp_gamma
    do ind_b = 0, 1
      sim_U_pf(:,ind_z,ind_g,ind_b) = regrid(gp_a,a_values,U_pf(:,ind_z,ind_g,ind_b),&
                                             sim_gp_a,sim_a_values)
    end do
    end do
  end do

  ! Initialise everybody with 0 assets, OLF, and not entitled
  assets = 1
  new_assets = -1

  LMstatus = 3
  new_LMstatus = -1

  entitled = 0
  new_entitled = -1

  ! Initialise macro variables to 0
  employed = 0.d0
  unemployed = 0.d0
  OLF = 0.d0
  tot_labincome = 0.d0
  tot_income = 0.d0
  tot_taxrev = 0.d0
  tot_bpaid = 0.d0
  tot_assets = 0.d0
  tot_z = 0.d0
  tot_q = 0.d0
  trans = 0.d0
  reps = 0

  do ind_p = 1, periods-1 ! CAREFUL with last period
  do ind_ag = 1, agents
    ! Initialise aux variables to 0
    income = 0.d0
    labincome = 0.d0
    taxrev = 0.d0
    bpaid = 0.d0

    ! Quantify shocks
    a = sim_a_values(assets(ind_ag))
    a_income = (1.d0+int_rate)*a + T ! CAREFUL: Lump-sum transfer already added here
    z = z_values(shock_z(ind_ag,ind_p))
    q = q_values(shock_q(ind_ag,ind_p))

    ! Employed
    if (LMstatus(ind_ag).eq.1) then
      ! Benefits - none
      ! Labor income
      labincome = wage*z*q
      ! Tax revenue
      taxrev = tau*labincome

      ! Assets next period
      new_assets(ind_ag) = sim_W_pf(assets(ind_ag),shock_z(ind_ag,ind_p),shock_q(ind_ag,ind_p))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*wage*z*q - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.lambda_e) then ! New offer
        ! Update value of q next period
        shock_q(ind_ag,ind_p+1) = max(shock_q(ind_ag,ind_p),shock_q(ind_ag,ind_p+1))
        ! Update entitlement next period
        new_entitled(ind_ag) = 0
        new_LMstatus(ind_ag) = choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                                        shock_g(ind_ag,ind_p+1),new_entitled(ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.lambda_e).and.&
              (shock_lm(ind_ag,ind_p).le.(1.d0-sigma))) then ! No fired, no new offer
        ! Update value of q next period
        shock_q(ind_ag,ind_p+1) = shock_q(ind_ag,ind_p+1)
        ! Update entitlement next period
        new_entitled(ind_ag) = 0
        new_LMstatus(ind_ag) = choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                                        shock_g(ind_ag,ind_p+1),new_entitled(ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.(1.d0-sigma)).and.&
              (shock_lm(ind_ag,ind_p).le.(1.d0-sigma+(lambda_u*sigma)))) then ! Fired, new offer
        ! Update entitlement next period
        new_entitled(ind_ag) = 1
        new_LMstatus(ind_ag) = choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                                        shock_g(ind_ag,ind_p+1),new_entitled(ind_ag))
      else ! Fired no new offer
        ! Update entitlement next period
        new_entitled(ind_ag) = 1
        new_LMstatus(ind_ag) = choice_J(ap,shock_z(ind_ag,ind_p+1),shock_g(ind_ag,ind_p+1),&
                                        new_entitled(ind_ag))
      end if

    ! Unemployed
    elseif (LMstatus(ind_ag).eq.2) then
      ! Benefits
      bpaid = real(entitled(ind_ag))*(benefits(z))
      ! Next period benefits
      if ((entitled(ind_ag).eq.1).and.(shock_mu(ind_ag,ind_p).ne.1)) then
        ! Agent is entitled and does NOT get hit by a mu shock
        new_entitled(ind_ag) = 1
      else
        new_entitled(ind_ag) = 0
      end if
      ! Labor income - None
      ! Tax revenue
      taxrev = tau*bpaid

      ! Assets next period
      new_assets(ind_ag) = sim_U_pf(assets(ind_ag),shock_z(ind_ag,ind_p),shock_g(ind_ag,ind_p),&
                                    entitled(ind_ag))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*bpaid - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.lambda_u) then ! Job offer received
        new_LMstatus(ind_ag) = choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                                        shock_g(ind_ag,ind_p+1),new_entitled(ind_ag))
      else ! NO Job offer received
        new_LMstatus(ind_ag) = choice_J(ap,shock_z(ind_ag,ind_p+1),shock_g(ind_ag,ind_p+1),&
                                        new_entitled(ind_ag))
      end if

    ! OLF
    elseif (LMstatus(ind_ag).eq.3) then
      ! Benefits - None
      new_entitled(ind_ag) = 0
      ! Labor income - None
      ! Tax revenue - None

      ! Assets next period
      new_assets(ind_ag) = sim_N_pf(assets(ind_ag),shock_z(ind_ag,ind_p))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.lambda_n) then ! Job offer received
        new_LMstatus(ind_ag) = choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                                        shock_g(ind_ag,ind_p+1),new_entitled(ind_ag))
      else ! NO Job offer received
        new_LMstatus(ind_ag) = choice_J(ap,shock_z(ind_ag,ind_p+1),shock_g(ind_ag,ind_p+1),&
                                        new_entitled(ind_ag))
      end if
    else
      print *, "Simulation: LMstatus not in range"
    end if

    ! Add to totals
    if (ind_p.ge.(periods/2)) then
      tot_income = tot_income + income
      tot_labincome = tot_labincome + labincome
      tot_taxrev = tot_taxrev + taxrev
      tot_bpaid = tot_bpaid + bpaid
      tot_assets = tot_assets + a
      trans(LMstatus(ind_ag),new_LMstatus(ind_ag)) = &
          trans(LMstatus(ind_ag),new_LMstatus(ind_ag)) + 1.d0
      ! Add to totals relevant for employed
      if (LMstatus(ind_ag).eq.1) then
        ! Register labor market status
        employed = employed + 1.d0
        tot_z = tot_z + z
        tot_q = tot_q + q
      elseif(LMstatus(ind_ag).eq.2) then
        ! Register labor market status
        unemployed = unemployed + 1.d0
      else
        ! Register labor market status
        OLF = OLF + 1.d0
      end if
    end if
  end do
    if (ind_p.ge.(periods/2)) then
      reps = reps + 1
    end if
    ! Update variables
    assets = new_assets
    LMstatus = new_LMstatus
    entitled = new_entitled
  end do

  ! Compute averages
  tot_income = tot_income/real(reps)
  tot_labincome = tot_labincome/real(reps)
  tot_taxrev = tot_taxrev/real(reps)
  tot_bpaid = tot_bpaid/real(reps)
  tot_assets = tot_assets/real(reps)
  tot_z = tot_z/real(reps)
  tot_q = tot_q/real(reps)
  employed = employed/real(reps)
  unemployed = unemployed/real(reps)
  OLF = OLF/real(reps)

  ! Compute results
  new_KLratio = tot_assets/tot_z
  new_average_z = tot_z/employed
  new_T = (tot_taxrev-tot_bpaid)/real(agents)
  Erate = employed/real(agents)
  Urate = unemployed/(employed+unemployed)
  Nrate = OLF/real(agents)
  transitions(1, :) = trans(1,:)/sum(trans(1,:))
  transitions(2, :) = trans(2,:)/sum(trans(2,:))
  transitions(3, :) = trans(3,:)/sum(trans(3,:))

CONTAINS
  !----- CHOICE V ---------------------------------------------------------------------------------
  integer function choice_V(val_a,ind_z,ind_q,ind_g,ind_b)
    implicit none
    integer :: ind_z, ind_q, ind_g, ind_b
    real(8) :: val_a
    real(8), dimension(1:3) :: aux_choice

    ! Interpolate value of being OLF
    aux_choice(3) = my_inter(a_values,N_vf(:,ind_z),gp_a,val_a)

    ! Interpolate value of being unemployed
    aux_choice(2) = my_inter(a_values,U_vf(:,ind_z,ind_g,ind_b),gp_a,val_a)

    ! Interpolate value of being employed
    aux_choice(1) = my_inter(a_values,W_vf(:,ind_z,ind_q),gp_a,val_a)

    choice_V = maxloc(aux_choice, dim=1)
  end function choice_V

  !----- CHOICE J ---------------------------------------------------------------------------------
  integer function choice_J(val_a,ind_z,ind_g,ind_b)
    implicit none
    integer :: ind_z, ind_g, ind_b
    real(8) :: val_a
    real(8) :: aux_U, aux_N

    ! Interpolate value of being OLF
    aux_N = my_inter(a_values,N_vf(:,ind_z),gp_a,val_a)

    ! Interpolate value of being unemployed
    aux_U = my_inter(a_values,U_vf(:,ind_z,ind_g,ind_b),gp_a,val_a)

    ! Compare values
    if (aux_U.ge.aux_N) then ! Unemployment is preferred
      choice_J = 2
    else
      choice_J = 3
    end if
  end function choice_J

  !----- REGRID A DECISION RULE -------------------------------------------------------------------
  function regrid(gp_old,grid_old,pf_old,gp_new,grid_new)
    ! Uses interpolation to reshape a policy function to a new grid
    use Utils
    implicit none

    integer :: gp_old, gp_new
    integer :: ind
    integer, dimension(gp_new) :: regrid, pf_new
    real(8) :: my_x, my_y
    real(8), dimension(gp_old) :: grid_old, pf_old
    real(8), dimension(gp_new) :: grid_new

    do ind = 1, gp_new
      ! Get the value from the new grid
      my_x = grid_new(ind)

      ! Interpolate using the old policy function
      my_y = my_inter(grid_old,pf_old,gp_old,my_x)

      ! Find closest position in the new grid for my_y
      pf_new(ind) = my_closest(grid_new,gp_new,my_y)
    end do

    regrid = pf_new
  end function regrid
end subroutine simulation
