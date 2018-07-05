subroutine simulation()
  use Globals
  use Utils

  implicit none

  integer, parameter :: sim_gp_a = 1000
  integer :: ind_ag, ind_p, ind_z, ind_q, ind_g, ind_b
  integer, dimension(agents) :: LMstatus, new_LMstatus, entitled, new_entitled, assets, new_assets
  integer, dimension(sim_gp_a,gp_z) :: sim_N_pf
  integer, dimension(sim_gp_a,gp_z,gp_q) :: sim_W_pf
  integer, dimension(sim_gp_a,gp_z,gp_gamma,0:1) :: sim_U_pf
  real(8) :: employed, unemployed, OLF, tot_labincome, tot_income, tot_taxrev, tot_bpaid, &
             tot_assets, tot_z, tot_q, income, a, a_income, labincome, taxrev, bpaid, z, q, ap
  real(8), dimension(sim_gp_a) :: sim_a_values
  real(8), dimension(3,3) :: trans

  ! Create grid of assets for the simulation (easier to determine if distribution is statinary)
  call loggrid(min_a, max_a, sim_gp_a, sim_a_values)

  ! Reshape decision rules with new grid
  do ind_z = 1, gp_z
    call regrid(gp_a,a_values,N_pf(:,ind_z),sim_gp_a,sim_a_values,sim_N_pf(:,ind_z))
    do ind_q = 1, gp_q
      call regrid(gp_a,a_values,W_pf(:,ind_z,ind_q),sim_gp_a,sim_a_values,sim_W_pf(:,ind_z,ind_q))
    end do
    do ind_g = 1, gp_gamma
    do ind_b = 0, 1
      call regrid(gp_a,a_values,U_pf(:,ind_z,ind_g,ind_b),sim_gp_a,sim_a_values,&
                  sim_U_pf(:,ind_z,ind_g,ind_b))
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

  ! ----------
  open(unit = 71, file = "employment_rate.txt")
  open(unit = 72, file = "unemployment_rate.txt")
  open(unit = 73, file = "share_OLF.txt")
  open(unit = 74, file = "transitions.txt")
  open(unit = 75, file = "sum_assets.txt")
  open(unit = 76, file = "sum_z.txt")
  open(unit = 77, file = "sum_income.txt")
  open(unit = 78, file = "sum_taxrev.txt")
  open(unit = 79, file = "sum_bpaid.txt")
  open(unit = 80, file = "sum_entitled.txt")
  open(unit = 81, file = "sum_q.txt")
  ! ----------

  do ind_p = 1, periods-1 ! CAREFUL with last period
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

  do ind_ag = 1, agents
    ! Initialise aux variables to 0
    income = 0.d0
    labincome = 0.d0
    taxrev = 0.d0
    bpaid = 0.d0

    ! Agent realisations
    a = sim_a_values(assets(ind_ag))
    a_income = (1.d0+int_rate)*a + T ! CAREFUL: Lump-sum transfer already added here
    z = z_values(shock_z(ind_ag,ind_p))
    q = q_values(shock_q(ind_ag,ind_p))

    ! Employed
    if (LMstatus(ind_ag).eq.1) then
      ! Register labor market status
      employed = employed + 1.d0
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
        call choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                      shock_g(ind_ag,ind_p+1),new_entitled(ind_ag),new_LMstatus(ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.lambda_e).and.&
              (shock_lm(ind_ag,ind_p).le.(1.d0-sigma))) then ! No fired, no new offer
        ! Update value of q next period
        shock_q(ind_ag,ind_p+1) = shock_q(ind_ag,ind_p)
        ! Update entitlement next period
        new_entitled(ind_ag) = 0
        call choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                      shock_g(ind_ag,ind_p+1),new_entitled(ind_ag),new_LMstatus(ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.(1.d0-sigma)).and.&
              (shock_lm(ind_ag,ind_p).le.(1.d0-sigma+(lambda_u*sigma)))) then ! Fired, new offer
        ! Update entitlement next period
        new_entitled(ind_ag) = 1
        call choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                      shock_g(ind_ag,ind_p+1),new_entitled(ind_ag),new_LMstatus(ind_ag))
      else ! Fired no new offer
        ! Update entitlement next period
        new_entitled(ind_ag) = 1
        call choice_J(ap,shock_z(ind_ag,ind_p+1),shock_g(ind_ag,ind_p+1),&
                      new_entitled(ind_ag),new_LMstatus(ind_ag))
      end if

    ! Unemployed
    elseif (LMstatus(ind_ag).eq.2) then
      ! Register labor market status
      unemployed = unemployed + 1.d0
      ! Benefits
      bpaid = real(entitled(ind_ag))*(benefits(z))
      ! Next period benefits
      if ((entitled(ind_ag).eq.1).and.(shock_mu(ind_ag,ind_p).eq.1)) then
        ! Agent is entitled and does NOT get by a mu shock
        new_entitled(ind_ag) = 1
      else
        new_entitled(ind_ag) = 0
      end if
      ! Labor income - None
      ! Tax revenue
      taxrev = tau*bpaid

      ! Assets next period
      new_assets(ind_ag) = sim_U_pf(assets(ind_ag),shock_z(ind_ag,ind_p),&
                          shock_g(ind_ag,ind_p),entitled(ind_ag))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*bpaid - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.lambda_u) then ! Job offer received
        call choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                      shock_g(ind_ag,ind_p+1),new_entitled(ind_ag),new_LMstatus(ind_ag))
      else ! NO Job offer received
        call choice_J(ap,shock_z(ind_ag,ind_p+1),shock_g(ind_ag,ind_p+1),&
                      new_entitled(ind_ag),new_LMstatus(ind_ag))
      end if

    ! OLF
    elseif (LMstatus(ind_ag).eq.3) then
      ! Register labor market status
      OLF = OLF + 1.d0
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
        call choice_V(ap,shock_z(ind_ag,ind_p+1),shock_q(ind_ag,ind_p+1),&
                      shock_g(ind_ag,ind_p+1),new_entitled(ind_ag),new_LMstatus(ind_ag))
      else ! NO Job offer received
        call choice_J(ap,shock_z(ind_ag,ind_p+1),shock_g(ind_ag,ind_p+1),&
                      new_entitled(ind_ag),new_LMstatus(ind_ag))
      end if
    else
      print *, "Simulation: LMstatus not in range"
    end if

    ! Add to totals
    tot_income = tot_income + income
    tot_labincome = tot_labincome + labincome
    tot_taxrev = tot_taxrev + taxrev
    tot_bpaid = tot_bpaid + bpaid
    tot_assets = tot_assets + a
    tot_z = tot_z + z
    tot_q = tot_q + q
    trans(LMstatus(ind_ag),new_LMstatus(ind_ag)) = &
        trans(LMstatus(ind_ag),new_LMstatus(ind_ag)) + 1.d0
  end do
    ! Compute errors
    ! Check convergence
    ! Update variables
    assets = new_assets
    LMstatus = new_LMstatus
    entitled = new_entitled

    ! ------------------
    write(71,*) employed/real(agents)
    write(72,*) unemployed/(employed+unemployed)
    write(73,*) OLF/real(agents)
    write(74,*) trans(1,:)/sum(trans(1,:)), trans(2,:)/sum(trans(2,:)), trans(3,:)/sum(trans(3,:))
    write(75,*) tot_assets
    write(76,*) tot_z
    write(77,*) tot_income
    write(78,*) tot_taxrev
    write(79,*) tot_bpaid
    write(80,*) sum(real(entitled))
    write(81,*) tot_q
  end do

  close(71)
  close(72)
  close(73)
  close(74)
  close(75)
  close(76)
  close(77)
  close(78)
  close(79)
  close(80)
  close(81)

  ! Print results
  ! print *, "Transition:"
  print '(3f7.4)', trans(1,:)/sum(trans(1,:))
  print '(3f7.4)', trans(2,:)/sum(trans(2,:))
  print '(3f7.4)', trans(3,:)/sum(trans(3,:))
  ! print *, "Employment rate:", employed/real(agents)
  ! print *, "Unemployment rate:", unemployed/(employed+unemployed)
  ! print *, "Share OLF:", OLF/real(agents)

end subroutine simulation

!===== CHOICE V ===================================================================================
subroutine choice_V(val_a,ind_z,ind_q,ind_g,ind_b,choice)
  use Globals
  use Utils
  implicit none
  integer, intent(in) :: ind_z, ind_q, ind_g, ind_b
  integer, intent(out) :: choice
  real(8), intent(in) :: val_a
  real(8), dimension(1:3) :: aux_choice

  ! Interpolate value of being OLF
  call my_inter(a_values,N_vf(:,ind_z),gp_a,val_a,aux_choice(3))

  ! Interpolate value of being unemployed
  call my_inter(a_values,U_vf(:,ind_z,ind_g,ind_b),gp_a,val_a,aux_choice(2))

  ! Interpolate value of being employed
  call my_inter(a_values,W_vf(:,ind_z,ind_q),gp_a,val_a,aux_choice(1))

  choice = maxloc(aux_choice, dim=1)
end subroutine choice_V

!===== CHOICE J ===================================================================================
subroutine choice_J(val_a,ind_z,ind_g,ind_b,choice)
  use Globals
  use Utils
  implicit none
  integer, intent(in) :: ind_z, ind_g, ind_b
  integer, intent(out) :: choice
  real(8), intent(in) :: val_a
  real(8) :: aux_U, aux_N

  ! Interpolate value of being OLF
  call my_inter(a_values,N_vf(:,ind_z),gp_a,val_a,aux_N)

  ! Interpolate value of being unemployed
  call my_inter(a_values,U_vf(:,ind_z,ind_g,ind_b),gp_a,val_a,aux_U)

  ! Compare values
  if (aux_U.ge.aux_N) then ! Unemployment is preferred
    choice = 2
  else
    choice = 3
  end if
end subroutine choice_J

!===== REGRID =====================================================================================
subroutine regrid(gp_old,grid_old,pf_old,gp_new,grid_new,pf_new)
  ! Uses interpolation to reshape a policy function to a new grid
  use Utils
  implicit none

  integer, intent(in) :: gp_old, gp_new
  integer :: ind
  integer, dimension(gp_new) :: pf_new
  real(8) :: my_x, my_y
  real(8), dimension(gp_old), intent(in) :: grid_old, pf_old
  real(8), dimension(gp_new), intent(in) :: grid_new

  do ind = 1, gp_new
    ! Get the value from the new grid
    my_x = grid_new(ind)

    ! Interpolate using the old policy function
    call my_inter(grid_old,pf_old,gp_old,my_x,my_y)

    ! Find closest position in the new grid for my_y
    pf_new(ind) = my_closest(grid_new,gp_new,my_y)
  end do
end subroutine regrid

!===== REALISE SHOCKS =============================================================================
subroutine realise_shocks(agents, periods, trans_matrix, gp_tm, shocks)
  ! Computes the realised shock implied by trans_matrix, and returns it in shocks

  implicit none

  integer, intent(in) :: agents, periods, gp_tm
  integer :: ind_ag, ind_p, ind_sh
  integer, dimension(agents, periods), intent(out) :: shocks
  real(8) :: aux_sum
  real(8), dimension(gp_tm) :: ss_dist, aux_vec
  real(8), dimension(gp_tm,gp_tm), intent(in) :: trans_matrix
  real(8), dimension(agents, periods) :: aux_shocks

  call random_number(aux_shocks)

  ! Find steady state distribution
  ss_dist = my_ss(trans_matrix, gp_tm)

  ! Simulate for all periods
  do ind_ag = 1, agents
  do ind_p = 1, periods
    ! Chose auxiliary vector for first period
    if (ind_p.eq.1) then
      aux_vec = ss_dist
    else
      aux_vec = trans_matrix(shocks(ind_ag,ind_p-1),:)
    end if
    ! Find out shock
    aux_sum = 0.d0
    do ind_sh = 1, gp_tm
      aux_sum  = aux_sum + aux_vec(ind_sh)
      if (aux_sum.ge.aux_shocks(ind_ag,ind_p)) exit
    end do
    shocks(ind_ag,ind_p) = ind_sh
  end do
  end do

CONTAINS
  function my_ss(tmatrix,gp)
    implicit none
    integer :: gp
    integer :: row, col, iter
    real(8) :: aux_sum
    real(8), dimension(gp) :: my_ss
    real(8), dimension(gp) :: dist, ndist
    real(8), dimension(gp,gp) :: tmatrix

    ! Initialise distribution
    dist = 1.d0/real(gp)

    do iter = 1, 10000
      ndist = 0.d0
      do col = 1, gp
        do row = 1, gp
          ndist(col) = ndist(col) + (dist(row)*tmatrix(row,col))
        end do
      end do
      aux_sum = sum(abs(ndist-dist))
      dist = ndist
      if (aux_sum.lt.1.0d-8) then
        exit
      end if
    end do

    my_ss = dist
  end function my_ss
end subroutine realise_shocks

!===== REALISE SIMPLE IID SHOCKS ==================================================================
subroutine realise_iid(agents, periods, param_value, shocks)
  ! Computes the realised shock of an iid process

  implicit none

  integer, intent(in) :: agents, periods
  integer :: ind_ag, ind_p
  integer, dimension(agents, periods), intent(out) :: shocks
  real(8), intent(in) :: param_value
  real(8), dimension(agents, periods) :: aux_shocks

  shocks = 0

  call random_number(aux_shocks)

  do ind_ag = 1, agents
  do ind_p = 1, periods
    if (aux_shocks(ind_ag,ind_p).le.param_value) then
      shocks(ind_ag,ind_p) = 1
    end if
  end do
  end do
end subroutine realise_iid

!===== RANDOM INTEGERS ============================================================================
subroutine random_integers(rows, columns, lb, ub, int_matrix)
  ! Returns a int_matrix(rows, columns) of random integer number numbers in [lb,ub]
  integer, intent(in) :: rows, columns, lb, ub
  integer, dimension(rows, columns), intent(out) :: int_matrix
  real(8), dimension(rows, columns) :: aux_mat

  call random_number(aux_mat)

  int_matrix = lb + floor(real(ub+1-lb)*aux_mat)
end subroutine random_integers
