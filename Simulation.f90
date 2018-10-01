!===== SIMULATION SINGLES =========================================================================
subroutine SimSingles(mysex, gen_output)
  use Globals
  use GlobalsSingles
  use Utils
  use csv_file

  implicit none

  logical, intent(in) :: gen_output
  integer, intent(in) :: mysex
  integer, parameter :: sim_gp_a = 1000
  integer :: ind_ag, ind_p, ind_z, ind_g, ind_b, reps
  integer, dimension(agents) :: LMstatus, new_LMstatus, entitled, new_entitled, assets, new_assets
  integer, dimension(agents,periods) :: myshock_z, myshock_g
  integer, dimension(sim_gp_a,gp_z) :: sim_N_pf, sim_W_pf
  integer, dimension(sim_gp_a,gp_z,gp_gamma,0:1) :: sim_U_pf
  real(8) :: employed, unemployed, OLF, tot_labincome, tot_income, tot_taxrev, tot_bpaid, &
             tot_assets, tot_z, income, a, a_income, labincome, taxrev, bpaid, z, ap, top_assets
  real(8), dimension(sim_gp_a) :: sim_a_values
  real(8), dimension(gp_z) :: aux_vec
  real(8), dimension(3,3) :: trans
  character(len=1024) :: aux_name

  ! Create grid of assets for the simulation (easier to determine if distribution is stationary)
  sim_a_values = loggrid(min_a, max_a, sim_gp_a)

  ! Reshape decision rules with new grid
  do ind_z = 1, gp_z
    sim_N_pf(:,ind_z) = regrid(gp_a,a_values,N_pf(:,ind_z),sim_gp_a,sim_a_values)
    sim_W_pf(:,ind_z) = regrid(gp_a,a_values,W_pf(:,ind_z),sim_gp_a,sim_a_values)
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

  ! Initialise myshocks
  myshock_z = -1
  myshock_g = -1

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
  trans = 0.d0
  reps = 0
  top_assets = 0.d0

  ! Open CSV file where to store simulation
  if (gen_output) then
    write (aux_name,"(A11,I1,A4)") "simulation_", mysex, ".csv"
    open(unit=14, file=aux_name,status='unknown')
    call csv_write(14,"ID",.false.)
    call csv_write(14,"Period",.false.)
    call csv_write(14,"Status today",.false.)
    call csv_write(14,"Status tomorrow",.false.)
    call csv_write(14,"Labour Income",.false.)
    call csv_write(14,"Taxes Paid",.false.)
    call csv_write(14,"Benefits Received",.false.)
    call csv_write(14,"Wealth",.true.)
  end if

  do ind_p = 1, periods-1 ! CAREFUL with last period
  do ind_ag = 1, agents
    ! Initialise aux variables to 0
    income = 0.d0
    labincome = 0.d0
    taxrev = 0.d0
    bpaid = 0.d0

    ! Assign shocks
    if (ind_p.eq.1) then
      myshock_z(ind_ag,ind_p) = realise(shock_z(ind_ag,ind_p),z_ssdist,gp_z)
      myshock_g(ind_ag,ind_p) = realise(shock_g(mysex,ind_ag,ind_p),gamma_trans, gp_gamma)
    end if
    aux_vec = z_trans(myshock_z(ind_ag,ind_p),:)
    myshock_z(ind_ag,ind_p+1) = realise(shock_z(ind_ag,ind_p+1),aux_vec,gp_z)
    myshock_g(ind_ag,ind_p+1) = realise(shock_g(mysex,ind_ag,ind_p+1),gamma_trans, gp_gamma)

    ! Quantify shocks
    a = sim_a_values(assets(ind_ag))
    a_income = (1.d0+int_rate)*a + T ! CAREFUL: Lump-sum transfer already added here
    z = z_values(myshock_z(ind_ag,ind_p))

    ! Employed
    if (LMstatus(ind_ag).eq.1) then
      ! Benefits - none
      ! Labor income
      labincome = wage*z
      ! Tax revenue
      taxrev = tau*labincome

      ! Assets next period
      new_assets(ind_ag) = sim_W_pf(assets(ind_ag),myshock_z(ind_ag,ind_p))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*wage*z - ap

      ! LM status next period
      if ((shock_lm(ind_ag,ind_p).le.(1.d0-sigma))) then ! No fired
        ! Update entitlement next period
        new_entitled(ind_ag) = 0
        new_LMstatus(ind_ag) = choice_V(ap,myshock_z(ind_ag,ind_p+1),&
                                        myshock_g(ind_ag,ind_p+1),new_entitled(ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.(1.d0-sigma)).and.&
              (shock_lm(ind_ag,ind_p).le.(1.d0-sigma+(lambda_u*sigma)))) then ! Fired, new offer
        ! Update entitlement next period
        new_entitled(ind_ag) = 1
        new_LMstatus(ind_ag) = choice_V(ap,myshock_z(ind_ag,ind_p+1),&
                                        myshock_g(ind_ag,ind_p+1),new_entitled(ind_ag))
      else ! Fired no new offer
        ! Update entitlement next period
        new_entitled(ind_ag) = 1
        new_LMstatus(ind_ag) = choice_J(ap,myshock_z(ind_ag,ind_p+1),myshock_g(ind_ag,ind_p+1),&
                                        new_entitled(ind_ag))
      end if

    ! Unemployed
    elseif (LMstatus(ind_ag).eq.2) then
      ! Benefits
      bpaid = real(entitled(ind_ag))*(benefits(z))
      ! Next period benefits
      if ((entitled(ind_ag).eq.1).and.(shock_mu(mysex,ind_ag,ind_p).gt.mu)) then
        ! Agent is entitled and does NOT get hit by a mu shock
        new_entitled(ind_ag) = 1
      else
        new_entitled(ind_ag) = 0
      end if
      ! Labor income - None
      ! Tax revenue
      taxrev = tau*bpaid

      ! Assets next period
      new_assets(ind_ag) = sim_U_pf(assets(ind_ag),myshock_z(ind_ag,ind_p),myshock_g(ind_ag,ind_p),&
                                    entitled(ind_ag))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*bpaid - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.lambda_u) then ! Job offer received
        new_LMstatus(ind_ag) = choice_V(ap,myshock_z(ind_ag,ind_p+1),&
                                        myshock_g(ind_ag,ind_p+1),new_entitled(ind_ag))
      else ! NO Job offer received
        new_LMstatus(ind_ag) = choice_J(ap,myshock_z(ind_ag,ind_p+1),myshock_g(ind_ag,ind_p+1),&
                                        new_entitled(ind_ag))
      end if

    ! OLF
    elseif (LMstatus(ind_ag).eq.3) then
      ! Benefits - None
      new_entitled(ind_ag) = 0
      ! Labor income - None
      ! Tax revenue - None

      ! Assets next period
      new_assets(ind_ag) = sim_N_pf(assets(ind_ag),myshock_z(ind_ag,ind_p))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.lambda_n) then ! Job offer received
        new_LMstatus(ind_ag) = choice_V(ap,myshock_z(ind_ag,ind_p+1),&
                                        myshock_g(ind_ag,ind_p+1),new_entitled(ind_ag))
      else ! NO Job offer received
        new_LMstatus(ind_ag) = choice_J(ap,myshock_z(ind_ag,ind_p+1),myshock_g(ind_ag,ind_p+1),&
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
      elseif(LMstatus(ind_ag).eq.2) then
        ! Register labor market status
        unemployed = unemployed + 1.d0
      else
        ! Register labor market status
        OLF = OLF + 1.d0
      end if
    end if

    ! Print to CSV
    if (gen_output) then
      if (ind_p.ge.(periods-printed)) then
        call csv_write(14,real(ind_ag),.false.)
        call csv_write(14,real(ind_p),.false.)
        call csv_write(14,real(LMstatus(ind_ag)),.false.)
        call csv_write(14,real(new_LMstatus(ind_ag)),.false.)
        call csv_write(14,labincome,.false.)
        call csv_write(14,taxrev,.false.)
        call csv_write(14,bpaid,.false.)
        call csv_write(14,a,.true.)
      end if
    end if

    ! Agents using top assets
    if (assets(ind_ag).eq.sim_gp_a) then
      top_assets = top_assets + 1.d0
    end if
  end do ! Agents
    if (ind_p.ge.(periods/2)) then
      reps = reps + 1
    end if
    ! Update variables
    assets = new_assets
    LMstatus = new_LMstatus
    entitled = new_entitled
  end do ! Periods

  ! Close csv file
  if (gen_output) then
    close(14)
  end if

  ! Compute averages
  tot_income = tot_income/real(reps)
  tot_labincome = tot_labincome/real(reps)
  tot_taxrev = tot_taxrev/real(reps)
  tot_bpaid = tot_bpaid/real(reps)
  tot_assets = tot_assets/real(reps)
  tot_z = tot_z/real(reps)
  employed = employed/real(reps)
  unemployed = unemployed/real(reps)
  OLF = OLF/real(reps)

  ! Compute results
  wealth(mysex) = tot_assets
  aux_KLratio(mysex) = tot_assets/tot_z
  aux_average_z(mysex) = tot_z/employed
  aux_T(mysex) = (tot_taxrev-tot_bpaid)/real(agents)
  Erate(single,mysex) = employed/real(agents)
  Urate(single,mysex) = unemployed/(employed+unemployed)
  Nrate(single,mysex) = OLF/real(agents)
  transitions(single,mysex,1, :) = trans(1,:)/sum(trans(1,:))
  transitions(single,mysex,2, :) = trans(2,:)/sum(trans(2,:))
  transitions(single,mysex,3, :) = trans(3,:)/sum(trans(3,:))

  ! Print share of agents using top assets
  print '(a,I2,a,f7.4)', " Share of single", mysex, " HH using top assets:",&
                                                              top_assets/(real(reps*agents))
  print *, ""


CONTAINS
  !----- CHOICE V ---------------------------------------------------------------------------------
  integer function choice_V(val_a,ind_z,ind_g,ind_b)
    implicit none
    integer :: ind_z, ind_g, ind_b
    real(8) :: val_a
    real(8), dimension(1:3) :: aux_choice

    ! Interpolate value of being OLF
    aux_choice(3) = my_inter(a_values,N_vf(:,ind_z),gp_a,val_a)

    ! Interpolate value of being unemployed
    aux_choice(2) = my_inter(a_values,U_vf(:,ind_z,ind_g,ind_b),gp_a,val_a)

    ! Interpolate value of being employed
    aux_choice(1) = my_inter(a_values,W_vf(:,ind_z),gp_a,val_a)

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
end subroutine SimSingles

!===== SIMULATION MARRIED =========================================================================
subroutine SimMarried(gen_output)
  use Globals
  use GlobalsMarried
  use Utils
  use csv_file

  implicit none

  logical, intent(in) :: gen_output
  integer, parameter :: sim_gp_a = 1000
  integer :: ind_ag, ind_p, ind_z, ind_g, ind_b, ind_g_f, ind_b_f, reps, mysex
  integer, dimension(agents) :: assets, new_assets
  integer, dimension(1:2,agents) :: LMstatus, new_LMstatus, entitled, new_entitled
  integer, dimension(agents,periods) :: myshock_z
  integer, dimension(1:2,agents,periods) :: myshock_g
  integer, dimension(1:9,1:2) :: mystates
  integer, dimension(sim_gp_a,gp_z2) :: sim_NN_pf, sim_WW_pf, sim_WN_pf, sim_NW_pf
  integer, dimension(sim_gp_a,gp_z2,gp_gamma,0:1) :: sim_WU_pf, sim_UW_pf, sim_NU_pf, sim_UN_pf
  integer, dimension(sim_gp_a,gp_z2,gp_gamma,gp_gamma,0:1,0:1) :: sim_UU_pf

  real(8) :: tot_labincome, tot_income, tot_taxrev, tot_bpaid, tot_assets, tot_z, income, a, &
             a_income, labincome, taxrev, bpaid, ap, top_assets
  real(8), dimension(1:2) :: z, employed, unemployed, OLF
  real(8), dimension(:), allocatable :: aux_vec
  real(8), dimension(sim_gp_a) :: sim_a_values
  real(8), dimension(1:9) :: aux_WW
  real(8), dimension(1:6) :: aux_WU, aux_WN, aux_UW, aux_NW
  real(8), dimension(1:4) :: aux_UU, aux_UN, aux_NU, aux_NN
  real(8), dimension(1:2,1:3,1:3) :: itrans
  real(8), dimension(1:3,1:3,1:3,1:3) :: jtrans

  ! Create an auxilary array to sort joint states
  ! WW
  mystates(1,male) = 1
  mystates(1,female) = 1
  ! WU
  mystates(2,male) = 1
  mystates(2,female) = 2
  ! WN
  mystates(3,male) = 1
  mystates(3,female) = 3
  ! UW
  mystates(4,male) = 2
  mystates(4,female) = 1
  ! UU
  mystates(5,male) = 2
  mystates(5,female) = 2
  ! UN
  mystates(6,male) = 2
  mystates(6,female) = 3
  ! NW
  mystates(7,male) = 3
  mystates(7,female) = 1
  ! NU
  mystates(8,male) = 3
  mystates(8,female) = 2
  ! NN
  mystates(9,male) = 3
  mystates(9,female) = 3

  ! Probabilities for WW household
  aux_WW = 0.d0
  ! Male keeps, female keeps
  aux_WW(1) = (1.d0-sigma(male))*(1.d0-sigma(female))
  ! Male keeps, female fired, finds new
  aux_WW(2) = aux_WW(1) + (1.d0-sigma(male))*sigma(female)*lambda_u(female)
  ! Male keeps, female fired, no new
  aux_WW(3) = aux_WW(2) + (1.d0-sigma(male))*sigma(female)*(1.d0-lambda_u(female))
  ! Male fired, find new, female keeps
  aux_WW(4) = aux_WW(3) + sigma(male)*lambda_u(male)*(1.d0-sigma(female))
  ! Male fired, finds new, female fired, finds new
  aux_WW(5) = aux_WW(4) + sigma(male)*lambda_u(male)*sigma(female)*lambda_u(female)
  ! Male fired, finds new, female fired, no new
  aux_WW(6) = aux_WW(5) + sigma(male)*lambda_u(male)*sigma(female)*(1.d0-lambda_u(female))
  ! Male fired, no new, female keeps
  aux_WW(7) = aux_WW(6) + sigma(male)*(1.d0-lambda_u(male))*(1.d0-sigma(female))
  ! Male fired, no new, female fired, finds new
  aux_WW(8) = aux_WW(7) + sigma(male)*(1.d0-lambda_u(male))*sigma(female)*lambda_u(female)
  ! Male fired, no new, female fired, no new
  aux_WW(9) = aux_WW(8) + sigma(male)*(1.d0-lambda_u(male))*sigma(female)*(1.d0-lambda_u(female))

  ! Probabilities for WU household
  aux_WU = 0.d0
  ! Male keeps, female finds
  aux_WU(1) = (1.d0-sigma(male))*lambda_u(female)
  ! Male keeps, female no finds
  aux_WU(2) = aux_WU(1) + (1.d0-sigma(male))*(1.d0-lambda_u(female))
  ! Male fired, finds new, female finds
  aux_WU(3) = aux_WU(2) + sigma(male)*lambda_u(male)*lambda_u(female)
  ! Male fired, finds new, female no finds
  aux_WU(4) = aux_WU(3) + sigma(male)*lambda_u(male)*(1.d0-lambda_u(female))
  ! Male fired, no new, female finds
  aux_WU(5) = aux_WU(4) + sigma(male)*(1.d0-lambda_u(male))*lambda_u(female)
  ! Male fired, no new, female no finds
  aux_WU(6) = aux_WU(5) + sigma(male)*(1.d0-lambda_u(male))*(1.d0-lambda_u(female))

  ! Probabilities for WN household
  aux_WN = 0.d0
  ! Male keeps, female finds
  aux_WN(1) = (1.d0-sigma(male))*lambda_n(female)
  ! Male keeps, female no finds
  aux_WN(2) = aux_WN(1) + (1.d0-sigma(male))*(1.d0-lambda_n(female))
  ! Male fired, finds new, female finds
  aux_WN(3) = aux_WN(2) + sigma(male)*lambda_u(male)*lambda_n(female)
  ! Male fired, finds new, female no finds
  aux_WN(4) = aux_WN(3) + sigma(male)*lambda_u(male)*(1.d0-lambda_n(female))
  ! Male fired, no new, female finds
  aux_WN(5) = aux_WN(4) + sigma(male)*(1.d0-lambda_u(male))*lambda_n(female)
  ! Male fired, no new, female no finds
  aux_WN(6) = aux_WN(5) + sigma(male)*(1.d0-lambda_u(male))*(1.d0-lambda_n(female))

  ! Probabilities for UW household
  aux_UW = 0.d0
  ! Female keeps, male finds
  aux_UW(1) = (1.d0-sigma(female))*lambda_u(male)
  ! Female keeps, male no finds
  aux_UW(2) = aux_UW(1) + (1.d0-sigma(female))*(1.d0-lambda_u(male))
  ! Female fired, finds new, male finds
  aux_UW(3) = aux_UW(2) + sigma(female)*lambda_u(female)*lambda_u(male)
  ! Female fired, finds new, male no finds
  aux_UW(4) = aux_UW(3) + sigma(female)*lambda_u(female)*(1.d0-lambda_u(male))
  ! Female fired, no new, male finds
  aux_UW(5) = aux_UW(4) + sigma(female)*(1.d0-lambda_u(female))*lambda_u(male)
  ! Female fired, no new, male no finds
  aux_UW(6) = aux_UW(5) + sigma(female)*(1.d0-lambda_u(female))*(1.d0-lambda_u(male))

  ! Probabilities for UU household
  aux_UU = 0.d0
  ! Male finds, female finds
  aux_UU(1) = lambda_u(male)*lambda_u(female)
  ! Male finds, female no finds
  aux_UU(2) = aux_UU(1) + lambda_u(male)*(1.d0-lambda_u(female))
  ! Male no finds, female no finds
  aux_UU(3) = aux_UU(2) + (1.d0-lambda_u(male))*lambda_u(female)
  ! Male no finds, female no finds
  aux_UU(4) = aux_UU(3) + (1.d0-lambda_u(male))*(1.d0-lambda_u(female))

  ! Probabilities for UN household
  aux_UN = 0.d0
  ! Male finds, female finds
  aux_UN(1) = lambda_u(male)*lambda_n(female)
  ! Male finds, female no finds
  aux_UN(2) = aux_UN(1) + lambda_u(male)*(1.d0-lambda_n(female))
  ! Male no finds, female no finds
  aux_UN(3) = aux_UN(2) + (1.d0-lambda_u(male))*lambda_n(female)
  ! Male no finds, female no finds
  aux_UN(4) = aux_UN(3) + (1.d0-lambda_u(male))*(1.d0-lambda_n(female))

  ! Probabilities for NW household
  aux_NW = 0.d0
  ! Female keeps, male finds
  aux_NW(1) = (1.d0-sigma(female))*lambda_n(male)
  ! Female keeps, male no finds
  aux_NW(2) = aux_NW(1) + (1.d0-sigma(female))*(1.d0-lambda_n(male))
  ! Female fired, finds new, male finds
  aux_NW(3) = aux_NW(2) + sigma(female)*lambda_u(female)*lambda_n(male)
  ! Female fired, finds new, male no finds
  aux_NW(4) = aux_NW(3) + sigma(female)*lambda_u(female)*(1.d0-lambda_n(male))
  ! Female fired, no new, male finds
  aux_NW(5) = aux_NW(4) + sigma(female)*(1.d0-lambda_u(female))*lambda_n(male)
  ! Female fired, no new, male no finds
  aux_NW(6) = aux_NW(5) + sigma(female)*(1.d0-lambda_u(female))*(1.d0-lambda_n(male))

  ! Probabilities for NU household
  aux_NU = 0.d0
  ! Male finds, female finds
  aux_NU(1) = lambda_n(male)*lambda_u(female)
  ! Male finds, female no finds
  aux_NU(2) = aux_NU(1) + lambda_n(male)*(1.d0-lambda_u(female))
  ! Male no finds, female no finds
  aux_NU(3) = aux_NU(2) + (1.d0-lambda_n(male))*lambda_u(female)
  ! Male no finds, female no finds
  aux_NU(4) = aux_NU(3) + (1.d0-lambda_n(male))*(1.d0-lambda_u(female))

  ! Probabilities for NN household
  aux_NN = 0.d0
  ! Male finds, female finds
  aux_NN(1) = lambda_n(male)*lambda_n(female)
  ! Male finds, female no finds
  aux_NN(2) = aux_NN(1) + lambda_n(male)*(1.d0-lambda_n(female))
  ! Male no finds, female no finds
  aux_NN(3) = aux_NN(2) + (1.d0-lambda_n(male))*lambda_n(female)
  ! Male no finds, female no finds
  aux_NN(4) = aux_NN(3) + (1.d0-lambda_n(male))*(1.d0-lambda_n(female))

  ! Create a bigger grid of assets for the simulation
  sim_a_values = loggrid(min_a, max_a, sim_gp_a)

  ! Reshape decision rules with new grid
  do ind_z = 1, gp_z2
    sim_NN_pf(:,ind_z) = regrid(gp_a,a_values,NN_pf(:,ind_z),sim_gp_a,sim_a_values)
    sim_NW_pf(:,ind_z) = regrid(gp_a,a_values,NW_pf(:,ind_z),sim_gp_a,sim_a_values)
    sim_WN_pf(:,ind_z) = regrid(gp_a,a_values,WN_pf(:,ind_z),sim_gp_a,sim_a_values)
    sim_WW_pf(:,ind_z) = regrid(gp_a,a_values,WW_pf(:,ind_z),sim_gp_a,sim_a_values)

    do ind_g = 1, gp_gamma
    do ind_b = 0, 1
      sim_UN_pf(:,ind_z,ind_g,ind_b) = regrid(gp_a,a_values,UN_pf(:,ind_z,ind_g,ind_b),&
                                              sim_gp_a,sim_a_values)
      sim_NU_pf(:,ind_z,ind_g,ind_b) = regrid(gp_a,a_values,NU_pf(:,ind_z,ind_g,ind_b),&
                                              sim_gp_a,sim_a_values)
      sim_UW_pf(:,ind_z,ind_g,ind_b) = regrid(gp_a,a_values,UW_pf(:,ind_z,ind_g,ind_b),&
                                              sim_gp_a,sim_a_values)
      sim_WU_pf(:,ind_z,ind_g,ind_b) = regrid(gp_a,a_values,WU_pf(:,ind_z,ind_g,ind_b),&
                                              sim_gp_a,sim_a_values)
      do ind_g_f = 1, gp_gamma
      do ind_b_f = 0, 1
        sim_UU_pf(:,ind_z,ind_g,ind_g_f,ind_b,ind_b_f) = regrid(gp_a,a_values,&
            UU_pf(:,ind_z,ind_g,ind_g_f,ind_b,ind_b_f),sim_gp_a,sim_a_values)
      end do
      end do
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

  ! Initialise myshocks
  myshock_z = -1
  myshock_g = -1

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
  itrans = 0.d0
  jtrans = 0.d0
  reps = 0
  top_assets = 0.d0

  ! Open CSV file where to store simulation
  if (gen_output) then
    open(unit=14, file='simulation_3.csv',status='unknown')
    call csv_write(14,"ID",.false.)
    call csv_write(14,"Period",.false.)
    call csv_write(14,"Status Male Today",.false.)
    call csv_write(14,"Status Male Tomorrow",.false.)
    call csv_write(14,"Status Female Today",.false.)
    call csv_write(14,"Status Female Tomorrow",.false.)
    call csv_write(14,"Labour Income",.false.)
    call csv_write(14,"Taxes Paid",.false.)
    call csv_write(14,"Benefits Received",.false.)
    call csv_write(14,"Wealth",.true.)
  end if

  do ind_p = 1, periods-1 ! CAREFUL with last period
  do ind_ag = 1, agents
    ! Initialise aux variables to 0
    income = 0.d0
    labincome = 0.d0
    taxrev = 0.d0
    bpaid = 0.d0

    ! Assign shocks
    allocate(aux_vec(gp_gamma))
    if (ind_p.eq.1) then
      myshock_z(ind_ag,ind_p) = realise(shock_z(ind_ag,ind_p),z_ssdist,gp_z)
      do mysex = 1, 2
        aux_vec = gamma_trans(mysex,:)
        myshock_g(mysex,ind_ag,ind_p) = realise(shock_g(mysex,ind_ag,ind_p),aux_vec, gp_gamma)
      end do
    end if
    do mysex = 1, 2
      aux_vec = gamma_trans(mysex,:)
      myshock_g(mysex,ind_ag,ind_p+1) = realise(shock_g(mysex,ind_ag,ind_p+1),aux_vec, gp_gamma)
    end do
    deallocate(aux_vec)
    allocate(aux_vec(gp_z2))
    aux_vec = z_trans(myshock_z(ind_ag,ind_p),:)
    myshock_z(ind_ag,ind_p+1) = realise(shock_z(ind_ag,ind_p+1),aux_vec,gp_z2)
    deallocate(aux_vec)

    ! Quantify shocks
    a = sim_a_values(assets(ind_ag))
    a_income = (1.d0+int_rate)*a + T ! CAREFUL: Lump-sum transfer already added here
    do mysex = 1, 2
      z(mysex) = z_values(mysex,myshock_z(ind_ag,ind_p))
    end do

    ! WW ------------------------------------------------------------------------------------------
    if ((LMstatus(male,ind_ag).eq.1).and.(LMstatus(female,ind_ag).eq.1)) then
      ! Benefits - none
      ! Labor income
      labincome = wage*(z(male)+z(female))
      ! Tax revenue
      taxrev = tau*labincome

      ! Assets next period
      new_assets(ind_ag) = sim_WW_pf(assets(ind_ag),myshock_z(ind_ag,ind_p))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*wage*(z(male)+z(female)) - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.aux_WW(1)) then
        ! Male keeps, female keeps
        new_entitled(male,ind_ag) = 0
        new_entitled(female,ind_ag) = 0
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WW(1)).and.(shock_lm(ind_ag,ind_p).le.aux_WW(2))) then
        ! Male keeps, female fired, finds new
        new_entitled(male,ind_ag) = 0
        new_entitled(female,ind_ag) = 1
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WW(2)).and.(shock_lm(ind_ag,ind_p).le.aux_WW(3))) then
        ! Male keeps, female fired, no new
        new_entitled(male,ind_ag) = 0
        new_entitled(female,ind_ag) = 1
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WW(3)).and.(shock_lm(ind_ag,ind_p).le.aux_WW(4))) then
        ! Male fired, find new, female keeps
        new_entitled(male,ind_ag) = 1
        new_entitled(female,ind_ag) = 0
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WW(4)).and.(shock_lm(ind_ag,ind_p).le.aux_WW(5))) then
        ! Male fired, finds new, female fired, finds new
        new_entitled(male,ind_ag) = 1
        new_entitled(female,ind_ag) = 1
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WW(5)).and.(shock_lm(ind_ag,ind_p).le.aux_WW(6))) then
        ! Male fired, finds new, female fired, no new
        new_entitled(male,ind_ag) = 1
        new_entitled(female,ind_ag) = 1
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WW(6)).and.(shock_lm(ind_ag,ind_p).le.aux_WW(7))) then
        ! Male fired, no new, female keeps
        new_entitled(male,ind_ag) = 1
        new_entitled(female,ind_ag) = 0
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WW(7)).and.(shock_lm(ind_ag,ind_p).le.aux_WW(8))) then
        ! Male fired, no new, female fired, finds new
        new_entitled(male,ind_ag) = 1
        new_entitled(female,ind_ag) = 1
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      else
        ! Male fired, no new, female fired, no new
        new_entitled(male,ind_ag) = 1
        new_entitled(female,ind_ag) = 1
        call choice_JJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      end if

    ! WU ------------------------------------------------------------------------------------------
    elseif ((LMstatus(male,ind_ag).eq.1).and.(LMstatus(female,ind_ag).eq.2)) then
      ! Benefits - check female
      bpaid = real(entitled(female,ind_ag))*(benefits(z(female)))
      if ((entitled(female,ind_ag).eq.1).and.(shock_mu(female,ind_ag,ind_p).gt.mu)) then
        ! Agent is entitled and does NOT get hit by a mu shock
        new_entitled(female,ind_ag) = 1
      else
        new_entitled(female,ind_ag) = 0
      end if
      ! Labor income
      labincome = wage*z(male)
      ! Tax revenue
      taxrev = tau*labincome + tau*bpaid

      ! Assets next period
      new_assets(ind_ag) = sim_WU_pf(assets(ind_ag),myshock_z(ind_ag,ind_p),&
                            myshock_g(female,ind_ag,ind_p),entitled(female,ind_ag))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*wage*z(male) + (1.d0-tau)*bpaid - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.aux_WU(1)) then
        ! Male keeps, female finds
        new_entitled(male,ind_ag) = 0
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WU(1)).and.(shock_lm(ind_ag,ind_p).le.aux_WU(2))) then
        ! Male keeps, female no finds
        new_entitled(male,ind_ag) = 0
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WU(2)).and.(shock_lm(ind_ag,ind_p).le.aux_WU(3))) then
        ! Male fired, finds new, female finds
        new_entitled(male,ind_ag) = 1
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WU(3)).and.(shock_lm(ind_ag,ind_p).le.aux_WU(4))) then
        ! Male fired, finds new, female no finds
        new_entitled(male,ind_ag) = 1
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WU(4)).and.(shock_lm(ind_ag,ind_p).le.aux_WU(5))) then
        ! Male fired, no new, female finds
        new_entitled(male,ind_ag) = 1
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      else
        ! Male fired, no new, female no finds
        new_entitled(male,ind_ag) = 1
        call choice_JJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      end if

    ! WN ------------------------------------------------------------------------------------------
    elseif ((LMstatus(male,ind_ag).eq.1).and.(LMstatus(female,ind_ag).eq.3)) then
      ! Benefits - none for female
      new_entitled(female,ind_ag) = 0
      ! Labor income
      labincome = wage*z(male)
      ! Tax revenue
      taxrev = tau*labincome

      ! Assets next period
      new_assets(ind_ag) = sim_WN_pf(assets(ind_ag),myshock_z(ind_ag,ind_p))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*wage*z(male) - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.aux_WN(1)) then
        ! Male keeps, female finds
        new_entitled(male,ind_ag) = 0
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WN(1)).and.(shock_lm(ind_ag,ind_p).le.aux_WN(2))) then
        ! Male keeps, female no finds
        new_entitled(male,ind_ag) = 0
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WN(2)).and.(shock_lm(ind_ag,ind_p).le.aux_WN(3))) then
        ! Male fired, finds new, female finds
        new_entitled(male,ind_ag) = 1
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WN(3)).and.(shock_lm(ind_ag,ind_p).le.aux_WN(4))) then
        ! Male fired, finds new, female no finds
        new_entitled(male,ind_ag) = 1
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_WN(4)).and.(shock_lm(ind_ag,ind_p).le.aux_WN(5))) then
        ! Male fired, no new, female finds
        new_entitled(male,ind_ag) = 1
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      else
        ! Male fired, no new, female no finds
        new_entitled(male,ind_ag) = 1
        call choice_JJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      end if

    ! UW ------------------------------------------------------------------------------------------
    elseif ((LMstatus(male,ind_ag).eq.2).and.(LMstatus(female,ind_ag).eq.1)) then
      ! Benefits - check male
      bpaid = real(entitled(male,ind_ag))*(benefits(z(male)))
      if ((entitled(male,ind_ag).eq.1).and.(shock_mu(male,ind_ag,ind_p).gt.mu)) then
        ! Agent is entitled and does NOT get hit by a mu shock
        new_entitled(male,ind_ag) = 1
      else
        new_entitled(male,ind_ag) = 0
      end if
      ! Labor income
      labincome = wage*z(female)
      ! Tax revenue
      taxrev = tau*labincome + tau*bpaid

      ! Assets next period
      new_assets(ind_ag) = sim_UW_pf(assets(ind_ag),myshock_z(ind_ag,ind_p),&
                            myshock_g(male,ind_ag,ind_p),entitled(male,ind_ag))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*wage*z(female) + (1.d0-tau)*bpaid - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.aux_UW(1)) then
        ! Female keeps, male finds
        new_entitled(female,ind_ag) = 0
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_UW(1)).and.(shock_lm(ind_ag,ind_p).le.aux_UW(2))) then
        ! Female keeps, male no finds
        new_entitled(female,ind_ag) = 0
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_UW(2)).and.(shock_lm(ind_ag,ind_p).le.aux_UW(3))) then
        ! Female fired, finds new, male finds
        new_entitled(female,ind_ag) = 1
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_UW(3)).and.(shock_lm(ind_ag,ind_p).le.aux_UW(4))) then
        ! Female fired, finds new, male no finds
        new_entitled(female,ind_ag) = 1
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_UW(4)).and.(shock_lm(ind_ag,ind_p).le.aux_UW(5))) then
        ! Female fired, no new, male finds
        new_entitled(female,ind_ag) = 1
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      else
        ! Female fired, no new, male no finds
        new_entitled(female,ind_ag) = 1
        call choice_JJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      end if

    ! UU ------------------------------------------------------------------------------------------
    elseif ((LMstatus(male,ind_ag).eq.2).and.(LMstatus(female,ind_ag).eq.2)) then
      ! Benefits - check both
      bpaid = (real(entitled(male,ind_ag))*(benefits(z(male)))) + &
              (real(entitled(female,ind_ag))*(benefits(z(female))))
      do mysex = 1, 2
        if ((entitled(mysex,ind_ag).eq.1).and.(shock_mu(mysex,ind_ag,ind_p).gt.mu)) then
          ! Agent is entitled and does NOT get hit by a mu shock
          new_entitled(mysex,ind_ag) = 1
        else
          new_entitled(mysex,ind_ag) = 0
        end if
      end do
      ! Labor income - None
      ! Tax revenue
      taxrev = tau*bpaid

      ! Assets next period
      new_assets(ind_ag) = sim_UU_pf(assets(ind_ag),myshock_z(ind_ag,ind_p),&
                            myshock_g(male,ind_ag,ind_p),myshock_g(female,ind_ag,ind_p),&
                            entitled(male,ind_ag),entitled(female,ind_ag))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*bpaid - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.aux_UU(1)) then
        ! Male finds, female finds
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_UU(1)).and.(shock_lm(ind_ag,ind_p).le.aux_UU(2))) then
        ! Male finds, female no finds
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_UU(2)).and.(shock_lm(ind_ag,ind_p).le.aux_UU(3))) then
        ! Male no finds, female no finds
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      else
        ! Male no finds, female no finds
        call choice_JJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      end if

    ! UN ------------------------------------------------------------------------------------------
    elseif ((LMstatus(male,ind_ag).eq.2).and.(LMstatus(female,ind_ag).eq.3)) then
      ! Benefits - None for female, check male
      new_entitled(female,ind_ag) = 0
      bpaid = real(entitled(male,ind_ag))*(benefits(z(male)))
      if ((entitled(male,ind_ag).eq.1).and.(shock_mu(male,ind_ag,ind_p).gt.mu)) then
        ! Agent is entitled and does NOT get hit by a mu shock
        new_entitled(male,ind_ag) = 1
      else
        new_entitled(male,ind_ag) = 0
      end if
      ! Labor income - None
      ! Tax revenue
      taxrev = tau*bpaid

      ! Assets next period
      new_assets(ind_ag) = sim_UN_pf(assets(ind_ag),myshock_z(ind_ag,ind_p),&
                            myshock_g(male,ind_ag,ind_p),entitled(male,ind_ag))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*bpaid - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.aux_UN(1)) then
        ! Male finds, female finds
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_UN(1)).and.(shock_lm(ind_ag,ind_p).le.aux_UN(2))) then
        ! Male finds, female no finds
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_UN(2)).and.(shock_lm(ind_ag,ind_p).le.aux_UN(3))) then
        ! Male no finds, female no finds
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      else
        ! Male no finds, female no finds
        call choice_JJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      end if

    ! NW ------------------------------------------------------------------------------------------
    elseif ((LMstatus(male,ind_ag).eq.3).and.(LMstatus(female,ind_ag).eq.1)) then
      ! Benefits - none for male
      new_entitled(male,ind_ag) = 0
      ! Labor income
      labincome = wage*z(female)
      ! Tax revenue
      taxrev = tau*bpaid

      ! Assets next period
      new_assets(ind_ag) = sim_NW_pf(assets(ind_ag),myshock_z(ind_ag,ind_p))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*wage*z(female) - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.aux_NW(1)) then
        ! Female keeps, male finds
        new_entitled(female,ind_ag) = 0
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_NW(1)).and.(shock_lm(ind_ag,ind_p).le.aux_NW(2))) then
        ! Female keeps, male no finds
        new_entitled(female,ind_ag) = 0
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_NW(2)).and.(shock_lm(ind_ag,ind_p).le.aux_NW(3))) then
        ! Female fired, finds new, male finds
        new_entitled(female,ind_ag) = 1
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_NW(3)).and.(shock_lm(ind_ag,ind_p).le.aux_NW(4))) then
        ! Female fired, finds new, male no finds
        new_entitled(female,ind_ag) = 1
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_NW(4)).and.(shock_lm(ind_ag,ind_p).le.aux_NW(5))) then
        ! Female fired, no new, male finds
        new_entitled(female,ind_ag) = 1
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      else
        ! Female fired, no new, male no finds
        new_entitled(female,ind_ag) = 1
        call choice_JJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      end if

    ! NU ------------------------------------------------------------------------------------------
    elseif ((LMstatus(male,ind_ag).eq.3).and.(LMstatus(female,ind_ag).eq.2)) then
      ! Benefits - None for male, check female
      new_entitled(male,ind_ag) = 0
      bpaid = real(entitled(female,ind_ag))*(benefits(z(female)))
      if ((entitled(female,ind_ag).eq.1).and.(shock_mu(female,ind_ag,ind_p).gt.mu)) then
        ! Agent is entitled and does NOT get hit by a mu shock
        new_entitled(female,ind_ag) = 1
      else
        new_entitled(female,ind_ag) = 0
      end if
      ! Labor income - None
      ! Tax revenue
      taxrev = tau*bpaid

      ! Assets next period
      new_assets(ind_ag) = sim_NU_pf(assets(ind_ag),myshock_z(ind_ag,ind_p),&
                            myshock_g(female,ind_ag,ind_p),entitled(female,ind_ag))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income + (1.d0-tau)*bpaid - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.aux_NU(1)) then
        ! Male finds, female finds
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_NU(1)).and.(shock_lm(ind_ag,ind_p).le.aux_NU(2))) then
        ! Male finds, female no finds
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_NU(2)).and.(shock_lm(ind_ag,ind_p).le.aux_NU(3))) then
        ! Male no finds, female no finds
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      else
        ! Male no finds, female no finds
        call choice_JJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      end if

    ! NN ------------------------------------------------------------------------------------------
    elseif ((LMstatus(male,ind_ag).eq.3).and.(LMstatus(female,ind_ag).eq.3)) then
      ! Benefits - None
      new_entitled(male,ind_ag) = 0
      new_entitled(female,ind_ag) = 0
      ! Labor income - None
      ! Tax revenue - None

      ! Assets next period
      new_assets(ind_ag) = sim_NN_pf(assets(ind_ag),myshock_z(ind_ag,ind_p))
      ap = sim_a_values(new_assets(ind_ag))

      ! Income/consumption
      income = a_income - ap

      ! LM status next period
      if (shock_lm(ind_ag,ind_p).le.aux_NN(1)) then
        ! Male finds, female finds
        call choice_VV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_NN(1)).and.(shock_lm(ind_ag,ind_p).le.aux_NN(2))) then
        ! Male finds, female no finds
        call choice_VJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      elseif ((shock_lm(ind_ag,ind_p).gt.aux_NN(2)).and.(shock_lm(ind_ag,ind_p).le.aux_NN(3))) then
        ! Male no finds, female no finds
        call choice_JV(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      else
        ! Male no finds, female no finds
        call choice_JJ(ap,myshock_z(ind_ag,ind_p+1),&
                      myshock_g(male,ind_ag,ind_p+1),myshock_g(female,ind_ag,ind_p+1),&
                      new_entitled(male,ind_ag),new_entitled(female,ind_ag), &
                      new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag))
      end if
    else
      print *, "Simulation: LMstatus not in range"
    end if ! LMstatus

    ! Add to totals
    if (ind_p.ge.(periods/2)) then
      tot_income = tot_income + income
      tot_labincome = tot_labincome + labincome
      tot_taxrev = tot_taxrev + taxrev
      tot_bpaid = tot_bpaid + bpaid
      tot_assets = tot_assets + a
      jtrans(LMstatus(male,ind_ag),LMstatus(female,ind_ag),&
            new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag)) = &
            jtrans(LMstatus(male,ind_ag),LMstatus(female,ind_ag),&
                  new_LMstatus(male,ind_ag),new_LMstatus(female,ind_ag)) + 1.d0
      ! Add to totals relevant for employed
      do mysex = 1, 2
        itrans(mysex,LMstatus(mysex,ind_ag),new_LMstatus(mysex,ind_ag)) = &
          itrans(mysex,LMstatus(mysex,ind_ag),new_LMstatus(mysex,ind_ag)) + 1.d0
        if (LMstatus(mysex,ind_ag).eq.1) then
          ! Register labor market status
          employed(mysex) = employed(mysex) + 1.d0
          tot_z = tot_z + z(mysex)
        elseif(LMstatus(mysex,ind_ag).eq.2) then
          ! Register labor market status
          unemployed(mysex) = unemployed(mysex) + 1.d0
        elseif (LMstatus(mysex,ind_ag).eq.3) then
          ! Register labor market status
          OLF(mysex) = OLF(mysex) + 1.d0
        end if
      end do
    end if

    ! Agents using top assets
    if (assets(ind_ag).eq.sim_gp_a) then
      top_assets = top_assets + 1.d0
    end if

    ! Print to CSV
    if (gen_output) then
      if (ind_p.ge.(periods-printed)) then
        call csv_write(14,real(ind_ag),.false.)
        call csv_write(14,real(ind_p),.false.)
        call csv_write(14,real(LMstatus(male,ind_ag)),.false.)
        call csv_write(14,real(new_LMstatus(male,ind_ag)),.false.)
        call csv_write(14,real(LMstatus(female,ind_ag)),.false.)
        call csv_write(14,real(new_LMstatus(female,ind_ag)),.false.)
        call csv_write(14,labincome,.false.)
        call csv_write(14,taxrev,.false.)
        call csv_write(14,bpaid,.false.)
        call csv_write(14,a,.true.)
      end if
    end if
  end do ! Agents
    ! Update variables
    assets = new_assets
    LMstatus = new_LMstatus
    entitled = new_entitled

    if (ind_p.ge.(periods/2)) then
      reps = reps + 1
    end if
  end do ! Periods

  ! Close csv file
  if (gen_output) then
    close(14)
  end if

  ! Compute averages
  tot_income = tot_income/real(reps)
  tot_labincome = tot_labincome/real(reps)
  tot_taxrev = tot_taxrev/real(reps)
  tot_bpaid = tot_bpaid/real(reps)
  tot_assets = tot_assets/real(reps)
  tot_z = tot_z/real(reps)
  do mysex = 1, 2
    employed(mysex) = employed(mysex)/real(reps)
    unemployed(mysex) = unemployed(mysex)/real(reps)
    OLF(mysex) = OLF(mysex)/real(reps)
  end do

  ! Compute results
  wealth(male+female) = tot_assets
  aux_KLratio(male+female) = tot_assets/tot_z
  aux_average_z(male+female) = tot_z/sum(employed)
  aux_T(male+female) = (tot_taxrev-tot_bpaid)/real(agents)
  do mysex = 1, 2
    Erate(married,mysex) = employed(mysex)/real(agents)
    Urate(married,mysex) = unemployed(mysex)/(employed(mysex)+unemployed(mysex))
    Nrate(married,mysex) = OLF(mysex)/real(agents)
    transitions(married,mysex,1, :) = itrans(mysex,1,:)/sum(itrans(mysex,1,:))
    transitions(married,mysex,2, :) = itrans(mysex,2,:)/sum(itrans(mysex,2,:))
    transitions(married,mysex,3, :) = itrans(mysex,3,:)/sum(itrans(mysex,3,:))
  end do

  ! Print share of agents using top assets
  print '(a,f7.4)', " Share of married HH using top assets:", top_assets/(real(reps*agents))
  print *, ""

CONTAINS

  !----- CHOICE VV --------------------------------------------------------------------------------
  subroutine choice_VV(val_a,my_z,my_gm,my_gf,my_bm,my_bf,choice_m,choice_f)
    implicit none
    integer, intent(in) :: my_z, my_gm, my_gf, my_bm, my_bf
    integer :: aux_max
    integer, intent(out) :: choice_m, choice_f
    real(8), intent(in) :: val_a
    real(8), dimension(1:9) :: aux_choice

    ! Interpolate value of WW
    aux_choice(1) = my_inter(a_values,WW_vf(:,my_z),gp_a,val_a)

    ! Interpolate value of WU
    aux_choice(2) = my_inter(a_values,WU_vf(:,my_z,my_gf,my_bf),gp_a,val_a)

    ! Interpolate value of WN
    aux_choice(3) = my_inter(a_values,WN_vf(:,my_z),gp_a,val_a)

    ! Interpolate value of UW
    aux_choice(4) = my_inter(a_values,UW_vf(:,my_z,my_gm,my_bm),gp_a,val_a)

    ! Interpolate value of UU
    aux_choice(5) = my_inter(a_values,UU_vf(:,my_z,my_gm,my_gf,my_bm,my_bf),gp_a,val_a)

    ! Interpolate value of UN
    aux_choice(6) = my_inter(a_values,UN_vf(:,my_z,my_gm,my_bm),gp_a,val_a)

    ! Interpolate value of NW
    aux_choice(7) = my_inter(a_values,NW_vf(:,my_z),gp_a,val_a)

    ! Interpolate value of NU
    aux_choice(8) = my_inter(a_values,NU_vf(:,my_z,my_gf,my_bf),gp_a,val_a)

    ! Interpolate value of NN
    aux_choice(9) = my_inter(a_values,NN_vf(:,my_z),gp_a,val_a)

    aux_max = maxloc(aux_choice, dim=1)

    choice_m = mystates(aux_max,male)
    choice_f = mystates(aux_max,female)
  end subroutine choice_VV

  !----- CHOICE VJ --------------------------------------------------------------------------------
  subroutine choice_VJ(val_a,my_z,my_gm,my_gf,my_bm,my_bf,choice_m,choice_f)
    implicit none
    integer, intent(in) :: my_z, my_gm, my_gf, my_bm, my_bf
    integer :: aux_max
    integer, intent(out) :: choice_m, choice_f
    real(8), intent(in) :: val_a
    real(8), dimension(1:9) :: aux_choice
    logical, dimension(1:9) :: aux_possible

    ! Initialise all values of mask to true
    aux_possible = .true.

    ! Interpolate value of WW
    aux_possible(1) = .false.

    ! Interpolate value of WU
    aux_choice(2) = my_inter(a_values,WU_vf(:,my_z,my_gf,my_bf),gp_a,val_a)

    ! Interpolate value of WN
    aux_choice(3) = my_inter(a_values,WN_vf(:,my_z),gp_a,val_a)

    ! Interpolate value of UW
    aux_possible(4) = .false.

    ! Interpolate value of UU
    aux_choice(5) = my_inter(a_values,UU_vf(:,my_z,my_gm,my_gf,my_bm,my_bf),gp_a,val_a)

    ! Interpolate value of UN
    aux_choice(6) = my_inter(a_values,UN_vf(:,my_z,my_gm,my_bm),gp_a,val_a)

    ! Interpolate value of NW
    aux_possible(7) = .false.

    ! Interpolate value of NU
    aux_choice(8) = my_inter(a_values,NU_vf(:,my_z,my_gf,my_bf),gp_a,val_a)

    ! Interpolate value of NN
    aux_choice(9) = my_inter(a_values,NN_vf(:,my_z),gp_a,val_a)

    aux_max = maxloc(aux_choice, dim=1, mask = aux_possible)

    choice_m = mystates(aux_max,male)
    choice_f = mystates(aux_max,female)
  end subroutine choice_VJ

  !----- CHOICE JV --------------------------------------------------------------------------------
  subroutine choice_JV(val_a,my_z,my_gm,my_gf,my_bm,my_bf,choice_m,choice_f)
    implicit none
    integer, intent(in) :: my_z, my_gm, my_gf, my_bm, my_bf
    integer :: aux_max
    integer, intent(out) :: choice_m, choice_f
    real(8), intent(in) :: val_a
    real(8), dimension(1:9) :: aux_choice
    logical, dimension(1:9) :: aux_possible

    ! Initialise all values of mask to true
    aux_possible = .true.

    ! Interpolate value of WW
    aux_possible(1) = .false.

    ! Interpolate value of WU
    aux_possible(2) = .false.

    ! Interpolate value of WN
    aux_possible(3) = .false.

    ! Interpolate value of UW
    aux_choice(4) = my_inter(a_values,UW_vf(:,my_z,my_gm,my_bm),gp_a,val_a)

    ! Interpolate value of UU
    aux_choice(5) = my_inter(a_values,UU_vf(:,my_z,my_gm,my_gf,my_bm,my_bf),gp_a,val_a)

    ! Interpolate value of UN
    aux_choice(6) = my_inter(a_values,UN_vf(:,my_z,my_gm,my_bm),gp_a,val_a)

    ! Interpolate value of NW
    aux_choice(7) = my_inter(a_values,NW_vf(:,my_z),gp_a,val_a)

    ! Interpolate value of NU
    aux_choice(8) = my_inter(a_values,NU_vf(:,my_z,my_gf,my_bf),gp_a,val_a)

    ! Interpolate value of NN
    aux_choice(9) = my_inter(a_values,NN_vf(:,my_z),gp_a,val_a)

    aux_max = maxloc(aux_choice, dim=1, mask = aux_possible)

    choice_m = mystates(aux_max,male)
    choice_f = mystates(aux_max,female)
  end subroutine choice_JV

  !----- CHOICE JV --------------------------------------------------------------------------------
  subroutine choice_JJ(val_a,my_z,my_gm,my_gf,my_bm,my_bf,choice_m,choice_f)
    implicit none
    integer, intent(in) :: my_z, my_gm, my_gf, my_bm, my_bf
    integer :: aux_max
    integer, intent(out) :: choice_m, choice_f
    real(8), intent(in) :: val_a
    real(8), dimension(1:9) :: aux_choice
    logical, dimension(1:9) :: aux_possible

    ! Initialise all values of mask to true
    aux_possible = .true.

    ! Interpolate value of WW
    aux_possible(1) = .false.

    ! Interpolate value of WU
    aux_possible(2) = .false.

    ! Interpolate value of WN
    aux_possible(3) = .false.

    ! Interpolate value of UW
    aux_possible(4) = .false.

    ! Interpolate value of UU
    aux_choice(5) = my_inter(a_values,UU_vf(:,my_z,my_gm,my_gf,my_bm,my_bf),gp_a,val_a)

    ! Interpolate value of UN
    aux_choice(6) = my_inter(a_values,UN_vf(:,my_z,my_gm,my_bm),gp_a,val_a)

    ! Interpolate value of NW
    aux_possible(7) = .false.

    ! Interpolate value of NU
    aux_choice(8) = my_inter(a_values,NU_vf(:,my_z,my_gf,my_bf),gp_a,val_a)

    ! Interpolate value of NN
    aux_choice(9) = my_inter(a_values,NN_vf(:,my_z),gp_a,val_a)

    aux_max = maxloc(aux_choice, dim=1, mask = aux_possible)

    choice_m = mystates(aux_max,male)
    choice_f = mystates(aux_max,female)
  end subroutine choice_JJ
end subroutine SimMarried
