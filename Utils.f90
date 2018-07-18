module Utils
  implicit none

CONTAINS
  !===== TIME IN HOURS, MINUTES, AND SECONDS ======================================================
  subroutine mytime(total)
    integer :: hours, minutes
    real(8), intent(in) :: total
    real(8) :: seconds

    hours = int(total / 3600.d0)
    minutes = int((total - 3600.d0*real(hours)) / 60.d0)
    seconds = total - 3600.d0*hours - 60.d0*minutes

    if (hours.ne.0) then
      print '(a4,i2,a7)',"  + ", hours, " hours."
    end if
    if (minutes.ne.0) then
      print '(a4,i2,a9)',"  + ", minutes, " minutes."
    end if
    if (seconds.gt.0.d0) then
      print '(a4,f9.6,a9)',"  + ", seconds, " seconds."
    end if
    print *, ""
  end subroutine mytime

  !===== LINSPACE =================================================================================
  function linspace(my_start, my_stop, n)
    implicit none
    integer :: n, i
    real(8) :: my_start, my_stop, step
    real(8), dimension(n) :: grid, linspace

    if (n.eq.1) then
      grid(n) = (my_start+my_stop)/2.d0
    elseif (n.ge.2) then
      grid(1) = my_start
      if (n.gt.2) then
        step = (my_stop-my_start)/(real(n-1))
        do i = 2, n-1
          grid(i) = grid(i-1) + step
        end do
      end if
      grid(n) = my_stop
    endif
    linspace = grid
  end function linspace

  !===== LOGGRID ==================================================================================
  function loggrid(my_start, my_stop, n)
    implicit none
    integer :: n, i
    real(8) :: my_start, my_stop, step
    real(8), dimension(n) :: grid, loggrid

    if (n.eq.1) then
      grid(n) = (my_start+my_stop)/2.d0
    elseif (n.ge.2) then
      grid(1) = my_start
      if (n.gt.2) then
        step=(log(my_stop+2.d0)-log(my_start+2.d0))/real(n-1)
        do i = 2, n-1
          grid(i) = exp(log(grid(i-1)+2.d0)+step)-2.d0
        end do
      end if
      grid(n) = my_stop
    endif
    loggrid = grid
  end function loggrid

  !===== TAUCHEN AR(1) DISCRETISATION =============================================================
  subroutine tauchen(rho, sigma, cover, gp, values, trans)
    ! the simple case of approximating first-order autoregressive process with Markov chain
    !
    ! y_t = rho * y_(t-1) + u_t
    !
    ! u_t is a Gaussian white noise process with standard deviation sigma.
    !
    ! cover determines the width of discretized state space, Tauchen uses m=3
    !
    ! gp is the number of possible states chosen to approximate
    ! the y_t process
    !
    ! trans is the transition matrix of the Markov chain
    !
    ! values is the discretized state space of y_t
    !
    ! Adapted from https://github.com/lucaguerrieri/
    implicit none
    integer, intent(in) :: gp
    integer :: j, k
    real(8), intent(in) :: rho, sigma, cover
    real(8) :: sd_y, ymin, ymax, w
    real(8), dimension(gp), intent(out) :: values
    real(8), dimension(gp,gp), intent(out) :: trans

    ! standard deviation of y_t
    sd_y = sqrt(sigma**2.d0/(1.d0-rho**2.d0))

    ymax = cover*sd_y   ! upper boundary of state space
    ymin = -ymax        ! lower boundary of state space
    w = (ymax-ymin)/real(gp-1) ! length of interval

    values = linspace(ymin, ymax, gp)

    ! Compute transition matrix
    do j = 1, gp
        do k = 2, gp-1
            trans(j,k) = &
            normcdf(values(k)-rho*values(j)+(values(k+1)-values(k))/2.d0,0.d0,sigma) &
            - normcdf(values(k)-rho*values(j)-(values(k)-values(k-1))/2.d0,0.d0,sigma)
        end do
        ! only subtract half the interval on the right
        trans(j,1) = normcdf(values(1)-rho*values(j)+w/2.d0,0.d0,sigma)
        ! only subtract half the interval on the left
        trans(j,gp) = 1.d0 - normcdf(values(gp)-rho*values(j)-w/2.d0,0.d0,sigma)
    end do
  CONTAINS
    real(8) function normcdf(x,mu,sigma)
      implicit none
      real(8), intent(in) :: x, mu, sigma

      normcdf = (1.d0+erf((x-mu)/sqrt(2.d0*sigma**2)))/2.d0
    end function normcdf
  end subroutine tauchen

  !===== GOLDEN SEARCH SECTION ====================================================================
  subroutine golden_method(f, a, b, x1, f1, mytol, mymaxit)
  ! Applies Golden-section search to search for the _maximum_ of a function in the interval (a, b)
  !
  ! https://en.wikipedia.org/wiki/Golden-section_search
  ! Adapted to Fortran90 from: https://github.com/QuantEcon
    integer, optional :: mymaxit
    integer :: maxit, it
    real(8), external :: f
    real(8), intent(in) :: a, b
    real(8), intent(out) :: x1, f1
    real(8), optional :: mytol
    real(8) :: tol, alpha1, alpha2, d, f2, x2, s

    ! Assign default value to maxit if not defined by user
    if (present(mymaxit)) then
        maxit = mymaxit
    else
        maxit = 1000
    end if

    ! Assign default value to tol if not defined by user
    if (present(mytol)) then
        tol = mytol
    else
        tol = 1.0d-6
    end if

    alpha1 = (3.d0 - sqrt(5.d0)) / 2.d0
    alpha2 = 1.d0 - alpha1
    d = b - a
    x1 = a + alpha1*d
    x2 = a + alpha2*d
    s = 1.d0
    f1 = f(x1)
    f2 = f(x2)
    d = alpha1*alpha2*d

    it = 0

    do while ((d.gt.tol).and.(it.lt.maxit))
        it = it + 1
        d = d*alpha2

        if (f2.gt.f1) then
            x1 = x2
            f1 = f2
            x2 = x1 + s*d
        else
            x2 = x1 - s*d
        end if

        s = sign(s, x2-x1)
        f2 = f(x2)
    end do

    if (it.ge.maxit) then
        print *, "Golden method: Maximum iterations exceeded"
    end if

    if (f2.gt.f1) then
        x1 = x2
        f1 = f2
    end if
  end subroutine golden_method

  !===== FINDS CLOSEST TWO NUMBERS IN A VECTOR ====================================================
  subroutine my_smin(svector,g_points,svalue,low_bnd,upp_bnd)
    implicit none
    integer, intent(in) :: g_points
    real(8), intent(in) :: svalue
    real(8), dimension(g_points), intent(in) :: svector
    integer, intent(out) :: low_bnd, upp_bnd
    integer :: ind_x

    low_bnd = g_points

    do ind_x = 2, g_points
        if (svalue.le.svector(ind_x)) then
            low_bnd = ind_x - 1
            exit
        endif
    enddo

    if (low_bnd.eq.g_points) then
        low_bnd = g_points - 1
    endif

    upp_bnd = low_bnd + 1
  end subroutine my_smin

  !===== FINDS CLOSEST POSITION OF A REAL IN A VECTOR OF REALS ====================================
  integer function my_closest(myvector,gp,myvalue)
    implicit none
    integer, intent(in) :: gp
    real(8), intent(in) :: myvalue
    real(8), dimension(gp), intent(in) :: myvector
    real(8), dimension(gp) :: aux

    aux = abs(myvector-myvalue)

    my_closest = minloc(aux, dim=1)
  end function my_closest

  !===== LINEAR INTERPOLATION =====================================================================
  real(8) function my_inter(xvector,yvector,gp_xy,x_inter)
    ! For each value in xvector there is an image in yvector
    ! This subroutine interpolate the value for x_inter that
    ! would have in y_vector
    implicit none
    integer :: gp_xy
    real(8), dimension(gp_xy) :: xvector, yvector
    real(8) :: x_inter
    integer :: x0, x1

    ! Find closest values in vector x
    call my_smin(xvector,gp_xy,x_inter,x0,x1)

    ! Linear interpolation
    my_inter = yvector(x0) + ( (yvector(x1)-yvector(x0)) * &
    ((x_inter-xvector(x0))/(xvector(x1)-xvector(x0))))
  end function my_inter

  !===== SET SEED =================================================================================
  subroutine setseed(my_seed)
    ! sets seed so random numbers
    ! are the same across model
    ! simulations
    !==========================!
    implicit none
    integer, optional ,intent(in) :: my_seed
    integer,allocatable :: seed(:)
    integer :: the_size,j !,lengthr

    call random_seed(size=the_size) ! how big is the intrisic seed?
    allocate(seed(the_size))        ! allocate space for seed
    do j=1,the_size                 ! create the seed
        seed(j)=abs(my_seed)+(j-1)
    enddo
    call random_seed(put=seed)      ! assign the seed
    deallocate(seed)                ! deallocate space

  end subroutine setseed

  !===== STEADY STATE MARKOV CHAIN ================================================================
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

  !===== REGRID A DECISION RULE ===================================================================
  function regrid(gp_old,grid_old,pf_old,gp_new,grid_new)
    ! Uses interpolation to reshape a policy function to a new grid
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

  !===== REALISE SHOCK ============================================================================
  integer function realise(shock, myvec, gp)
    implicit none
    integer :: gp, ind_sh
    real(8) :: shock, aux_sum
    real(8), dimension(gp) :: myvec

    aux_sum = 0.d0
    do ind_sh = 1, gp
      aux_sum  = aux_sum + myvec(ind_sh)
      if (aux_sum.ge.shock) exit
    end do
    realise = ind_sh
  end function realise

  !===== PRINT MATRIX TO SCREEN ===================================================================
  subroutine write_array(a,decimals)
    real(8) :: a(:,:)
    integer             :: decimals, i, j
    character(len=100)  :: colnum
    character(len=100)  :: fmt,fmt2

    write(colnum,*) ubound(a,2)
    write(fmt,*) decimals
    write(fmt2,*) 'f12.'//trim(fmt)

    do i = lbound(a,1), ubound(a,1)
      write(*, '('//trim(colnum)//trim(fmt2)//')' ) (a(i, j), j = lbound(a,2), ubound(a,2) )
    end do
    return
  end subroutine write_array

  !===== PRINT VECTOR TO SCREEN ===================================================================
  subroutine write_vect(a,decimals)
    real(8) :: a(:)
    integer :: decimals, i
    character(len=100) :: fmt,fmt2

    write(fmt,*) decimals
    write(fmt2,*) 'f12.'//trim(fmt)

    do i = lbound(a,1), ubound(a,1)
        write(*, '('//'1'//trim(fmt2)//')' ) (a(i) )
    end do
    return
  end subroutine write_vect

  !===== DISCRETIZE VAR ===========================================================================
  subroutine discretize2vars(A0x,vex,nbar,ntune,pn,yn)
    ! Taken from Gospodinov and Lkhagvasuren (2013)
    ! SEE: https://sites.google.com/site/dlkhagva/var_mmm

    real(8), intent(in)     :: A0x(2,2), vex(2,2)
    integer, intent(in)     :: nbar, ntune
    real(8), intent(out)    :: pn(nbar*nbar,nbar*nbar),yn(2,nbar*nbar)

    real(8)                 :: pmat(2,nbar,nbar,nbar)
    real(8)                 :: z0(nbar),pz0(nbar,nbar),pzz0(nbar,nbar),y1(nbar),y2(nbar)
    real(8)                 :: mu,vact,r,vactx,vactz,rz
    real(8)                 :: p,v1,v1x,px,pz
    real(8), dimension(2,2) :: A0new, vynew, vyold, venew
    integer                 :: nx,n,n1,n2
    integer                 :: na,nb,dummy_exceed,nax,nbx,dummy_exceedx,naz,nbz
    integer                 :: i,j,k,ix,ixx,iz
    real(8), allocatable    :: B(:,:), bvectemp(:),dif1(:)

    integer                 :: i1,i2,i3,i4,ix1,ix2

    nx = ntune+1
    n  = nbar
    n1 = n
    n2 = n

    allocate(B(nx,6),bvectemp(nx),dif1(nx))
    B  = 999.d0

    call rouwenhurst(0.d0,0.d0,1.d0,n,z0,pz0)

    y1 = z0
    y2 = z0

    call var_norm(A0x,vex,A0new,vynew,vyold,venew)

    !print*, "A0x"
    !call write_array(A0x,4)
    !print*, "vex"
    !call write_array(vex,4)
    !print*, "A0new mat"
    !call write_array(A0new,4)
    !print*, "vyold mat"
    !call write_array(vyold,4)
    !print*, "vynew mat"
    !call write_array(venew,4)
    !print*, " "

    do i = 1,n
    do j = 1,n
    do k = 1,2

        mu      = A0new(k,1)*y1(i)+A0new(k,2)*y2(j)
        vact    = venew(k,k)
        r       = sqrt(1.d0-vact)

        call rouwenhurst(r,0.d0,1.d0*sqrt(1.d0-r**2.d0),n,z0,pz0)
        call cal_mu_fast(mu,vact,n,z0,v1,p,na,nb,dummy_exceed)

        if (nx<2)   then
            if (na==nb) then
                pmat(k,i,j,:)=pz0(na,:)
            else
                pmat(k,i,j,:)=p*pz0(na,:)+(1-p)*pz0(nb,:)
            endif
        else
            if (na==nb) then
                pmat(k,i,j,:)=pz0(na,:)
            else
                ixx=0
                do ix=1,nx
                    vactx=max(1e-14,vact*(1.d0-(ix-1.d0)/(nx-1.d0)))
                    call cal_mu_fast(mu,vactx,n,z0,v1x,px,nax,nbx,dummy_exceedx)
                    if (abs(dummy_exceedx)<0.5) then
                        ixx = ixx+1
                        B(ixx,:) = (/ v1x, px, 1.d0*nax, 1.d0*nbx, 1.d0*dummy_exceedx, vactx /)
                    endif
                enddo

                if (ixx<1)  then
                    pmat(k,i,j,:) = p*pz0(na,:)+(1-p)*pz0(nb,:)
                else
                    bvectemp = B(:,1)-vact
                    dif1 = abs(bvectemp)
                    iz = minloc(dif1,dim=1)

                    pz = B(iz,2);
                    naz = int(B(iz,3))
                    nbz = int(B(iz,4))
                    vactz = B(iz,6)

                    rz=sqrt(1-vactz)
                    call rouwenhurst(rz,0.d0,1.d0*sqrt(1.d0-rz**2.d0),n,z0,pzz0)
                    pmat(k,i,j,:) = pz*pzz0(naz,:)+(1-pz)*pzz0(nbz,:)
                endif
            endif

        endif   ! end of (nx<2) conditional

    enddo   ! k loop
    enddo   ! j loop
    enddo   ! i loop

    ! CREATE TRANSITION MATRIX
    ix2 = 0
    do i1 = 1,n
    do i2 = 1,n
        ix2 = ix2 + 1
        do i3 = 1,n
        do i4 = 1,n
            ix1 = (i3-1)*n+i4
            Pn(ix1,ix2) = pmat(1,i1,i2,i3)*pmat(2,i1,i2,i4)
        enddo
        enddo
    enddo
    enddo

    do i = 1,n*n
        Pn(:,i) = Pn(:,i) / sum(Pn(:,i))
    enddo

    Pn = transpose(Pn)

    ! CREATE MATRIX WITH DISCRETE VALUES
    ix = 0
    do i1 = 1,n
    do i2 = 1,n
        ix = ix + 1
        Yn(1,ix) = y1(i1)
        Yn(2,ix) = y2(i2)
    enddo
    enddo

    Yn(1,:) = Yn(1,:)*sqrt(vyold(1,1))
    Yn(2,:) = Yn(2,:)*sqrt(vyold(2,2))
  end subroutine discretize2vars

  !===== NORM =====================================================================================
  subroutine var_norm(A,ve,Anew,vynew,vyold,venew)
    implicit none
    real(8), intent(in)    :: A(:,:), ve(:,:)
    real(8), intent(out)   :: Anew(:,:),vynew(:,:),vyold(:,:),venew(:,:)
    real(8), allocatable   :: v(:,:),v0(:,:)
    integer                :: nn,i,j
    real(8)                :: dif

    dif = 100.0
    nn  = ubound(A,1)

    allocate(v0(nn,nn))
    v0 = 0.d0

    do while (dif>1e-12)
        v=matmul(A,matmul(v0,transpose(A)))+ve
        dif=maxval(v-v0)
        v0=v
    enddo

    vyold=v0

    do i=1,nn
        venew(i,i)=ve(i,i)/vyold(i,i)
        do j=1,nn
            Anew(i,j)=A(i,j)*sqrt(vyold(j,j))/sqrt(vyold(i,i))
        enddo
    enddo

    do i=1,nn
    do j=1,nn
        vynew(i,j) = vyold(i,j)/(sqrt(vyold(i,i))*sqrt(vyold(j,j)) )
    enddo
    enddo
  end subroutine var_norm

  subroutine cal_mu_fast(mu,v0,n,z,v1,p,na,nb,dummy_exceed)
    real(8), intent(in)     :: mu,v0
    integer, intent(in)     :: n
    real(8), intent(in)     :: z(:)
    real(8), intent(out)    :: v1,p
    integer, intent(out)    :: na,nb,dummy_exceed
    real(8), allocatable    :: zm(:)

    allocate(zm(n))

    zm = z*sqrt(1.d0-v0)

    if ( mu>=zm(n) )    then
        dummy_exceed=1
        na=n
        nb=n
        p=0.d00
        v1=v0
    elseif ( mu<=zm(1) )    then
        dummy_exceed=-1
        na=1
        nb=1
        p=1.d0
        v1=v0
    else
        dummy_exceed=0
        na=1+floor((mu-zm(1))/(zm(2)-zm(1)))
        nb=na+1

        p=(zm(nb)-mu)/(zm(nb)-zm(na))

        v1=v0+p*(1-p)*(zm(nb)-zm(na))**2.d0
    endif
  end subroutine cal_mu_fast

  !===== ROUWENHORST ==============================================================================
  subroutine rouwenhurst(rho, mu_eps, sigma_eps, n, zvect, pmat)
    ! discretizes an ar(1) process, with persistence parameter
    ! 'rho', mean 'mu_eps' and standard deviation 'sigma_eps'
    ! stores results in zvect(n) and pmat(n,n)

    implicit none

    real(8), intent(in):: rho, mu_eps, sigma_eps   !, coverage
    integer, intent(in):: n
    real(8), intent(out):: zvect(n)
    real(8), intent(out):: pmat(n,n)

    real(8) mu_z, sigma_z, q, eps
    real(8), allocatable, dimension(:,:):: p1, p2
    integer status, i, j

    mu_z = mu_eps/(1-rho)
    sigma_z = sigma_eps/sqrt(1-rho**2.d0)

    q = (rho+1)/2
    eps = sqrt(dble(n-1)) * sigma_z

    if (n == 1) then
        pmat = 1.0d0
        zvect = mu_z
        return
    else if (n == 2) then
        pmat = reshape((/q, 1-q, 1-q, q/),(/2,2/))
        zvect = (/mu_z-eps,mu_z+eps/)
        return
    end if

    allocate(p1(2,2),stat=status)
    p1 = reshape((/q, 1-q, 1-q, q/),(/2,2/))

    do i=2,n-1
        allocate(p2(i+1,i+1),stat=status)
        p2 = q * reshape( (/  (/(p1(:,j),0.0d0 ,j=1,i)/) ,  (/(0.0d0,j=1,i+1)/)    /), (/i+1,i+1/) ) + &
             (1-q) * reshape( (/  (/(0.0d0,j=1,i+1)/), (/ (p1(:,j),0.0d0 ,j=1,i)/)   /) ,   (/i+1,i+1/) ) + &
             (1-q) * reshape( (/  (/ (0.0d0,p1(:,j) ,j=1,i) /) ,  (/(0.0d0,j=1,i+1)/)  /), (/i+1,i+1/) ) + &
             q * reshape( (/ (/(0.0d0,j=1,i+1)/), (/(0.0d0,p1(:,j) ,j=1,i)/)   /) ,   (/i+1,i+1/) )

        p2(2:i,:) = p2(2:i,:)/2

        deallocate(p1,stat=status)

        if (i==n-1) then
            pmat = p2
        else
            allocate(p1(i+1,i+1), stat=status)
            p1 = p2
        end if

        deallocate(p2,stat=status)
    end do

    zvect = (/ (mu_z-eps + (2.d0*eps)*i/(n-1),i=0,n-1) /)
  end subroutine rouwenhurst
end module Utils
