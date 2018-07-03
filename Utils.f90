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
  subroutine linspace(my_start, my_stop, n, grid)
    implicit none
    integer, intent(in) :: n
    integer :: i
    real(8), intent(in):: my_start, my_stop
    real(8) :: step
    real(8), dimension(n), intent(out) :: grid

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
  end subroutine linspace

  !===== LOGGRID ==================================================================================
  subroutine loggrid(my_start, my_stop, n, grid)
    implicit none
    integer, intent(in) :: n
    integer :: i
    real(8), intent(in) :: my_start, my_stop
    real(8) ::step
    real(8), dimension(n), intent(out) :: grid

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
  end subroutine loggrid

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

    call linspace(ymin, ymax, gp, values)

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
  subroutine my_inter(xvector,yvector,gp_xy,x_inter,y_inter)
    ! For each value in xvector there is an image in yvector
    ! This subroutine interpolate the value for x_inter that
    ! would have in y_vector
    implicit none
    integer, intent(in) :: gp_xy
    real(8), dimension(gp_xy), intent(in) :: xvector, yvector
    real(8), intent(in) :: x_inter
    real(8), intent(out) :: y_inter
    integer :: x0, x1

    ! Find closest values in vector x
    call my_smin(xvector,gp_xy,x_inter,x0,x1)

    ! Linear interpolation
    y_inter = yvector(x0) + ( (yvector(x1)-yvector(x0)) * &
    ((x_inter-xvector(x0))/(xvector(x1)-xvector(x0))))
  end subroutine my_inter
end module Utils
