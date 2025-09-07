! Gaussian_Copula_Meteorology.f90
! Self-contained Gaussian Copula demo for meteorology
! Variables:
!  - Temperature T (°C): Normal(mu_T, sd_T)
!  - Relative Humidity RH (%): Logit-Normal mapped to [0,100]
!  - Wind speed W (m/s): Weibull(k, lambda)
!
! Dependence via Gaussian Copula in Z-space:
!   1) Draw z_raw ~ N(0, I)
!   2) z = L * z_raw   where L = chol(CorrZ)
!   3) u = Phi(z)      (componentwise)
!   4) x_i = F_i^{-1}(u_i) using chosen marginals
!
! Prints:
!  - Estimated P(T >= hot_thr  &  RH <= dry_thr  &  W >= wind_thr)
!  - Sample Z-space correlation (should approximate CorrZ)

program copula_meteo
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  ! ----- kinds & constants -----
  integer, parameter :: dp = real64
  real(dp), parameter :: PI = acos(-1.0_dp)
  real(dp), parameter :: EPS_JITTER = 1.0e-12_dp

  ! ----- dimensions -----
  integer, parameter :: p = 3    ! [T, RH, W]

  ! ----- variables -----
  integer :: i, n, r, c, count
  real(dp) :: CorrZ(p,p), L(p,p)
  real(dp) :: zraw(p), z(p), u(p)
  real(dp) :: T_val, RH_val, W_val
  real(dp) :: hot_thr, dry_thr, wind_thr, prob_est

  ! streaming moments for Z (diagnostics)
  real(dp) :: meanZ(p), oldMeanZ(p)
  real(dp) :: M2Z(p,p), covZ(p,p), corrZhat(p,p)

  ! ----- marginal params (edit for your site) -----
  ! Temperature ~ Normal
  real(dp), parameter :: mu_T = 30.0_dp, sd_T = 6.0_dp

  ! RH (%) ~ Logit-Normal over [0,100]
  ! Choose logit-scale mean/sd so median ~ 35% and reasonable spread
  real(dp), parameter :: RH_med = 0.35_dp
  real(dp), parameter :: mu_L = log(RH_med/(1.0_dp-RH_med))
  real(dp), parameter :: sd_L = 0.8_dp

  ! Wind (m/s) ~ Weibull(k, lambda)
  real(dp), parameter :: k_w = 2.0_dp, lambda_w = 8.0_dp

  ! ----- copula correlation in Z-space -----
  ! Adjust to your dependence. Must be symmetric, SPD.
  CorrZ = reshape([ 1.0_dp, -0.65_dp,  0.30_dp,  &
                   -0.65_dp,  1.0_dp, -0.20_dp,  &
                    0.30_dp, -0.20_dp,  1.0_dp ], [p,p])

  ! ----- simulation settings / thresholds -----
  n        = 100000
  hot_thr  = 35.0_dp
  dry_thr  = 20.0_dp
  wind_thr = 10.0_dp

  ! ----- factorize CorrZ with tiny jitter for robustness -----
  if (.not. cholesky_spd(CorrZ + EPS_JITTER*eye(p), L)) then
     print *, "ERROR: CorrZ is not SPD even after jitter. Shrink correlations and retry."
     stop 1
  end if

  call seed_with_time()

  prob_est = 0.0_dp
  meanZ    = 0.0_dp
  M2Z      = 0.0_dp
  count    = 0

  do i = 1, n
     call randn_vec(p, zraw)    ! z_raw ~ N(0, I)
     z = matmul(L, zraw)        ! correlate

     ! uniforms in (0,1) from standard normal CDF
     u(1) = phi(z(1))
     u(2) = phi(z(2))
     u(3) = phi(z(3))

     ! marginals
     T_val  = q_normal(mu_T, sd_T, u(1))
     RH_val = 100.0_dp * logistic( mu_L + sd_L * invnorm(u(2)) )
     W_val  = q_weibull(k_w, lambda_w, u(3))

     ! compound event
     if (T_val >= hot_thr .and. RH_val <= dry_thr .and. W_val >= wind_thr) then
        prob_est = prob_est + 1.0_dp
     end if

     ! streaming moments for Z (matrix Welford)
     count    = count + 1
     oldMeanZ = meanZ
     meanZ    = oldMeanZ + (z - oldMeanZ) / real(count, dp)
     do r = 1, p
        do c = 1, p
           M2Z(r,c) = M2Z(r,c) + (z(r) - meanZ(r)) * (z(c) - oldMeanZ(c))
        end do
     end do
  end do

  prob_est = prob_est / real(n, dp)

  ! sample corr in Z-space (should match CorrZ up to Monte Carlo error)
  if (count > 1) then
     covZ = M2Z / real(count - 1, dp)
  else
     covZ = 0.0_dp
  end if

  do r = 1, p
     do c = 1, p
        corrZhat(r,c) = covZ(r,c) / sqrt( max(covZ(r,r),1e-30_dp) * max(covZ(c,c),1e-30_dp) )
     end do
  end do

  write(*,'(a,f12.6)') "Estimated P( T>=hot & RH<=dry & W>=wind ) = ", prob_est
  print *, "Sample Z-space correlation (should approximate CorrZ):"
  call print_mat(corrZhat)

contains
  !===================== Helpers =====================

  pure function eye(n) result(Id)
    ! Identity matrix (n x n)
    integer, intent(in) :: n
    real(dp) :: Id(n,n)
    integer :: ii
    Id = 0.0_dp
    do ii = 1, n
       Id(ii, ii) = 1.0_dp
    end do
  end function eye

  logical function cholesky_spd(A, L)
    ! Cholesky factorization A = L*L^T for symmetric positive-definite A
    real(dp), intent(in)  :: A(:,:)
    real(dp), intent(out) :: L(size(A,1), size(A,2))
    integer :: i, j, k, n
    real(dp) :: sumv
    n = size(A,1)
    L = 0.0_dp
    cholesky_spd = .true.
    do i = 1, n
       do j = 1, i
          sumv = A(i,j)
          do k = 1, j-1
             sumv = sumv - L(i,k)*L(j,k)
          end do
          if (i == j) then
             if (sumv <= 0.0_dp) then
                cholesky_spd = .false.
                return
             end if
             L(i,j) = sqrt(sumv)
          else
             L(i,j) = sumv / L(j,j)
          end if
       end do
    end do
  end function cholesky_spd

  subroutine randn_vec(m, z)
    ! Box–Muller vector of length m (uses pairs)
    integer, intent(in) :: m
    real(dp), intent(out) :: z(m)
    integer :: idx
    real(dp) :: u1, u2, radius, theta
    idx = 1
    do
       call random_number(u1)
       call random_number(u2)
       if (u1 <= 1.0e-16_dp) cycle
       radius = sqrt(-2.0_dp*log(u1))
       theta  = 2.0_dp*PI*u2
       z(idx) = radius * cos(theta)
       if (idx + 1 <= m) then
          z(idx+1) = radius * sin(theta)
          idx = idx + 2
       else
          exit
       end if
       if (idx > m) exit
    end do
  end subroutine randn_vec

  subroutine seed_with_time()
    ! Reasonable portable seeding from system time
    integer :: t(8), nseed, i
    integer, allocatable :: seed(:)
    call date_and_time(values=t)
    call random_seed(size=nseed)
    allocate(seed(nseed))
    do i = 1, nseed
       seed(i) = max(1, modulo( t(8) + 1000*t(7) + 100000*t(6) + 37*i, 2147483647 ))
    end do
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine seed_with_time

  ! ---- Gaussian CDF and quantile ----
  pure real(dp) function phi(x) result(pu)
    real(dp), intent(in) :: x
    pu = 0.5_dp * erfc( -x / sqrt(2.0_dp) )
  end function phi

  pure real(dp) function invnorm(u) result(z)
    ! Acklam’s rational approximation for Phi^{-1}(u), u in (0,1)
    real(dp), intent(in) :: u
    real(dp) :: x, q, t
    real(dp), parameter :: a1=-3.969683028665376d+01, a2= 2.209460984245205d+02, a3=-2.759285104469687d+02, &
                           a4= 1.383577518672690d+02, a5=-3.066479806614716d+01, a6= 2.506628277459239d+00
    real(dp), parameter :: b1=-5.447609879822406d+01, b2= 1.615858368580409d+02, b3=-1.556989798598866d+02, &
                           b4= 6.680131188771972d+01, b5=-1.328068155288572d+01
    real(dp), parameter :: c1=-7.784894002430293d-03, c2=-3.223964580411365d-01, c3=-2.400758277161838d+00, &
                           c4=-2.549732539343734d+00, c5= 4.374664141464968d+00, c6= 2.938163982698783d+00
    real(dp), parameter :: d1= 7.784695709041462d-03, d2= 3.224671290700398d-01, d3= 2.445134137142996d+00, &
                           d4= 3.754408661907416d+00
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow

    x = max( min(u, 1.0_dp-1.0e-16_dp), 1.0e-16_dp )
    if (x < plow) then
       q = sqrt(-2.0_dp*log(x))
       z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
       z = -z
    else if (x <= phigh) then
      q = x - 0.5_dp
      t = q*q
      z = (((((a1*t+a2)*t+a3)*t+a4)*t+a5)*t+a6)*q / (((((b1*t+b2)*t+b3)*t+b4)*t+b5)*t+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-x))
      z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
  end function invnorm

  ! ---- marginal quantiles ----
  pure real(dp) function q_normal(mu, sd, u) result(x)
    real(dp), intent(in) :: mu, sd, u
    x = mu + sd * invnorm(u)
  end function q_normal

  pure real(dp) function logistic(t) result(y)
    real(dp), intent(in) :: t
    y = 1.0_dp / (1.0_dp + exp(-t))
  end function logistic

  pure real(dp) function q_weibull(k, lambda, u) result(x)
    real(dp), intent(in) :: k, lambda, u
    real(dp) :: uu
    uu = max( min(u, 1.0_dp-1.0e-16_dp), 1.0e-16_dp )
    x  = lambda * (-log(1.0_dp - uu))**(1.0_dp / k)
  end function q_weibull

  subroutine print_mat(A)
    real(dp), intent(in) :: A(:,:)
    integer :: irow, jcol, nrow, ncol
    nrow = size(A,1); ncol = size(A,2)
    do irow = 1, nrow
       write(*,'(100(f12.6,1x))') (A(irow,jcol), jcol=1,ncol)
    end do
  end subroutine print_mat

end program copula_meteo

