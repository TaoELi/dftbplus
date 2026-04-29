!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains Types and subroutine to build up and query a Slater-Koster table where the integrals are
!> specified on an equidistant grid.
module dftbp_dftb_slakoeqgrid
  use dftbp_common_accuracy, only : distFudge, distFudgeOld, dp
  use dftbp_io_message, only : error
  use dftbp_math_interpolation, only : freeCubicSpline, poly5ToZero, polyInterUniform
  implicit none

  private
  public :: TSlakoEqGrid, init
  public :: getSKIntegrals, getSKIntegralDerivs, getNIntegrals, getCutoff
  public :: skEqGridOld, skEqGridNew


  !> Represents an equally spaced Slater-Koster grid
  type TSlakoEqGrid
    private
    integer :: nGrid
    integer :: nInteg
    real(dp) :: dist
    real(dp), allocatable :: skTab(:,:)
    integer :: skIntMethod
    logical :: tInit = .false.
  end type TSlakoEqGrid


  !> Initialises SlakoEqGrid.
  interface init
    module procedure SlakoEqGrid_init
  end interface init


  !> Returns the integrals for a given distance.
  interface getSKIntegrals
    module procedure SlakoEqGrid_getSKIntegrals
  end interface getSKIntegrals

  !> Returns the integrals and their radial derivatives for a given distance.
  interface getSKIntegralDerivs
    module procedure SlakoEqGrid_getSKIntegralDerivs
  end interface getSKIntegralDerivs


  !> Returns the number of integrals the table contains
  interface getNIntegrals
    module procedure SlakoEqGrid_getNIntegrals
  end interface getNIntegrals


  !> Returns the cutoff of the interaction.
  interface getCutoff
    module procedure SlakoEqGrid_getCutoff
  end interface getCutoff

  ! Interpolation methods

  !> Historical method
  integer, parameter :: skEqGridOld = 1

  !> Current method
  integer, parameter :: skEqGridNew = 2

  ! Nr. of grid points to use for the polynomial interpolation

  !> Historical choice
  integer, parameter :: nInterOld_ = 3

  !> Present choice
  integer, parameter :: nInterNew_ = 8

  ! Nr. of grid points on the right of the interpolated point.

  ! For an odd number of intervals, the number of right points should be bigger than the number of
  ! left points, to remain compatible with the old code.


  !> Value nRightInterOld: floor(real(nInterOld_, dp) / 2.0_dp + 0.6_dp)
  integer, parameter :: nRightInterOld_ = 2

  !> Value nRightInterNew: floor(real(nInterNew_, dp) / 2.0_dp + 0.6_dp)
  integer, parameter :: nRightInterNew_ = 4


  !> Displacement for deriving interpolated polynomials
  real(dp), parameter :: deltaR_ = 1e-5_dp

contains


  !> Initialises SlakoEqGrid.
  subroutine SlakoEqGrid_init(this, dist, table, skIntMethod)

    !> SlakoEqGrid instance.
    type(TSlakoEqGrid), intent(out) :: this

    !> Distance between the grid points.
    real(dp), intent(in) :: dist

    !> Slater-Koster table (first entry belongs to first grid point)
    real(dp), intent(in) :: table(:,:)

    !> Method for the interpolation between the entries.
    integer, intent(in) :: skintMethod

    @:ASSERT(.not. this%tInit)
    @:ASSERT(dist >= 0.0_dp)
    @:ASSERT(skIntMethod == skEqGridOld .or. skIntMethod == skEqGridNew)

    this%dist = dist
    this%nGrid = size(table, dim=1)
    this%nInteg = size(table, dim=2)
    allocate(this%skTab(this%nGrid, this%nInteg))
    this%skTab(:,:) = table(:,:)
    this%skIntMethod = skIntMethod
    this%tInit = .true.

  end subroutine SlakoEqGrid_init


  !> Returns the integrals for a given distance.
  subroutine SlakoEqGrid_getSKIntegrals(this, sk, dist)

    !> SlakoEqGrid instance.
    type(TSlakoEqGrid), intent(in) :: this

    !> Contains the interpolated integrals on exit
    real(dp), intent(out) :: sk(:)

    !> Distance for which the integrals should be interpolated.
    real(dp), intent(in) :: dist

    @:ASSERT(this%tInit)
    @:ASSERT(size(sk) >= this%nInteg)
    @:ASSERT(dist >= 0.0_dp)

    if (this%skIntMethod == skEqGridOld) then
      call SlakoEqGrid_interOld_(this, sk, dist)
    else
      call SlakoEqGrid_interNew_(this, sk, dist)
    end if

  end subroutine SlakoEqGrid_getSKIntegrals


  !> Returns the integrals and their radial derivatives for a given distance.
  subroutine SlakoEqGrid_getSKIntegralDerivs(this, sk, dSkDr, dist)

    !> SlakoEqGrid instance.
    type(TSlakoEqGrid), intent(in) :: this

    !> Contains the interpolated integrals on exit.
    real(dp), intent(out) :: sk(:)

    !> Contains the radial derivative of the interpolated integrals on exit.
    real(dp), intent(out) :: dSkDr(:)

    !> Distance for which the integrals should be interpolated.
    real(dp), intent(in) :: dist

    @:ASSERT(this%tInit)
    @:ASSERT(size(sk) >= this%nInteg)
    @:ASSERT(size(dSkDr) >= this%nInteg)
    @:ASSERT(dist >= 0.0_dp)

    if (this%skIntMethod == skEqGridOld) then
      call SlakoEqGrid_interOldWithDeriv_(this, sk, dSkDr, dist)
    else
      call SlakoEqGrid_interNewWithDeriv_(this, sk, dSkDr, dist)
    end if

  end subroutine SlakoEqGrid_getSKIntegralDerivs


  !> Returns the number of intgrals the table contains
  function SlakoEqGrid_getNIntegrals(this) result(nInt)

    !> SlakoEqGrid instance.
    type(TSlakoEqGrid), intent(in) :: this

    !> Number of integrals.
    integer :: nInt

    nInt = this%nInteg

  end function SlakoEqGrid_getNIntegrals


  !> Returns the cutoff of the interaction.
  function SlakoEqGrid_getCutoff(this) result(cutoff)

    !>  SlakoEqGrid instance.
    type(TSlakoEqGrid), intent(in) :: this

    !> Grid cutoff
    real(dp) :: cutoff

    cutoff = real(this%nGrid, dp) * this%dist
    if (this%skIntMethod == skEqGridOld) then
      cutoff = cutoff + distFudgeOld
    else
      cutoff = cutoff + distFudge
    end if

  end function SlakoEqGrid_getCutoff


  !> Inter- and extrapolation for SK-tables, new method.
  subroutine SlakoEqGrid_interNew_(this, dd, rr)

    !> SlakoEqGrid table on equiv. grid
    type(TSlakoEqGrid), intent(in) :: this

    !> Output table of interpolated values.
    real(dp), intent(out) :: dd(:)

    !> Distance between two atoms of interest
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterNew_), ya(nInterNew_), yb(this%nInteg,nInterNew_), y1, y1p, y1pp
    real(dp) :: incr, dr, rMax, y0(this%nInteg), y2(this%nInteg)
    integer :: leng, ind, iLast
    integer :: ii

    real(dp), parameter :: invdistFudge = -1.0_dp / distFudge

    leng = this%nGrid
    incr = this%dist
    rMax = real(leng, dp) * incr + distFudge
    ind = floor(rr / incr)

    !! Consistency check, does the SK-table contain enough entries?
    if (leng < nInterNew_ + 1) then
      call error("SlakoEqGrid: Not enough points in the SK-table for &
          &interpolation!")
    end if

    dd(:) = 0.0_dp
    if (rr >= rMax) then
      !! Beyond last grid point + distFudge => no interaction
      dd(:) = 0.0_dp
    elseif (ind < leng) then
      !! Closer to origin than last grid point => polynomial fit
      iLast = min(leng, ind + nRightInterNew_)
      iLast = max(iLast, nInterNew_)
      do ii = 1, nInterNew_
        xa(ii) = real(iLast - nInterNew_ + ii, dp) * incr
      end do
      yb = transpose(this%skTab(iLast-nInterNew_+1:iLast,:this%nInteg))
      dd(:this%nInteg) = polyInterUniform(xa, yb, rr)
    else
      !! Beyond the grid => extrapolation with polynomial of 5th order
      dr = rr - rMax
      iLast = leng
      do ii = 1, nInterNew_
        xa(ii) = real(iLast - nInterNew_ + ii, dp) * incr
      end do
      yb = transpose(this%skTab(iLast-nInterNew_+1:iLast,:this%nInteg))
      y0 = polyInterUniform(xa, yb, xa(nInterNew_) - deltaR_)
      y2 = polyInterUniform(xa, yb, xa(nInterNew_) + deltaR_)
      do ii = 1, this%nInteg
        ya(:) = this%skTab(iLast-nInterNew_+1:iLast, ii)
        y1 = ya(nInterNew_)
        y1p = (y2(ii) - y0(ii)) / (2.0_dp * deltaR_)
        y1pp = (y2(ii) + y0(ii) - 2.0_dp * y1) / (deltaR_ * deltaR_)
        dd(ii) = poly5ToZero(y1, y1p, y1pp, dr, -1.0_dp * distFudge, invDistFudge)
      end do
    end if

  end subroutine SlakoEqGrid_interNew_


  !> Inter- and extrapolation for SK-tables, new method, with radial derivative.
  subroutine SlakoEqGrid_interNewWithDeriv_(this, dd, dddr, rr)

    !> SlakoEqGrid table on equiv. grid.
    type(TSlakoEqGrid), intent(in) :: this

    !> Output table of interpolated values.
    real(dp), intent(out) :: dd(:)

    !> Output table of radial derivatives.
    real(dp), intent(out) :: dddr(:)

    !> Distance between two atoms of interest.
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterNew_), yb(this%nInteg,nInterNew_)
    real(dp) :: incr, dr, rMax, y0(this%nInteg), y1p(this%nInteg), y1pp(this%nInteg)
    integer :: leng, ind, iLast
    integer :: ii

    real(dp), parameter :: invdistFudge = -1.0_dp / distFudge

    leng = this%nGrid
    incr = this%dist
    rMax = real(leng, dp) * incr + distFudge
    ind = floor(rr / incr)

    if (leng < nInterNew_ + 1) then
      call error("SlakoEqGrid: Not enough points in the SK-table for &
          &interpolation!")
    end if

    dd(:) = 0.0_dp
    dddr(:) = 0.0_dp
    if (rr >= rMax) then
      dd(:) = 0.0_dp
      dddr(:) = 0.0_dp
    elseif (ind < leng) then
      iLast = min(leng, ind + nRightInterNew_)
      iLast = max(iLast, nInterNew_)
      do ii = 1, nInterNew_
        xa(ii) = real(iLast - nInterNew_ + ii, dp) * incr
      end do
      yb = transpose(this%skTab(iLast-nInterNew_+1:iLast,:this%nInteg))
      call polyInterUniformDerivs_(xa, yb, rr, dd(:this%nInteg), dddr(:this%nInteg))
    else
      dr = rr - rMax
      iLast = leng
      do ii = 1, nInterNew_
        xa(ii) = real(iLast - nInterNew_ + ii, dp) * incr
      end do
      yb = transpose(this%skTab(iLast-nInterNew_+1:iLast,:this%nInteg))
      call polyInterUniformDerivs_(xa, yb, xa(nInterNew_), y0, y1p, y1pp)
      do ii = 1, this%nInteg
        dd(ii) = poly5ToZero(y0(ii), y1p(ii), y1pp(ii), dr, -1.0_dp * distFudge,&
            & invDistFudge)
        dddr(ii) = poly5ToZeroDeriv_(y0(ii), y1p(ii), y1pp(ii), dr,&
            & -1.0_dp * distFudge, invDistFudge)
      end do
    end if

  end subroutine SlakoEqGrid_interNewWithDeriv_


  !> Inter- and extra-polation for SK-tables equivalent to the old DFTB code.
  subroutine SlakoEqGrid_interOld_(this, dd, rr)

    !> Data structure for SK interpolation
    type(TSlakoEqGrid), intent(in) :: this

    !> Output table of interpolated values.
    real(dp), intent(out) :: dd(:)

    !> Distance between two atoms of interest
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterOld_), yb(this%nInteg,nInterOld_),y0, y1, y2, y1p, y1pp
    real(dp) :: incr, dr
    integer :: leng, ind, mInd, iLast
    integer :: ii
    real(dp) :: r1, r2
    real(dp) :: invdistFudge

    leng = this%nGrid
    incr = this%dist
    mInd = leng + floor(distFudgeOld/incr)
    ind = floor(rr / incr)

    invdistFudge = -1.0_dp / (real(mInd - leng -1, dp) * incr)

    !! Consistency check, does the SK-table contain enough entries?
    if (leng < nInterOld_ + 1) then
      call error("skspar: Not enough points in the SK-table for interpolation!")
    end if

    dd(:) = 0.0_dp
    if (ind < leng-1) then
      !! Distance closer than penultimate grid point => polynomial fit
      iLast = min(leng, ind + nRightInterOld_)
      iLast = max(iLast, nInterOld_)
      do ii = 1, nInterOld_
        xa(ii) = real(iLast - nInterOld_ + ii, dp) * incr
      end do
      yb = transpose(this%skTab(iLast-nInterOld_+1:iLast,:this%nInteg))
      dd(:this%nInteg) = polyInterUniform(xa, yb, rr)
    elseif (ind < leng) then
      !! Distance between penultimate and last grid point => free cubic spline
      dr = rr - real(leng - 1, dp) * incr
      do ii = 1, this%nInteg
        y0 = this%skTab(leng-2, ii)
        y1 = this%skTab(leng-1, ii)
        y2 = this%skTab(leng, ii)
        y1p = (y2 - y0) / (2.0_dp * incr)
        y1pp = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, y1p, y1pp, incr, y2, dr, dd(ii))
      end do
    elseif (ind < mInd - 1) then
      !! Extrapolation
      dr = rr - real(mInd - 1, dp) * incr
      do ii = 1, this%nInteg
        y0 = this%skTab(leng-2, ii)
        y1 = this%skTab(leng-1, ii)
        y2 = this%skTab(leng, ii)
        r1 = (y2 - y0) / (2.0_dp * incr)
        r2 = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, r1, r2, incr, y2, incr, yp=y1p, ypp=y1pp)
        dd(ii) = poly5ToZero(y2, y1p, y1pp, dr,&
            & -1.0_dp * real(mInd - leng -1, dp)*incr, invdistFudge)
      end do
    else
      !! Dist. greater than tabulated sk range + distFudge => no interaction
      dd(:) = 0.0_dp
    end if

  end subroutine SlakoEqGrid_interOld_


  !> Inter- and extra-polation for SK-tables equivalent to the old DFTB code, with derivative.
  subroutine SlakoEqGrid_interOldWithDeriv_(this, dd, dddr, rr)

    !> Data structure for SK interpolation.
    type(TSlakoEqGrid), intent(in) :: this

    !> Output table of interpolated values.
    real(dp), intent(out) :: dd(:)

    !> Output table of radial derivatives.
    real(dp), intent(out) :: dddr(:)

    !> Distance between two atoms of interest.
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterOld_), yb(this%nInteg,nInterOld_), y0, y1, y2, y1p, y1pp, yppTmp
    real(dp) :: incr, dr
    integer :: leng, ind, mInd, iLast
    integer :: ii
    real(dp) :: r1, r2
    real(dp) :: invdistFudge

    leng = this%nGrid
    incr = this%dist
    mInd = leng + floor(distFudgeOld/incr)
    ind = floor(rr / incr)

    invdistFudge = -1.0_dp / (real(mInd - leng -1, dp) * incr)

    if (leng < nInterOld_ + 1) then
      call error("skspar: Not enough points in the SK-table for interpolation!")
    end if

    dd(:) = 0.0_dp
    dddr(:) = 0.0_dp
    if (ind < leng-1) then
      iLast = min(leng, ind + nRightInterOld_)
      iLast = max(iLast, nInterOld_)
      do ii = 1, nInterOld_
        xa(ii) = real(iLast - nInterOld_ + ii, dp) * incr
      end do
      yb = transpose(this%skTab(iLast-nInterOld_+1:iLast,:this%nInteg))
      call polyInterUniformDerivs_(xa, yb, rr, dd(:this%nInteg), dddr(:this%nInteg))
    elseif (ind < leng) then
      dr = rr - real(leng - 1, dp) * incr
      do ii = 1, this%nInteg
        y0 = this%skTab(leng-2, ii)
        y1 = this%skTab(leng-1, ii)
        y2 = this%skTab(leng, ii)
        y1p = (y2 - y0) / (2.0_dp * incr)
        y1pp = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, y1p, y1pp, incr, y2, dr, dd(ii), yp=dddr(ii), ypp=yppTmp)
      end do
    elseif (ind < mInd - 1) then
      dr = rr - real(mInd - 1, dp) * incr
      do ii = 1, this%nInteg
        y0 = this%skTab(leng-2, ii)
        y1 = this%skTab(leng-1, ii)
        y2 = this%skTab(leng, ii)
        r1 = (y2 - y0) / (2.0_dp * incr)
        r2 = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, r1, r2, incr, y2, incr, yp=y1p, ypp=y1pp)
        dd(ii) = poly5ToZero(y2, y1p, y1pp, dr,&
            & -1.0_dp * real(mInd - leng -1, dp)*incr, invdistFudge)
        dddr(ii) = poly5ToZeroDeriv_(y2, y1p, y1pp, dr,&
            & -1.0_dp * real(mInd - leng -1, dp)*incr, invdistFudge)
      end do
    else
      dd(:) = 0.0_dp
      dddr(:) = 0.0_dp
    end if

  end subroutine SlakoEqGrid_interOldWithDeriv_


  !> Polynomial interpolation through uniformly spaced points with first and second derivatives.
  subroutine polyInterUniformDerivs_(xp, yp, xx, yy, yyp, yypp)

    !> x-coordinates of the fit points.
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points.
    real(dp), intent(in) :: yp(:,:)

    !> The point where the polynomial should be evaluated.
    real(dp), intent(in) :: xx

    !> The polynomial value.
    real(dp), intent(out) :: yy(:)

    !> The first derivative.
    real(dp), intent(out) :: yyp(:)

    !> The second derivative.
    real(dp), intent(out), optional :: yypp(:)

    integer :: iFit, jFit, kFit, lFit, nFit
    real(dp) :: basis, basisPrime, basisSecond, term

    nFit = size(xp)
    @:ASSERT(size(yp, dim=2) == nFit)
    @:ASSERT(size(yy) >= size(yp, dim=1))
    @:ASSERT(size(yyp) >= size(yp, dim=1))
    if (present(yypp)) then
      @:ASSERT(size(yypp) >= size(yp, dim=1))
    end if

    yy(:) = 0.0_dp
    yyp(:) = 0.0_dp
    if (present(yypp)) then
      yypp(:) = 0.0_dp
    end if

    do iFit = 1, nFit
      basis = 1.0_dp
      do jFit = 1, nFit
        if (jFit == iFit) cycle
        basis = basis * (xx - xp(jFit)) / (xp(iFit) - xp(jFit))
      end do
      yy(:size(yp, dim=1)) = yy(:size(yp, dim=1)) + yp(:, iFit) * basis

      basisPrime = 0.0_dp
      do jFit = 1, nFit
        if (jFit == iFit) cycle
        term = 1.0_dp / (xp(iFit) - xp(jFit))
        do kFit = 1, nFit
          if (kFit == iFit .or. kFit == jFit) cycle
          term = term * (xx - xp(kFit)) / (xp(iFit) - xp(kFit))
        end do
        basisPrime = basisPrime + term
      end do
      yyp(:size(yp, dim=1)) = yyp(:size(yp, dim=1)) + yp(:, iFit) * basisPrime

      if (present(yypp)) then
        basisSecond = 0.0_dp
        do jFit = 1, nFit
          if (jFit == iFit) cycle
          do kFit = 1, nFit
            if (kFit == iFit .or. kFit == jFit) cycle
            term = 1.0_dp / ((xp(iFit) - xp(jFit)) * (xp(iFit) - xp(kFit)))
            do lFit = 1, nFit
              if (lFit == iFit .or. lFit == jFit .or. lFit == kFit) cycle
              term = term * (xx - xp(lFit)) / (xp(iFit) - xp(lFit))
            end do
            basisSecond = basisSecond + term
          end do
        end do
        yypp(:size(yp, dim=1)) = yypp(:size(yp, dim=1)) + yp(:, iFit) * basisSecond
      end if
    end do

  end subroutine polyInterUniformDerivs_


  !> First derivative of the fifth-order cutoff polynomial used by poly5ToZero().
  pure function poly5ToZeroDeriv_(y0, y0p, y0pp, xx, dx, invdx) result(yyp)

    !> Value of the polynomial at x = dx.
    real(dp), intent(in) :: y0

    !> Value of the first derivative at x = dx.
    real(dp), intent(in) :: y0p

    !> Value of the second derivative at x = dx.
    real(dp), intent(in) :: y0pp

    !> The point where the polynomial should be differentiated.
    real(dp), intent(in) :: xx

    !> The point where the polynomial matches the provided values.
    real(dp), intent(in) :: dx

    !> Reciprocal of dx.
    real(dp), intent(in) :: invdx

    !> First derivative of the polynomial at xx.
    real(dp) :: yyp

    real(dp) :: dx1, dx2, dd, ee, ff, xr

    dx1 = y0p * dx
    dx2 = y0pp * dx * dx
    dd =  10.0_dp * y0 - 4.0_dp * dx1 + 0.5_dp * dx2
    ee = -15.0_dp * y0 + 7.0_dp * dx1 - 1.0_dp * dx2
    ff =   6.0_dp * y0 - 3.0_dp * dx1 + 0.5_dp * dx2
    xr = xx * invdx
    yyp = (5.0_dp * ff * xr**4 + 4.0_dp * ee * xr**3 + 3.0_dp * dd * xr**2) * invdx

  end function poly5ToZeroDeriv_

end module dftbp_dftb_slakoeqgrid
