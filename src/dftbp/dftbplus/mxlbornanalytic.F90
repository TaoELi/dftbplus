!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#:if WITH_SOCKETS

!> Analytic Slater-Koster matrix derivatives for Born effective-charge response.
module dftbp_dftbplus_mxlbornanalytic
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_schedule, only : assembleChunks, distributeRangeInChunks
  use dftbp_dftb_slakocont, only : getMIntegrals, getSKIntegralDerivs, TSlakoCont
  use dftbp_io_message, only : error
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: getFirstDerivAnalytic


  !> First-order automatic-differentiation scalar with derivatives wrt x, y and z.
  type :: TDual3
    real(dp) :: val = 0.0_dp
    real(dp) :: der(3) = 0.0_dp
  end type TDual3


  interface assignment(=)
    module procedure assignRealToDual
  end interface assignment(=)

  interface operator(+)
    module procedure addDualDual, addDualReal, addRealDual, addDualInteger, addIntegerDual
  end interface operator(+)

  interface operator(-)
    module procedure subDualDual, subDualReal, subRealDual, subDualInteger, subIntegerDual, negDual
  end interface operator(-)

  interface operator(*)
    module procedure mulDualDual, mulDualReal, mulRealDual, mulDualInteger, mulIntegerDual
  end interface operator(*)

  interface operator(/)
    module procedure divDualDual, divDualReal, divRealDual, divDualInteger, divIntegerDual
  end interface operator(/)

  interface operator(**)
    module procedure powDualInteger, powDualReal
  end interface operator(**)

  interface sqrt
    module procedure sqrtDual
  end interface sqrt

  interface getFirstDerivAnalytic
    module procedure getFirstDerivAnalyticMatrix
    module procedure getFirstDerivAnalyticBlock
  end interface getFirstDerivAnalytic


  !> Maximal angular momentum, for which rotations are present.
  integer, parameter :: mAngRot_ = 3

contains


  !> Assign a real scalar to a dual scalar with zero derivative.
  elemental subroutine assignRealToDual(lhs, rhs)

    type(TDual3), intent(out) :: lhs
    real(dp), intent(in) :: rhs

    lhs%val = rhs
    lhs%der(:) = 0.0_dp

  end subroutine assignRealToDual


  function makeDual(val, der) result(out)

    real(dp), intent(in) :: val
    real(dp), intent(in) :: der(3)
    type(TDual3) :: out

    out%val = val
    out%der(:) = der(:)

  end function makeDual


  elemental function addDualDual(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs, rhs
    type(TDual3) :: out
    out%val = lhs%val + rhs%val
    out%der(:) = lhs%der(:) + rhs%der(:)
  end function addDualDual


  elemental function addDualReal(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    real(dp), intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs%val + rhs
    out%der(:) = lhs%der(:)
  end function addDualReal


  elemental function addRealDual(lhs, rhs) result(out)
    real(dp), intent(in) :: lhs
    type(TDual3), intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs + rhs%val
    out%der(:) = rhs%der(:)
  end function addRealDual


  elemental function addDualInteger(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    integer, intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs%val + real(rhs, dp)
    out%der(:) = lhs%der(:)
  end function addDualInteger


  elemental function addIntegerDual(lhs, rhs) result(out)
    integer, intent(in) :: lhs
    type(TDual3), intent(in) :: rhs
    type(TDual3) :: out
    out%val = real(lhs, dp) + rhs%val
    out%der(:) = rhs%der(:)
  end function addIntegerDual


  elemental function subDualDual(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs, rhs
    type(TDual3) :: out
    out%val = lhs%val - rhs%val
    out%der(:) = lhs%der(:) - rhs%der(:)
  end function subDualDual


  elemental function subDualReal(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    real(dp), intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs%val - rhs
    out%der(:) = lhs%der(:)
  end function subDualReal


  elemental function subRealDual(lhs, rhs) result(out)
    real(dp), intent(in) :: lhs
    type(TDual3), intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs - rhs%val
    out%der(:) = -rhs%der(:)
  end function subRealDual


  elemental function subDualInteger(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    integer, intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs%val - real(rhs, dp)
    out%der(:) = lhs%der(:)
  end function subDualInteger


  elemental function subIntegerDual(lhs, rhs) result(out)
    integer, intent(in) :: lhs
    type(TDual3), intent(in) :: rhs
    type(TDual3) :: out
    out%val = real(lhs, dp) - rhs%val
    out%der(:) = -rhs%der(:)
  end function subIntegerDual


  elemental function negDual(arg) result(out)
    type(TDual3), intent(in) :: arg
    type(TDual3) :: out
    out%val = -arg%val
    out%der(:) = -arg%der(:)
  end function negDual


  elemental function mulDualDual(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs, rhs
    type(TDual3) :: out
    out%val = lhs%val * rhs%val
    out%der(:) = lhs%der(:) * rhs%val + rhs%der(:) * lhs%val
  end function mulDualDual


  elemental function mulDualReal(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    real(dp), intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs%val * rhs
    out%der(:) = lhs%der(:) * rhs
  end function mulDualReal


  elemental function mulRealDual(lhs, rhs) result(out)
    real(dp), intent(in) :: lhs
    type(TDual3), intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs * rhs%val
    out%der(:) = lhs * rhs%der(:)
  end function mulRealDual


  elemental function mulDualInteger(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    integer, intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs%val * real(rhs, dp)
    out%der(:) = lhs%der(:) * real(rhs, dp)
  end function mulDualInteger


  elemental function mulIntegerDual(lhs, rhs) result(out)
    integer, intent(in) :: lhs
    type(TDual3), intent(in) :: rhs
    type(TDual3) :: out
    out%val = real(lhs, dp) * rhs%val
    out%der(:) = real(lhs, dp) * rhs%der(:)
  end function mulIntegerDual


  elemental function divDualDual(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs, rhs
    type(TDual3) :: out
    out%val = lhs%val / rhs%val
    out%der(:) = (lhs%der(:) * rhs%val - lhs%val * rhs%der(:)) / rhs%val**2
  end function divDualDual


  elemental function divDualReal(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    real(dp), intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs%val / rhs
    out%der(:) = lhs%der(:) / rhs
  end function divDualReal


  elemental function divRealDual(lhs, rhs) result(out)
    real(dp), intent(in) :: lhs
    type(TDual3), intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs / rhs%val
    out%der(:) = -lhs * rhs%der(:) / rhs%val**2
  end function divRealDual


  elemental function divDualInteger(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    integer, intent(in) :: rhs
    type(TDual3) :: out
    out%val = lhs%val / real(rhs, dp)
    out%der(:) = lhs%der(:) / real(rhs, dp)
  end function divDualInteger


  elemental function divIntegerDual(lhs, rhs) result(out)
    integer, intent(in) :: lhs
    type(TDual3), intent(in) :: rhs
    type(TDual3) :: out
    out%val = real(lhs, dp) / rhs%val
    out%der(:) = -real(lhs, dp) * rhs%der(:) / rhs%val**2
  end function divIntegerDual


  elemental function powDualInteger(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    integer, intent(in) :: rhs
    type(TDual3) :: out
    integer :: ii

    if (rhs == 0) then
      out = 1.0_dp
    elseif (rhs < 0) then
      out = 1.0_dp / powDualInteger(lhs, -rhs)
    else
      out = lhs
      do ii = 2, rhs
        out = out * lhs
      end do
    end if

  end function powDualInteger


  elemental function powDualReal(lhs, rhs) result(out)
    type(TDual3), intent(in) :: lhs
    real(dp), intent(in) :: rhs
    type(TDual3) :: out
    integer :: iPower

    iPower = nint(rhs)
    if (abs(rhs - real(iPower, dp)) <= 10.0_dp * epsilon(1.0_dp)) then
      out = lhs ** iPower
    else
      out%val = lhs%val ** rhs
      out%der(:) = rhs * lhs%val ** (rhs - 1.0_dp) * lhs%der(:)
    end if

  end function powDualReal


  elemental function sqrtDual(arg) result(out)
    type(TDual3), intent(in) :: arg
    type(TDual3) :: out
    real(dp) :: root

    root = sqrt(arg%val)
    out%val = root
    out%der(:) = 0.5_dp * arg%der(:) / root

  end function sqrtDual


  !> Calculates the analytic first derivative of all of H0 or S wrt iAt position.
  subroutine getFirstDerivAnalyticMatrix(deriv, env, skCont, coords, species, iAt, orb,&
      & nNeighbourSK, iNeighbours, iPair)

    !> Derivative of H0 or S matrix, with respect to x,y,z.
    real(dp), intent(out) :: deriv(:,:)

    !> Computational environment settings.
    type(TEnvironment), intent(in) :: env

    !> Container for the SK integrals.
    type(TSlakoCont), intent(in) :: skCont

    !> List of all coordinates, including possible periodic images of atoms.
    real(dp), intent(in) :: coords(:,:)

    !> Chemical species of each atom.
    integer, intent(in) :: species(:)

    !> Atom to differentiate wrt.
    integer, intent(in) :: iAt

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours for each central cell atom.
    integer, intent(in) :: nNeighbourSK(:)

    !> List of neighbours for each central cell.
    integer, intent(in) :: iNeighbours(0:,:)

    !> Starting position of atom-neighbor interaction in the sparse matrix.
    integer, intent(in) :: iPair(0:,:)

    real(dp) :: blockDeriv(orb%mOrb, orb%mOrb, 3)
    real(dp) :: sign
    integer :: iAt1, iAt2, iAtFirst, iAtLast, iNeigh1, iSp1, iSp2, ind, nOrb1, nOrb2
    integer :: iCart

    deriv(:,:) = 0.0_dp

    call distributeRangeInChunks(env, 1, size(nNeighbourSK), iAtFirst, iAtLast)

    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh1 = 1, nNeighbourSK(iAt1)
        iAt2 = iNeighbours(iNeigh1, iAt1)
        if (iAt /= iAt1 .and. iAt /= iAt2) cycle
        iSp2 = species(iAt2)
        nOrb2 = orb%nOrbSpecies(iSp2)
        ind = iPair(iNeigh1, iAt1)
        sign = 1.0_dp
        if (iAt == iAt1) sign = -1.0_dp
        call getPairDerivative_(blockDeriv, skCont, coords, species, iAt1, iAt2, orb, sign)
        do iCart = 1, 3
          deriv(ind + 1 : ind + nOrb2 * nOrb1, iCart) =&
              & reshape(blockDeriv(1:nOrb2, 1:nOrb1, iCart), [nOrb2 * nOrb1])
        end do
      end do
    end do

    call assembleChunks(env, deriv)

  end subroutine getFirstDerivAnalyticMatrix


  !> Calculates the analytic first derivative of a diatomic H0 or S block wrt atomI.
  subroutine getFirstDerivAnalyticBlock(deriv, skCont, coords, species, atomI, atomJ, orb)

    !> Derivative of matrix block, with respect to x,y,z.
    real(dp), intent(out) :: deriv(:,:,:)

    !> Container for SK integrals.
    type(TSlakoCont), intent(in) :: skCont

    !> List of all coordinates.
    real(dp), intent(in) :: coords(:,:)

    !> Chemical species of each atom.
    integer, intent(in) :: species(:)

    !> The first atom in the diatomic block.
    integer, intent(in) :: atomI

    !> The second atom in the diatomic block.
    integer, intent(in) :: atomJ

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    call getPairDerivative_(deriv, skCont, coords, species, atomI, atomJ, orb, -1.0_dp)

  end subroutine getFirstDerivAnalyticBlock


  !> Analytic derivative of the rotated two-centre block wrt sign * d/d(R_J - R_I).
  subroutine getPairDerivative_(deriv, skCont, coords, species, atomI, atomJ, orb, sign)

    real(dp), intent(out) :: deriv(:,:,:)
    type(TSlakoCont), intent(in) :: skCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: atomI, atomJ
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: sign

    type(TDual3) :: vect(3), dist, ll, mm, nn
    type(TDual3) :: skDual(getMIntegrals(skCont))
    type(TDual3) :: blockDual(orb%mOrb, orb%mOrb)
    real(dp) :: vectVal(3), skVal(getMIntegrals(skCont)), dSkDr(getMIntegrals(skCont))
    real(dp) :: zeroDer(3)
    integer :: iCart, iOrb, jOrb, iInt, iSp1, iSp2

    deriv(:,:,:) = 0.0_dp
    zeroDer(:) = 0.0_dp

    vectVal(:) = coords(:, atomJ) - coords(:, atomI)
    do iCart = 1, 3
      zeroDer(:) = 0.0_dp
      zeroDer(iCart) = 1.0_dp
      vect(iCart) = makeDual(vectVal(iCart), zeroDer)
    end do

    dist = sqrt(vect(1) * vect(1) + vect(2) * vect(2) + vect(3) * vect(3))
    if (dist%val <= epsilon(1.0_dp)) then
      call error("Analytic Slater-Koster derivative is undefined for zero-distance atom pairs")
    end if

    ll = vect(1) / dist
    mm = vect(2) / dist
    nn = vect(3) / dist

    iSp1 = species(atomI)
    iSp2 = species(atomJ)
    skVal(:) = 0.0_dp
    dSkDr(:) = 0.0_dp
    call getSKIntegralDerivs(skCont, skVal, dSkDr, dist%val, iSp1, iSp2)
    do iInt = 1, size(skDual)
      skDual(iInt) = makeDual(skVal(iInt), dSkDr(iInt) * dist%der(:))
    end do

    call rotateH0Dual(blockDual, skDual, ll, mm, nn, iSp1, iSp2, orb)

    do iCart = 1, 3
      do iOrb = 1, size(deriv, dim=2)
        do jOrb = 1, size(deriv, dim=1)
          deriv(jOrb, iOrb, iCart) = sign * blockDual(jOrb, iOrb)%der(iCart)
        end do
      end do
    end do

  end subroutine getPairDerivative_


  !> Driver for making the non-SCC hhamiltonian or overlap matrices for a given diatomic block
  !> Caveat: Only angular momenta up to f are currently allowed
  subroutine rotateH0Dual(hh, skIntegs, ll, mm, nn, iSp1, iSp2, orb)

    !> The rectangular matrix containing the resulting diatomic matrix elements
    type(TDual3), intent(out) :: hh(:,:)

    !> Slater-Koster table for dimer of species i-j
    type(TDual3), intent(in), target :: skIntegs(:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Chemical species of atom i
    integer, intent(in) :: iSp1

    !> Chemical species of atom j
    integer, intent(in) :: iSp2

    !> Information about the orbitals of chemical species in the system.
    type(TOrbitals), intent(in) :: orb

    integer :: iCol, iRow, ind, iSh1, iSh2
    integer :: ang1, ang2, nOrb1, nOrb2
    type(TDual3), pointer :: pSK(:)
    type(TDual3) :: tmpH(2*mAngRot_+1,2*mAngRot_+1)

    @:ASSERT(maxval(orb%angShell) <=  mAngRot_)
    @:ASSERT(all(shape(hh) >= (/ orb%nOrbSpecies(iSp1), orb%nOrbSpecies(iSp2) /)))

    hh(:,:) = 0.0_dp
    ind = 1
    iCol = 1
    do iSh1 = 1, orb%nShell(iSp1)
      ang1 = orb%angShell(iSh1, iSp1)
      nOrb1 = 2 * ang1 + 1
      iRow = 1
      do iSh2 = 1, orb%nShell(iSp2)
        ang2 = orb%angShell(iSh2, iSp2)
        nOrb2 = 2 * ang2 + 1
        @:ASSERT(size(skIntegs) >= ind + min(ang1,ang2))
        pSK => skIntegs(ind:ind+min(ang1,ang2))
        select case (ang1)
        case (0)
          select case (ang2)
          case (0)
            call ss(tmpH,pSK)
          case (1)
            call sp(tmpH,ll,mm,nn,pSK)
          case (2)
            call sd(tmpH,ll,mm,nn,pSK)
          case (3)
            call sf(tmpH,ll,mm,nn,pSK)
          end select
        case (1)
          select case (ang2)
          case (0)
            call sp(tmpH,ll,mm,nn,pSK)
          case (1)
            call pp(tmpH,ll,mm,nn,pSK)
          case (2)
            call pd(tmpH,ll,mm,nn,pSK)
          case (3)
            call pf(tmpH,ll,mm,nn,pSK)
          end select
        case (2)
          select case (ang2)
          case(0)
            call sd(tmpH,ll,mm,nn,pSK)
          case(1)
            call pd(tmpH,ll,mm,nn,pSK)
          case (2)
            call dd(tmpH,ll,mm,nn,pSK)
          case (3)
            call df(tmpH,ll,mm,nn,pSK)
          end select
        case (3)
          select case (ang2)
          case(0)
            call sf(tmpH,ll,mm,nn,pSK)
          case(1)
            call pf(tmpH,ll,mm,nn,pSK)
          case(2)
            call df(tmpH,ll,mm,nn,pSK)
          case(3)
            call ff(tmpH,ll,mm,nn,pSK)
          end select
        end select

        if (ang1 <= ang2) then
          hh(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1) = tmpH(1:nOrb2,1:nOrb1)
        else
          hh(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1) = (-1.0_dp)**(ang1+ang2) &
              &* transpose(tmpH(1:nOrb1,1:nOrb2))
        end if
        ind = ind + min(ang1,ang2) + 1
        iRow = iRow + nOrb2
      end do
      iCol = iCol + nOrb1
    end do

  end subroutine rotateH0Dual


  !> Rotation routine for interaction of an s orbital with an s orbital
  subroutine ss(hh, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 1)
    @:ASSERT(all(shape(hh) >= (/ 1, 1 /)))

    hh(1,1) = sk(1)

  end subroutine ss


  !> Rotation routine for interaction of an s orbital with a p orbital
  subroutine sp(hh, ll, mm, nn, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 1)
    @:ASSERT(all(shape(hh) >= (/ 3, 1 /)))

    hh(1,1) = mm*sk(1)
    hh(2,1) = nn*sk(1)
    hh(3,1) = ll*sk(1)

  end subroutine sp


  !> Rotation routine for interaction of an s orbital with a d orbital
  subroutine sd(hh, ll, mm, nn, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 1)
    @:ASSERT(all(shape(hh) >= (/ 5, 1 /)))

    hh(1,1) = ll*mm*sqrt(3.0_dp)*sk(1)
    hh(2,1) = mm*sqrt(3.0_dp)*nn*sk(1)
    hh(3,1) = (3.0_dp/ 2.0_dp*nn**2- 1.0_dp/ 2.0_dp)*sk(1)
    hh(4,1) = ll*sqrt(3.0_dp)*nn*sk(1)
    hh(5,1) = (2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*sk(1)/2.0_dp

  end subroutine sd


  !> Rotation routine for interaction of an s orbital with an f orbital
  subroutine sf(hh, ll, mm, nn, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 1)
    @:ASSERT(all(shape(hh) >= (/ 1, 7 /)))

    hh(1,1) = sqrt(2.0_dp)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)&
        &*sk(1)/ 4.0_dp
    hh(2,1) = ll*mm*sqrt(15.0_dp)*nn*sk(1)
    hh(3,1) = sqrt(2.0_dp)*mm*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)/ 4.0_dp
    hh(4,1) = (nn*(5.0_dp*nn**2-3.0_dp)*sk(1))/ 2.0_dp
    hh(5,1) = sqrt(2.0_dp)*ll*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)/4.0_dp
    hh(6,1) = (2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(15.0_dp)*nn*sk(1)/ 2.0_dp
    hh(7,1) = sqrt(2.0_dp)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp

  end subroutine sf


  !> Rotation routine for interaction of a p orbital with a p orbital
  subroutine pp(hh, ll, mm, nn, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 2)
    @:ASSERT(all(shape(hh) >= (/ 3, 3 /)))

    hh(1,1) = (1.0_dp-nn**2-ll**2)*sk(1)+(nn**2+ll**2)*sk(2)
    hh(2,1) = nn*mm*sk(1)-nn*mm*sk(2)
    hh(3,1) = ll*mm*sk(1)-ll*mm*sk(2)
    hh(1,2) = hh(2,1)
    hh(2,2) = nn**2*sk(1)+(1.0_dp-nn**2)*sk(2)
    hh(3,2) = nn*ll*sk(1)-nn*ll*sk(2)
    hh(1,3) = hh(3,1)
    hh(2,3) = hh(3,2)
    hh(3,3) = ll**2*sk(1)+(1.0_dp-ll**2)*sk(2)

  end subroutine pp


  !> Rotation routine for interaction of a p orbital with a d orbital
  subroutine pd(hh, ll, mm, nn, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 2)
    @:ASSERT(all(shape(hh) >= (/ 3, 5 /)))

    hh(1,1) = -(-1.0_dp+nn**2+ll**2)*ll*sqrt(3.0_dp)&
        &*sk(1)+((2.0_dp*nn**2+2.0_dp*ll**2-1.0_dp)*ll*sk(2))
    hh(2,1) = -(-1.0_dp+nn**2+ll**2)*sqrt(3.0_dp)*nn*&
        &sk(1)+((2.0_dp*nn**2+2.0_dp*ll**2-1.0_dp)*nn*sk(2))
    hh(3,1) = mm*(3.0_dp*nn**2-1.0_dp)*sk(1)/2.0_dp&
        &-sqrt(3.0_dp)*(nn**2)*mm*sk(2)
    hh(4,1) = mm*ll*sqrt(3.0_dp)*nn*sk(1)-2.0_dp*ll*mm*nn*sk(2)
    hh(5,1) = mm*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*sk(1)/2.0_dp&
        &-(nn**2+2.0_dp*ll**2)*mm*sk(2)
    hh(1,2) = ll*mm*nn*sqrt(3.0_dp)*sk(1)-2.0_dp*nn*ll*mm*sk(2)
    hh(2,2) = mm*(nn**2)*sqrt(3.0_dp)*sk(1)&
        &-(2.0_dp*nn**2-1.0_dp)*mm*sk(2)
    hh(3,2) = (nn*(3.0_dp*nn**2-1.0_dp)*sk(1))/2.0_dp&
        &-nn*sqrt(3.0_dp)*(-1.0_dp+nn**2)*sk(2)
    hh(4,2) = ll*nn**2*sqrt(3.0_dp)*sk(1)-(2.0_dp*nn**2-1.0_dp)*ll*&
        &sk(2)
    hh(5,2) = (2.0_dp*ll**2-1.0_dp+nn**2)*nn*sqrt(3.0_dp)*sk(1)/2.0_dp&
        &-(nn*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(2))
    hh(1,3) = (ll**2)*mm*sqrt(3.0_dp)*sk(1)&
        &-(2.0_dp*ll**2-1.0_dp)*mm*sk(2)
    hh(2,3) = ll*mm*sqrt(3.0_dp)*nn*sk(1)-2.0_dp*mm*ll*nn*sk(2)
    hh(3,3) = (ll*(3.0_dp*nn**2-1.0_dp)*sk(1))/2.0_dp&
        &-sqrt(3.0_dp)*(nn**2)*ll*sk(2)
    hh(4,3) = ll**2*sqrt(3.0_dp)*nn*sk(1)-(2.0_dp*ll**2-1.0_dp)*nn&
        &*sk(2)
    hh(5,3) = ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*sk(1)/2.0_dp&
        &-((nn**2-2.0_dp+2.0_dp*ll**2)*ll*sk(2))

  end subroutine pd


  !> Rotation routine for interaction of a p orbital with an f orbital
  subroutine pf(hh, ll, mm, nn, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 2)
    @:ASSERT(all(shape(hh) >= (/ 3, 7 /)))

    hh(1,1) = -(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(2.0_dp)&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp+sqrt(15.0_dp)&
        &*(nn**4-nn**2+5.0_dp*nn**2*ll**2-3.0_dp*ll**2+4.0_dp*ll**4)&
        &*sk(2)/4.0_dp
    hh(2,1) = -(-1.0_dp+nn**2+ll**2)*ll*sqrt(15.0_dp)*nn*sk(1)&
        &+(3.0_dp*nn**2+3.0_dp*ll**2-2.0_dp)*ll*nn*sqrt(10.0_dp)*sk(2)&
        &/2.0_dp
    hh(3,1) = -(-1.0_dp+nn**2+ll**2)*sqrt(2.0_dp)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)/4.0_dp+(15.0_dp/4.0_dp&
        &*(nn**4)+ 15.0_dp/ 4.0_dp*(ll**2)*(nn**2)- 11.0_dp&
        &/4.0_dp*(nn**2)-(ll**2)/ 4.0_dp)*sk(2)
    hh(4,1) = mm*nn*(5.0_dp*nn**2-3.0_dp)*sk(1)/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*nn*mm&
        &*sk(2)/4.0_dp
    hh(5,1) = mm*ll*sqrt(2.0_dp)*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*sk(1)&
        &/4.0_dp-(15.0_dp*nn**2-1.0_dp)*ll*mm*sk(2)/4.0_dp
    hh(6,1) = mm*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(15.0_dp)*nn*&
        &sk(1)/2.0_dp-(3.0_dp*nn**2+6.0_dp*ll**2-1.0_dp)*nn*mm&
        &*sqrt(10.0_dp)*sk(2)/4.0_dp
    hh(7,1) = mm*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(2.0_dp)&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-ll*mm*sqrt(15.0_dp)&
        &*(3.0_dp*nn**2+4.0_dp*ll**2-1.0_dp)*sk(2)/4.0_dp
    hh(1,2) = sqrt(2.0_dp)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(15.0_dp)*nn*sk(2)/4.0_dp
    hh(2,2) = ll*mm*(nn**2)*sqrt(15.0_dp)*sk(1)-(3.0_dp*nn**2-1.0_dp)*ll*mm&
        &*sqrt(10.0_dp)*sk(2)/2.0_dp
    hh(3,2) = sqrt(2.0_dp)*mm*nn*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*sk(1)&
        &/4.0_dp-(15.0_dp*nn**2-11.0_dp)*nn*mm*sk(2)/4.0_dp
    hh(4,2) = (nn**2*(5.0_dp*nn**2-3.0_dp)*sk(1))/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*(-1.0_dp+nn**2)&
        &*sk(2)/4.0_dp
    hh(5,2) = sqrt(2.0_dp)*ll*nn*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*sk(1)&
        &/4.0_dp-(15.0_dp*nn**2-11.0_dp)*nn*ll*sk(2)/4.0_dp
    hh(6,2) = (2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)*sqrt(15.0_dp)*sk(1)&
        &/2.0_dp-(3.0_dp*nn**2-1.0_dp)*(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(10.0_dp)*sk(2)/4.0_dp
    hh(7,2) = sqrt(2.0_dp)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*sqrt(15.0_dp)*nn*sk(2)/4.0_dp
    hh(1,3) = ll*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(2.0_dp)*sqrt(5.0_dp)&
        &*sk(1)/4.0_dp-ll*mm*sqrt(15.0_dp)*(nn**2+4.0_dp*ll**2-3.0_dp)&
        &*sk(2)/4.0_dp
    hh(2,3) = (ll**2)*mm*sqrt(15.0_dp)*nn*sk(1)-(3.0_dp*ll**2-1.0_dp)*nn*mm&
        &*sqrt(10.0_dp)*sk(2)/2.0_dp
    hh(3,3) = ll*mm*sqrt(2.0_dp)*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)/4.0_dp-(15.0_dp*nn**2-1.0_dp)*ll*mm*sk(2)/4.0_dp
    hh(4,3) = (ll*nn*(5.0_dp*nn**2-3.0_dp)*sk(1))/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*nn*ll*sk(2)&
        &/4.0_dp
    hh(5,3) = ll**2*sqrt(2.0_dp)*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)/4.0_dp+(-15.0_dp/4.0_dp*ll**2*(nn**2)+5.0_dp&
        &/4.0_dp*(nn**2)- 1.0_dp/4.0_dp+ll**2/4.0_dp)*sk(2)
    hh(6,3) = ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(15.0_dp)*nn*sk(1)&
        &/2.0_dp-(3.0_dp*nn**2+6.0_dp*ll**2-5.0_dp)*nn*ll*sqrt(10.0_dp)&
        &*sk(2)/4.0_dp
    hh(7,3) = (ll**2)*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(2.0_dp)&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-sqrt(15.0_dp)*(3.0_dp*ll**2*&
        &nn**2-nn**2-5.0_dp*ll**2+4.0_dp*ll**4+1.0_dp)*sk(2)/4.0_dp

  end subroutine pf


  !> Rotation routine for interaction of a d orbital with a d orbital
  subroutine dd(hh, ll, mm, nn, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 3)
    @:ASSERT(all(shape(hh) >= (/ 5, 5 /)))

    hh(1,1) = -3.0_dp*ll**2*(-1.0_dp+nn**2+ll**2)*sk(1)&
        &+(4.0_dp*ll**2*nn**2-nn**2+4.0_dp*ll**4-4.0_dp*ll**2+1.0_dp)*sk(2)&
        &+(-ll**2*nn**2+nn**2+ll**2-ll**4)*sk(3)
    hh(2,1) = -3.0_dp*ll*(-1.0_dp+nn**2+ll**2)*nn*sk(1)&
        &+(4.0_dp*nn**2+4.0_dp*ll**2-3.0_dp)*nn*ll*sk(2)&
        &-ll*(nn**2+ll**2)*nn*sk(3)
    hh(3,1) = ll*mm*sqrt(3.0_dp)*(3.0_dp*nn**2-1.0_dp)*sk(1)/2.0_dp&
        &-2.0_dp*sqrt(3.0_dp)*mm*ll*(nn**2)*sk(2)&
        &+ll*mm*(nn**2+1.0_dp)*sqrt(3.0_dp)*sk(3)/2.0_dp
    hh(4,1) = 3.0_dp*(ll**2)*mm*nn*sk(1)-(4.0_dp*ll**2-1.0_dp)*nn*mm&
        &*sk(2)+mm*(-1.0_dp+ll**2)*nn*sk(3)
    hh(5,1) = 3.0_dp/2.0_dp*mm*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(1)&
        &-2.0_dp*mm*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(2)&
        &+mm*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(3)/2.0_dp
    hh(1,2) = hh(2,1)
    hh(2,2) = -3.0_dp*(-1.0_dp+nn**2+ll**2)*nn**2*sk(1)&
        &+(4.0_dp*nn**4-4.0_dp*nn**2+4.0_dp*ll**2*nn**2+1.0_dp-ll**2)*sk(2)&
        &-(-1.0_dp+nn)*(nn**3+nn**2+ll**2*nn+ll**2)*sk(3)
    hh(3,2) = mm*sqrt(3.0_dp)*nn*(3.0_dp*nn**2-1.0_dp)*sk(1)/2.0_dp&
        &-nn*sqrt(3.0_dp)*(2.0_dp*nn**2-1.0_dp)*mm*sk(2)&
        &+(-1.0_dp+nn**2)*sqrt(3.0_dp)*nn*mm*sk(3)/2.0_dp
    hh(4,2) = 3.0_dp*mm*ll*(nn**2)*sk(1)-(4.0_dp*nn**2-1.0_dp)*mm*ll&
        &*sk(2)+ll*mm*(-1.0_dp+nn**2)*sk(3)
    hh(5,2) = 3.0_dp/ 2.0_dp*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*nn*sk(1)&
        &-(2.0_dp*nn**2-1.0_dp+4.0_dp*ll**2)*nn*mm*sk(2)&
        &+mm*(nn**2+2.0_dp*ll**2+1.0_dp)*nn*sk(3)/2.0_dp
    hh(1,3) = hh(3,1)
    hh(2,3) = hh(3,2)
    hh(3,3) = ((3.0_dp*nn**2-1.0_dp)**2*sk(1))/4.0_dp&
        &-(3.0_dp*(-1.0_dp+nn**2)*nn**2*sk(2))&
        &+3.0_dp/4.0_dp*((-1.0_dp+nn**2)**2)*sk(3)
    hh(4,3) = ll*(3.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*nn*sk(1)/2.0_dp&
        &-(2.0_dp*nn**2-1.0_dp)*ll*nn*sqrt(3.0_dp)*sk(2)&
        &+nn*ll*sqrt(3.0_dp)*(-1.0_dp+nn**2)*sk(3)/2.0_dp
    hh(5,3) = (2.0_dp*ll**2-1.0_dp+nn**2)*(3.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)&
        &*sk(1)/4.0_dp-(2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)&
        &*sqrt(3.0_dp)*sk(2)+sqrt(3.0_dp)*(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*(nn**2+1.0_dp)*sk(3)/4.0_dp
    hh(1,4) = hh(4,1)
    hh(2,4) = hh(4,2)
    hh(3,4) = hh(4,3)
    hh(4,4) = 3.0_dp*ll**2*nn**2*sk(1)+(-4.0_dp*ll**2*nn**2+nn**2+ll**2)&
        &*sk(2)+(-1.0_dp+nn)*(-nn+ll**2*nn-1.0_dp+ll**2)*sk(3)
    hh(5,4) = 3.0_dp/2.0_dp*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*nn*sk(1)&
        &-((2.0_dp*nn**2-3.0_dp+4.0_dp*ll**2)*nn*ll*sk(2))&
        &+(ll*(nn**2-3.0_dp+2.0_dp*ll**2)*nn*sk(3))/2.0_dp
    hh(1,5) = hh(5,1)
    hh(2,5) = hh(5,2)
    hh(3,5) = hh(5,3)
    hh(4,5) = hh(5,4)
    hh(5,5) = 3.0_dp/4.0_dp*((2.0_dp*ll**2-1.0_dp+nn**2)**2)*sk(1)&
        &+((-nn**4+nn**2-4.0_dp*ll**2*nn**2-4.0_dp*ll**4&
        &+4.0_dp*ll**2)*sk(2))+((nn**4)/4.0_dp&
        &+(ll**2*nn**2)+(nn**2)/2.0_dp+1.0_dp/4.0_dp-(ll**2)+(ll**4))*sk(3)

  end subroutine dd


  !> Rotation routine for interaction of a d orbital with an f orbital
  subroutine df(hh, ll, mm, nn, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 3)
    @:ASSERT(all(shape(hh) >= (/ 5, 7 /)))

    hh(1,1) = -ll*(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(6.0_dp)*sqrt(5.0_dp)*sk(1)/4.0_dp&
        &+sqrt(15.0_dp)*ll*(2.0_dp*nn**4-5.0_dp*nn**2+10.0_dp*ll**2*nn**2&
        &+3.0_dp-10.0_dp*ll**2+8.0_dp*ll**4)*sk(2)/4.0_dp&
        &-sqrt(6.0_dp)*ll*(nn**4+5.0_dp*ll**2*nn**2-4.0_dp*nn**2+1.0_dp&
        &+4.0_dp*ll**4-5.0_dp*ll**2)*sk(3)/4.0_dp
    hh(2,1) = -3.0_dp*(ll**2)*(-1.0_dp+nn**2+ll**2)*sqrt(5.0_dp)*nn*sk(1)&
        &+(6.0_dp*ll**2*nn**2-nn**2+1.0_dp+6.0_dp*ll**4-6.0_dp*ll**2)&
        &*sqrt(10.0_dp)*nn*sk(2)/2.0_dp-(nn*(3.0_dp*ll**2*nn**2-2.0_dp*nn**2&
        &+1.0_dp-3.0_dp*ll**2+3.0_dp*ll**4)*sk(3))
    hh(3,1) = -3.0_dp/4.0_dp*ll*(-1.0_dp+nn**2+ll**2)*sqrt(2.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)+((30.0_dp*nn**4+30.0_dp*ll**2*nn**2&
        &-27.0_dp*nn**2-2.0_dp*ll**2+1.0_dp)*ll*sk(2))/4.0_dp&
        &-ll*sqrt(10.0_dp)*(3.0_dp*nn**4+3.0_dp*ll**2*nn**2+ll**2-1.0_dp)&
        &*sk(3)/4.0_dp
    hh(4,1) = ll*mm*sqrt(3.0_dp)*nn*(5.0_dp*nn**2-3.0_dp)*sk(1)/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*nn*ll*mm*sqrt(6.0_dp)*sk(2)/2.0_dp&
        &+ll*mm*(nn**2+1.0_dp)*sqrt(15.0_dp)*nn*sk(3)/2.0_dp
    hh(5,1) = 3.0_dp/4.0_dp*(ll**2)*mm*sqrt(2.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-(30.0_dp*ll**2*nn**2-5.0_dp*nn**2-2.0_dp*ll**2+1.0_dp)*mm&
        &*sk(2)/4.0_dp+mm*sqrt(10.0_dp)*(3.0_dp*ll**2*nn**2-2.0_dp*nn**2&
        &+ll**2)*sk(3)/4.0_dp
    hh(6,1) = 3.0_dp/2.0_dp*ll*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)*nn&
        &*sk(1)-3.0_dp/2.0_dp*nn*sqrt(10.0_dp)*mm*ll&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(2)+3.0_dp/2.0_dp*ll*mm&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*nn*sk(3)
    hh(7,1) = (ll**2)*mm*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(6.0_dp)&
        &*sqrt(5.0_dp)*sk(1)/ 4.0_dp-sqrt(15.0_dp)*mm*&
        &(6.0_dp*ll**2*nn**2-nn**2+1.0_dp+8.0_dp*ll**4-6.0_dp*ll**2)&
        &*sk(2)/4.0_dp+sqrt(6.0_dp)*mm*(3.0_dp*ll**2*nn**2-2.0_dp*nn**2&
        &+4.0_dp*ll**4-3.0_dp*ll**2)*sk(3)/4.0_dp
    hh(1,2) = -(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(6.0_dp)*nn*sqrt(5.0_dp)*sk(1)/4.0_dp&
        &+sqrt(15.0_dp)*nn*(2.0_dp*nn**4-3.0_dp*nn**2+10.0_dp*ll**2*nn**2&
        &+1.0_dp+8.0_dp*ll**4-8.0_dp*ll**2)*sk(2)/4.0_dp&
        &-sqrt(6.0_dp)*nn*(nn**4+5.0_dp*ll**2*nn**2-1.0_dp+4.0_dp*ll**4-ll**2)&
        &*sk(3)/4.0_dp
    hh(2,2) = -3.0_dp*(-1.0_dp+nn**2+ll**2)*ll*(nn**2)*sqrt(5.0_dp)*sk(1)&
        &+(6.0_dp*nn**4-6.0_dp*nn**2+6.0_dp*ll**2*nn**2+1.0_dp-ll**2)*ll&
        &*sqrt(10.0_dp)*sk(2)/2.0_dp-(ll*(3.0_dp*nn**4+3.0_dp*ll**2*nn**2&
        &-3.0_dp*nn**2+1.0_dp-2.0_dp*ll**2)*sk(3))
    hh(3,2) = -3.0_dp/4.0_dp*(-1.0_dp+nn**2+ll**2)*sqrt(2.0_dp)*nn&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)+((30.0_dp*nn**4-37.0_dp*nn**2&
        &+30.0_dp*ll**2*nn**2+11.0_dp-12.0_dp*ll**2)*nn*sk(2))&
        &/4.0_dp-(-1.0_dp+nn)*sqrt(10.0_dp)*(3.0_dp*nn**3+3.0_dp*nn**2-nn&
        &+3.0_dp*ll**2*nn-1.0_dp+3.0_dp*ll**2)*nn*sk(3)/4.0_dp
    hh(4,2) = mm*sqrt(3.0_dp)*(nn**2)*(5.0_dp*nn**2-3.0_dp)*sk(1)&
        &/2.0_dp-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*&
        &(2.0_dp*nn**2-1.0_dp)*mm*sk(2)/4.0_dp+(-1.0_dp+nn**2)&
        &*sqrt(15.0_dp)*(nn**2)*mm*sk(3)/2.0_dp
    hh(5,2) = 3.0_dp/4.0_dp*mm*ll*sqrt(2.0_dp)*nn*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-3.0_dp/2.0_dp*(5.0_dp*nn**2-2.0_dp)*mm*ll*nn&
        &*sk(2)+3.0_dp/4.0_dp*nn*sqrt(10.0_dp)*ll*mm*(-1.0_dp+nn**2)&
        &*sk(3)
    hh(6,2) = 3.0_dp/2.0_dp*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)&
        &*sqrt(5.0_dp)*sk(1)-(6.0_dp*nn**4+12.0_dp*ll**2*nn**2-5.0_dp*nn**2&
        &+1.0_dp-2.0_dp*ll**2)*mm*sqrt(10.0_dp)*sk(2)/4.0_dp&
        &+mm*(3.0_dp*nn**4-nn**2+6.0_dp*ll**2*nn**2-4.0_dp*ll**2)*sk(3)&
        &/2.0_dp
    hh(7,2) = mm*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(6.0_dp)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-mm*ll*sqrt(15.0_dp)*nn&
        &*(3.0_dp*nn**2-2.0_dp+4.0_dp*ll**2)*sk(2)/2.0_dp&
        &+sqrt(6.0_dp)*mm*ll*nn*(3.0_dp*nn**2+1.0_dp+4.0_dp*ll**2)&
        &*sk(3)/4.0_dp
    hh(1,3) = sqrt(2.0_dp)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*(3.0_dp*nn**2-1.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp&
        &-3.0_dp/4.0_dp*(nn**2)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)&
        &*sk(2)+3.0_dp/8.0_dp*sqrt(2.0_dp)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*(nn**2+1.0_dp)*sk(3)
    hh(2,3) = ll*mm*(3.0_dp*nn**2-1.0_dp)*sqrt(15.0_dp)*nn*sk(1)/2.0_dp&
        &-(3.0_dp*nn**2-1.0_dp)*ll*nn*mm*sqrt(30.0_dp)*sk(2)/2.0_dp&
        &+sqrt(3.0_dp)*ll*mm*nn*(3.0_dp*nn**2-1.0_dp)*sk(3)/2.0_dp
    hh(3,3) = sqrt(2.0_dp)*mm*(3.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp-(15.0_dp*nn**2-11.0_dp)&
        &*mm*(nn**2)*sqrt(3.0_dp)*sk(2)/4.0_dp&
        &+(3.0_dp*nn**3+3.0_dp*nn**2-nn-1.0_dp)*(-1.0_dp+nn)*sqrt(5.0_dp)&
        &*sqrt(2.0_dp)*mm*sqrt(3.0_dp)*sk(3)/8.0_dp
    hh(4,3) = ((3.0_dp*nn**2-1.0_dp)*nn*(5.0_dp*nn**2-3.0_dp)*sk(1))&
        &/4.0_dp-3.0_dp/4.0_dp*(5.0_dp*nn**2-1.0_dp)*(-1.0_dp+nn**2)*nn&
        &*sqrt(2.0_dp)*sk(2)+3.0_dp/4.0_dp*((-1.0_dp+nn**2)**2)&
        &*sqrt(5.0_dp)*nn*sk(3)
    hh(5,3) = sqrt(2.0_dp)*ll*(3.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*&
        &(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp-(15.0_dp*nn**2-11.0_dp)&
        &*ll*(nn**2)*sqrt(3.0_dp)*sk(2)/4.0_dp+(3.0_dp*nn**3+3.0_dp*nn**2&
        &-nn-1.0_dp)*(-1.0_dp+nn)*sqrt(5.0_dp)*sqrt(2.0_dp)*ll*sqrt(3.0_dp)&
        &*sk(3)/8.0_dp
    hh(6,3) = (2.0_dp*ll**2-1.0_dp+nn**2)*(3.0_dp*nn**2-1.0_dp)*sqrt(15.0_dp)&
        &*nn*sk(1)/4.0_dp-(3.0_dp*nn**2-1.0_dp)*(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*nn*sqrt(30.0_dp)*sk(2)/4.0_dp+sqrt(3.0_dp)&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*nn*(3.0_dp*nn**2-1.0_dp)*sk(3)/&
        &4.0_dp
    hh(7,3) = sqrt(2.0_dp)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*(3.0_dp*nn**2-1.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp&
        &-3.0_dp/4.0_dp*(nn**2)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*sqrt(5.0_dp)*sk(2)+ 3.0_dp/8.0_dp*sqrt(2.0_dp)*ll&
        &*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*(nn**2+1.0_dp)*sk(3)
    hh(1,4) = ll*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(6.0_dp)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-ll*mm*sqrt(15.0_dp)*nn&
        &*(nn**2-2.0_dp+4.0_dp*ll**2)*sk(2)/2.0_dp&
        &+sqrt(6.0_dp)*mm*ll*nn*(nn**2+4.0_dp*ll**2-5.0_dp)*sk(3)/4.0_dp
    hh(2,4) = 3.0_dp*(ll**2)*mm*(nn**2)*sqrt(5.0_dp)*sk(1)&
        &-(6.0_dp*ll**2.0_dp*nn**2-nn**2-ll**2)*mm*sqrt(10.0_dp)*sk(2)&
        &/2.0_dp+mm*(-2.0_dp*nn**2+3.0_dp*ll**2.0_dp*nn**2+1.0_dp-2.0_dp*ll**2)&
        &*sk(3)
    hh(3,4) = 3.0_dp/4.0_dp*ll*mm*sqrt(2.0_dp)*nn*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-3.0_dp/2.0_dp*(5.0_dp*nn**2-2.0_dp)*ll*mm*nn*&
        &sk(2)+3.0_dp/4.0_dp*nn*sqrt(10.0_dp)*mm*ll*(-1.0_dp+nn**2)&
        &*sk(3)
    hh(4,4) = ll*sqrt(3.0_dp)*(nn**2)*(5.0_dp*nn**2-3.0_dp)*sk(1)/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*(2.0_dp*nn**2-1.0_dp)&
        &*ll*sk(2)/4.0_dp+(-1.0_dp+nn**2)*sqrt(15.0_dp)*(nn**2)*ll&
        &*sk(3)/2.0_dp
    hh(5,4) = 3.0_dp/ 4.0_dp*ll**2.0_dp*sqrt(2.0_dp)*nn*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-(30.0_dp*ll**2.0_dp*(nn**2)-(5.0_dp*nn**2)&
        &-12.0_dp*ll**2+1.0_dp)*nn*sk(2)/4.0_dp+(-1.0_dp+nn)&
        &*sqrt(10.0_dp)*(-(2.0_dp*nn)+3.0_dp*ll**2.0_dp*nn-2.0_dp+3.0_dp*ll**2)&
        &*nn*sk(3)/4.0_dp
    hh(6,4) = 3.0_dp/2.0_dp*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)&
        &*sqrt(5.0_dp)&
        &*sk(1)-(6.0_dp*nn**4+12.0_dp*ll**2.0_dp*nn**2-9.0_dp*nn**2&
        &+1.0_dp-2.0_dp*ll**2)*ll*sqrt(10.0_dp)*sk(2)/4.0_dp&
        &+(ll*(3.0_dp*nn**4-9.0_dp*nn**2+6.0_dp*ll**2.0_dp*nn**2+4.0_dp&
        &-4.0_dp*ll**2)*sk(3))/2.0_dp
    hh(7,4) = (ll**2)*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(6.0_dp)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-sqrt(15.0_dp)*nn&
        &*(6.0_dp*ll**2.0_dp*nn**2-nn**2+1.0_dp+8.0_dp*ll**4-8.0_dp*ll**2)&
        &*sk(2)/4.0_dp+sqrt(6.0_dp)*nn*(3.0_dp*ll**2.0_dp*nn**2&
        &-2.0_dp*nn**2-7.0_dp*ll**2+2.0_dp+4.0_dp*ll**4)*sk(3)/4.0_dp
    hh(1,5) = (2.0_dp*ll**2-1.0_dp+nn**2)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(6.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp-sqrt(15.0_dp)*mm&
        &*(nn**4-nn**2+6.0_dp*ll**2.0_dp*nn**2+8.0_dp*ll**4-6.0_dp*ll**2)&
        &*sk(2)/4.0_dp+sqrt(6.0_dp)*mm*(nn**4+6.0_dp*ll**2.0_dp*nn**2&
        &+2.0_dp*nn**2+1.0_dp+8.0_dp*ll**4-6.0_dp*ll**2)*sk(3)/8.0_dp
    hh(2,5) = 3.0_dp/2.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*mm*sqrt(5.0_dp)*nn&
        &*sk(1)-3.0_dp/2.0_dp*nn*sqrt(10.0_dp)*(2.0_dp*ll**2-1.0_dp&
        &+nn**2)*ll*mm*sk(2)+ 3.0_dp/2.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*mm&
        &*nn*sk(3)
    hh(3,5) = 3.0_dp/8.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*mm*sqrt(2.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)-(15.0_dp*nn**4+30.0_dp*ll**2*nn**2&
        &-11.0_dp*nn**2-2.0_dp*ll**2)*mm*sk(2)/4.0_dp&
        &+mm*sqrt(10.0_dp)*(3.0_dp*nn**4+2.0_dp*nn**2+6.0_dp*ll**2*nn**2-1.0_dp&
        &+2.0_dp*ll**2)*sk(3)/8.0_dp
    hh(4,5) = (2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*nn*(5.0_dp*nn**2-3.0_dp)&
        &*sk(1)/4.0_dp-(5.0_dp*nn**2-1.0_dp)*nn*(2.0_dp*ll**2-1.0_dp&
        &+nn**2)*sqrt(6.0_dp)*sk(2)/4.0_dp+(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*(nn**2+1.0_dp)*sqrt(15.0_dp)*nn*sk(3)/4.0_dp
    hh(5,5) = 3.0_dp/8.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*sqrt(2.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)-((15.0_dp*nn**4+30.0_dp*ll**2*nn**2&
        &-21.0_dp*nn**2+2.0_dp-2.0_dp*ll**2)*ll*sk(2))/4.0_dp&
        &+ll*sqrt(10.0_dp)*(3.0_dp*nn**4+6.0_dp*ll**2*nn**2-6.0_dp*nn**2&
        &+2.0_dp*ll**2-1.0_dp)*sk(3)/8.0_dp
    hh(6,5) = 3.0_dp/4.0_dp*((2.0_dp*ll**2-1.0_dp+nn**2)**2)*sqrt(5.0_dp)*nn*&
        &sk(1)-(3.0_dp*nn**4+12.0_dp*ll**2*nn**2-4.0_dp*nn**2+12.0_dp*ll**4&
        &+1.0_dp-12.0_dp*ll**2)*sqrt(10.0_dp)*nn*sk(2)/4.0_dp&
        &+(nn*(3.0_dp*nn**4+12.0_dp*ll**2*nn**2+2.0_dp*nn**2-12.0_dp*ll**2&
        &-1.0_dp+12.0_dp*ll**4)*sk(3))/4.0_dp
    hh(7,5) = (2.0_dp*ll**2-1.0_dp+nn**2)*ll*(4.0_dp*ll**2-3.0_dp&
        &+3.0_dp*nn**2)*sqrt(6.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp&
        &-sqrt(15.0_dp)*ll*(3.0_dp*nn**4+10.0_dp*ll**2*nn**2-5.0_dp*nn**2&
        &-10.0_dp*ll**2+8.0_dp*ll**4+2.0_dp)*sk(2)/4.0_dp+sqrt(6.0_dp)*ll&
        &*(3.0_dp*nn**4+10.0_dp*ll**2*nn**2-2.0_dp*nn**2+3.0_dp+8.0_dp*ll**4&
        &-10.0_dp*ll**2)*sk(3)/8.0_dp

  end subroutine df


  !> Rotation routine for interaction of an f orbital with an f orbital
  subroutine ff(hh, ll, mm, nn, sk)

    !> Dimeric block to put the results in to
    type(TDual3), intent(inout) :: hh(:,:)

    !> Directional cosine ll
    type(TDual3), intent(in) :: ll

    !> Directional cosine mm
    type(TDual3), intent(in) :: mm

    !> Directional cosine nn
    type(TDual3), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    type(TDual3), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 4)
    @:ASSERT(all(shape(hh) >= (/ 7, 7 /)))

    hh(1,1) = - 5.0_dp/ 8.0_dp*(-1.0_dp+nn**2+ll**2)*((4.0_dp*ll**2-1.0_dp&
        &+nn**2)**2)*sk(1)+(15.0_dp/16.0_dp*(nn**6)- 15.0_dp&
        &/8.0_dp*(nn**4)+135.0_dp/ 16.0_dp*(nn**4)*(ll**2)-135.0_dp/8.0_dp&
        &*(ll**2)*(nn**2)+15.0_dp/16.0_dp*(nn**2)+45.0_dp/2.0_dp*(nn**2)&
        &*(ll**4)+135.0_dp/16.0_dp*(ll**2)-45.0_dp/2.0_dp*(ll**4)&
        &+(15.0_dp*ll**6))*sk(2)+(- 3.0_dp/ 8.0_dp*(nn**6)&
        &- 27.0_dp/8.0_dp*(nn**4)&
        &*(ll**2)-3.0_dp/8.0_dp*(nn**4)+3.0_dp/8.0_dp*(nn**2)-(9.0_dp*nn**2&
        &*ll**4)+27.0_dp/4.0_dp*(ll**2)*(nn**2)+3.0_dp/8.0_dp-(6.0_dp*ll**6)&
        &-27.0_dp/8.0_dp*(ll**2)+(9.0_dp*ll**4))*sk(3)&
        &+((nn**6)/16.0_dp+3.0_dp/8.0_dp*(nn**4)+9.0_dp/16.0_dp*(nn**4)*(ll**2)&
        &-9.0_dp/8.0_dp*(ll**2)*(nn**2)+9.0_dp/16.0_dp*(nn**2)&
        &+3.0_dp/2.0_dp*(nn**2)*(ll**4)+9.0_dp/16.0_dp*(ll**2)-3.0_dp/2.0_dp&
        &*(ll**4)+(ll**6))*sk(4)
    hh(2,1) = - 5.0_dp/4.0_dp*(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp&
        &+nn**2)*ll*sqrt(6.0_dp)*nn*sk(1)+5.0_dp/8.0_dp*sqrt(6.0_dp)*nn*ll&
        &*(3.0_dp*nn**4+15.0_dp*ll**2*nn**2-7.0_dp*nn**2+4.0_dp-15.0_dp*ll**2&
        &+12.0_dp*ll**4)*sk(2)-sqrt(6.0_dp)*nn*ll*(3.0_dp*nn**4&
        &+15.0_dp*ll**2*nn**2-10.0_dp*nn**2+5.0_dp-15.0_dp*ll**2+12.0_dp*ll**4)&
        &*sk(3)/4.0_dp+ll*sqrt(6.0_dp)*(nn**4+5.0_dp*ll**2*nn**2&
        &-5.0_dp*nn**2+4.0_dp*ll**4-5.0_dp*ll**2)*nn*sk(4)/8.0_dp
    hh(3,1) = -(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)&
        &*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp&
        &+sqrt(15.0_dp)*(15.0_dp*nn**6-26.0_dp*nn**4+75.0_dp*ll**2*nn**4&
        &-70.0_dp*ll**2*nn**2+11.0_dp*nn**2+60.0_dp*ll**4*nn**2-4.0_dp*ll**4&
        &+3.0_dp*ll**2)*sk(2)/16.0_dp-sqrt(15.0_dp)*(3.0_dp*nn**6-nn**4&
        &+15.0_dp*ll**2*nn**4-3.0_dp*nn**2-2.0_dp*ll**2*nn**2+12.0_dp&
        &*ll**4*nn**2+1.0_dp+4.0_dp*ll**4-5.0_dp*ll**2)*sk(3)/8.0_dp&
        &+sqrt(15.0_dp)*(nn**6+2.0_dp*nn**4+5.0_dp*ll**2*nn**4-3.0_dp*nn**2&
        &+6.0_dp*ll**2*nn**2+4.0_dp*ll**4*nn**2-3.0_dp*ll**2+4.0_dp*ll**4)&
        &*sk(4)/16.0_dp
    hh(4,1) = mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(2.0_dp)&
        &*sqrt(5.0_dp)*nn*(5.0_dp*nn**2-3.0_dp)*sk(1)/8.0_dp-3.0_dp&
        &/16.0_dp*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)*nn*sqrt(2.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(2)+3.0_dp/8.0_dp*sqrt(10.0_dp)*mm&
        &*(4.0_dp*ll**2-1.0_dp+nn**2)*(nn**2+1.0_dp)*nn*sk(3)&
        &-sqrt(5.0_dp)*sqrt(2.0_dp)*(nn**2+3.0_dp)*nn*(4.0_dp*ll**2-1.0_dp&
        &+nn**2)*mm*sk(4)/16.0_dp
    hh(5,1) = mm*(4.0_dp*ll**2-1.0_dp+nn**2)*ll*sqrt(5.0_dp)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp-mm*sqrt(15.0_dp)*ll&
        &*(15.0_dp*nn**4+60.0_dp*ll**2*nn**2-26.0_dp*nn**2+3.0_dp-4.0_dp*ll**2)&
        &*sk(2)/ 16.0_dp+mm*sqrt(15.0_dp)*ll*(3.0_dp*nn**4&
        &+12.0_dp*ll**2*nn**2-10.0_dp*nn**2+4.0_dp*ll**2-1.0_dp)*sk(3)&
        &/8.0_dp-mm*sqrt(15.0_dp)*ll*(nn**4+4.0_dp*ll**2*nn**2-6.0_dp*nn**2&
        &+4.0_dp*ll**2-3.0_dp)*sk(4)/16.0_dp
    hh(6,1) = 5.0_dp/8.0_dp*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*(2.0_dp*ll**2&
        &-1.0_dp+nn**2)*sqrt(6.0_dp)*nn*sk(1)-5.0_dp/16.0_dp*sqrt(6.0_dp)*mm*nn&
        &*(3.0_dp*nn**4+18.0_dp*ll**2*nn**2-4.0_dp*nn**2+24.0_dp*ll**4+1.0_dp&
        &-18.0_dp*ll**2)*sk(2)+sqrt(6.0_dp)*mm*nn*(3.0_dp*nn**4&
        &+18.0_dp*ll**2*nn**2+2.0_dp*nn**2+24.0_dp*ll**4-18.0_dp*ll**2-1.0_dp)&
        &*sk(3)/8.0_dp-mm*sqrt(6.0_dp)*(nn**4+6.0_dp*ll**2*nn**2&
        &+4.0_dp*nn**2+3.0_dp-6.0_dp*ll**2+8.0_dp*ll**4)*nn*sk(4)/16.0_dp
    hh(7,1) = 5.0_dp/8.0_dp*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*mm&
        &*(4.0_dp*ll**2-1.0_dp+nn**2)*sk(1)-15.0_dp/16.0_dp*ll&
        &*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sk(2)+3.0_dp/8.0_dp*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*mm&
        &*(4.0_dp*ll**2-1.0_dp+nn**2)*sk(3)-ll*(4.0_dp*ll**2-3.0_dp&
        &+3.0_dp*nn**2)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sk(4)/16.0_dp
    hh(1,2) = hh(2,1)
    hh(2,2) = -(15.0_dp*ll**2*(-1.0_dp+nn**2+ll**2)*nn**2*sk(1))&
        &+(45.0_dp/2.0_dp*(ll**2)*(nn**4)-5.0_dp/2.0_dp*(nn**4)&
        &-(25.0_dp*ll**2*nn**2)+45.0_dp/2.0_dp*(ll**4)*(nn**2)+5.0_dp/2.0_dp&
        &*(nn**2)+5.0_dp/2.0_dp*(ll**2)-5.0_dp/2.0_dp*(ll**4))*sk(2)&
        &+((-9.0_dp*ll**2*nn**4+4.0_dp*nn**4+13.0_dp*ll**2*nn**2-9.0_dp*ll**4&
        &*nn**2-4.0_dp*nn**2-4.0_dp*ll**2+4.0_dp*ll**4+1.0_dp)*sk(3))&
        &+(3.0_dp/2.0_dp*(ll**2)*(nn**4)-3.0_dp/2.0_dp*(nn**4)+3.0_dp/2.0_dp&
        &*(nn**2)+3.0_dp/2.0_dp*(ll**4)*(nn**2)-(3.0_dp*ll**2*nn**2)&
        &+3.0_dp/2.0_dp*(ll**2)-3.0_dp/2.0_dp*(ll**4))*sk(4)
    hh(3,2) = - 3.0_dp/4.0_dp*ll*(-1.0_dp+nn**2+ll**2)*sqrt(10.0_dp)*nn&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)+(45.0_dp*nn**4-53.0_dp*nn**2&
        &+45.0_dp*ll**2*nn**2+12.0_dp-13.0_dp*ll**2)*nn*ll*sqrt(10.0_dp)&
        &*sk(2)/8.0_dp-sqrt(10.0_dp)*ll*(9.0_dp*nn**4+9.0_dp*ll**2*nn**2&
        &-10.0_dp*nn**2+3.0_dp-5.0_dp*ll**2)*nn*sk(3)/4.0_dp&
        &+3.0_dp/8.0_dp*ll*sqrt(2.0_dp)*(-1.0_dp+nn)*sqrt(5.0_dp)&
        &*(nn**3+nn**2+ll**2*nn+ll**2)*nn*sk(4)
    hh(4,2) = ll*mm*sqrt(15.0_dp)*(nn**2)*(5.0_dp*nn**2-3.0_dp)*sk(1)&
        &/2.0_dp-(5.0_dp*nn**2-1.0_dp)*(3.0_dp*nn**2-1.0_dp)*ll*mm&
        &*sqrt(15.0_dp)*sk(2)/4.0_dp+ll*mm*(nn**2)*(3.0_dp*nn**2-1.0_dp)&
        &*sqrt(15.0_dp)*sk(3)/2.0_dp-ll*mm*sqrt(3.0_dp)*(nn**3+nn**2+nn+1.0_dp)&
        &*sqrt(5.0_dp)*(-1.0_dp+nn)*sk(4)/4.0_dp
    hh(5,2) = 3.0_dp/4.0_dp*(ll**2)*mm*sqrt(10.0_dp)*nn*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-(45.0_dp*ll**2*nn**2-5.0_dp*nn**2-13.0_dp*ll**2+1.0_dp)*nn&
        &*mm*sqrt(10.0_dp)*sk(2)/8.0_dp+sqrt(10.0_dp)*mm*(-4.0_dp*nn**2&
        &+9.0_dp*ll**2*nn**2+2.0_dp-5.0_dp*ll**2)*nn*sk(3)/4.0_dp&
        &-3.0_dp/8.0_dp*mm*sqrt(2.0_dp)*(-1.0_dp+nn)*sqrt(5.0_dp)&
        &*(-nn+ll**2*nn-1.0_dp+ll**2)*nn*sk(4)
    hh(6,2) = 15.0_dp/2.0_dp*ll*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)&
        &*sk(1)-5.0_dp/4.0_dp*(9.0_dp*nn**2-1.0_dp)*(2.0_dp*ll**2&
        &-1.0_dp+nn**2)*ll*mm*sk(2)+(2.0_dp*ll**2-1.0_dp+nn**2)*ll*mm&
        &*(9.0_dp*nn**2-4.0_dp)*sk(3)/2.0_dp-3.0_dp/4.0_dp&
        &*(-1.0_dp+nn**2)*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*mm*sk(4)
    hh(7,2) = 5.0_dp/4.0_dp*(ll**2)*mm*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*sqrt(6.0_dp)*nn*sk(1)-5.0_dp/8.0_dp*sqrt(6.0_dp)*mm*nn&
        &*(9.0_dp*nn**2*ll**2-nn**2+1.0_dp+12.0_dp*ll**4-9.0_dp*ll**2)&
        &*sk(2)+sqrt(6.0_dp)*mm*nn*(9.0_dp*nn**2*ll**2-4.0_dp*nn**2&
        &+12.0_dp*ll**4-9.0_dp*ll**2+2.0_dp)*sk(3)/4.0_dp&
        &-mm*sqrt(6.0_dp)*(3.0_dp*nn**2*ll**2-3.0_dp*nn**2-1.0_dp-3.0_dp*ll**2&
        &+4.0_dp*ll**4)*nn*sk(4)/8.0_dp
    hh(1,3) = hh(3,1)
    hh(2,3) = hh(3,2)
    hh(3,3) = -3.0_dp/8.0_dp*(-1.0_dp+nn**2+ll**2)*((5.0_dp*nn**2-1.0_dp)**2)&
        &*sk(1)+(225.0_dp/16.0_dp*(nn**6)-165.0_dp/8.0_dp*(nn**4)&
        &+225.0_dp/16.0_dp*(ll**2)*(nn**4)+121.0_dp/16.0_dp*(nn**2)-65.0_dp&
        &/8.0_dp*(nn**2)*(ll**2)+(ll**2)/16.0_dp)*sk(2)-5.0_dp/8.0_dp&
        &*(-1.0_dp+nn)*(9.0_dp*nn**5+9.0_dp*nn**4-6.0_dp*nn**3+9.0_dp*ll**2&
        &*nn**3-6.0_dp*nn**2+9.0_dp*nn**2*ll**2+nn-ll**2*nn+1.0_dp-ll**2)*sk(3)&
        &+15.0_dp/16.0_dp*(-nn**2+nn**4+nn**2*ll**2-ll**2)*(-1.0_dp+nn**2)&
        &*sk(4)
    hh(4,3) = mm*sqrt(2.0_dp)*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*nn&
        &*(5.0_dp*nn**2-3.0_dp)*sk(1)/8.0_dp-(5.0_dp*nn**2-1.0_dp)&
        &*sqrt(3.0_dp)*sqrt(2.0_dp)*(15.0_dp*nn**2-11.0_dp)*nn*mm&
        &*sk(2)/ 16.0_dp+5.0_dp/8.0_dp*(-1.0_dp+nn**2)*nn*sqrt(3.0_dp)&
        &*(3.0_dp*nn**2-1.0_dp)*sqrt(2.0_dp)*mm*sk(3)-5.0_dp/16.0_dp&
        &*(-2.0_dp*nn**2+1.0_dp+nn**4)*sqrt(2.0_dp)*nn*sqrt(3.0_dp)*mm&
        &*sk(4)
    hh(5,3) = 3.0_dp/8.0_dp*mm*ll*((5.0_dp*nn**2-1.0_dp)**2)*sk(1)&
        &-(225.0_dp*nn**4-130.0_dp*nn**2+1.0_dp)*mm*ll*sk(2)/16.0_dp&
        &+5.0_dp/8.0_dp*(9.0_dp*nn**3+9.0_dp*nn**2-nn-1.0_dp)*ll*(-1.0_dp+nn)&
        &*mm*sk(3)- 15.0_dp/16.0_dp*((-1.0_dp+nn**2)**2)*mm*ll*sk(4)
    hh(6,3) = 3.0_dp/8.0_dp*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*(5.0_dp*nn**2&
        &-1.0_dp)*sqrt(10.0_dp)*nn*sk(1)-(45.0_dp*nn**4+90.0_dp*nn**2*ll**2&
        &-48.0_dp*nn**2+11.0_dp-26.0_dp*ll**2)*nn*mm*sqrt(10.0_dp)&
        &*sk(2)/16.0_dp+sqrt(10.0_dp)*mm*(9.0_dp*nn**4-6.0_dp*nn**2&
        &+18.0_dp*nn**2*ll**2-10.0_dp*ll**2+1.0_dp)*nn*sk(3)/8.0_dp&
        &-3.0_dp/16.0_dp*mm*(-1.0_dp+nn)*sqrt(5.0_dp)*sqrt(2.0_dp)*(nn**3+nn**2&
        &+2.0_dp*nn*ll**2+nn+2.0_dp*ll**2+1.0_dp)*nn*sk(4)
    hh(7,3) = mm*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp-mm*ll*sqrt(15.0_dp)&
        &*(45.0_dp*nn**4+60.0_dp*nn**2*ll**2-38.0_dp*nn**2+1.0_dp-4.0_dp*ll**2)&
        &*sk(2)/16.0_dp+mm*ll*sqrt(15.0_dp)*(9.0_dp*nn**4+2.0_dp*nn**2&
        &+12.0_dp*nn**2*ll**2+4.0_dp*ll**2-3.0_dp)*sk(3)/8.0_dp-mm&
        &*sqrt(15.0_dp)*ll*(3.0_dp*nn**4+4.0_dp*nn**2*ll**2+6.0_dp*nn**2-1.0_dp&
        &+4.0_dp*ll**2)*sk(4)/ 16.0_dp
    hh(1,4) = hh(4,1)
    hh(2,4) = hh(4,2)
    hh(3,4) = hh(4,3)
    hh(4,4) = (nn**2*(5.0_dp*nn**2-3.0_dp)**2*sk(1))/4.0_dp-3.0_dp&
        &/8.0_dp*(-1.0_dp+nn**2)*((5.0_dp*nn**2-1.0_dp)**2)*sk(2)&
        &+15.0_dp/4.0_dp*(nn**2)*((-1.0_dp+nn**2)**2)*sk(3)-5.0_dp&
        &/8.0_dp*(-1.0_dp+nn**2)*(-2.0_dp*nn**2+1.0_dp+nn**4)*sk(4)
    hh(5,4) = sqrt(2.0_dp)*ll*nn*(5.0_dp*nn**2-3.0_dp)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp-(15.0_dp*nn**2-11.0_dp)&
        &*nn*ll*(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)&
        &*sk(2)/16.0_dp+5.0_dp/8.0_dp*(3.0_dp*nn**3&
        &+3.0_dp*nn**2-nn-1.0_dp)*(-1.0_dp+nn)*sqrt(2.0_dp)*ll*nn*sqrt(3.0_dp)&
        &*sk(3)-5.0_dp/16.0_dp*nn*sqrt(3.0_dp)*ll*((-1.0_dp+nn**2)**2)&
        &*sqrt(2.0_dp)*sk(4)
    hh(6,4) = (2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)*(5.0_dp*nn**2-3.0_dp)&
        &*sqrt(15.0_dp)*sk(1)/4.0_dp-(3.0_dp*nn**2-1.0_dp)&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*(5.0_dp*nn**2-1.0_dp)*sqrt(15.0_dp)&
        &*sk(2)/8.0_dp+sqrt(15.0_dp)*(nn**2)*(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*(3.0_dp*nn**2-1.0_dp)*sk(3)/4.0_dp-sqrt(5.0_dp)&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*(-1.0_dp+nn**4)*sk(4)&
        &/8.0_dp
    hh(7,4) = sqrt(2.0_dp)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*nn&
        &*(5.0_dp*nn**2-3.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp-3.0_dp&
        &/16.0_dp*sqrt(2.0_dp)*(5.0_dp*nn**2-1.0_dp)*ll*(4.0_dp*ll**2-3.0_dp&
        &+3.0_dp*nn**2)*sqrt(5.0_dp)*nn*sk(2)+3.0_dp/8.0_dp&
        &*sqrt(10.0_dp)*nn*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*(nn**2+1.0_dp)&
        &*sk(3)-sqrt(2.0_dp)*sqrt(5.0_dp)*ll*(4.0_dp*ll**2-3.0_dp&
        &+3.0_dp*nn**2)*nn*(nn**2+3.0_dp)*sk(4)/16.0_dp
    hh(1,5) = hh(5,1)
    hh(2,5) = hh(5,2)
    hh(3,5) = hh(5,3)
    hh(4,5) = hh(5,4)
    hh(5,5) = 3.0_dp/8.0_dp*(ll**2)*((5.0_dp*nn**2-1.0_dp)**2)*sk(1)&
        &+(-225.0_dp/ 16.0_dp*(ll**2)*(nn**4)+25.0_dp/16.0_dp*(nn**4)&
        &+65.0_dp/8.0_dp*(ll**2)*(nn**2)-5.0_dp/8.0_dp*(nn**2)+1.0_dp/16.0_dp&
        &-(ll**2)/16.0_dp)*sk(2)+5.0_dp/8.0_dp*(-1.0_dp+nn)&
        &*(9.0_dp*ll**2*nn**3-4.0_dp*nn**3+9.0_dp*ll**2*nn**2&
        &-4.0_dp*nn**2-ll**2*nn&
        &-ll**2)*sk(3)- 15.0_dp/16.0_dp*(ll**2*nn**2+1.0_dp-nn**2-ll**2)&
        &*(-1.0_dp+nn**2)*sk(4)
    hh(6,5) = 3.0_dp/ 8.0_dp*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*&
        &(5.0_dp*nn**2-1.0_dp)*sqrt(10.0_dp)*nn*sk(1)-(45.0_dp*nn**4&
        &+90.0_dp*ll**2*nn**2-68.0_dp*nn**2+15.0_dp-26.0_dp*ll**2)*nn*ll&
        &*sqrt(10.0_dp)*sk(2)/ 16.0_dp+sqrt(10.0_dp)*ll&
        &*(9.0_dp*nn**4-22.0_dp*nn**2+18.0_dp*ll**2*nn**2+9.0_dp-10.0_dp*ll**2)&
        &*nn*sk(3)/ 8.0_dp-3.0_dp/ 16.0_dp*ll*(-1.0_dp+nn)*sqrt(5.0_dp)&
        &*sqrt(2.0_dp)*(nn**3+nn**2+2.0_dp*ll**2*nn-3.0_dp*nn+2.0_dp*ll**2&
        &-3.0_dp)*nn*sk(4)
    hh(7,5) = (ll**2)*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp&
        &-sqrt(15.0_dp)*(-5.0_dp*nn**4+45.0_dp*ll**2*nn**4-58.0_dp*ll**2*nn**2&
        &+6.0_dp*nn**2+60.0_dp*ll**4*nn**2-1.0_dp+5.0_dp*ll**2-4.0_dp*ll**4)&
        &*sk(2)/16.0_dp+sqrt(15.0_dp)*(-4.0_dp*nn**4+9.0_dp*ll**2*nn**4&
        &-14.0_dp*ll**2*nn**2+12.0_dp*ll**4*nn**2+4.0_dp*nn**2-3.0_dp*ll**2&
        &+4.0_dp*ll**4)*sk(3)/8.0_dp-sqrt(15.0_dp)&
        &*(3.0_dp*ll**2.0_dp*nn**4-3.0_dp*nn**4+2*nn**2+4.0_dp*ll**4*nn**2&
        &-6.0_dp*ll**2*nn**2+1.0_dp-5.0_dp*ll**2+4.0_dp*ll**4)*sk(4)&
        &/16.0_dp
    hh(1,6) = hh(6,1)
    hh(2,6) = hh(6,2)
    hh(3,6) = hh(6,3)
    hh(4,6) = hh(6,4)
    hh(5,6) = hh(6,5)
    hh(6,6) = 15.0_dp/4.0_dp*((2.0_dp*ll**2-1.0_dp+nn**2)**2)*(nn**2)&
        &*sk(1)+(- 45.0_dp/8.0_dp*(nn**6)-45.0_dp/2.0_dp*(ll**2)*(nn**4)&
        &+75.0_dp/8.0_dp*(nn**4)+(25.0_dp*ll**2*nn**2)-35.0_dp/8.0_dp*(nn**2)&
        &-45.0_dp/2.0_dp*(ll**4)*(nn**2)-5.0_dp/2.0_dp*(ll**2)+5.0_dp/8.0_dp&
        &+5.0_dp/2.0_dp*(ll**4))*sk(2)+(9.0_dp/4.0_dp*(nn**6)&
        &+(9.0_dp*ll**2*nn**4)-3.0_dp/2.0_dp*(nn**4)-(13.0_dp*ll**2*nn**2)&
        &+(nn**2)/4.0_dp+(9.0_dp*ll**4*nn**2)-(4.0_dp*ll**4)+(4.0_dp*ll**2))&
        &*sk(3)+(- 3.0_dp/8.0_dp*(nn**6)-3.0_dp/2.0_dp*(ll**2)*(nn**4)&
        &-3.0_dp/8.0_dp*(nn**4)+ 3.0_dp/8.0_dp*(nn**2)+(3.0_dp*ll**2*nn**2)&
        &-3.0_dp/2.0_dp*(ll**4)*(nn**2)+3.0_dp/2.0_dp*(ll**4)-3.0_dp/2.0_dp&
        &*(ll**2)+3.0_dp/8.0_dp)*sk(4)
    hh(7,6) = 5.0_dp/8.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*(4.0_dp*ll**2&
        &-3.0_dp+3.0_dp*nn**2)*sqrt(6.0_dp)*nn*sk(1)-5.0_dp/ 16.0_dp&
        &*sqrt(6.0_dp)*ll*nn*(9.0_dp*nn**4+30.0_dp*ll**2*nn**2-16.0_dp*nn**2&
        &-30.0_dp*ll**2&
        &+24.0_dp*ll**4+7.0_dp)*sk(2)+sqrt(6.0_dp)*ll*nn*(9.0_dp*nn**4&
        &+30.0_dp*ll**2*nn**2-10.0_dp*nn**2+24.0_dp*ll**4-30.0_dp*ll**2+5.0_dp)&
        &*sk(3)/8.0_dp-sqrt(6.0_dp)*ll*(3.0_dp*nn**4+10.0_dp*ll**2*nn**2&
        &+8.0_dp*ll**4+5.0_dp-10.0_dp*ll**2)*nn*sk(4)/16.0_dp
    hh(1,7) = hh(7,1)
    hh(2,7) = hh(7,2)
    hh(3,7) = hh(7,3)
    hh(4,7) = hh(7,4)
    hh(5,7) = hh(7,5)
    hh(6,7) = hh(7,6)
    hh(7,7) = 5.0_dp/8.0_dp*(ll**2)*((4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)**2)&
        &*sk(1)+(-135.0_dp/16.0_dp*(ll**2)*(nn**4)+15.0_dp/16.0_dp&
        &*(nn**4)-45.0_dp/2.0_dp*(ll**4)*(nn**2)-15.0_dp/8.0_dp*(nn**2)&
        &+135.0_dp/8.0_dp*(ll**2)*(nn**2)+45.0_dp/2.0_dp*(ll**4)&
        &+15.0_dp/16.0_dp-135.0_dp/16.0_dp*(ll**2)-(15.0_dp*ll**6))*sk(2)&
        &+(27.0_dp/8.0_dp*(ll**2)*(nn**4)-3.0_dp/2.0_dp*(nn**4)+3.0_dp/2.0_dp&
        &*(nn**2)+(9.0_dp*ll**4*nn**2)-27.0_dp/4.0_dp*(ll**2)*(nn**2)+27.0_dp&
        &/8.0_dp*(ll**2)-(9.0_dp*ll**4)+(6.0_dp*ll**6))*sk(3)&
        &+(9.0_dp/16.0_dp*(nn**4)-9.0_dp/ 16.0_dp*(ll**2)*(nn**4)-3.0_dp/2.0_dp&
        &*(ll**4)*(nn**2)+3.0_dp/8.0_dp*(nn**2)+9.0_dp/8.0_dp*(ll**2)*(nn**2)&
        &+3.0_dp/2.0_dp*(ll**4)+1.0_dp/16.0_dp-9.0_dp/16.0_dp*(ll**2)-(ll**6))&
        &*sk(4)

  end subroutine ff

end module dftbp_dftbplus_mxlbornanalytic

#:endif
