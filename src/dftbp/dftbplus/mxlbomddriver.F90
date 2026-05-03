!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Driver-level helpers for MaxwellLink-coupled Born-Oppenheimer MD.
module dftbp_dftbplus_mxlbomddriver
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_dftbplus_initprogram, only : TDftbPlusMain
  use dftbp_dftbplus_outputfiles, only : autotestTag, bornChargesOut, resultsTag
  use dftbp_md_mdintegrator, only : state
  implicit none

  private
  public :: applyMxlBomdFieldDerivs, getMxlBomdSoluteDipole, receiveMxlBomdField
  public :: sendMxlBomdSource


contains


  !> Receive and store the MaxwellLink electric field for the current MD step.
  subroutine receiveMxlBomdField(this, env, tStopDriver)

    !> Global variables.
    type(TDftbPlusMain), intent(inout) :: this

    !> Environment settings.
    type(TEnvironment), intent(inout) :: env

    !> Whether the driver should stop.
    logical, intent(inout) :: tStopDriver

    real(dp) :: field(3)

    call this%mxlBomd%receiveField(env, this%deltaT, field, tStopDriver)
    if (tStopDriver) then
      return
    end if
    field(:) = this%eFieldScaling%scaledExtEField(field)
    call this%mxlBomd%setField(field)

  end subroutine receiveMxlBomdField


  !> Update MaxwellLink Born charges when needed and add field-force derivatives.
  subroutine applyMxlBomdFieldDerivs(this, env, iGeoStep, errStatus)

    !> Global variables.
    type(TDftbPlusMain), intent(inout) :: this

    !> Environment settings.
    type(TEnvironment), intent(inout) :: env

    !> Geometry step.
    integer, intent(in) :: iGeoStep

    !> Error status.
    type(TStatus), intent(inout) :: errStatus

    if (this%mxlBomd%needsBornUpdate(iGeoStep)) then
      call updateMxlBomdBornCharges(this, env, errStatus)
      if (errStatus%hasError()) then
        return
      end if
    end if
    call this%mxlBomd%addFieldDerivs(this%qOutput, this%q0, this%derivs)

  end subroutine applyMxlBomdFieldDerivs


  !> Return the current scaled solute dipole used for MaxwellLink source diagnostics.
  subroutine getMxlBomdSoluteDipole(this, soluteDipole)

    !> Global variables.
    type(TDftbPlusMain), intent(inout) :: this

    !> Scaled solute dipole in atomic units.
    real(dp), intent(out) :: soluteDipole(3)

    call this%mxlBomd%getDipole(this%qOutput, this%q0, this%coord0,&
        & this%iAtInCentralRegion, this%dipoleMoment(:, this%deltaDftb%iFinal), soluteDipole)
    soluteDipole(:) = this%eFieldScaling%scaledSoluteDipole(soluteDipole)

  end subroutine getMxlBomdSoluteDipole


  !> Send MaxwellLink source data after the velocity-Verlet geometry update.
  subroutine sendMxlBomdSource(this, env, iGeoStep, soluteDipole, tStopDriver, soluteDipoleNext)

    !> Global variables.
    type(TDftbPlusMain), intent(inout) :: this

    !> Environment settings.
    type(TEnvironment), intent(inout) :: env

    !> Geometry step.
    integer, intent(in) :: iGeoStep

    !> Scaled solute dipole at the force evaluation time.
    real(dp), intent(in) :: soluteDipole(3)

    !> Whether the driver should stop.
    logical, intent(inout) :: tStopDriver

    !> Optional scaled solute dipole at the endpoint geometry.
    real(dp), intent(in), optional :: soluteDipoleNext(3)

    character(4 * lc) :: extraJson
    logical :: tMxlStop
    real(dp) :: dipole(3), dipoleMiddle(3), source(3)
    real(dp) :: energyKin, energyMerminKin
    real(dp), allocatable :: veloHalf(:,:)

    energyMerminKin = this%dftbEnergy(this%deltaDftb%iFinal)%EMerminKin
    energyKin = this%dftbEnergy(this%deltaDftb%iFinal)%Ekin

    if (present(soluteDipoleNext)) then
      call this%mxlBomd%getFiniteDifferenceSource(this%deltaT, soluteDipole, soluteDipoleNext,&
          & source)
    else
      allocate(veloHalf(3, this%nMovedAtom))
      call state(this%pMdIntegrator, velocities=veloHalf)
      call this%mxlBomd%getSource(this%indMovedAtom, veloHalf, source)
      deallocate(veloHalf)
      source(:) = this%eFieldScaling%scaledSoluteDipole(source)
    end if

    call this%mxlBomd%getDipoles(this%deltaT, soluteDipole, source, dipole, dipoleMiddle)
    call this%mxlBomd%buildExtraJson((real(iGeoStep, dp) + 0.5_dp) * this%deltaT,&
        & energyMerminKin, energyKin, dipole, dipoleMiddle, extraJson)
    call this%mxlBomd%sendSource(env, energyMerminKin, source, trim(extraJson), tMxlStop)
    tStopDriver = tStopDriver .or. tMxlStop

  end subroutine sendMxlBomdSource


  !> Recompute Born effective charges for MaxwellLink-coupled BOMD.
  subroutine updateMxlBomdBornCharges(this, env, errStatus)

    !> Global variables.
    type(TDftbPlusMain), intent(inout) :: this

    !> Environment settings.
    type(TEnvironment), intent(inout) :: env

    !> Error status.
    type(TStatus), intent(out) :: errStatus

    real(dp), allocatable :: bornCharges(:,:,:)
    integer :: maxPerturbIter
    real(dp) :: perturbSccTol

    call this%mxlBomd%getPerturbSettings(maxPerturbIter, perturbSccTol)
    call this%response%dxAtom(env, this%parallelKS, this%filling, this%eigen, this%eigVecsReal,&
        & this%eigvecsCplx, this%rhoPrim, this%potential, this%qOutput, this%q0,&
        & this%ints%hamiltonian, this%ints%overlap, this%skHamCont, this%skOverCont,&
        & this%mxlBomdNonSccDeriv, this%orb, this%nAtom, this%species, this%speciesName,&
        & this%neighbourList, this%nNeighbourSK, this%denseDesc, this%iSparseStart,&
        & this%img2CentCell, this%coord, this%scc, maxPerturbIter, perturbSccTol,&
        & this%nMixElements, this%nIneqOrb, this%iEqOrbitals, this%tempElec, this%Ef,&
        & this%tFixEf, this%spinW, this%thirdOrd, this%dftbU, this%iEqBlockDftbu,&
        & this%onSiteElements, this%iEqBlockOnSite, this%hybridXc, this%nNeighbourCam,&
        & this%chrgMixerReal, .false., this%taggedWriter, .false., autotestTag, .false.,&
        & resultsTag, .false., this%fdDetailedOut%unit, this%kPoint, this%kWeight,&
        & this%iCellVec, this%cellVec, this%tPeriodic, this%tHelical, this%tMulliken, errStatus,&
        & bornChargesOut=bornCharges, tWriteResponseOutput=.false.)
    if (errStatus%hasError()) then
      return
    end if

    call this%mxlBomd%updateBornChargesFromResponse(bornCharges)

  end subroutine updateMxlBomdBornCharges

end module dftbp_dftbplus_mxlbomddriver
