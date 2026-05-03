!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> MaxwellLink coupling helpers for Born-Oppenheimer molecular dynamics.
!!
!! This module reuses dftbp_io_mxlsocket and only contains the BOMD-specific pieces:
!! converting an electric field into nuclear forces through dmu/dR and returning source currents
!! consistent with the velocity-Verlet time step.
module dftbp_md_mxlbomd
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_environment, only : TEnvironment
  use dftbp_io_message, only : error
  use dftbp_io_mxlcommon, only : buildMxlExtraJson, checkMxlInit, makeMxlSocketCommInput,&
      & TMxlSocketInput
  use dftbp_io_mxlsocket, only : MxlSocketComm, MxlSocketComm_init, MxlSocketCommInp
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_bcast
#:endif
  implicit none

  private

  public :: TMxlBomd, TMxlBomd_init
  public :: TMxlBomdInput, mxlBomdDerivTypes


  !> Namespace for available approximations to dmu/dR.
  type :: TMxlBomdDerivTypes

    !> Use fixed user-supplied partial charges.
    integer :: fixedCharges

    !> Use current DFTB Mulliken-like net atomic charges.
    integer :: mullikenCharges

    !> Recompute Born effective charges with DFTB perturbation theory during MD.
    integer :: bornChargesOnTheFly

  end type TMxlBomdDerivTypes


  !> Available dmu/dR approximations.
  type(TMxlBomdDerivTypes), parameter :: mxlBomdDerivTypes =&
      & TMxlBomdDerivTypes(1, 2, 3)


  !> Input data for MaxwellLink-coupled velocity-Verlet BOMD.
  type :: TMxlBomdInput

    !> Shared MaxwellLink socket settings.
    type(TMxlSocketInput) :: socket

    !> Approximation used for dmu/dR.
    integer :: derivType = mxlBomdDerivTypes%mullikenCharges

    !> User supplied charges for fixed-charge coupling.
    real(dp), allocatable :: fixedCharges(:)

    !> Refresh interval for on-the-fly Born effective charges.
    integer :: bornUpdateEvery = 1

    !> Degeneracy tolerance for on-the-fly coordinate perturbation.
    real(dp) :: perturbDegenTol = 1.0e-9_dp

    !> SCC tolerance for on-the-fly coordinate perturbation.
    real(dp) :: perturbSccTol = 1.0e-5_dp

    !> Maximum perturbation SCC iterations for on-the-fly Born charges.
    integer :: maxPerturbIter = 100

  end type TMxlBomdInput


  !> Runtime state for MaxwellLink-coupled velocity-Verlet BOMD.
  type :: TMxlBomd
    private

    !> Input options.
    type(TMxlBomdInput) :: input

    !> Socket communicator. Only used on the global lead process.
    type(MxlSocketComm) :: socket

    !> Whether this object has been initialized.
    logical :: tInitialized = .false.

    !> Whether the socket has been connected.
    logical :: tConnected = .false.

    !> Whether source data is available for MaxwellLink.
    logical :: tHaveResult = .false.

    !> Current electric field used by BOMD, in atomic units.
    real(dp) :: field(3) = 0.0_dp

    !> Initial dipole offset for ResetDipole.
    real(dp) :: initialDipole(3) = 0.0_dp

    !> Whether initialDipole has been set.
    logical :: tInitialDipoleSet = .false.

    !> Current Born effective charges, shape (dipole component, coordinate component, atom).
    real(dp), allocatable :: bornCharges(:,:,:)

  contains

    !> Receive MaxwellLink electric field for the current MD step.
    procedure :: receiveField

    !> Send source-current data back to MaxwellLink.
    procedure :: sendSource

    !> Set the current, scaled electric field used for force/source contractions.
    procedure :: setField

    !> Add MaxwellLink field force contribution to energy derivatives.
    procedure :: addFieldDerivs

    !> Compute fixed-charge dmu/dt from velocity-Verlet half-step velocities.
    procedure :: getSource

    !> Return whether source should be computed from finite-difference dipoles.
    procedure :: usesFiniteDifferenceSource

    !> Compute dmu/dt from two endpoint dipoles.
    procedure :: getFiniteDifferenceSource

    !> Compute the dipole consistent with the selected BOMD coupling model.
    procedure :: getDipole

    !> Prepare dipole diagnostics and midpoint estimate.
    procedure :: getDipoles

    !> Build JSON metadata to return to MaxwellLink.
    procedure :: buildExtraJson

    !> Update currently stored Born effective charges from response-theory output.
    procedure :: updateBornChargesFromResponse

    !> Return whether on-the-fly Born charges are enabled.
    procedure :: usesBornChargesOnTheFly

    !> Return whether the Born charges should be refreshed on this geometry step.
    procedure :: needsBornUpdate

    !> Return perturbation settings for on-the-fly Born charges.
    procedure :: getPerturbSettings

    !> Shut down the socket.
    procedure :: shutdown

  end type TMxlBomd


contains


  !> Initialize MaxwellLink BOMD state.
  subroutine TMxlBomd_init(this, input, nAtom)

    !> Instance.
    type(TMxlBomd), intent(out) :: this

    !> Input data.
    type(TMxlBomdInput), intent(in) :: input

    !> Number of central-cell atoms.
    integer, intent(in) :: nAtom

    this%input = input
    this%field(:) = 0.0_dp
    this%initialDipole(:) = 0.0_dp
    this%tInitialDipoleSet = .false.
    this%tHaveResult = .false.

    if (this%input%bornUpdateEvery < 1) then
      call error("MaxwellLinkSocket BornUpdateEvery must be positive")
    end if

    select case (this%input%derivType)
    case (mxlBomdDerivTypes%fixedCharges)
      if (.not. allocated(this%input%fixedCharges)) then
        call error("MaxwellLinkSocket FixedCharges requires Charges")
      end if
      if (size(this%input%fixedCharges) /= nAtom) then
        call error("MaxwellLinkSocket Charges length must match the number of atoms")
      end if

    case (mxlBomdDerivTypes%mullikenCharges)
      continue

    case (mxlBomdDerivTypes%bornChargesOnTheFly)
      continue

    case default
      call error("Internal error: unknown MaxwellLinkSocket DipoleDerivative mode")
    end select

    this%tInitialized = .true.

  end subroutine TMxlBomd_init


  !> Receive field data from MaxwellLink.
  subroutine receiveField(this, env, deltaT, field, tStop)

    !> Instance.
    class(TMxlBomd), intent(inout) :: this

    !> Environment settings.
    type(TEnvironment), intent(inout) :: env

    !> MD time step in atomic units.
    real(dp), intent(in) :: deltaT

    !> Electric field in atomic units.
    real(dp), intent(out) :: field(3)

    !> Whether MaxwellLink requested termination.
    logical, intent(out) :: tStop

    type(MxlSocketCommInp) :: socketInput
    character(lc) :: msg
    logical :: tReceivedInit
    real(dp) :: initDt
    integer :: receivedMoleculeId

    if (.not. this%tInitialized) then
      call error("MaxwellLink BOMD state used before initialization")
    end if

  #:if WITH_MPI
    if (env%mpi%tGlobalLead) then
  #:endif
      if (.not. this%tConnected) then
        call makeMxlSocketCommInput(this%input%socket, socketInput)
        call MxlSocketComm_init(this%socket, socketInput)
        this%tConnected = .true.
      end if

      call this%socket%receiveField(this%tHaveResult, field, tStop, tReceivedInit)
      if (tReceivedInit) then
        initDt = this%socket%getInitDt()
        receivedMoleculeId = this%socket%getMoleculeId()
      else
        initDt = -1.0_dp
        receivedMoleculeId = -1
      end if
  #:if WITH_MPI
    else
      field(:) = 0.0_dp
      tStop = .false.
      tReceivedInit = .false.
      initDt = -1.0_dp
      receivedMoleculeId = -1
    end if
    call mpifx_bcast(env%mpi%globalComm, field)
    call mpifx_bcast(env%mpi%globalComm, tStop)
    call mpifx_bcast(env%mpi%globalComm, tReceivedInit)
    call mpifx_bcast(env%mpi%globalComm, initDt)
    call mpifx_bcast(env%mpi%globalComm, receivedMoleculeId)
  #:endif

    if (tReceivedInit) then
      call checkMxlInit(initDt, receivedMoleculeId, deltaT, this%input%socket%moleculeId,&
          & "VelocityVerlet", msg)
      if (len_trim(msg) > 0) then
        call error(trim(msg))
      end if
    end if

  end subroutine receiveField


  !> Send source data to MaxwellLink.
  subroutine sendSource(this, env, energy, source, extraJson, tStop)

    !> Instance.
    class(TMxlBomd), intent(inout) :: this

    !> Environment settings.
    type(TEnvironment), intent(inout) :: env

    !> Total molecular energy in Hartree.
    real(dp), intent(in) :: energy

    !> Source current dmu/dt in atomic units.
    real(dp), intent(in) :: source(3)

    !> Additional JSON payload.
    character(*), intent(in) :: extraJson

    !> Whether MaxwellLink requested termination.
    logical, intent(out) :: tStop

  #:if WITH_MPI
    if (env%mpi%tGlobalLead) then
  #:endif
      call this%socket%sendSource(energy, source, trim(extraJson), tStop)
  #:if WITH_MPI
    else
      tStop = .false.
    end if
    call mpifx_bcast(env%mpi%globalComm, tStop)
  #:endif
    this%tHaveResult = .false.

  end subroutine sendSource


  !> Set current, scaled electric field.
  subroutine setField(this, field)

    !> Instance.
    class(TMxlBomd), intent(inout) :: this

    !> Electric field in atomic units.
    real(dp), intent(in) :: field(3)

    this%field(:) = field(:)

  end subroutine setField


  !> Add MaxwellLink external-field forces to energy derivatives.
  subroutine addFieldDerivs(this, qOutput, q0, derivs)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Electron populations.
    real(dp), intent(in) :: qOutput(:,:,:)

    !> Reference populations.
    real(dp), intent(in) :: q0(:,:,:)

    !> Energy derivatives with respect to coordinates.
    real(dp), intent(inout) :: derivs(:,:)

    real(dp) :: force(3)
    integer :: iAtom, nAtom

    nAtom = size(derivs, dim=2)

    do iAtom = 1, nAtom
      call getAtomFieldForce(this, qOutput, q0, iAtom, force)
      derivs(:, iAtom) = derivs(:, iAtom) - force(:)
    end do

  end subroutine addFieldDerivs


  !> Return fixed-charge dmu/dt from velocity-Verlet half-step velocities.
  subroutine getSource(this, indMovedAtom, movedVeloHalf, source)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Moved atom indices.
    integer, intent(in) :: indMovedAtom(:)

    !> Half-step velocities of moved atoms.
    real(dp), intent(in) :: movedVeloHalf(:,:)

    !> Source current dmu/dt.
    real(dp), intent(out) :: source(3)

    integer :: iMoved, iAtom

    source(:) = 0.0_dp

    if (this%input%derivType /= mxlBomdDerivTypes%fixedCharges) then
      call error("MaxwellLink half-step velocity source is only valid for FixedCharges")
    end if

    do iMoved = 1, size(indMovedAtom)
      iAtom = indMovedAtom(iMoved)
      source(:) = source(:) + this%input%fixedCharges(iAtom) * movedVeloHalf(:, iMoved)
    end do

  end subroutine getSource


  !> Return whether the selected coupling needs finite-difference dipole current.
  function usesFiniteDifferenceSource(this) result(tFiniteDiff)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Whether to compute the source from endpoint dipoles.
    logical :: tFiniteDiff

    tFiniteDiff = this%input%derivType /= mxlBomdDerivTypes%fixedCharges

  end function usesFiniteDifferenceSource


  !> Return dmu/dt from a finite difference of endpoint dipoles.
  subroutine getFiniteDifferenceSource(this, deltaT, dipole0, dipole1, source)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Time step in atomic units.
    real(dp), intent(in) :: deltaT

    !> Dipole at the force evaluation time.
    real(dp), intent(in) :: dipole0(3)

    !> Dipole one full time step after the force evaluation time.
    real(dp), intent(in) :: dipole1(3)

    !> Source current dmu/dt.
    real(dp), intent(out) :: source(3)

    if (deltaT <= 0.0_dp) then
      call error("MaxwellLink BOMD finite-difference source requires a positive time step")
    end if

    source(:) = (dipole1(:) - dipole0(:)) / deltaT

  end subroutine getFiniteDifferenceSource


  !> Return the dipole consistent with the selected coupling model.
  subroutine getDipole(this, qOutput, q0, coord, iAtInCentralRegion, dftbDipole, dipole)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Electron populations.
    real(dp), intent(in) :: qOutput(:,:,:)

    !> Reference populations.
    real(dp), intent(in) :: q0(:,:,:)

    !> Atomic coordinates in atomic units.
    real(dp), intent(in) :: coord(:,:)

    !> Atoms included in the molecular dipole.
    integer, intent(in) :: iAtInCentralRegion(:)

    !> DFTB+ molecular dipole in atomic units.
    real(dp), intent(in) :: dftbDipole(3)

    !> Coupling dipole in atomic units.
    real(dp), intent(out) :: dipole(3)

    real(dp) :: charge
    integer :: ii, iAtom

    dipole(:) = 0.0_dp

    select case (this%input%derivType)
    case (mxlBomdDerivTypes%fixedCharges)
      do ii = 1, size(iAtInCentralRegion)
        iAtom = iAtInCentralRegion(ii)
        dipole(:) = dipole(:) + this%input%fixedCharges(iAtom) * coord(:, iAtom)
      end do

    case (mxlBomdDerivTypes%mullikenCharges)
      do ii = 1, size(iAtInCentralRegion)
        iAtom = iAtInCentralRegion(ii)
        charge = sum(q0(:, iAtom, 1) - qOutput(:, iAtom, 1))
        dipole(:) = dipole(:) + charge * coord(:, iAtom)
      end do

    case (mxlBomdDerivTypes%bornChargesOnTheFly)
      dipole(:) = dftbDipole(:)
    end select

  end subroutine getDipole


  !> Return current and midpoint dipoles for diagnostics.
  subroutine getDipoles(this, deltaT, dipole, source, dipoleCurrent, dipoleMiddle)

    !> Instance.
    class(TMxlBomd), intent(inout) :: this

    !> Time step in atomic units.
    real(dp), intent(in) :: deltaT

    !> Current dipole in atomic units.
    real(dp), intent(in) :: dipole(3)

    !> Source current in atomic units.
    real(dp), intent(in) :: source(3)

    !> Dipole after reset-dipole adjustment.
    real(dp), intent(out) :: dipoleCurrent(3)

    !> Midpoint dipole estimate.
    real(dp), intent(out) :: dipoleMiddle(3)

    if (this%input%socket%resetDipole .and. .not. this%tInitialDipoleSet) then
      this%initialDipole(:) = dipole(:)
      this%tInitialDipoleSet = .true.
    end if

    dipoleCurrent(:) = dipole(:)
    if (this%input%socket%resetDipole) then
      dipoleCurrent(:) = dipoleCurrent(:) - this%initialDipole(:)
    end if
    dipoleMiddle(:) = dipoleCurrent(:) + 0.5_dp * deltaT * source(:)

  end subroutine getDipoles


  !> Build JSON metadata returned with source data.
  subroutine buildExtraJson(this, time, energy, energyKin, dipole, dipoleMiddle, extraJson)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Elapsed simulation time in atomic units.
    real(dp), intent(in) :: time

    !> Total energy in Hartree.
    real(dp), intent(in) :: energy

    !> Nuclear kinetic energy in Hartree.
    real(dp), intent(in) :: energyKin

    !> Dipole at the force evaluation time.
    real(dp), intent(in) :: dipole(3)

    !> Dipole half a time step after the force evaluation time.
    real(dp), intent(in) :: dipoleMiddle(3)

    !> JSON metadata.
    character(*), intent(out) :: extraJson

    call buildMxlExtraJson(time, energy, energyKin, dipoleMiddle, dipole, extraJson)

  end subroutine buildExtraJson


  !> Store current Born effective charges from response theory.
  subroutine updateBornChargesFromResponse(this, bornCharges)

    !> Instance.
    class(TMxlBomd), intent(inout) :: this

    !> Born effective charges in response-theory force/electron-dipole convention.
    real(dp), intent(in) :: bornCharges(:,:,:)

    if (allocated(this%bornCharges)) then
      deallocate(this%bornCharges)
    end if
    allocate(this%bornCharges(size(bornCharges, dim=1), size(bornCharges, dim=2),&
        & size(bornCharges, dim=3)))
    ! Response theory returns the force/electron-dipole convention. MaxwellLink contracts
    ! derivatives of the molecular dipole returned by getDipoleMoment.
    this%bornCharges(:,:,:) = -bornCharges(:,:,:)

  end subroutine updateBornChargesFromResponse


  !> Return whether on-the-fly Born charges are enabled.
  function usesBornChargesOnTheFly(this) result(tEnabled)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Whether enabled.
    logical :: tEnabled

    tEnabled = this%input%derivType == mxlBomdDerivTypes%bornChargesOnTheFly

  end function usesBornChargesOnTheFly


  !> Return whether Born charges should be refreshed on this geometry step.
  function needsBornUpdate(this, iGeoStep) result(tUpdate)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Geometry step.
    integer, intent(in) :: iGeoStep

    !> Whether to update.
    logical :: tUpdate

    tUpdate = this%usesBornChargesOnTheFly() .and.&
        & (.not. allocated(this%bornCharges) .or. mod(iGeoStep, this%input%bornUpdateEvery) == 0)

  end function needsBornUpdate


  !> Return perturbation settings for on-the-fly Born charges.
  subroutine getPerturbSettings(this, maxPerturbIter, perturbSccTol)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Maximum perturbation SCC iterations.
    integer, intent(out) :: maxPerturbIter

    !> Perturbation SCC tolerance.
    real(dp), intent(out) :: perturbSccTol

    maxPerturbIter = this%input%maxPerturbIter
    perturbSccTol = this%input%perturbSccTol

  end subroutine getPerturbSettings


  !> Shut down socket connection.
  subroutine shutdown(this, env)

    !> Instance.
    class(TMxlBomd), intent(inout) :: this

    !> Environment settings.
    type(TEnvironment), intent(in) :: env

  #:if WITH_MPI
    if (env%mpi%tGlobalLead) then
  #:endif
      if (this%tConnected) then
        call this%socket%shutdown()
      end if
  #:if WITH_MPI
    end if
  #:endif
    this%tConnected = .false.

  end subroutine shutdown


  !> Return field force on one atom.
  subroutine getAtomFieldForce(this, qOutput, q0, iAtom, force)

    !> Instance.
    class(TMxlBomd), intent(in) :: this

    !> Electron populations.
    real(dp), intent(in) :: qOutput(:,:,:)

    !> Reference populations.
    real(dp), intent(in) :: q0(:,:,:)

    !> Atom index.
    integer, intent(in) :: iAtom

    !> Field force.
    real(dp), intent(out) :: force(3)

    real(dp) :: charge
    integer :: iCart

    select case (this%input%derivType)
    case (mxlBomdDerivTypes%fixedCharges)
      force(:) = this%input%fixedCharges(iAtom) * this%field(:)

    case (mxlBomdDerivTypes%mullikenCharges)
      charge = sum(q0(:, iAtom, 1) - qOutput(:, iAtom, 1))
      force(:) = charge * this%field(:)

    case (mxlBomdDerivTypes%bornChargesOnTheFly)
      if (.not. allocated(this%bornCharges)) then
        call error("MaxwellLink BOMD force requested before Born charges are available")
      end if
      do iCart = 1, 3
        force(iCart) = dot_product(this%field(:), this%bornCharges(:, iCart, iAtom))
      end do

    end select

  end subroutine getAtomFieldForce

end module dftbp_md_mxlbomd
