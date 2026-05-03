!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Shared MaxwellLink input and diagnostic helpers.
module dftbp_io_mxlcommon
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_io_mxlsocket, only : MxlSocketCommInp
  implicit none

  private
  public :: TMxlSocketInput
  public :: buildMxlExtraJson, checkMxlInit, makeMxlSocketCommInput


  !> Shared socket address and metadata checks for MaxwellLink coupling.
  type :: TMxlSocketInput

    !> Host name for TCP, or socket path for UNIX sockets.
    character(lc) :: host = 'localhost'

    !> TCP port. Values below 1 select UNIX sockets.
    integer :: port = 31415

    !> MaxwellLink communication verbosity.
    integer :: verbosity = 0

    !> Optional molecule id expected from MaxwellLink. Negative values disable checking.
    integer :: moleculeId = -1

    !> Whether to subtract the initial dipole in diagnostics.
    logical :: resetDipole = .false.

  end type TMxlSocketInput


contains


  !> Convert shared MaxwellLink settings into the socket communicator input type.
  subroutine makeMxlSocketCommInput(input, socketInput)

    !> Shared MaxwellLink settings.
    type(TMxlSocketInput), intent(in) :: input

    !> Socket communicator input.
    type(MxlSocketCommInp), intent(out) :: socketInput

    socketInput%host = trim(input%host)
    socketInput%port = input%port
    socketInput%verbosity = input%verbosity

  end subroutine makeMxlSocketCommInput


  !> Validate MaxwellLink INIT metadata against the DFTB+ integrator settings.
  subroutine checkMxlInit(initDt, receivedMoleculeId, expectedDt, expectedMoleculeId, context,&
      & message)

    !> MaxwellLink INIT dt_au value, or a negative value when not supplied.
    real(dp), intent(in) :: initDt

    !> MaxwellLink INIT molecule id.
    integer, intent(in) :: receivedMoleculeId

    !> DFTB+ integrator time step in atomic units.
    real(dp), intent(in) :: expectedDt

    !> User-requested molecule id. Negative values disable checking.
    integer, intent(in) :: expectedMoleculeId

    !> Context string for diagnostics.
    character(*), intent(in) :: context

    !> Empty if valid; otherwise contains an error message.
    character(lc), intent(out) :: message

    message = ''

    if (initDt > 0.0_dp .and.&
        & abs(initDt - expectedDt) > 1.0e-10_dp * max(1.0_dp, abs(expectedDt))) then
      write(message, '(A)') 'MaxwellLink INIT dt_au does not match ' // trim(context)&
          & // ' TimeStep'
      return
    end if

    if (expectedMoleculeId >= 0 .and. receivedMoleculeId /= expectedMoleculeId) then
      write(message, '(A)') 'MaxwellLink INIT molecule id does not match ' // trim(context)&
          & // ' MaxwellLinkSocket MoleculeId'
      return
    end if

  end subroutine checkMxlInit


  !> Build JSON metadata returned with MaxwellLink source data.
  subroutine buildMxlExtraJson(time, energy, energyKin, dipole, dipoleMiddle, extraJson)

    !> Elapsed simulation time in atomic units.
    real(dp), intent(in) :: time

    !> Total energy in Hartree.
    real(dp), intent(in) :: energy

    !> Nuclear kinetic energy in Hartree.
    real(dp), intent(in) :: energyKin

    !> Dipole reported through the primary diagnostic keys.
    real(dp), intent(in) :: dipole(3)

    !> Dipole reported through the midpoint diagnostic keys.
    real(dp), intent(in) :: dipoleMiddle(3)

    !> JSON metadata.
    character(*), intent(out) :: extraJson

    write(extraJson, '(A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,&
        & A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A)')&
        & '{"time_au":', time,&
        & ',"mux_au":', dipole(1),&
        & ',"muy_au":', dipole(2),&
        & ',"muz_au":', dipole(3),&
        & ',"mux_m_au":', dipoleMiddle(1),&
        & ',"muy_m_au":', dipoleMiddle(2),&
        & ',"muz_m_au":', dipoleMiddle(3),&
        & ',"energy_au":', energy,&
        & ',"energy_kin_au":', energyKin,&
        & ',"energy_pot_au":', energy - energyKin,&
        & '}'

  end subroutine buildMxlExtraJson

end module dftbp_io_mxlcommon
