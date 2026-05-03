!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Socket communication with MaxwellLink molecular-driver hubs.
!!
!! MaxwellLink uses an i-PI-compatible packet layout. In the EM coupling convention used here,
!! POSDATA/FIELDDATA carries a single 3-vector electric field in atomic units, and FORCEREADY
!! returns a single 3-vector source current, also in atomic units.
module dftbp_io_mxlsocket
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_extlibs_fsockets, only : close_socket, connect_inet_socket, connect_unix_socket,&
      & readbuffer, writebuffer
  use dftbp_io_logger, only : LogWriter
  use dftbp_io_message, only : error, warning
  implicit none

  private
  public :: MxlSocketCommInp
  public :: MxlSocketComm, MxlSocketComm_init


  !> Input for initialising MxlSocketComm.
  type :: MxlSocketCommInp

    !> Host name for TCP, or socket path for UNIX sockets.
    character(:), allocatable :: host

    !> Port to connect to if using TCP; if less than 1, host is treated as a UNIX socket path.
    integer :: port

    !> Verbosity level of communication logging.
    integer :: verbosity

  end type MxlSocketCommInp


  !> Communicator for MaxwellLink socket coupling.
  type :: MxlSocketComm
    private

    !> Used to log messages.
    type(LogWriter) :: logger

    !> Socket number.
    integer :: socket

    !> Socket has been connected.
    logical :: tInit = .false.

    !> INIT packet has been received.
    logical :: tHandshake = .false.

    !> Molecule id assigned by MaxwellLink.
    integer :: moleculeId = -1

    !> Time step requested by MaxwellLink INIT payload, if provided.
    real(dp) :: initDt = -1.0_dp

    !> Raw INIT JSON payload.
    character(:), allocatable :: initJson

  contains

    !> Receive a field packet, handling STATUS/INIT/STOP packets on the way.
    procedure :: receiveField

    !> Send a source-current packet, handling STATUS/STOP packets on the way.
    procedure :: sendSource

    !> Shut the socket down.
    procedure :: shutdown

    !> Return the molecule id assigned in the INIT packet.
    procedure :: getMoleculeId

    !> Return the dt_au value from the INIT packet, or a negative value if absent.
    procedure :: getInitDt

  end type MxlSocketComm


  !> Length of MaxwellLink/i-PI message headers.
  integer, parameter :: MXL_MSGLEN = 12


contains


  !> Construct MxlSocketComm instance.
  subroutine MxlSocketComm_init(this, input)

    !> Instance.
    type(MxlSocketComm), intent(out) :: this

    !> Input data.
    type(MxlSocketCommInp), intent(in) :: input

    character(lc) :: msg
    logical :: tUnix

    if (.not. allocated(input%host)) then
      call error("MaxwellLink socket host/path was not set")
    end if

    this%logger = LogWriter(input%verbosity)

    tUnix = input%port < 1
    if (tUnix) then
      call this%logger%write('mxlSocketCreate: Establishing UNIX socket connection to '&
          & // trim(input%host), 1)
      call connect_unix_socket(this%socket, input%host)
    else
      call this%logger%write('mxlSocketCreate: Establishing TCP connection', 1)
      call this%logger%write('Host: ' // trim(input%host), 1)
      write(msg, '(A,I0)') 'Port: ', input%port
      call this%logger%write(msg, 1)
      call connect_inet_socket(this%socket, input%host, input%port)
    end if

    this%tInit = .true.
    call this%logger%write('mxlSocketCreate: ...Done', 1)

  end subroutine MxlSocketComm_init


  !> Receive field data from MaxwellLink.
  subroutine receiveField(this, tHaveResult, efield, tStop, tReceivedInit)

    !> Instance.
    class(MxlSocketComm), intent(inout) :: this

    !> Whether DFTB+ already has source data ready.
    logical, intent(in) :: tHaveResult

    !> Electric field in atomic units.
    real(dp), intent(out) :: efield(3)

    !> MaxwellLink requested termination.
    logical, intent(out) :: tStop

    !> INIT packet was consumed during this call.
    logical, intent(out) :: tReceivedInit

    character(len=MXL_MSGLEN) :: buffer, header

    if (.not. this%tInit) then
      call error("MaxwellLink socket receive requested before initialisation")
    end if

    tStop = .false.
    tReceivedInit = .false.
    efield(:) = 0.0_dp

    listen: do
      call readbuffer(this%socket, header)
      call this%logger%write('mxlsocket%receiveField: read from socket: '&
          & // trim(header), 3)

      select case (trim(header))
      case ('STATUS')
        if (.not. this%tHandshake) then
          buffer = 'NEEDINIT'
        else if (tHaveResult) then
          buffer = 'HAVEDATA'
        else
          buffer = 'READY'
        end if
        call writebuffer(this%socket, buffer)
        call this%logger%write('mxlsocket%receiveField: write to socket: '&
            & // trim(buffer), 3)

      case ('INIT')
        call receiveInit(this)
        tReceivedInit = .true.

      case ('POSDATA', 'FIELDDATA')
        call receiveFieldData(this, efield)
        exit listen

      case ('STOP', 'EXIT')
        buffer = 'BYE'
        call writebuffer(this%socket, buffer)
        call warning("mxlsocket%receiveField: MaxwellLink requested stop.")
        tStop = .true.
        exit listen

      case default
        call error("mxlsocket%receiveField: unexpected MaxwellLink message '"&
            & // trim(header) // "'")
      end select
    end do listen

  end subroutine receiveField


  !> Send source data to MaxwellLink.
  subroutine sendSource(this, energy, source, extraJson, tStop)

    !> Instance.
    class(MxlSocketComm), intent(inout) :: this

    !> Total energy in Hartree.
    real(dp), intent(in) :: energy

    !> Source-current vector in atomic units.
    real(dp), intent(in) :: source(3)

    !> Additional JSON payload.
    character(*), intent(in) :: extraJson

    !> MaxwellLink requested termination.
    logical, intent(out) :: tStop

    character(len=MXL_MSGLEN) :: buffer, header

    if (.not. this%tInit) then
      call error("MaxwellLink socket send requested before initialisation")
    end if

    tStop = .false.

    listen: do
      call readbuffer(this%socket, header)
      call this%logger%write('mxlsocket%sendSource: read from socket: '&
          & // trim(header), 3)

      select case (trim(header))
      case ('STATUS')
        buffer = 'HAVEDATA'
        call writebuffer(this%socket, buffer)
        call this%logger%write('mxlsocket%sendSource: write to socket: HAVEDATA', 3)

      case ('GETFORCE', 'GETSOURCE')
        call sendSourceReady(this, energy, source, extraJson)
        exit listen

      case ('STOP', 'EXIT')
        buffer = 'BYE'
        call writebuffer(this%socket, buffer)
        call warning("mxlsocket%sendSource: MaxwellLink requested stop.")
        tStop = .true.
        exit listen

      case default
        call error("mxlsocket%sendSource: unexpected MaxwellLink message '"&
            & // trim(header) // "'")
      end select
    end do listen

  end subroutine sendSource


  !> Receive INIT metadata.
  subroutine receiveInit(this)

    !> Instance.
    class(MxlSocketComm), intent(inout) :: this

    integer :: nChar
    logical :: tFoundDt

    call readbuffer(this%socket, this%moleculeId)
    call readbuffer(this%socket, nChar)

    if (nChar < 0) then
      call error("mxlsocket%receiveInit: negative INIT payload length")
    end if

    if (allocated(this%initJson)) then
      deallocate(this%initJson)
    end if
    if (nChar > 0) then
      allocate(character(len=nChar) :: this%initJson)
      call readbuffer(this%socket, this%initJson)
    else
      this%initJson = '{}'
    end if

    call extractJsonReal(this%initJson, 'dt_au', this%initDt, tFoundDt)
    if (.not. tFoundDt) then
      this%initDt = -1.0_dp
    end if

    this%tHandshake = .true.

  end subroutine receiveInit


  !> Receive field data after POSDATA/FIELDDATA header has been consumed.
  subroutine receiveFieldData(this, efield)

    !> Instance.
    class(MxlSocketComm), intent(inout) :: this

    !> Electric field in atomic units.
    real(dp), intent(out) :: efield(3)

    integer :: nAtom
    real(dp) :: cell(9), invCell(9)
    real(dp) :: fieldData(3)

    call readbuffer(this%socket, cell)
    call readbuffer(this%socket, invCell)
    call readbuffer(this%socket, nAtom)

    if (nAtom /= 1) then
      call error("mxlsocket%receiveFieldData: expected exactly one field vector")
    end if

    call readbuffer(this%socket, fieldData)
    efield(:) = fieldData(:)

  end subroutine receiveFieldData


  !> Send FORCEREADY packet with source-current data.
  subroutine sendSourceReady(this, energy, source, extraJson)

    !> Instance.
    class(MxlSocketComm), intent(inout) :: this

    !> Total energy in Hartree.
    real(dp), intent(in) :: energy

    !> Source-current vector in atomic units.
    real(dp), intent(in) :: source(3)

    !> Additional JSON payload.
    character(*), intent(in) :: extraJson

    character(len=MXL_MSGLEN) :: buffer
    integer :: nChar, nAtom
    real(dp) :: virial(9)

    buffer = 'FORCEREADY'
    nAtom = 1
    virial(:) = 0.0_dp

    call writebuffer(this%socket, buffer)
    call writebuffer(this%socket, energy)
    call writebuffer(this%socket, nAtom)
    call writebuffer(this%socket, source)
    call writebuffer(this%socket, virial)

    nChar = len_trim(extraJson)
    call writebuffer(this%socket, nChar)
    if (nChar > 0) then
      call writebuffer(this%socket, extraJson(:nChar))
    end if

    call this%logger%write('mxlsocket%sendSource: write to socket: FORCEREADY', 3)

  end subroutine sendSourceReady


  !> Extract a real scalar from a small JSON object without depending on a JSON parser.
  subroutine extractJsonReal(json, key, value, tFound)

    !> JSON text.
    character(*), intent(in) :: json

    !> Key to search for.
    character(*), intent(in) :: key

    !> Parsed value.
    real(dp), intent(out) :: value

    !> Whether the value was found and parsed.
    logical, intent(out) :: tFound

    character(:), allocatable :: token
    integer :: first, last, pos, relColon, relEnd, stat

    tFound = .false.
    value = 0.0_dp

    pos = index(json, '"' // trim(key) // '"')
    if (pos == 0) then
      pos = index(json, trim(key))
    end if
    if (pos == 0) then
      return
    end if

    relColon = index(json(pos:), ':')
    if (relColon == 0) then
      return
    end if

    first = pos + relColon
    do while (first <= len(json))
      if (json(first:first) /= ' ' .and. json(first:first) /= '"') then
        exit
      end if
      first = first + 1
    end do
    if (first > len(json)) then
      return
    end if

    relEnd = scan(json(first:), ',}')
    if (relEnd == 0) then
      last = len(json)
    else
      last = first + relEnd - 2
    end if
    if (last < first) then
      return
    end if

    token = json(first:last)
    read(token, *, iostat=stat) value
    tFound = stat == 0

  end subroutine extractJsonReal


  !> Shuts down the socket.
  subroutine shutdown(this)

    !> Instance.
    class(MxlSocketComm), intent(inout) :: this

    if (this%tInit) then
      call close_socket(this%socket)
      this%tInit = .false.
      this%tHandshake = .false.
    end if

  end subroutine shutdown


  !> Return assigned molecule id.
  function getMoleculeId(this) result(moleculeId)

    !> Instance.
    class(MxlSocketComm), intent(in) :: this

    !> Molecule id.
    integer :: moleculeId

    moleculeId = this%moleculeId

  end function getMoleculeId


  !> Return INIT dt_au value.
  function getInitDt(this) result(initDt)

    !> Instance.
    class(MxlSocketComm), intent(in) :: this

    !> Time step in atomic units.
    real(dp) :: initDt

    initDt = this%initDt

  end function getInitDt

end module dftbp_io_mxlsocket
