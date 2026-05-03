!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> MaxwellLink socket driver for real-time TD-DFTB dynamics.
submodule (dftbp_timedep_timeprop) dftbp_timedep_mxlrtdynamics
  use dftbp_io_mxlcommon, only : buildMxlExtraJson, checkMxlInit, makeMxlSocketCommInput
  use dftbp_io_mxlsocket, only : MxlSocketComm, MxlSocketComm_init, MxlSocketCommInp
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_bcast
#:endif
  implicit none

contains

  module procedure runMxlSocketDynamics

    type(MxlSocketComm) :: mxlSocket
    type(MxlSocketCommInp) :: mxlInput
    character(lc) :: mxlExtra, mxlMessage
    logical :: mxlHaveResult, mxlReceivedInit, mxlStop
    real(dp) :: mxlDipoleCurrent(3), mxlDipoleEnd(3), mxlDipoleInitial(3)
    real(dp) :: mxlDipoleStart(3)
    real(dp) :: mxlEnergy, mxlEnergyEnd, mxlEnergyStart, mxlField(3), mxlInitDt
    integer :: iStep, mxlReceivedMoleculeId
    real(dp) :: mxlSource(3), mxlTime, timeElec

    call makeMxlSocketCommInput(this%mxlSocketInput, mxlInput)
  #:if WITH_MPI
    if (env%mpi%tGlobalLead) then
  #:endif
    call MxlSocketComm_init(mxlSocket, mxlInput)
  #:if WITH_MPI
    end if
  #:endif

    mxlHaveResult = .false.
    mxlDipoleInitial(:) = 0.0_dp
    if (this%mxlSocketInput%resetDipole) then
      call getMxlDipole(this, mxlDipoleInitial)
    end if

    do iStep = 1, this%nSteps

    #:if WITH_MPI
      if (env%mpi%tGlobalLead) then
    #:endif
      call mxlSocket%receiveField(mxlHaveResult, mxlField, mxlStop, mxlReceivedInit)
      if (mxlReceivedInit) then
        mxlInitDt = mxlSocket%getInitDt()
        mxlReceivedMoleculeId = mxlSocket%getMoleculeId()
      else
        mxlInitDt = -1.0_dp
        mxlReceivedMoleculeId = -1
      end if
    #:if WITH_MPI
      else
        mxlField(:) = 0.0_dp
        mxlStop = .false.
        mxlReceivedInit = .false.
        mxlInitDt = -1.0_dp
        mxlReceivedMoleculeId = -1
      end if
      call mpifx_bcast(env%mpi%globalComm, mxlField)
      call mpifx_bcast(env%mpi%globalComm, mxlStop)
      call mpifx_bcast(env%mpi%globalComm, mxlReceivedInit)
      call mpifx_bcast(env%mpi%globalComm, mxlInitDt)
      call mpifx_bcast(env%mpi%globalComm, mxlReceivedMoleculeId)
    #:endif
      if (mxlReceivedInit) then
        call checkMxlInit(mxlInitDt, mxlReceivedMoleculeId, this%dt,&
            & this%mxlSocketInput%moleculeId, "ElectronDynamics", mxlMessage)
        if (len_trim(mxlMessage) > 0) then
          @:RAISE_ERROR(errStatus, -1, trim(mxlMessage))
        end if
      end if
      if (mxlStop) then
        exit
      end if

      call runOneMxlTdStep(this, boundaryCond, iStep, coord, orb, neighbourList, nNeighbourSK,&
          & symNeighbourList, nNeighbourCamSym, iSquare, iSparseStart, img2CentCell, skHamCont,&
          & skOverCont, ints, env, coordAll, q0, referenceN0, spinW, tDualSpinOrbit, xi,&
          & thirdOrd, dftbU, onSiteElements, refExtPot, solvation, eFieldScaling, hybridXc,&
          & repulsive, iAtInCentralRegion, tFixEf, Ef, electronicSolver, qDepExtPot, mxlField,&
          & mxlDipoleInitial, mxlDipoleStart, mxlDipoleEnd, mxlEnergyStart, mxlEnergyEnd,&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)

      mxlTime = this%time
      mxlEnergy = 0.5_dp * (mxlEnergyStart + mxlEnergyEnd)
      mxlDipoleCurrent(:) = 0.5_dp * (mxlDipoleStart(:) + mxlDipoleEnd(:))
      mxlSource(:) = (mxlDipoleEnd(:) - mxlDipoleStart(:)) / this%dt

      call buildMxlExtraJson(mxlTime, mxlEnergy, this%energyKin, mxlDipoleCurrent,&
          & mxlDipoleCurrent, mxlExtra)

      mxlHaveResult = .true.
    #:if WITH_MPI
      if (env%mpi%tGlobalLead) then
    #:endif
      call mxlSocket%sendSource(mxlEnergy, mxlSource, trim(mxlExtra), mxlStop)
    #:if WITH_MPI
      else
        mxlStop = .false.
      end if
      call mpifx_bcast(env%mpi%globalComm, mxlStop)
    #:endif
      mxlHaveResult = .false.
      this%tdFieldIsSet = .false.
      if (mxlStop) then
        exit
      end if

      if (mod(iStep, max(this%nSteps / 10, 1)) == 0) then
        call loopTime%stop()
        timeElec  = loopTime%getWallClockTime()
        write(stdOut, "(A,2x,I6,2(2x,A,F10.6))") 'Step ', iStep, 'elapsed loop time: ',&
            & timeElec, 'average time per loop ', timeElec / (iStep + 1)
      end if

    end do

  #:if WITH_MPI
    if (env%mpi%tGlobalLead) then
  #:endif
    call mxlSocket%shutdown()
  #:if WITH_MPI
    end if
  #:endif

  end procedure runMxlSocketDynamics

end submodule dftbp_timedep_mxlrtdynamics
