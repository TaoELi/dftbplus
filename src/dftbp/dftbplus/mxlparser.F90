!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Parser helpers for MaxwellLink socket coupling.
module dftbp_dftbplus_mxlparser
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_unitconversion, only : energyUnits
  use dftbp_extlibs_xmlf90, only : char, fnode, getNodeName, string
  use dftbp_io_charmanip, only : unquote
  use dftbp_io_hsdutils, only : detailedError, getChild, getChildValue
  use dftbp_io_hsdutils2, only : convertUnitHsd
  use dftbp_io_mxlcommon, only : TMxlSocketInput
  use dftbp_md_mxlbomd, only : TMxlBomdInput, mxlBomdDerivTypes
  use dftbp_timedep_timeprop, only : TElecDynamicsInp
  use dftbp_type_linkedlist, only : asArray, destruct, init, len, TListReal
  implicit none

  private
  public :: parseMaxwellLinkBomdInput, parseMaxwellLinkElecInput


contains


  !> Read MaxwellLink socket address and shared options.
  subroutine parseMaxwellLinkSocketAddress(mxlNode, input)

    !> MaxwellLinkSocket node.
    type(fnode), pointer, intent(in) :: mxlNode

    !> Shared MaxwellLink socket settings.
    type(TMxlSocketInput), intent(out) :: input

    type(fnode), pointer :: child, fileNode, hostNode
    type(string) :: buffer, buffer2
    character(lc) :: mxlFile

    call getChild(mxlNode, "File", child=fileNode, requested=.false.)
    call getChild(mxlNode, "Host", child=hostNode, requested=.false.)
    if (associated(fileNode) .and. associated(hostNode)) then
      call detailedError(mxlNode, "Either Host or File, but not both, must be set for&
          & MaxwellLinkSocket")
    end if

    if (associated(fileNode)) then
      call getChildValue(fileNode, "", buffer2)
      mxlFile = unquote(char(buffer2))
      if (len_trim(mxlFile) == 0) then
        call detailedError(fileNode, "MaxwellLinkSocket File must not be empty")
      end if
      if (mxlFile(1:1) == "/") then
        input%host = trim(mxlFile)
      else
        call getChildValue(mxlNode, "Prefix", buffer, "/tmp/socketmxl_")
        input%host = trim(unquote(char(buffer))) // trim(mxlFile)
      end if
      input%port = 0
    else
      if (associated(hostNode)) then
        call getChildValue(hostNode, "", buffer2)
        input%host = unquote(char(buffer2))
      else
        input%host = "localhost"
      end if
      call getChildValue(mxlNode, "Port", input%port, 31415, child=child)
      if (input%port <= 0) then
        call detailedError(child, "Invalid MaxwellLinkSocket port number")
      end if
    end if

    call getChildValue(mxlNode, "Verbosity", input%verbosity, 0)
    call getChildValue(mxlNode, "MoleculeId", input%moleculeId, -1)
    call getChildValue(mxlNode, "ResetDipole", input%resetDipole, .false.)

  end subroutine parseMaxwellLinkSocketAddress


  !> Read MaxwellLink input for real-time TDDFTB/Ehrenfest dynamics.
  subroutine parseMaxwellLinkElecInput(mxlNode, input)

    !> MaxwellLinkSocket node.
    type(fnode), pointer, intent(in) :: mxlNode

    !> Electronic dynamics input.
    type(TElecDynamicsInp), intent(inout) :: input

    input%tMxlSocket = .true.
    call parseMaxwellLinkSocketAddress(mxlNode, input%mxlSocketInput)

  end subroutine parseMaxwellLinkElecInput


  !> Read MaxwellLink input for Born-Oppenheimer molecular dynamics.
  subroutine parseMaxwellLinkBomdInput(mxlNode, input)

    !> MaxwellLinkSocket node.
    type(fnode), pointer, intent(in) :: mxlNode

    !> MaxwellLink BOMD input.
    type(TMxlBomdInput), intent(inout) :: input

    type(fnode), pointer :: child, child2, value
    type(string) :: buffer, modifier
    type(TListReal) :: mxlCharges

    call parseMaxwellLinkSocketAddress(mxlNode, input%socket)

    call getChildValue(mxlNode, "DipoleDerivative", value, "MullikenCharges", child=child)
    call getNodeName(value, buffer)
    select case (char(buffer))
    case ("fixedcharges")
      input%derivType = mxlBomdDerivTypes%fixedCharges
      call init(mxlCharges)
      call getChildValue(mxlNode, "Charges", mxlCharges, child=child2)
      if (len(mxlCharges) < 1) then
        call detailedError(child2, "MaxwellLinkSocket Charges must not be empty")
      end if
      allocate(input%fixedCharges(len(mxlCharges)))
      call asArray(mxlCharges, input%fixedCharges)
      call destruct(mxlCharges)

    case ("mullikencharges")
      input%derivType = mxlBomdDerivTypes%mullikenCharges

    case ("bornchargesfile")
      call detailedError(child, "MaxwellLinkSocket DipoleDerivative = BornChargesFile has&
          & been removed; use BornChargesOnTheFly")

    case ("bornchargesonthefly")
      input%derivType = mxlBomdDerivTypes%bornChargesOnTheFly
      call getChildValue(mxlNode, "BornUpdateEvery", input%bornUpdateEvery, 1, child=child2)
      if (input%bornUpdateEvery < 1) then
        call detailedError(child2, "MaxwellLinkSocket BornUpdateEvery must be positive")
      end if
      call getChildValue(mxlNode, "MaxPerturbIter", input%maxPerturbIter, 100)
      call getChildValue(mxlNode, "PerturbSccTol", input%perturbSccTol, 1.0E-5_dp)
      call getChildValue(mxlNode, "PerturbDegenTol", input%perturbDegenTol, 1.0E-9_dp,&
          & modifier=modifier, child=child2)
      call convertUnitHsd(char(modifier), energyUnits, child2, input%perturbDegenTol)
      if (input%perturbDegenTol < epsilon(0.0_dp)) then
        call detailedError(child2, "MaxwellLinkSocket PerturbDegenTol must be above machine&
            & epsilon")
      end if

    case default
      call detailedError(child, "Invalid MaxwellLinkSocket DipoleDerivative '"&
          & // char(buffer) // "'")
    end select

  end subroutine parseMaxwellLinkBomdInput

end module dftbp_dftbplus_mxlparser
