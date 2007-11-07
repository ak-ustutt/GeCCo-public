*----------------------------------------------------------------------*
      subroutine oneprop_ao(ffdao,dens,
     &     env_type,orb_info)
*----------------------------------------------------------------------*
*     contract 1-particle density matrix (on ffdao) with available
*     integrals from the environment
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'

      type(filinf), intent(inout) ::
     &     ffdao
      type(operator), intent(in) ::
     &     dens
      character, intent(in) ::
     &     env_type*(*)
      type(orbinf), intent(in) ::
     &     orb_info

      select case(trim(env_type))
      case ('dalton','DALTON')
        call oneprop_ao_dalton(ffdao,dens,orb_info)
      case ('intern','INTERN')
        call quit(1,'oneprop_ao','type INTERN not implemented')
      case ('aces2','ACES2')
        call quit(1,'oneprop_ao','type ACES2 not implemented')
      case ('tmole','TMOLE')
        call quit(1,'oneprop_ao','type TMOLE not implemented')
      case default
        call quit(1,'oneprop_ao','unknown type '//trim(env_type))
      end select

      return
      end

