*----------------------------------------------------------------------*
      subroutine import_cmo(ffcmo,cmo_type,env_type,orb_info)
*----------------------------------------------------------------------*
*     import MO-coefficients from environment
*     cmo_type:
*       on entry -- what we want
*       on exit  -- we we get
*       <0 -- anything (on entry), nothing found (on exit)
*        1 -- SAO basis
*        2 -- AO  basis
*        3 -- CAO basis
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'

      integer, intent(inout) ::
     &     cmo_type
      character, intent(in) ::
     &     env_type*(*)
      type(filinf), intent(inout) ::
     &     ffcmo
      type(orbinf), intent(in) ::
     &     orb_info

      select case(trim(env_type))
      case ('dalton','DALTON')
        call import_cmo_dalton(ffcmo,orb_info)
        cmo_type = 1  ! DALTON provides SAO basis only
      case ('intern','INTERN')
        call quit(1,'import_cmo','type INTERN not implemented')
      case ('aces2','ACES2')
        call quit(1,'import_cmo','type ACES2 not implemented')
      case ('tmole','TMOLE')
        call quit(1,'import_cmo','type TMOLE not implemented')
      case default
        call quit(1,'import_cmo','unknown type '//trim(env_type))
      end select

      end
