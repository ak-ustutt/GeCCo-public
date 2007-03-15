      subroutine read_env(env_type,orb_info)

      implicit none
      include 'stdunit.h'
      include 'def_orbinf.h'

      character, intent(in) ::
     &     env_type*(*)
      type(orbinf), intent(out) ::
     &     orb_info

      if (iprlvl.ge.1)
     &     write(luout,*) 'Reading data from environment ....'
      if (iprlvl.ge.2)
     &     write(luout,*) 'Environment type: ',
     &     trim(env_type)

      select case(trim(env_type))
      case ('dalton','DALTON')
        call read_env_dalton(orb_info)
      case ('intern','INTERN')
        call quit(1,'read_env','type INTERN not implemented')
      case ('aces2','ACES2')
        call quit(1,'read_env','type ACES2 not implemented')
      case ('tmole','TMOLE')
        call quit(1,'read_env','type TMOLE not implemented')
      case default
        call quit(1,'read_env','unknown type '//trim(env_type))
      end select

      return
      end
