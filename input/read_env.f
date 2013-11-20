      subroutine read_env(env_type,orb_info)

      implicit none
      include 'stdunit.h'
      include 'def_orbinf.h'

      character, intent(in) ::
     &     env_type*(*)
      type(orbinf), intent(out) ::
     &     orb_info

      if (iprlvl.ge.1)
     &     write(lulog,*) 'Reading data from environment ....'
      if (iprlvl.ge.2)
     &     write(lulog,*) 'Environment type: ',
     &     trim(env_type)

      select case(trim(env_type))
      case ('dalton','DALTON','dalton_special','DALTON_SPECIAL')
        call read_env_dalton(orb_info)
      case ('dalton64','DALTON64')
        call read_env_dalton64(orb_info)
      case ('gamess','GAMESS')
        call read_env_gamess(orb_info)
      case ('molpro_dump','MOLPRO_DUMP')
        call read_env_molpro_dump(orb_info)
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
