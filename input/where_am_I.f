*----------------------------------------------------------------------*
      subroutine where_am_I(l_infile,name_infile,env_type)
*----------------------------------------------------------------------*
*     find out, which type of environment to use
*     very initial version
*----------------------------------------------------------------------*

      implicit none

      logical, intent(inout) ::
     &     l_infile
      character, intent(inout) ::
     &     name_infile*(*), env_type*(*)

      logical ::
     &     l_exist, ok

      logical, external ::
     &     is_dalton64

      ! we start with DALTON
      ! we need (at least): SIRIFC, MOTWOINT
      inquire(file='SIRIFC',exist=l_exist)
      ok = l_exist

      if (ok) inquire(file='MO_G',exist=l_exist)

      if (ok.and.l_exist) then
        env_type='DALTON_SPECIAL'
      else
        if (ok) inquire(file='MOTWOINT',exist=l_exist)
        ok = ok.and.l_exist
        if (ok) then
          if (is_dalton64()) then
            env_type='DALTON64 '
          else
            env_type='DALTON '
          end if
        end if
      end if

      ! try GAMESS
      if (.not.ok) then
        inquire(file='DICTNRY',exist=l_exist)
        ok = l_exist
        if (ok) then
          inquire(file='MOINTS',exist=l_exist)
          ok = l_exist
          if (ok) env_type='GAMESS '
        end if
      end if

      ! try new MOLPRO interface file
      if (.not.ok) then
        inquire(file='mpro_gecco_ifc.dat',exist=l_exist)
        ok = l_exist
        if (ok) env_type='MOLPRO_IFC'
      end if

      ! try MOLPRO fci interface
      if (.not.ok) then
        inquire(file='FCIDUMP',exist=l_exist)
        ok = l_exist
        if (ok) env_type='MOLPRO_DUMP'
      end if

      ! future: try other possibilities here ...
      if (.not.ok) call quit(0,'where_am_I',
     &     'did not find proper environment')

      ! no input file given: try standard name
      if (.not.l_infile) then
        inquire(file='gecco.inp',exist=l_exist)
        if (.not.l_exist) call quit(0,'where_am_I',
     &       'did not find any input file')
        name_infile = 'gecco.inp'
      end if

      return
      end
