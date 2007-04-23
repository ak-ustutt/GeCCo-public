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

      ! we start with DALTON
      ! we need (at least): SIRIFC, MOTWOINT
      inquire(file='SIRIFC',exist=l_exist)
      ok = l_exist
c      if (ok) inquire(file='AOONEINT',exist=l_exist)
c      ok = ok.and.l_exist
      if (ok) inquire(file='MOTWOINT',exist=l_exist)
      ok = ok.and.l_exist

      if (ok) env_type='DALTON '

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
