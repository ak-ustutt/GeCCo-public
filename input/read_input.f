      subroutine read_input(ffinput)

      use keyword_trees, only:create_keyword_trees
      use parse_input,only : inp_parse
      implicit none
      include 'stdunit.h'
      include 'def_filinf.h'
      integer,parameter::
     &     file_loc_len=28
      character(len=file_loc_len),parameter ::
     &     rel_file_loc="/data/keyword_registry3.xml"
c does not work, why?
c      include 'ifc_fileutil.h'
 
      type(filinf), intent(inout) ::
     &     ffinput

c work around for problem with interface file
      logical, external :: file_exists

      write(lulog,*) 'Reading input file ....'
      write(lulog,*) 'Input file: ',
     &     trim(ffinput%name)
      if (iprlvl.ge.2.and.luout.ne.lulog)
     &     write(luout,*) 'Input file: ',
     &     trim(ffinput%name)

      call create_keyword_trees(get_keyword_file())

      if (.not.file_exists(ffinput))
     &     call quit(0,'read_input','file does not exist - "'//
     &     trim(ffinput%name)//'"')

      call file_open(ffinput)

      call list_file(lulog,ffinput%unit)
      if (luout.ne.lulog) call list_file(luout,ffinput%unit)

      call inp_parse(ffinput%unit)

      call file_close_keep(ffinput)

      return
      contains

*----------------------------------------------------------------------*
!>    returns name of the keyword_file
!!    reads an environment variable and appends a specified relative path
!!    variable and path are named in par_keyreg.h
*----------------------------------------------------------------------*
      function get_keyword_file() result(file_name)
*----------------------------------------------------------------------*
      implicit none
      include "par_keyreg.h"
      integer,parameter::
     &     ntest= 00
      character(len=16),parameter ::
     &     i_am="get_keyword_file"
      character(len=256)::
     &     path_name,file_name
      integer ::
     &     length,file_loc_len

      file_loc_len=len(keyreg_file_loc)

      call get_environment_variable( keyreg_variable, value=path_name,
     &     length = length)

      if (length.EQ.0)
     &     call quit(0,i_am,
     &     "Please, set the "//keyreg_variable
     &     //" environment variable.")

      if (length .gt.(256-file_loc_len)  )
     &     call quit(0,i_am,
     &     "GECCO_DIR to long, cannot set keyword_registry")

      file_name=trim(path_name)//keyreg_file_loc

      end function
      end
