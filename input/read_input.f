      subroutine read_input(ffinput)

      use keyword_trees, only:reg_import
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

      call reg_import(get_keyword_file())

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
!!    returns name of the keyword_file
*----------------------------------------------------------------------*
      function get_keyword_file() result(file_name)
*----------------------------------------------------------------------*
      implicit none
      integer,parameter::
     &     ntest= 00
      character(len=16),parameter ::
     &     i_am="get_keyword_file"
      character(len=256)::
     &     path_name,file_name
      integer ::
     &     len


      call get_environment_variable( "GECCO_DIR", value=path_name,
     &     length = len)

      if (len.EQ.0)
     &     call quit(0,i_am,
     &     "Please, set the GECCO_DIR environment variable.")
      if (len .gt.(256-file_loc_len)  )
     &     call quit(0,i_am,
     &     "GECCO_DIR to long, cannot set keyword_registry")
      file_name=trim(path_name)//rel_file_loc

      end function
      end
