      subroutine read_input(ffinput)

      use parse_input
      implicit none
      include 'stdunit.h'
      include 'def_filinf.h'
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

      call set_keywords()

      if (.not.file_exists(ffinput))
     &     call quit(0,'read_input','file does not exist - "'//
     &     trim(ffinput%name)//'"')

      call file_open(ffinput)

      call list_file(lulog,ffinput%unit)
      if (luout.ne.lulog) call list_file(luout,ffinput%unit)

      call keyword_parse(ffinput%unit)

      call file_close_keep(ffinput)

      return

      end
