*----------------------------------------------------------------------*
      subroutine file_init(fhand,name,type,reclen)
*----------------------------------------------------------------------*
*     a little auxiliary routine to init the hande fhand
*     reclen: in real(8) words (system specific factors are handled
*             automatically in the inner routines)
*     type: see def_filinf.h
*----------------------------------------------------------------------*
      implicit none

      include 'def_filinf.h'

      type(filinf), intent(inout) ::
     &     fhand
      integer, intent(in) ::
     &     type, reclen
      character, intent(in) ::
     &     name*(*)

      integer ::
     &     len
      
      len = len_trim(name)
      if (len.gt.maxfilnam)
     &     call quit(1,'file_init','filename too long')
      if (type.lt.1.or.type.gt.4)
     &     call quit(1,'file_init','illegal file type')
      if (type.eq.1.and.reclen.le.0)
     &     call quit(1,'file_init','zero or negative record length')

      fhand%name(1:maxfilnam) = ' '
      fhand%name(1:len) = name
      fhand%unit = -1
      fhand%type = type
      if (type.eq.1) then
        fhand%reclen = reclen
      else
        fhand%reclen = -1
      end if
      fhand%buffered = .false.  ! default behaviour
      fhand%recoff_superfile = 0 ! unused
      fhand%length_of_record = 0 ! super-records, init to zero
      fhand%current_record = 1   ! current super-record settings:
      fhand%active_records = (/1,1/)   

      return

      end
