
*----------------------------------------------------------------------*
      subroutine file_open(fhand)
*----------------------------------------------------------------------*
*     open file with handle fhand (unit number is set)
*----------------------------------------------------------------------*

      implicit none

      include 'def_filinf.h'
      include 'ioparam.h'

      integer, external ::
     &     iopen_nud, iopen_nus, iopen_nfs, iopen_nustr

      type(filinf), intent(inout) ::
     &     fhand


      if (fhand%buffered) return

      if (fhand%unit.gt.0)
     &     call quit(1,'file_open',
     &     'file is already open: '//trim(fhand%name))

      if (fhand%type.eq.ftyp_da_unf) then
        fhand%unit = iopen_nud(fhand%name,fhand%reclen)
      else if (fhand%type.eq.ftyp_sq_unf) then
        fhand%unit = iopen_nus(fhand%name)
      else if (fhand%type.eq.ftyp_sq_frm) then
        fhand%unit = iopen_nfs(fhand%name)
      else if (fhand%type.eq.ftyp_st_unf) then
        fhand%unit = iopen_nustr(fhand%name)
      else
        call quit(1,'file_open','illegal file type')
      end if

      return
      end
