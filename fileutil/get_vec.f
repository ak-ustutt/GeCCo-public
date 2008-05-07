*----------------------------------------------------------------------*
      subroutine get_vec(fhand,buf,idxst,idxnd)
*----------------------------------------------------------------------*
      implicit none

      include 'def_filinf.h'

      type(filinf), intent(in) ::
     &     fhand
      real(8), intent(out) ::
     &     buf(*)
      integer, intent(in) ::
     &     idxst, idxnd

      integer ::
     &     un, lenr, irecst, irecnd, ioffrec1, idxrecl, nread, irec,
     &     idxbfst, idxbfnd, idx
      real(8) ::
     &     xdum
*----------------------------------------------------------------------*

      ! do not use structure elements directly
      un = fhand%unit
      lenr = fhand%reclen

      ! first and last record to read from
      irecst = (idxst-1)/lenr+1
      irecnd = (idxnd-1)/lenr+1

      ! offset in first record
      ioffrec1 = idxst-1 - (irecst-1)*lenr
      ! last index in last record
      idxrecl = idxnd - (irecnd-1)*lenr

      if (irecst.eq.irecnd) then
        ! special case -- only one record to read:
        nread = idxrecl-ioffrec1
        read(un,rec=irecst) (xdum,idx=1,ioffrec1),buf(1:nread)
      else
        ! first record
        nread = lenr-ioffrec1
        read(un,rec=irecst) (xdum,idx=1,ioffrec1),buf(1:nread)
        idxbfst = nread+1
        ! 2nd to (last-1)st record
        do irec = irecst+1, irecnd-1
          idxbfnd = idxbfst+lenr-1
          read(un,rec=irec) buf(idxbfst:idxbfnd)
          idxbfst = idxbfnd+1
        end do
        ! last record
        idxbfnd = idxbfst+idxrecl-1
        read(un,rec=irecnd) buf(idxbfst:idxbfnd)
      end if

      return
      end
*----------------------------------------------------------------------*
