*----------------------------------------------------------------------*
      subroutine put_vec(fhand,buf,idxst,idxnd)
*----------------------------------------------------------------------*
*     write vector elements idxst to idxnd (in buf) to file fhand
*----------------------------------------------------------------------*
      implicit none

      include 'def_filinf.h'
      include 'stdunit.h'

      type(filinf), intent(in) ::
     &     fhand
      real(8), intent(in) ::
     &     buf(*)
      integer, intent(in) ::
     &     idxst, idxnd

      integer ::
     &     un, lenr, irecst, irecnd, ioffrec1, idxrecl, nwrite, irec,
     &     idxbfst, idxbfnd, idx
      real(8) ::
     &     xdum
      real(8), allocatable ::
     &     recbuf(:)
*----------------------------------------------------------------------*

      ! do not use structure elements directly
      un = fhand%unit
      lenr = fhand%reclen

c dbg
      write(luout,*)'put_vec unit',un
c dbg

      ! first and last record to read from
      irecst = (idxst-1)/lenr+1
      irecnd = (idxnd-1)/lenr+1

      ! offset in first record
      ioffrec1 = idxst-1 - (irecst-1)*lenr
      ! last index in last record
      idxrecl = idxnd - (irecnd-1)*lenr

      ! unless we hit the block-size exactly,
      ! we will have to handle the "borders":
      if (ioffrec1.gt.0.or.idxrecl.lt.lenr)
     &     allocate(recbuf(lenr))

      if (irecst.eq.irecnd) then
        ! special case -- only one record to read:
        if (ioffrec1.gt.0.or.idxrecl.lt.lenr) then
          nwrite = idxrecl-ioffrec1
          ! we have to read in the record first
          read(un,rec=irecst,err=10) recbuf(1:lenr)
          goto 20
          ! error exception for uninitialized record
 10       recbuf(1:ioffrec1) = 0d0
          recbuf(idxrecl+1:lenr) = 0d0
 20       recbuf(ioffrec1+1:ioffrec1+nwrite) = buf(1:nwrite)
          write(un,rec=irecst) recbuf(1:lenr)
        else
          ! exactly one record to write
          write(un,rec=irecst) buf(1:lenr)
        end if
      else
        ! first record
        if (ioffrec1.eq.0) then
          write(un,rec=irecst) buf(1:lenr)
          idxbfst = lenr+1
        else
          ! handle "border"
          read(un,rec=irecst,err=11) recbuf(1:ioffrec1)
          goto 21
          ! handle uninitilized record
 11       recbuf(1:ioffrec1) = 0d0          
 21       nwrite = lenr-ioffrec1
          recbuf(ioffrec1+1:ioffrec1+nwrite) = buf(1:nwrite)
          write(un,rec=irecst) recbuf(1:lenr)
          idxbfst = nwrite+1
        end if

        ! 2nd to (last-1)st record
        do irec = irecst+1, irecnd-1
          idxbfnd = idxbfst+lenr-1
          write(un,rec=irec) buf(idxbfst:idxbfnd)
          idxbfst = idxbfnd+1
        end do

        ! last record
        idxbfnd = idxbfst+idxrecl-1
        if (idxrecl.eq.lenr) then
          write(un,rec=irecnd) buf(idxbfst:idxbfnd)
        else
          ! handle "border"
          read(un,rec=irecnd,err=12) recbuf(1:lenr)
          goto 22
          ! handle uninitialized record
 12       recbuf(idxrecl+1:lenr) = 0d0
 22       recbuf(1:idxrecl) = buf(idxbfst:idxbfnd)
          write(un,rec=irecnd) recbuf(1:lenr)
        end if
      end if

      if (ioffrec1.gt.0.or.idxrecl.lt.lenr)
     &     deallocate(recbuf)

      return
      end
*----------------------------------------------------------------------*
