************************************************************************
      subroutine put_vec_a(ffvec,vec,lenvec,irecst)
************************************************************************
*     alternative version of put_vec
************************************************************************
      
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer,parameter::
     &     ntest=00
      type(filinf), intent(in) ::
     &     ffvec

      integer, intent(in) ::
     &     lenvec, irecst
      real(8), intent(in) ::
     &     vec(lenvec)

      integer ::
     &     luvec, lblk, nrec, len_last, len, irec, irecnd, idxst, idxnd

      if(ntest.ge.20)then
         write (lulog,'(1X,"depositing on file ",A20)')ffvec%name
      end if
      lblk = ffvec%reclen
      luvec = ffvec%unit
      nrec = lenvec/lblk
      len_last = mod(lenvec,lblk)
      if (len_last.gt.0) nrec = nrec+1
      if (len_last.eq.0) len_last = lblk

      len = lblk
      idxst = 1
      irecnd = irecst-1 + nrec
      do irec = irecst, irecnd
        if (irec.eq.irecnd) len = len_last
        idxnd = idxst-1 + len
        write(luvec,rec=irec) vec(idxst:idxnd)
        idxst = idxnd + 1
      end do

      return
      end
