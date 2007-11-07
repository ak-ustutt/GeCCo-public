*----------------------------------------------------------------------*
      subroutine da_listvec(ffvec,idxvec,lenvec,imod,xbuf,lenbuf)
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      type(filinf), intent(in) ::
     &     ffvec
      integer, intent(in) ::
     &     idxvec, lenvec, lenbuf, imod
      real(8), intent(inout) ::
     &     xbuf(lenbuf)

      integer ::
     &     luvec, lblk,
     &     nrecs, len_rest, nrecbuf, nbatch, nrec_lbat, len_lbat,
     &     irecst, lenbat, ibatch
      real(8) ::
     &     xnrm

      luvec  = ffvec%unit
      lblk = ffvec%reclen
      if (lenbuf.lt.min(lenvec,lblk)) then
        write(luout,*) 'Insufficient buffer size!'
        write(luout,*) ' buffer length = ',lenbuf,' blocksize = ',lblk
        call quit(1,'da_listvec','Insufficient buffer size!')
      end if

      ! number of records per vector
      nrecs = lenvec/lblk
      len_rest = mod(lenvec,lblk)
      if (len_rest.gt.0) nrecs = nrecs+1

      ! number of records that fit in buffer:
      nrecbuf = max(1,lenbuf/lblk)
      lenbat = nrecbuf*lblk
      
      ! number of batches to process
      nbatch = nrecs/nrecbuf
      nrec_lbat = mod(nrecs,nrecbuf)
      len_lbat = len_rest
      if (nrec_lbat.gt.0) then
        if (nrec_lbat.eq.1.and.lenbat+len_rest.le.lenbuf) then
          len_lbat = lenbat+len_rest
          nrec_lbat = nrecbuf+1
        else
          nbatch = nbatch+1
          len_lbat = (nrec_lbat-1)*lblk+len_rest
        end if
      end if

      irecst = (idxvec-1)*nrecs+1
      do ibatch = 1, nbatch
        if (ibatch.eq.nbatch) lenbat = len_lbat

        call get_vec_a(ffvec,xbuf,lenbat,irecst)

        call da_print_vec(xbuf,lenbat,irecst,imod,lblk)

        irecst = irecst + nrecbuf

      end do

      return

      contains
*----------------------------------------------------------------------*
      subroutine da_print_vec(vec,lenvec,irecst,imod,lblk)
*----------------------------------------------------------------------*
      
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     lenvec, irecst, lblk, imod
      real(8), intent(in) ::
     &     vec(lenvec)

      integer ::
     &     nrec, len_last, len, irec, irecnd, idxst, idxnd
      real(8) ::
     &     xnrm

      real(8), external ::
     &     ddot

      nrec = lenvec/lblk
      len_last = mod(lenvec,lblk)
      if (len_last.gt.0) nrec = nrec+1
      if (len_last.eq.0) len_last = lblk

      len = lblk
      idxst = 1
      irecnd = irecst-1 + nrec
      do irec = irecst, irecnd
        if (irec.eq.irecnd) len = len_last
c        idxnd = idxst-1 + len

        xnrm = sqrt(ddot(len,vec(idxst),1,vec(idxst),1))
        write(luout,'(x,"rec=",i10,2x,"norm=",g20.10)') irec,xnrm
        if (imod.gt.0) then
          call wrtmat2(vec(idxst),1,len,1,len)
        end if

        idxst = idxst + len
      end do

      return
      end subroutine

      end
