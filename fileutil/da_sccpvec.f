*----------------------------------------------------------------------*
      subroutine da_sccpvec(ffvecin,idxvecin,
     &                      ffvecout,idxvecout,
     &                      xfac,lenvec,xbuf,lenbuf)
*----------------------------------------------------------------------*
*
*     copy and/or scale vector on disc
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00
      
      type(filinf), intent(in) ::
     &     ffvecin, ffvecout
      integer, intent(in) ::
     &     idxvecin, idxvecout, lenvec, lenbuf
      real(8), intent(in) ::
     &     xfac
      real(8), intent(inout) ::
     &     xbuf(lenbuf)

      integer ::
     &     luvecin, luvecout, lblk,
     &     nrecs, len_rest, nrecbuf, nbatch, nrec_lbat, len_lbat,
     &     irecstin, irecstout, lenbat, ibatch

      if (ntest.ne.0) then
        write(lulog,*) '----------------------'
        write(lulog,*) ' info from da_sccpvec'
        write(lulog,*) '----------------------'
      end if
      if (ntest.ge.10) then
        write(lulog,*) 'ffvecin, idxvecin:  ',
     &       trim(ffvecin%name),idxvecin
        write(lulog,*) 'ffvecout,idxvecout: ',
     &       trim(ffvecout%name),idxvecout
        write(lulog,*) 'xfac: ',xfac
      end if
      if (ntest.ge.100) then
        write(lulog,*) 'lenvec,lenbuf: ',lenvec,lenbuf
      end if

      luvecin  = ffvecin%unit
      luvecout = ffvecout%unit
      lblk = ffvecin%reclen
      if (lblk.ne.ffvecout%reclen) then
        write(lulog,*) 'Incompatible block-lengthes: ',
     &       lblk,ffvecout%reclen
        call quit(1,'da_sccpvec','Incompatible block-lengthes')
      end if

      if (lenbuf.lt.min(lenvec,lblk)) then
        write(lulog,*) 'Insufficient buffer size!'
        write(lulog,*) ' buffer length = ',lenbuf,' blocksize = ',lblk
        call quit(1,'da_sccpvec','Insufficient buffer size!')
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
        if (nbatch.gt.1.and.  !we will need at least one batch
     &      nrec_lbat.eq.1.and.lenbat+len_rest.le.lenbuf) then
          len_lbat = lenbat+len_rest
          nrec_lbat = nrecbuf+1
        else
          nbatch = nbatch+1
          len_lbat = (nrec_lbat-1)*lblk+len_rest
        end if
      end if
      
      irecstin = (idxvecin-1)*nrecs+1
      irecstout = (idxvecout-1)*nrecs+1
      do ibatch = 1, nbatch
        if (ibatch.eq.nbatch) lenbat = len_lbat

        call get_vec_a(ffvecin,xbuf,lenbat,irecstin)

        if (xfac.ne.1d0)
     &       xbuf(1:lenbat) = xfac*xbuf(1:lenbat)

        call put_vec_a(ffvecout,xbuf,lenbat,irecstout)

        irecstin = irecstin + nrecbuf
        irecstout = irecstout + nrecbuf

      end do

      if (ntest.ge.10) then
        write(lulog,*) 'scale/copy needed ',nbatch,' batches'
      end if
      
      return
      end
