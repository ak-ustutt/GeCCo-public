
*----------------------------------------------------------------------*
      real(8) function da_fndmnx(ffvec,idxvec,minmax,
     &                           lenvec,xbuf,lenbuf)
*----------------------------------------------------------------------*
*
*                                 minmax:
*     result = min(abs(vec(i)))     1
*     result = max(abs(vec(i)))     2
*     result = min(vec(i))         -1
*     result = max(vec(i))         -2
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00
      
      type(filinf), intent(in) ::
     &     ffvec
      integer, intent(in) ::
     &     idxvec, minmax,
     &     lenvec, lenbuf
      real(8), intent(inout) ::
     &     xbuf(lenbuf)

      integer ::
     &     lblk, luvec,
     &     nrecs, len_rest, nrecbuf, nbatch, nrec_lbat, len_lbat,
     &     irecst, lenbat, ibatch
      real(8) ::
     &     res

      real(8), external ::
     &     fndmnx

      if (minmax.lt.-2.or.minmax.eq.0.or.minmax.gt.2) then
        write(luout,*) 'da_fndmnx: illegal value for minmax: ',minmax
        call quit(1,'da_fndmnx','Illegal value for minmax')
      end if

      luvec = ffvec%unit
      lblk = ffvec%reclen
      if (lenbuf.lt.min(lenvec,lblk)) then
        write(luout,*) 'Insufficient buffer size!'
        write(luout,*) ' buffer length = ',lenbuf,' blocksize = ',lblk
        call quit(1,'da_fndmnx','Insufficient buffer size!')
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

      if (ntest.ge.100) then
        write(luout,*) 'nrecs,nrecbuf,nbatch: ',nrecs,nrecbuf,nbatch
      end if
      
      irecst = (idxvec-1)*nrecs+1
      res = 0d0
      do ibatch = 1, nbatch
        if (ibatch.eq.nbatch) lenbat = len_lbat
        call get_vec_a(ffvec,xbuf,lenbat,irecst)

        if (ibatch.eq.1.and.minmax.gt.0) res = abs(xbuf(1))
        if (ibatch.eq.1.and.minmax.lt.0) res = xbuf(1)
        if (abs(minmax).eq.1) res = min(res,fndmnx(xbuf,lenbat,minmax))
        if (abs(minmax).eq.2) res = max(res,fndmnx(xbuf,lenbat,minmax))

        irecst = irecst + nrecbuf

      end do

      da_fndmnx = res

      if (ntest.ge.100) then
        write(luout,*) 'da_fndmnx needed ',nbatch,' batches'
      end if
      
      return
      end
