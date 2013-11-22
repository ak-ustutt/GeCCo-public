*----------------------------------------------------------------------*
      subroutine da_vecsum(ffvecr,idxvecr,
     &                     ffvec1,idxvec1,xfac1,
     &                     ffvec2,idxvec2,xfac2,
     &                     lenvec,xbuf1,xbuf2,lenbuf)
*----------------------------------------------------------------------*
*
*     vecr = xfac1*vec1 + xfac2*vec2
*
*     files/start-records may be identical
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 0
      
      type(filinf) ::
     &     ffvecr, ffvec1, ffvec2
      integer, intent(in) ::
     &     idxvecr,
     &     idxvec1,
     &     idxvec2,lenvec, lenbuf
      real(8), intent(in) ::
     &     xfac1, xfac2
      real(8), intent(inout) ::
     &     xbuf1(lenbuf), xbuf2(lenbuf)

      integer ::
     &     lblk, luvec1, luvec2, luvecr,
     &     nrecs, len_rest, nrecbuf, nbatch, nrec_lbat, len_lbat,
     &     irecst1, irecst2, irecstr, lenbat, ibatch

      if (ntest.ne.0) then
        write(lulog,*) '---------------------'
        write(lulog,*) ' info from da_vecsum'
        write(lulog,*) '---------------------'
      end if
      if (ntest.ge.10) then
        write(lulog,*) 'ffvecr, idxvecr:  ',trim(ffvecr%name), idxvecr
        write(lulog,*) 'ffvec1, idxvec1:  ',trim(ffvec1%name), idxvec1
        write(lulog,*) 'ffvec2, idxvec2:  ',trim(ffvec2%name), idxvec2
        write(lulog,*) 'xfac1, xfac2: ',xfac1,xfac2
      end if
      if (ntest.ge.100) then
        write(lulog,*) 'lenvec,lblk,lenbuf: ',lenvec,lenbuf
      end if

      luvecr = ffvecr%unit
      luvec1 = ffvec1%unit
      luvec2 = ffvec2%unit
      lblk = ffvecr%reclen
      if (lblk.ne.ffvec1%reclen.or.
     &    lblk.ne.ffvec2%reclen) then
        write(lulog,*) 'Incompatible block-lengthes: ',
     &       lblk,ffvec1%reclen,ffvec2%reclen
        call quit(1,'da_vecsum','Incompatible block-lengthes')
      end if

      if (lenbuf.lt.min(lenvec,lblk)) then
        write(lulog,*) 'Insufficient buffer size!'
        write(lulog,*) ' buffer length = ',lenbuf,' blocksize = ',lblk
        call quit(1,'da_vecsum','Insufficient buffer size!')
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
      
      irecst1 = (idxvec1-1)*nrecs+1
      irecst2 = (idxvec2-1)*nrecs+1
      irecstr = (idxvecr-1)*nrecs+1
      do ibatch = 1, nbatch
        if (ibatch.eq.nbatch) lenbat = len_lbat

        if (xfac1.ne.0d0)
     &       call get_vec_a(ffvec1,xbuf1,lenbat,irecst1)
        if (xfac2.ne.0d0)
     &       call get_vec_a(ffvec2,xbuf2,lenbat,irecst2)

        if (xfac1.eq.0d0.and.xfac2.eq.0d0) then
          xbuf1(1:lenbat) = 0d0
        else if (xfac1.ne.0d0.and.xfac2.ne.0d0) then
          xbuf1(1:lenbat) = xfac1*xbuf1(1:lenbat) +
     &                         xfac2*xbuf2(1:lenbat)
        else if (xfac1.eq.0d0) then
          xbuf2(1:lenbat) = xfac2*xbuf2(1:lenbat)
        else if (xfac2.eq.0d0) then
          xbuf1(1:lenbat) = xfac1*xbuf1(1:lenbat)
        end if

        if (xfac1.eq.0d0.and.xfac2.ne.0d0) then
          call put_vec_a(ffvecr,xbuf2,lenbat,irecstr)
        else
          call put_vec_a(ffvecr,xbuf1,lenbat,irecstr)
        end if

        irecst1 = irecst1 + nrecbuf
        irecst2 = irecst2 + nrecbuf
        irecstr = irecstr + nrecbuf

      end do

      if (ntest.ge.100) then
        write(lulog,*) 'vecsum needed ',nbatch,' batches'
      end if
      
      return
      end
