*----------------------------------------------------------------------*
      subroutine da_diavec(ffvecr,idxvecr,ist0_r,xfac1,
     &                     ffvec1,idxvec1,ist0_1,xfac2,
     &                     ffvec2,idxvec2,ist0_2,
     &                     xshift,xpot,
     &                     lenvec,xbuf1,xbuf2,lenbuf)
*----------------------------------------------------------------------*
*
*       vecr(i) = xfac1*vecr(i)+xfac2*vec1(i)*[vec2(i)+xshift]**xpot
*
*     files/start-records may be identical
*
*     special handling for xpot = -2,-1,0,1,2
*     caution: be sure to pass double-precision numbers as it else
*     may happen that 1.0 != 1.0d0
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 0
      
      type(filinf), intent(in) ::
     &     ffvecr, ffvec1, ffvec2
      integer, intent(in) ::
     &     idxvecr,ist0_r,
     &     idxvec1,ist0_1,
     &     idxvec2,ist0_2,lenvec, lenbuf
      real(8), intent(in) ::
     &     xfac1, xfac2, xshift, xpot
      real(8), intent(inout) ::
     &     xbuf1(lenbuf), xbuf2(lenbuf)

      integer ::
     &     luvecr, luvec1, luvec2, lblk,
     &     nrecs, len_rest, nrecbuf, nbatch, nrec_lbat, len_lbat,
     &     irecst1, irecst2, irecstr, lenbat, ibatch

      luvecr = ffvecr%unit
      luvec1 = ffvec1%unit
      luvec2 = ffvec2%unit
      lblk = ffvecr%reclen
      if (lblk.ne.ffvec1%reclen.or.
     &    lblk.ne.ffvec2%reclen) then
        write(luout,*) 'Incompatible block-lengthes: ',
     &       lblk,ffvec1%reclen,ffvec2%reclen
        call quit(1,'da_diavec','Incompatible block-lengthes')
      end if

      if (lenbuf.lt.min(lenvec,lblk)) then
        write(luout,*) 'Insufficient buffer size!'
        write(luout,*) ' buffer length = ',lenbuf,' blocksize = ',lblk
        call quit(1,'da_diavec','Insufficient buffer size!')
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
      
      irecst1 = (idxvec1-1)*nrecs+ist0_1
      irecst2 = (idxvec2-1)*nrecs+ist0_2
      irecstr = (idxvecr-1)*nrecs+ist0_r
      do ibatch = 1, nbatch
        if (ibatch.eq.nbatch) lenbat = len_lbat

        if (xfac2.eq.0d0) then
          xbuf1(1:lenbat) = 0d0
        else
          call get_vec_a(ffvec1,xbuf1,lenbat,irecst1)
          if (xpot.ne.0d0) then
            call get_vec_a(ffvec2,xbuf2,lenbat,irecst2)

            if (xshift.ne.0d0)
     &           xbuf2(1:lenbat) = xbuf2(1:lenbat) + xshift

            if (xpot.eq.-2d0) then
              xbuf1(1:lenbat) = xfac2*xbuf1(1:lenbat)/
     &             (xbuf2(1:lenbat)*xbuf2(1:lenbat))
            else if (xpot.eq.-1d0) then
              xbuf1(1:lenbat) = xfac2*xbuf1(1:lenbat)/xbuf2(1:lenbat)
            else if (xpot.eq.1d0) then 
              xbuf1(1:lenbat) = xfac2*xbuf1(1:lenbat)*xbuf2(1:lenbat)
            else if (xpot.eq.2d0) then 
              xbuf1(1:lenbat) = xfac2*xbuf1(1:lenbat)
     &             *xbuf2(1:lenbat)*xbuf2(1:lenbat)
            else
              xbuf1(1:lenbat) = xfac2*xbuf1(1:lenbat)
     &             *xbuf2(1:lenbat)**xpot
            end if

          else if (xfac2.ne.1d0) then
            xbuf1(1:lenbat) = xfac2*xbuf1(1:lenbat)
          end if
        end if

        if (xfac1.ne.0d0) then
          call get_vec_a(ffvecr,xbuf2,lenbat,irecstr)
          xbuf1(1:lenbat) = xbuf1(1:lenbat)+xfac1*xbuf2(1:lenbat)
        end if

        call put_vec_a(ffvecr,xbuf1,lenbat,irecstr)

        irecst1 = irecst1 + nrecbuf
        irecst2 = irecst2 + nrecbuf
        irecstr = irecstr + nrecbuf

      end do

      if (ntest.ge.100) then
        write(luout,*) 'diavec needed ',nbatch,' batches'
      end if
      
      return
      end
