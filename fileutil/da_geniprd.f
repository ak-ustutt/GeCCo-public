*----------------------------------------------------------------------*
      real(8) function da_geniprd(ffvec1,idxvec1,ffvec2,idxvec2,
     &                           ffvec3,idxvec3,xshift,xpot,
     &                           lenvec,xbuf1,xbuf2,lenbuf)
*----------------------------------------------------------------------*
*
*     result = vec1_i (vec3_i+shift)**xpot vec2_i
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
     &     ffvec1, ffvec2, ffvec3
      integer, intent(in) ::
     &     idxvec1,
     &     idxvec2,
     &     idxvec3,lenvec, lenbuf
      real(8), intent(in) ::
     &     xshift, xpot
      real(8), intent(inout) ::
     &     xbuf1(lenbuf), xbuf2(lenbuf)

      integer ::
     &     luvec1, luvec2, luvec3, lblk,
     &     nrecs, len_rest, nrecbuf, nbatch, nrec_lbat, len_lbat,
     &     irecst1, irecst2, irecst3, lenbat, ibatch, ii
      real(8) ::
     &     res

      luvec1 = ffvec1%unit
      luvec2 = ffvec2%unit
      luvec3 = ffvec3%unit
      lblk = ffvec1%reclen
      if (lblk.ne.ffvec2%reclen.or.
     &    lblk.ne.ffvec3%reclen) then
        write(luout,*) 'Incompatible block-lengthes: ',
     &       lblk,ffvec2%reclen,ffvec3%reclen
        call quit(1,'da_geniprd','Incompatible block-lengthes')
      end if

      if (lenbuf.lt.min(lenvec,lblk)) then
        write(luout,*) 'Insufficient buffer size!'
        write(luout,*) ' buffer length = ',lenbuf,' blocksize = ',lblk
        stop 'Insufficient buffer size!'
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
      
      irecst1 = (idxvec1-1)*nrecs+1
      irecst2 = (idxvec2-1)*nrecs+1
      irecst3 = (idxvec3-1)*nrecs+1
      res = 0d0
      do ibatch = 1, nbatch
        if (ibatch.eq.nbatch) lenbat = len_lbat

        call get_vec_a(ffvec1,xbuf1,lenbat,irecst1)

        call get_vec_a(ffvec3,xbuf2,lenbat,irecst3)

        if (xshift.ne.0d0)
     &       xbuf2(1:lenbat) = xbuf2(1:lenbat) + xshift

        if (xpot.eq.-2d0) then
          xbuf2(1:lenbat) = xbuf1(1:lenbat)/
     &         (xbuf2(1:lenbat)*xbuf2(1:lenbat))
        else if (xpot.eq.-1d0) then
          xbuf2(1:lenbat) = xbuf1(1:lenbat)/xbuf2(1:lenbat)
        else if (xpot.eq.0d0) then
          xbuf2(1:lenbat) = xbuf1(1:lenbat)
        else if (xpot.eq.1d0) then
          xbuf2(1:lenbat) = xbuf1(1:lenbat)*xbuf2(1:lenbat)
        else if (xpot.eq.2d0) then
          xbuf2(1:lenbat) = xbuf1(1:lenbat)
     &         *xbuf2(1:lenbat)*xbuf2(1:lenbat)
        else
          xbuf2(1:lenbat) = xbuf1(1:lenbat)*xbuf2(1:lenbat)**xpot
        end if

        if (luvec1.ne.luvec2.or.irecst1.ne.irecst2)
     &       call get_vec_a(ffvec2,xbuf1,lenbat,irecst2)
                
        do ii = 1, lenbat
          res = res + xbuf1(ii)*xbuf2(ii)
        end do

        irecst1 = irecst1 + nrecbuf
        irecst2 = irecst2 + nrecbuf
        irecst3 = irecst3 + nrecbuf

      end do

      da_geniprd = res

      if (ntest.ge.100) then
        write(luout,*) 'da_geniprd needed ',nbatch,' batches'
      end if
      
      return
      end
