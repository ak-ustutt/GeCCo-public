*----------------------------------------------------------------------*
      real(8) function da_ddot(ffvec1,idxvec1,ist0_1,
     &                         ffvec2,idxvec2,ist0_2,
     &                         lenvec,xbuf1,xbuf2,lenbuf)
*----------------------------------------------------------------------*
*
*     result = vec1^T*vec2
*
*     files/start-records may be identical
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00

      type(filinf), intent(in) ::
     &     ffvec1, ffvec2      
      integer, intent(in) ::
     &     idxvec1,ist0_1,
     &     idxvec2,ist0_2,lenvec,lenbuf
      real(8), intent(inout) ::
     &     xbuf1(lenbuf), xbuf2(lenbuf)

      integer ::
     &     lblk, luvec1, luvec2,
     &     nrecs, len_rest, nrecbuf, nbatch, nrec_lbat, len_lbat,
     &     irecst1, irecst2, lenbat, ibatch, ii
      real(8) ::
     &     res

      luvec1 = ffvec1%unit
      luvec2 = ffvec2%unit
      lblk = ffvec1%reclen
      if (lblk.ne.ffvec2%reclen) then
        write(luout,*) 'Incompatible block-lengthes: ',
     &       lblk,ffvec2%reclen
        call quit(1,'da_ddot','Incompatible block-lengthes')
      end if

      if (lenbuf.lt.min(lenvec,lblk)) then
        write(luout,*) 'Insufficient buffer size!'
        write(luout,*) ' buffer length = ',lenbuf,' blocksize = ',lblk
        call quit(1,'da_ddot','Insufficient buffer size!')
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
      
      irecst1 = (idxvec1-1)*nrecs+ist0_1
      irecst2 = (idxvec2-1)*nrecs+ist0_2
      res = 0d0
      do ibatch = 1, nbatch
        if (ibatch.eq.nbatch) lenbat = len_lbat
        call get_vec_a(ffvec1,xbuf1,lenbat,irecst1)
        if (luvec1.ne.luvec2.or.irecst1.ne.irecst2) then
          call get_vec_a(ffvec2,xbuf2,lenbat,irecst2)
          
          do ii = 1, lenbat
            res = res + xbuf1(ii)*xbuf2(ii)
          end do

        else

          do ii = 1, lenbat
            res = res + xbuf1(ii)*xbuf1(ii)
          end do

        end if

        if (ntest.ge.100) then
          write(luout,'(x,a,4i5,e20.10)') 'batch/st1/st2/len/res: ',
     &         ibatch,irecst1,irecst2,lenbat,res
        end if

        irecst1 = irecst1 + nrecbuf
        irecst2 = irecst2 + nrecbuf

      end do

      da_ddot = res

      if (ntest.ge.100) then
        write(luout,*) 'da_ddot needed ',nbatch,' batches'
      end if
      
      return
      end
