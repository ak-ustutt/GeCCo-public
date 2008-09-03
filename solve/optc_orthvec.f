*----------------------------------------------------------------------*
      subroutine optc_orthvec(nadd,
     &     ff_sbsp,iord_sbsp,ndim_sbsp,mxsbsp,zero_vec,
     &     ffnew,nnew,nopt,
     &     nwfpar,nincore,xbuf1,xbuf2,xbuf3,lenbuf)
*----------------------------------------------------------------------*
*     orthogonalize and add nnew new vectors (on ff_new) to subspace 
*     (on ff_sbsp); actual number of new vectors returned on nadd
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_file_array.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     nadd
      integer, intent(in) ::
     &     nnew, mxsbsp, lenbuf, nincore, nopt, nwfpar(nopt)
      integer, intent(inout) ::
     &     ndim_sbsp, iord_sbsp(mxsbsp)
      type(file_array), intent(in) ::
     &     ff_sbsp(nopt)
      type(filinf), intent(in) ::
     &     ffnew(nopt)
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)
      logical ::
     &     zero_vec(ndim_sbsp)

      integer ::
     &     nold, inew, jnew, iold, ii, jj, idx, ivec, irec, iopt
      real(8) ::
     &     xdum
      logical ::
     &     previous_zero

      integer ::
     &     iordnew(nnew)
      real(8), allocatable ::
     &     smat(:,:), xmat(:,:), scr(:)

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     ddot, da_ddot

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'optc_orthvec')
        write(luout,*) 'ndim, nnew, nopt: ',ndim_sbsp, nnew, nopt
      end if

      ! look for zero-vectors (in initial iteration); they should
      ! reside at the end of the initial subspace
      nold = 0
      previous_zero = .false.
      do ivec = 1, ndim_sbsp
        if (zero_vec(ivec)) then
          previous_zero = .true.
        else
          if (previous_zero) then
            write(luout,*) 'zero_vec: ',zero_vec(1:ndim_sbsp)
            call quit(1,'optc_orthvec',
     &           'zero vectors should reside at the end of the'//
     &           ' initial subspace')
          end if
          nold = nold+1
        end if
      end do

      ! calculate overlap matrix
      allocate(
     &     smat(nold+nnew,nold+nnew),
     &     xmat(nold+nnew,nold+nnew),
     &     scr(nold+nnew))

      if (ntest.ge.100)
     &     write(luout,*) 'nold, nnew: ',nold, nnew

      smat(1:nold+nnew,1:nold+nnew) = 0d0
      do iold = 1, nold
        smat(iold,iold) = 1d0
      end do

      do iopt = 1, nopt

      ! a) <new|old>
      if (nold.gt.0.and.nold.ge.nnew) then
        do inew = 1, nnew
          
          ii = inew+nold
          if (nincore.ge.2) then
            call vec_from_da(ffnew(iopt),inew,xbuf1,nwfpar(iopt))
          end if

          do iold = 1, nold
            
            jj = iord_sbsp(iold)
            if (nincore.ge.2) then
              call vec_from_da(ff_sbsp(iopt)%fhand,iold,
     &                                             xbuf2,nwfpar(iopt))
              smat(jj,ii) = smat(jj,ii) +
     &             ddot(nwfpar(iopt),xbuf1,1,xbuf2,1)
              smat(ii,jj) = smat(jj,ii)
            else
              smat(jj,ii) = smat(jj,ii) +
     &             da_ddot(ff_sbsp(iopt)%fhand,iold,ffnew(iopt),inew,
     &                     nwfpar(iopt),xbuf1,xbuf2,lenbuf)
              smat(ii,jj) = smat(jj,ii)
            end if
            
          end do

        end do

      else if (nold.gt.0.and.nnew.gt.nold) then

        do iold = 1, nold
          
          ii = iord_sbsp(iold)
          if (nincore.ge.2) then
            call vec_from_da(ff_sbsp(iopt)%fhand,iold,
     &                                            xbuf1,nwfpar(iopt))
          end if

          do inew = 1, nnew
            
            jj = inew+nold
            if (nincore.ge.2) then
              call vec_from_da(ffnew(iopt),inew,xbuf2,nwfpar(iopt))
              smat(jj,ii) = smat(jj,ii) +
     &             ddot(nwfpar(iopt),xbuf1,1,xbuf2,1)
              smat(ii,jj) = smat(jj,ii)
            else
              smat(jj,ii) = smat(jj,ii) +
     &             da_ddot(ff_sbsp(iopt)%fhand,iold,ffnew(iopt),inew,
     &             nwfpar(iopt),xbuf1,xbuf2,lenbuf)
            end if
            
          end do

        end do

      end if

      ! b) <new|new>
      do inew = 1, nnew
        if (nincore.ge.2) then
          call vec_from_da(ffnew(iopt),inew,xbuf1,nwfpar(iopt))
          smat(nold+inew,nold+inew) = smat(nold+inew,nold+inew) +
     &         ddot(nwfpar(iopt),xbuf1,1,xbuf1,1)
c dbg
c          print *,'|new|: ',sqrt(smat(nold+inew,nold+inew))
c dbg
        else
          smat(nold+inew,nold+inew) = smat(nold+inew,nold+inew) +
     &         da_ddot(ffnew(iopt),inew,ffnew(iopt),inew,
     &         nwfpar(iopt),xbuf1,xbuf2,lenbuf)
        end if

        do jnew = inew+1, nnew
          if (nincore.ge.2) then
            call vec_from_da(ffnew(iopt),jnew,xbuf2,nwfpar(iopt))
            smat(nold+inew,nold+jnew) = smat(nold+inew,nold+jnew) +
     &           ddot(nwfpar(iopt),xbuf1,1,xbuf2,1)
            smat(nold+jnew,nold+inew) = smat(nold+inew,nold+jnew)
          else
            smat(nold+inew,nold+jnew) = smat(nold+inew,nold+jnew) +
     &           da_ddot(ffnew,inew,ffnew,jnew,
     &           nwfpar(iopt),xbuf1,xbuf2,lenbuf)
            smat(nold+jnew,nold+inew) = smat(nold+inew,nold+jnew)
          end if
        end do

      end do

      if (ntest.ge.100) then
        write(luout,*) 'the overlap matrix after iopt = ',iopt
        call wrtmat2(smat,nold+nnew,nold+nnew,nold+nnew,nold+nnew)
      end if

      end do ! iopt

c dbg
c     allocate larger scr!!
c      print *,'eigenvalues of smat:'
c      xmat = smat
c      call rs(nold+nnew,nold+nnew,xmat,scr,0,
c     &        scr(2*(nold+nnew)+1),scr(3*(nold+nnew)+1),
c     &     scr(4*(nold+nnew)+1),idx)
c      print *,scr(1:nold+nnew)
c dbg

      ! call modified Gram-Schmidt:
      call mgs(xmat,smat,nold+nnew,scr)

      if (ntest.ge.100) then
        write(luout,*) 'the reorth. matrix: '
        call wrtmat2(xmat,nold+nnew,nold+nnew,nold+nnew,nold+nnew)
      end if

      ! assemble new vectors

      ! dummy iord for ffnew:
      do inew = 1, nnew        
        iordnew(inew) = inew
      end do

      
      do iopt = 1, nopt

      idx = nold
c dbg
c      print *,'iopt = ',iopt,'  idx = ',idx
c dbg
      do inew = 1, nnew
        ! check whether vector is nonzero
c dbg
c        print *,'test = ',ddot(nold+nnew,xmat(1,nold+inew),1,
c     &                     xmat(1,nold+inew),1),nold,nnew,inew
c dbg
        if (ddot(nold+nnew,xmat(1,nold+inew),1,
     &                     xmat(1,nold+inew),1).lt.1d-10)
     &       cycle
        idx = idx+1
c dbg
c        print *,'idx, ndim_sbsp: ',idx, ndim_sbsp
c dbg
        if (idx.le.ndim_sbsp) then
          ! overwrite previous (zero) record
          irec = ioptc_get_sbsp_rec(idx,iord_sbsp,ndim_sbsp,mxsbsp)
        else
          ! get new record
          irec = ioptc_get_sbsp_rec(0,iord_sbsp,ndim_sbsp,mxsbsp)
        end if
c dbg
c        print *,'idx = ',idx,'   irec = ',irec
c dbg
        
        if (ntest.ge.100)
     &       write(luout,*) 'next new record: ',irec

        ! add contributions from new space
        call optc_expand_vec(xmat(nold+1,nold+inew),nnew,xdum,.false.,
     &                       ff_sbsp(iopt)%fhand,irec,0d0,
     &                                            ffnew(iopt),iordnew,
     &                       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)
        ! add contributions from old space
        call optc_expand_vec(xmat(1,nold+inew),nold,xdum,.false.,
     &                      ff_sbsp(iopt)%fhand,irec,1d0,
     &                                  ff_sbsp(iopt)%fhand,iord_sbsp,
     &                       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)

      end do

      end do ! iopt

      nadd = idx-nold

      deallocate(smat,xmat,scr)

      return
      end
