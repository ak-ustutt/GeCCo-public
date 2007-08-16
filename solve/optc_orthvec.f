*----------------------------------------------------------------------*
      subroutine optc_orthvec(nadd,
     &     ff_sbsp,iord_sbsp,ndim_sbsp,mxsbsp,zero_vec,
     &     ffnew,nnew,
     &     nwfpar,nincore,xbuf1,xbuf2,xbuf3,lenbuf)
*----------------------------------------------------------------------*
*     orthogonalize and add nnew new vectors (on ff_new) to subspace 
*     (on ff_sbsp); actual number of new vectors returned on nadd
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(out) ::
     &     nadd
      integer, intent(in) ::
     &     nnew, mxsbsp, nwfpar, lenbuf, nincore
      integer, intent(inout) ::
     &     ndim_sbsp, iord_sbsp(mxsbsp)
      type(filinf), intent(in) ::
     &     ff_sbsp, ffnew
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)
      logical ::
     &     zero_vec(ndim_sbsp)

      integer ::
     &     nold, inew, jnew, iold, ii, jj, idx, ivec, irec
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
        write(luout,*) 'ndim, nnew: ',ndim_sbsp, nnew
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

      ! a) <new|old>
      if (nold.gt.0.and.nold.ge.nnew) then
        do inew = 1, nnew
          
          ii = inew+nold
          if (nincore.ge.2) then
            call vec_from_da(ffnew,inew,xbuf1,nwfpar)
          end if

          do iold = 1, nold
            
            jj = iord_sbsp(iold)
            if (nincore.ge.2) then
              call vec_from_da(ff_sbsp,iold,xbuf2,nwfpar)
              smat(jj,ii) = ddot(nwfpar,xbuf1,1,xbuf2,1)
              smat(ii,jj) = smat(jj,ii)
            else
              smat(jj,ii) = da_ddot(ff_sbsp,iold,ffnew,inew,
     &             nwfpar,xbuf1,xbuf2,lenbuf)
              smat(ii,jj) = smat(jj,ii)
            end if
            
          end do

        end do

      else if (nold.gt.0.and.nnew.gt.nold) then

        do iold = 1, nold
          
          ii = iord_sbsp(iold)
          if (nincore.ge.2) then
            call vec_from_da(ff_sbsp,iold,xbuf1,nwfpar)
          end if

          do inew = 1, nnew
            
            jj = inew+nold
            if (nincore.ge.2) then
              call vec_from_da(ffnew,inew,xbuf2,nwfpar)
              smat(jj,ii) = ddot(nwfpar,xbuf1,1,xbuf2,1)
              smat(ii,jj) = smat(jj,ii)
            else
              smat(jj,ii) = da_ddot(ff_sbsp,iold,ffnew,inew,
     &             nwfpar,xbuf1,xbuf2,lenbuf)
            end if
            
          end do

        end do

      end if

      ! b) <new|new>
      do inew = 1, nnew
        if (nincore.ge.2) then
          call vec_from_da(ffnew,inew,xbuf1,nwfpar)
          smat(nold+inew,nold+inew) = ddot(nwfpar,xbuf1,1,xbuf1,1)
        else
          smat(nold+inew,nold+inew) = da_ddot(ffnew,inew,ffnew,inew,
     &           nwfpar,xbuf1,xbuf2,lenbuf)
        end if

        do jnew = inew+1, nnew
          if (nincore.ge.2) then
            call vec_from_da(ffnew,jnew,xbuf2,nwfpar)
            smat(nold+inew,nold+jnew) = ddot(nwfpar,xbuf1,1,xbuf2,1)
            smat(nold+jnew,nold+inew) = smat(inew,jnew)
          else
            smat(nold+inew,nold+jnew) = da_ddot(ffnew,inew,ffnew,jnew,
     &           nwfpar,xbuf1,xbuf2,lenbuf)
            smat(nold+jnew,nold+inew) = smat(nold+inew,nold+jnew)
          end if
        end do

      end do

      if (ntest.ge.100) then
        write(luout,*) 'the overlap matrix:'
        call wrtmat2(smat,nold+nnew,nold+nnew,nold+nnew,nold+nnew)
      end if

      ! call modified Gram-Schmidt:
      call mgs(xmat,smat,nold+nnew,scr)

      ! assemble new vectors

      ! dummy iord for ffnew:
      do inew = 1, nnew        
        iordnew(inew) = inew
      end do

      idx = nold
      do inew = 1, nnew
        ! check whether vector is nonzero
        if (ddot(nold+nnew,xmat(1,nold+inew),1,
     &                     xmat(1,nold+inew),1).lt.1d-10)
     &       cycle
        idx = idx+1
        if (idx.lt.ndim_sbsp) then
          ! overwrite previous (zero) record
          irec = ioptc_get_sbsp_rec(idx,iord_sbsp,ndim_sbsp,mxsbsp)
        else
          ! get new record
          irec = ioptc_get_sbsp_rec(0,iord_sbsp,ndim_sbsp,mxsbsp)
        end if
        
        if (ntest.ge.100)
     &       write(luout,*) 'next new record: ',irec

        ! add contributions from new space
        call optc_expand_vec(xmat(nold+1,nold+inew),nnew,xdum,.false.,
     &                       ff_sbsp,irec,0d0,ffnew,iordnew,
     &                       nincore,nwfpar,lenbuf,xbuf1,xbuf2)
        ! add contributions from old space
        call optc_expand_vec(xmat(1,nold+inew),nold,xdum,.false.,
     &                       ff_sbsp,irec,1d0,ff_sbsp,iord_sbsp,
     &                       nincore,nwfpar,lenbuf,xbuf1,xbuf2)

      end do

      nadd = idx-nold

      deallocate(smat,xmat,scr)

      return
      end
