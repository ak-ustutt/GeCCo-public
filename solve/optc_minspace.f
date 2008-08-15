*----------------------------------------------------------------------*
      subroutine optc_minspace(
     &     iord_vsbsp,ffvsbsp,
     &     iord_rsbsp,ffrsbsp,
     &     iord_ssbsp,ffssbsp,use_s,
     &     vred,gred,mred,sred,nred,nroot,nrhs,mxdim,nopt,
     &     ffscr,ioffscr,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
*----------------------------------------------------------------------*
*     reduce current subspace to minimal space, spanning exactly the
*     currently best solutions
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_file_array.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(inout) ::
     &     nred, nrhs, mxdim, nopt,
     &     iord_vsbsp(*), iord_rsbsp(*), iord_ssbsp(*)
      logical, intent(in) ::
     &     use_s(nopt)
      type(file_array), intent(inout) ::
     &     ffvsbsp(nopt), ffrsbsp(nopt), ffssbsp(nopt)
      type(filinf), intent(inout) ::
     &     ffscr
      integer, intent(in) ::
     &     nroot, nincore, nwfpar(nopt), lenbuf, ioffscr
      real(8), intent(inout) ::
     &     vred(mxdim,nroot), gred(mxdim,nrhs),
     &     mred(mxdim,nred), sred(mxdim,nred)
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

      logical ::
     &     update_s
      integer ::
     &     iroot, idxvec, iopt
      real(8) ::
     &     xnrm
      real(8) ::
     &     smat(nroot,nroot), mscr(nred*nroot), vorth(nred*nroot),
     &     vscr(nroot)

      real(8), external ::
     &     dnrm2

      update_s = .false.
      do iopt = 1, nopt
        update_s = update_s.or.use_s(iopt)
      end do

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'optc_minspace')
        write(luout,*) 'nred, nroot: ',nred,nroot
        write(luout,*) 'vred on entry:'
        call wrtmat2(vred,nred,nroot,mxdim,nroot)
        if (nrhs.gt.0) then
          write(luout,*) 'gred on entry:'
          call wrtmat2(gred,nred,nroot,mxdim,nroot)
        end if
        write(luout,*) 'mred on entry:'
        call wrtmat2(mred,nred,nred,mxdim,nred)
        if (update_s) then
          write(luout,*) 'sred on entry:'
          call wrtmat2(sred,nred,nred,mxdim,nred)
        end if
      end if

      if (nroot.gt.nred)
     &     call quit(1,'optc_minspace','strange: nroot.gt.nred!')

      ! get overlap matrix of vectors S_ij = <v_i|v_j>
      call dgemm('T','N',nroot,nroot,nred,
     &           1d0,vred,mxdim,
     &               vred,mxdim,
     &           0d0,smat,nroot)

      if (ntest.ge.100) then
        write(luout,*) 'overlap matrix: '
        call wrtmat2(smat,nroot,nroot,nroot,nroot)
      end if

      ! obtain orthonormalization matrix
      call mgs(mscr,smat,nroot,vscr)

      ! obtain trafo: reduced space sol. ->  orth reduced space sol.
      ! |v_i^{orth}> = |v_j> t(j,i)
      call dgemm('N','N',nred,nroot,nroot,
     &           1d0,vred,mxdim,
     &               mscr,nroot,
     &           0d0,vorth,nred)

      if (ntest.ge.100) then
        write(luout,*) 'orthogonalized vred: '
        call wrtmat2(vorth,nred,nroot,nred,nroot)
      end if

      ! check norm of new vectors, lin. dependency occurred?
      do iroot = 1, nroot
        xnrm = dnrm2(nred,vorth,1)
c dbg
c        print *,'iroot, norm: ',iroot,xnrm
c dbg
        if (xnrm.lt.1d-4) then
          write(luout,*) 'root, norm: ',iroot, xnrm
          write(luout,*) 'looks like linear dependency ... oh no!'
          call quit(1,'optc_minspace',
     &         'linear dependency should not occur at this place!')
        end if
      end do

      ! trafo: old reduced space -> new reduced space
      ! |v_i> = |b_old><b_old|v_i> = |b_new><b_new|b_old><v_old|v_i>
      ! update solution vectors ...
      call dgemm('T','N',nroot,nroot,nred,
     &           1d0,vorth,nred,
     &               vred,mxdim,
     &           0d0,smat,nroot) ! use smat as scratch
      do iroot = 1, nroot
        vred(1:nroot,iroot) = smat(1:nroot,iroot)
      end do

      if (nrhs.gt.0) then
        ! ... update RHS vectors ...
        call dgemm('T','N',nroot,nrhs,nred,
     &             1d0,vorth,nred,
     &                 gred,mxdim,
     &             0d0,smat,nroot) ! use smat as scratch
        do iroot = 1, nrhs
          gred(1:nroot,iroot) = smat(1:nroot,iroot)
        end do
      end if

      if (ntest.ge.100) then
        write(luout,*) 'new vred: '
        call wrtmat2(vred,nroot,nroot,mxdim,nroot)
        if (nrhs.gt.0) then
          write(luout,*) 'new gred: '
          call wrtmat2(gred,nroot,nroot,mxdim,nroot)
        end if
      end if

      ! ... and subspace matrix in reduced space ...
      call dgemm('T','N',nroot,nred,nred,
     &           1d0,vorth,nred,
     &               mred,mxdim,
     &           0d0,mscr,nroot)
      call dgemm('N','N',nroot,nroot,nred,
     &           1d0,mscr,nroot,
     &               vorth,nred,
     &           0d0,mred,mxdim)

      if (ntest.ge.100) then
        write(luout,*) 'new mred: '
        call wrtmat2(mred,nroot,nroot,mxdim,nroot)
      end if

      ! do the same for metric:
      if (update_s) then
        call dgemm('T','N',nroot,nred,nred,
     &           1d0,vorth,nred,
     &               sred,mxdim,
     &           0d0,mscr,nroot)
        call dgemm('N','N',nroot,nroot,nred,
     &           1d0,mscr,nroot,
     &               vorth,nred,
     &           0d0,sred,mxdim)

        if (ntest.ge.100) then
          write(luout,*) 'new sred: '
          call wrtmat2(sred,nroot,nroot,mxdim,nroot)
        end if

      end if

      do iopt = 1, nopt

      ! update trial vectors and MV-products in full space
      ! assemble trial vectors:
      idxvec = 1
      do iroot = 1, nroot
        call optc_expand_vec(vorth(idxvec),nred,xnrm,.false.,
     &       ffscr,ioffscr+iroot,0d0,ffvsbsp(iopt)%fhand,iord_vsbsp,
     &       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)
        idxvec = idxvec+nred
      end do
      ! copy to ffvsbsp and reset iord_vsbsp:
      do iroot = 1, nroot
        iord_vsbsp(iroot) = iroot
        call da_sccpvec(ffscr,ioffscr+iroot,
     &                  ffvsbsp(iopt)%fhand,iroot,
     &                  1d0,nwfpar(iopt),xbuf1,lenbuf)
      end do
c dbg
c      print *,'ioffscr = ',ioffscr
c      print *,'contents on new ffscr file:'
c      call da_listvec(ffscr,1,nwfpar,0,xbuf1,lenbuf)
c      call da_listvec(ffscr,2,nwfpar,0,xbuf1,lenbuf)
c      print *,'contents on new ffvsbsp file:'
c      call da_listvec(ffvsbsp,1,nwfpar,0,xbuf1,lenbuf)
c dbg

      ! assemble MV products:
      idxvec = 1
      do iroot = 1, nroot
        call optc_expand_vec(vorth(idxvec),nred,xnrm,.false.,
     &       ffscr,ioffscr+iroot,0d0,ffrsbsp(iopt)%fhand,iord_rsbsp,
     &       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)
        idxvec = idxvec+nred
      end do
      ! copy to ffrsbsp and reset iord_rsbsp:
      do iroot = 1, nroot
        iord_rsbsp(iroot) = iroot
        call da_sccpvec(ffscr,ioffscr+iroot,
     &                  ffrsbsp(iopt)%fhand,iroot,
     &                  1d0,nwfpar(iopt),xbuf1,lenbuf)
      end do
c dbg
c      print *,'contents on new ffrsbsp file:'
c      call da_listvec(ffrsbsp(iopt)%fhand,1,nwfpar,0,xbuf1,lenbuf)
c dbg

      if (use_s(iopt)) then
        ! assemble metric-vector products:
        idxvec = 1
        do iroot = 1, nroot
          call optc_expand_vec(vorth(idxvec),nred,xnrm,.false.,
     &       ffscr,ioffscr+iroot,0d0,ffssbsp(iopt)%fhand,iord_ssbsp,
     &       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)
          idxvec = idxvec+nred
        end do
        ! copy to ffssbsp and reset iord_ssbsp:
        do iroot = 1, nroot
          iord_ssbsp(iroot) = iroot
          call da_sccpvec(ffscr,ioffscr+iroot,
     &                  ffssbsp(iopt)%fhand,iroot,
     &                  1d0,nwfpar(iopt),xbuf1,lenbuf)
        end do
c dbg
c        print *,'contents on new ffssbsp file:'
c        call da_listvec(ffssbsp(iopt)%fhand,1,nwfpar,0,xbuf1,lenbuf)
c dbg
      end if
      
      end do ! iopt

      ! reset nred:
      nred = nroot

      return
      end
