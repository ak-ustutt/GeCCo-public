*----------------------------------------------------------------------*
      subroutine optc_minspace(
     &     iord_vsbsp,ffvsbsp,iord_rsbsp,ffrsbsp,
     &     vred,gred,mred,nred,nroot,
     &     ffscr,ioffscr,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
*----------------------------------------------------------------------*
*     reduce current subspace to minimal space, spanning exactly the
*     currently best solutions
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(inout) ::
     &     nred, iord_vsbsp(*), iord_rsbsp(*)
      type(filinf), intent(inout) ::
     &     ffvsbsp, ffrsbsp, ffscr
      integer, intent(in) ::
     &     nroot, nincore, nwfpar, lenbuf, ioffscr
      real(8), intent(inout) ::
     &     vred(nred*nroot), gred(nred*nroot), mred(nred*nred)
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

      integer ::
     &     iroot, idxvec
      real(8) ::
     &     xnrm
      real(8) ::
     &     smat(nroot*nroot), mscr(nred*nroot), vorth(nred*nroot),
     &     vscr(nroot)

      real(8), external ::
     &     dnrm2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'optc_minspace')
        write(luout,*) 'nred, nroot: ',nred,nroot
        write(luout,*) 'vred on entry:'
        call wrtmat2(vred,nred,nroot,nred,nroot)
        write(luout,*) 'gred on entry:'
        call wrtmat2(gred,nred,nroot,nred,nroot)
        write(luout,*) 'mred on entry:'
        call wrtmat2(mred,nred,nred,nred,nred)
      end if

      if (nroot.gt.nred)
     &     call quit(1,'optc_minspace','strange: nroot.gt.nred!')

      ! get overlap matrix of vectors S_ij = <v_i|v_j>
      call dgemm('T','N',nroot,nroot,nred,
     &           1d0,vred,nred,
     &               vred,nred,
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
     &           1d0,vred,nred,
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
        print *,'iroot, norm: ',iroot,xnrm
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
     &               vred,nred,
     &           0d0,smat,nroot) ! use smat as scratch
      vred(1:nroot*nroot) = smat(1:nroot*nroot)

      ! ... update RHS vectors ...
      call dgemm('T','N',nroot,nroot,nred,
     &           1d0,vorth,nred,
     &               gred,nred,
     &           0d0,smat,nroot) ! use smat as scratch
      gred(1:nroot*nroot) = smat(1:nroot*nroot)

      if (ntest.ge.100) then
        write(luout,*) 'new vred: '
        call wrtmat2(vred,nroot,nroot,nroot,nroot)
        write(luout,*) 'new gred: '
        call wrtmat2(gred,nroot,nroot,nroot,nroot)
      end if

      ! ... and subspace matrix in reduced space ...
      call dgemm('T','N',nroot,nred,nred,
     &           1d0,vorth,nred,
     &               mred,nred,
     &           0d0,mscr,nroot)
      call dgemm('N','N',nroot,nroot,nred,
     &           1d0,mscr,nroot,
     &               vorth,nred,
     &           0d0,mred,nroot)

      if (ntest.ge.100) then
        write(luout,*) 'new mred: '
        call wrtmat2(mred,nroot,nroot,nroot,nroot)
      end if

      ! update trial vectors and MV-products in full space
      ! assemble trial vectors:
      idxvec = 1
      do iroot = 1, nroot
        call optc_expand_vec(vorth(idxvec),nred,xnrm,.false.,
     &       ffscr,ioffscr+iroot,0d0,ffvsbsp,iord_vsbsp,
     &       nincore,nwfpar,lenbuf,xbuf1,xbuf2)
        idxvec = idxvec+nred
      end do
      ! copy to ffvsbsp and reset iord_vsbsp:
      do iroot = 1, nroot
        iord_vsbsp(iroot) = iroot
        call da_sccpvec(ffscr,ioffscr+iroot,
     &                  ffvsbsp,iroot,
     &                  1d0,nwfpar,xbuf1,lenbuf)
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
     &       ffscr,ioffscr+iroot,0d0,ffrsbsp,iord_rsbsp,
     &       nincore,nwfpar,lenbuf,xbuf1,xbuf2)
        idxvec = idxvec+nred
      end do
      ! copy to ffrsbsp and reset iord_rsbsp:
      do iroot = 1, nroot
        iord_rsbsp(iroot) = iroot
        call da_sccpvec(ffscr,ioffscr+iroot,
     &                  ffrsbsp,iroot,
     &                  1d0,nwfpar,xbuf1,lenbuf)
      end do
c dbg
      print *,'contents on new ffrsbsp file:'
      call da_listvec(ffrsbsp,1,nwfpar,0,xbuf1,lenbuf)
c dbg
      
      ! reset nred:
      nred = nroot

      return
      end
