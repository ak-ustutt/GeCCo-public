*----------------------------------------------------------------------*
      subroutine optc_diis_sbsp_add(ndim_rsbsp,ndim_vsbsp,mxdim_sbsp,
     &     get_new_rec,
     &     iord_vsbsp,iord_rsbsp,
     &     ffamp,ffgrd,ffdia,ff_rsbsp,ff_vsbsp,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2)
*----------------------------------------------------------------------*
*
*     add (preconditioned) gradient and sum of current vector and
*      minus the prec. gradient to the respective subspace files
*
*     incore.ge.2: on output is the precond. gradient on xbuf1
*                      and vector - precond. gradient on xbuf2
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     get_new_rec
      type(filinf), intent(in) ::
     &     ffamp,ffgrd,ffdia,ff_rsbsp,ff_vsbsp
      integer, intent(inout) ::
     &     ndim_rsbsp,ndim_vsbsp,iord_vsbsp,iord_rsbsp
      integer, intent(in) ::
     &     mxdim_sbsp,
     &     nincore, nwfpar, lenbuf
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*)

      integer ::
     &     irecr, irecv, inum

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     dnrm2

      ! get record-number for new vector in subspace
      inum = ndim_rsbsp   ! get last entry
      if (get_new_rec) inum = 0 ! request new entry
      irecr = ioptc_get_sbsp_rec(inum,iord_rsbsp,ndim_rsbsp,mxdim_sbsp)
      inum = ndim_vsbsp   ! get last entry
      if (get_new_rec) inum = 0 ! request new entry
      irecv = ioptc_get_sbsp_rec(inum,iord_vsbsp,ndim_vsbsp,mxdim_sbsp)

      if (ntest.ge.100) then
        write(luout,*) 'added records: ',irecv, irecr
        write(luout,*) 'nwfpar: ',nwfpar
      end if

      if (nincore.ge.2) then

        call vec_from_da(ffgrd,1,xbuf1,nwfpar)
        call vec_from_da(ffdia,1,xbuf2,nwfpar)
c dbg
c        print *,'g,d:',dnrm2(nwfpar,xbuf1,1),dnrm2(nwfpar,xbuf2,1)
c dbg

        ! prelim. w/o damping
        xbuf1(1:nwfpar) = xbuf1(1:nwfpar)/xbuf2(1:nwfpar)
c dbg
c        print *,'g/d:',dnrm2(nwfpar,xbuf1,1)
c dbg

        call vec_from_da(ffamp,1,xbuf2,nwfpar)
c dbg
c        print *,'t:',dnrm2(nwfpar,xbuf2,1)
c dbg

        xbuf2(1:nwfpar) = xbuf2(1:nwfpar) - xbuf1(1:nwfpar)

        call vec_to_da(ff_rsbsp,irecr,xbuf1,nwfpar)
        call vec_to_da(ff_vsbsp,irecv,xbuf2,nwfpar)

      else

        ! prelim. w/o damping
        ! add D^-1|gradient(n)> to subspace
        call da_diavec(ff_rsbsp,irecr,0d0,
     &       ffgrd,1,1d0,
     &       ffdia,1,0d0,-1d0,
     &       nwfpar,xbuf1,xbuf2,lenbuf)

        ! add |vec(n)> - D^-1|gradient(n)> to subspace
        call da_vecsum(ff_vsbsp,irecv,
     &       ff_rsbsp,irecr,-1d0,
     &       ffamp,1,1d0,
     &       nwfpar,xbuf1,xbuf2,lenbuf)

      end if

      return
      end
