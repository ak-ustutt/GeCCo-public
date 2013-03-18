*----------------------------------------------------------------------*
      subroutine optc_sbsp_add(ndim_sbsp,mxdim_sbsp,iord_sbsp,
     &     ffvec,irec,ffdia,divide,ff_sbsp,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2)
*----------------------------------------------------------------------*
*
* add vector (or precond. vector) to subspace file
* no need to care for incore/out-of-core
*     
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00

      type(filinf) ::
     &     ffvec,ffdia,ff_sbsp
      logical, intent(in) ::
     &     divide
      integer, intent(inout) ::
     &     ndim_sbsp,iord_sbsp(mxdim_sbsp)
      integer, intent(in) ::
     &     mxdim_sbsp,irec,
     &     nincore, nwfpar, lenbuf
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*)

      integer ::
     &     isbsprec

      integer, external ::
     &     ioptc_get_sbsp_rec

      isbsprec = ioptc_get_sbsp_rec(0,iord_sbsp,ndim_sbsp,mxdim_sbsp)

      if (divide) then

        call da_diavec(ff_sbsp,isbsprec,1,0d0,
     &       ffvec,irec,1,1d0,
     &       ffdia,1,1,0d0,-1d0,
     &       nwfpar,xbuf1,xbuf2,lenbuf)

      else

        call da_sccpvec(ffvec,irec,
     &               ff_sbsp,isbsprec,1d0,
     &               nwfpar,xbuf1,lenbuf)  

      end if

      return
      end
