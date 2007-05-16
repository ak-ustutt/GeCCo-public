*----------------------------------------------------------------------*
      subroutine da_matvec(ffvecr,idxvecr,xfac,
     &                     ffvecm,idxvecm,xvec,nvec,thrsh,
     &                     lenvec,xbuf1,xbuf2,lenbuf)
*----------------------------------------------------------------------*
*
*     vecr = xfac*vecr + sum_i xvec_i * vecm_i
*
*      for all (|xvec_i| > thrsh)
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 0
      
      type(filinf), intent(in) ::
     &     ffvecr, ffvecm
      integer, intent(in) ::
     &     idxvecr,
     &     idxvecm, nvec,
     &     lenvec, lenbuf
      real(8), intent(in) ::
     &     xfac, xvec(nvec), thrsh
      real(8), intent(inout) ::
     &     xbuf1(lenbuf), xbuf2(lenbuf)

      logical ::
     &     first
      integer ::
     &     ivec
      real(8) ::
     &     xxfac

      first = .true.
      do ivec = 1, nvec
        if (abs(xvec(ivec)).le.thrsh) cycle
        
        if (first.and.abs(xfac).lt.thrsh) then
          call da_sccpvec(ffvecm,idxvecm-1+ivec,ffvecr,idxvecr,
     &         xvec(ivec),lenvec,xbuf1,lenbuf)
        else
          xxfac = 1d0
          if (first) xxfac = xfac

          call da_vecsum(ffvecr,idxvecr,
     &                 ffvecr,idxvecr,xxfac,
     &                 ffvecm,idxvecm-1+ivec,xvec(ivec),
     &                 lenvec,xbuf1,xbuf2,lenbuf)
        end if
        first = .false.

      end do

      return
      end
