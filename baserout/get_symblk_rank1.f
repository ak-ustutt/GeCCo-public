*----------------------------------------------------------------------*
      subroutine get_symblk_rank1(xblk,xfull,sym,psign,
     &                            ndimt,ndim,nsym)
*----------------------------------------------------------------------*
*     extract symmetry block from rank 1 matrix given in full form
*     xfull is in upper triangular form
*     psign = +1 for symmetric matrix
*     psign = -1 for anti-symmetric matrix
*     adapted from jeppe olsen
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     sym, nsym, ndim(nsym), ndimt
      real(8), intent(in) ::
     &     xfull(*), psign
      real(8), intent(out) ::
     &     xblk(*)

      integer ::
     &     irsym, icsym, isym, nrow, ncol,
     &     irow, icol, ibrow, ibcol, ioff,
     &     icabs, irabs, icrmin, icrmax, form
      real(8) ::
     &     fac

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'get_symblk_rank1')
        write(lulog,*) 'symmetry of block: ',sym
        write(lulog,*) 'input matrix:'
        call prtrlt(xfull,ndimt)
      end if

      ioff = 0
      do irsym = 1, nsym
        icsym = multd2h(sym,irsym)

        nrow = ndim(irsym)
        ncol = ndim(icsym)

        ibrow = 1
        do isym = 1, irsym-1
          ibrow = ibrow + ndim(isym)
        end do
        ibcol = 1
        do isym = 1, icsym-1
          ibcol = ibcol + ndim(isym)
        end do

        do irow = 1, nrow
          irabs = ibrow + irow-1          
          do icol = 1, ncol
            icabs = ibcol + icol-1
            icrmax = max(icabs,irabs)
            icrmin = min(icabs,irabs)

            fac = 1d0
            if (irabs.lt.icabs) fac = psign
            
c dbg
c            print *,'->',ioff + (icol-1)*nrow+irow,
c     &           icrmax*(icrmax-1)/2+icrmin
c dbg
            xblk(ioff + (icol-1)*nrow+irow) =
     &           fac*xfull(icrmax*(icrmax-1)/2+icrmin)

          end do
        end do
        ioff = ioff + nrow*ncol

      end do

      if (ntest.ge.100) then
        write(lulog,*) 'output matrix:'
        form = sym
        if (sym.eq.1) form=0
        call wr_blkmat(xblk,ndim,ndim,nsym,form)
      end if

      return
      end
