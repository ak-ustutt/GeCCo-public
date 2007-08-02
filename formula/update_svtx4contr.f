*----------------------------------------------------------------------*
      subroutine update_svtx4contr(contr)
*----------------------------------------------------------------------*
*     for a given contr%svertex array:
*     - check correct index sequence (i.e. first appearance in ascending
*       order
*     - set actual number of super vertices (contr%nsupvtx)
*     - set up the array contr%joined
*     all arrays are assumed to be allocated properly (by resize_contr)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      type(contraction), intent(inout) ::
     &     contr

      integer ::
     &     ivtx, jvtx, idx_ori, idx, npvtx

      integer ::
     &     svertex_new(contr%nvtx)

      ! ensure canonical order:
      idx = 0
      iloop: do ivtx = 1, contr%nvtx
        idx_ori = contr%svertex(ivtx)
        do jvtx = 1, ivtx-1
          if (idx_ori.eq.contr%svertex(jvtx)) then
            svertex_new(ivtx) = svertex_new(jvtx)
            cycle iloop
          end if
        end do
        idx = idx+1
        svertex_new(ivtx) = idx
      end do iloop
      ! set correct number of super-vertices
      contr%nsupvtx = idx

      contr%svertex(1:contr%nvtx) = svertex_new(1:contr%nvtx)

      contr%joined(0,1:contr%nsupvtx) = 0
      do ivtx = 1, contr%nvtx
        npvtx = contr%joined(0,contr%svertex(ivtx))
        npvtx = npvtx+1
        contr%joined(0,contr%svertex(ivtx)) = npvtx
        contr%joined(npvtx,contr%svertex(ivtx)) = ivtx
      end do

      return
      end
