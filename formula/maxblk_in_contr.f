*----------------------------------------------------------------------*
      integer function maxblk_in_contr(contr)
*----------------------------------------------------------------------*
*     return the maximum block index occurring in contr
*----------------------------------------------------------------------*
      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer ::
     &     idxmax, ivtx

      type(cntr_vtx), pointer ::
     &     vtx(:)

      vtx => contr%vertex
      idxmax = contr%iblk_res
      do ivtx = 1, contr%nvtx
        idxmax = max(idxmax,vtx(ivtx)%iblk_op)
      end do

      maxblk_in_contr = idxmax

      return
      end
