*----------------------------------------------------------------------*
      integer function maxop_in_contr(contr)
*----------------------------------------------------------------------*
*     return the maximum operator index occurring in contr
*----------------------------------------------------------------------*
      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer ::
     &     idxmax, ivtx

      idxmax = contr%idx_res
      do ivtx = 1, contr%nvtx
        idxmax = max(idxmax,contr%vertex(ivtx)%idx_op)
      end do

      maxop_in_contr = idxmax

      return
      end
