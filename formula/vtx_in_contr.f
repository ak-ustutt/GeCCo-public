*----------------------------------------------------------------------*
      integer function vtx_in_contr(idxop,contr)
*----------------------------------------------------------------------*
*     return the first place of operator idxop in contr
*----------------------------------------------------------------------*
      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     idxop

      integer ::
     &     jvtx, kvtx

      kvtx = 0
      do jvtx = 1, contr%nvtx
        if (contr%vertex(jvtx)%idx_op.eq.idxop) then
          kvtx = jvtx
          exit
        end if
      end do

      vtx_in_contr = kvtx

      return
      end
