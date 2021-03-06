*----------------------------------------------------------------------*
      integer function vtx_in_contr(idxop,adj,vtxst,contr)
*----------------------------------------------------------------------*
*     return the first place of operator idxop in contr
*----------------------------------------------------------------------*
      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     idxop, vtxst
      logical, intent(in) ::
     &     adj

      integer ::
     &     jvtx, kvtx

      kvtx = 0
      do jvtx = vtxst, contr%nvtx
        if (contr%vertex(jvtx)%idx_op.eq.idxop .and.
     &     (contr%vertex(jvtx)%dagger.eqv.adj)) then
          kvtx = jvtx
          exit
        end if
      end do

      vtx_in_contr = kvtx

      return
      end
