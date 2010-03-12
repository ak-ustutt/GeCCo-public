*----------------------------------------------------------------------*
      subroutine joinmap4contr(vtxmap,contr,
     &                         iop_intm,ipos_intm,
     &                         svmap_intm,nvtx_intm,njoined_intm)
*----------------------------------------------------------------------*
*     make a map: 
*     for each vertex in expanded contraction --
*       positive number == vertex in original contraction contr
*       negative number == vertex of expanded intermediate
*     the intermediate is either identified by iop_intm (if >=0) or
*     the positions are given by ipos_intm(1:njoined)
*     svmap_intm: cf. svmap4contr()
*     nvtx_intm: number of vertices of expanded intermediate
*     vtxmap must have the size  contr%nvtx-njoined_intm+nvtx_intm
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     nvtx_intm, njoined_intm,
     &     iop_intm, ipos_intm(njoined_intm), svmap_intm(nvtx_intm)
      integer, intent(out) ::
     &     vtxmap(*)

      integer ::
     &     ivtx_total, ijoin, ivtx_intm, ivtx, jvtx

      ivtx_total = 0
      ijoin = 1
      ivtx_intm = 1
      do ivtx = 1, contr%nvtx
        if ((iop_intm.ge.0.and.
     &                contr%vertex(ivtx)%idx_op.ne.iop_intm).or.
     &      (iop_intm.lt.0.and.
     &       ipos_intm(svmap_intm(ivtx_intm)).ne.ivtx)) then
          ! keep this vertex
          ivtx_total = ivtx_total+1
          vtxmap(ivtx_total) = ivtx
        else
          ! insert vertices of intermediate definition, which
          ! contribution to current primitive vertex of intermediate
          do jvtx = ivtx_intm, nvtx_intm
            if (svmap_intm(jvtx).le.ijoin) then
              ivtx_total = ivtx_total+1
              vtxmap(ivtx_total) = -jvtx
            else
              ivtx_intm = jvtx
              ijoin = ijoin+1
              exit
            end if
          end do
        end if
      end do
      
      return
      end
