*----------------------------------------------------------------------*
      subroutine proto2contr(contr,proto,ol_map,occ_vtx,op_info)
*----------------------------------------------------------------------*
*     input: a proto contraction in its final state (i.e. all
*       arcs are set such that all vertices are connected, but still
*       it contains open line vertices which must be removed and any
*       arc to such a vertex must be transformed to an xarc)
*       an array telling us which of the proto-vertices corresponds
*       to which open line vertex
*     output: the actual contraction correctly set up with xarcs
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(out) ::
     &     contr
      type(contraction), intent(in) ::
     &     proto
      integer, intent(in) ::
     &     ol_map(proto%nvtx), occ_vtx(ngastp,2,proto%nvtx)
      type(operator_info) ::
     &     op_info

      logical ::
     &     is_xarc, is_xarc_dag
      integer ::
     &     nvtx_p, narc_p, nvtx, narc, nxarc, ivtx, ivtx_p, 
     &     iarc_p, iarc, ixarc, jxarc, jvtx, ihelp
      integer ::
     &     vtx_map(proto%nvtx)
      integer, pointer ::
     &     svmap(:)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'proto2contr')
        write(luout,*) 'on input:'
        call prt_contr3(luout,proto,occ_vtx)
      end if

      nvtx_p = proto%nvtx
      narc_p = proto%narc

      ! count number of vertices, arcs, xarcs
      ! and set up a mapping vtx_map(ivtx_p) -> new_vtx/result_vtx
      nvtx = 0
      ivtx = 0
      do ivtx_p = 1, nvtx_p
        if (ol_map(ivtx_p).ne.0) then
          vtx_map(ivtx_p) = abs(ol_map(ivtx_p))
        else
          ivtx = ivtx+1
          vtx_map(ivtx_p) = ivtx
        end if
      end do
      nvtx = ivtx

      ! count number of arcs, xarcs
      narc = narc_p
      nxarc = 0
      do iarc_p = 1, narc_p
        if (ol_map(proto%arc(iarc_p)%link(1)).ne.0.or.
     &      ol_map(proto%arc(iarc_p)%link(2)).ne.0) then
          narc = narc-1
          nxarc = nxarc+1
        end if
      end do

      call resize_contr(contr,nvtx,narc,nxarc,0)

      contr%nvtx = nvtx
      contr%narc = narc
      contr%nxarc = nxarc
      contr%nfac = 0
      contr%fac  = proto%fac
      contr%idx_res = proto%idx_res
      contr%iblk_res = proto%iblk_res
      contr%dagger = proto%dagger

      ! set vertices and svertex array
      do ivtx_p = 1, nvtx_p
        if (ol_map(ivtx_p).ne.0) cycle
        ivtx = vtx_map(ivtx_p)
        contr%vertex(ivtx) = proto%vertex(ivtx_p)
        contr%svertex(ivtx) = proto%svertex(ivtx_p)
      end do

      call update_svtx4contr(contr)

      ! set arcs and xarcs
      iarc = 0
      ixarc = 0
      iarc_p_loop: do iarc_p = 1, narc_p
        ivtx = proto%arc(iarc_p)%link(1)
        jvtx = proto%arc(iarc_p)%link(2)
        is_xarc = ol_map(jvtx).ne.0
        is_xarc_dag = ol_map(ivtx).ne.0
        if (is_xarc.and.is_xarc_dag)
     &       call quit(1,'proto2contr','this must not happen')
        if (is_xarc.or.is_xarc_dag) then
          if (is_xarc) then
            ivtx = vtx_map(ivtx)
            jvtx = vtx_map(jvtx)
          end if
          if (is_xarc_dag) then 
            ihelp = vtx_map(jvtx)
            jvtx = vtx_map(ivtx)
            ivtx = ihelp
          end if
          ! scan for entry with same ivtx, jvtx
          do jxarc = 1, ixarc
            if (contr%xarc(jxarc)%link(1).eq.ivtx.and.
     &          contr%xarc(jxarc)%link(2).eq.jvtx) then
              if (is_xarc_dag) then
                contr%xarc(jxarc)%occ_cnt =
     &               contr%xarc(jxarc)%occ_cnt +
     &               iocc_dagger(proto%arc(iarc_p)%occ_cnt)
              else
                contr%xarc(jxarc)%occ_cnt =
     &               contr%xarc(jxarc)%occ_cnt +
     &               proto%arc(iarc_p)%occ_cnt
              end if
              cycle iarc_p_loop
            end if
          end do
          ixarc = ixarc+1
          contr%xarc(ixarc)%link(1) = ivtx
          contr%xarc(ixarc)%link(2) = jvtx
          if (is_xarc_dag) then
            contr%xarc(ixarc)%occ_cnt =
     &           iocc_dagger(proto%arc(iarc_p)%occ_cnt)
          else
            contr%xarc(ixarc)%occ_cnt =
     &           proto%arc(iarc_p)%occ_cnt
          end if
        else
          iarc = iarc+1
          contr%arc(iarc)%link(1) = vtx_map(proto%arc(iarc_p)%link(1))
          contr%arc(iarc)%link(2) = vtx_map(proto%arc(iarc_p)%link(2))
          contr%arc(iarc)%occ_cnt = proto%arc(iarc_p)%occ_cnt
        end if
      end do iarc_p_loop
      
      ! set actual number of xarcs:
      contr%nxarc = ixarc

      if (ntest.ge.100) then
        write(luout,*) 'on exit:'
        call prt_contr2(luout,contr,op_info)
      end if

      return
      end
