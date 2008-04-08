*----------------------------------------------------------------------*
      subroutine transpose_contr(contr,op_info)
*----------------------------------------------------------------------*
*     transpose a contraction
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction) ::
     &     contr
      type(operator_info) ::
     &     op_info

      integer ::
     &     nvtx, narc, nxarc, ivtx, jvtx, ihelp, isuper, idx,
     &     idxst, idxnd, njoined, njoined_res, iarc, ifac
      logical ::
     &     skip, reo

      type(operator), pointer ::
     &     op
      type(cntr_vtx) ::
     &     vhelp
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(cntr_arc), pointer ::
     &     arc(:), xarc(:)
      integer, pointer ::
     &     svertex(:), occ_vtx(:,:,:), vtx_reo(:), occ(:,:,:)
      logical, pointer ::
     &     fix_vtx(:)
      
      integer, external ::
     &     iblk_occ
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'transpose_contr')
        write(luout,*) 'on input:'
        call prt_contr2(luout,contr,op_info)
      end if

      nvtx = contr%nvtx
      vertex => contr%vertex
      svertex => contr%svertex
      narc = contr%narc
      arc => contr%arc
      nxarc = contr%nxarc
      xarc => contr%xarc
      
      ! set transposition info for result
      contr%dagger = .not.contr%dagger

      ! process vertices
      ! reverse sequence
      do ivtx = 1, nvtx/2
        vhelp = vertex(ivtx)
        vertex(ivtx) = vertex(nvtx-ivtx+1)
        vertex(nvtx-ivtx+1) = vhelp
      end do
      do ivtx = 1, nvtx/2
        ihelp = svertex(ivtx)
        svertex(ivtx) = svertex(nvtx-ivtx+1)
        svertex(nvtx-ivtx+1) = ihelp
      end do

      ! update transposition info or block of operator
      do ivtx = 1, nvtx
        if (vertex(ivtx)%dagger) then
          ! if the vertex was transpose, we know that iblk
          ! refers to the un-transpose operator, so we merely
          ! change the dagger label
          vertex(ivtx)%dagger = .false.
        else
          ! handle super-vertex only once
          isuper = svertex(ivtx)          
          skip = .false.
          do jvtx = 1, ivtx-1 ! previously visited?
            skip = skip.or.(isuper.eq.svertex(jvtx))
          end do
          if (skip) cycle
          ! find out whether the transposed block exists:
          op => op_info%op_arr(vertex(ivtx)%idx_op)%op
          njoined = op%njoined
          idxst = vertex(ivtx)%iblk_op
          idxnd = idxst+njoined-1
          occ => op%ihpvca_occ(1:ngastp,1:2,idxst:idxnd)
c dbg
          call wrt_occ_n(6,occ,op%njoined)
c dbg
          idx = iblk_occ(occ,.true.,op)
c dbg
          print *,'idx = ',idx
c dbg
          if (idx.gt.0) then
            idx = (idx-1)*njoined+1
            do jvtx = ivtx, nvtx
              if (svertex(jvtx).ne.isuper) cycle
              vertex(jvtx)%iblk_op = idx
              idx = idx+1
            end do
          else
            vertex(ivtx)%dagger = .true.
          end if
        end if
      end do

      call update_svtx4contr(contr)

      ! update arcs:
      do iarc = 1, narc
        ihelp =             nvtx-arc(iarc)%link(1)+1
        arc(iarc)%link(1) = nvtx-arc(iarc)%link(2)+1
        arc(iarc)%link(2) = ihelp
        ! no adjoint here !
        arc(iarc)%occ_cnt = arc(iarc)%occ_cnt
c        arc(iarc)%occ_cnt = iocc_dagger(arc(iarc)%occ_cnt)
      end do

      if (nxarc.gt.0) then
        op => op_info%op_arr(contr%idx_res)%op
        njoined_res = op%njoined
      end if

      ! update xarcs:
      do iarc = 1, nxarc
        xarc(iarc)%link(1) = nvtx-xarc(iarc)%link(1)+1
        xarc(iarc)%link(2) = njoined_res-xarc(iarc)%link(2)+1
        xarc(iarc)%occ_cnt = iocc_dagger(xarc(iarc)%occ_cnt)
      end do
c dbg
        call prt_contr2(luout,contr,op_info)
c dbg      
      allocate(occ_vtx(ngastp,2,nvtx),vtx_reo(nvtx),fix_vtx(nvtx))
      fix_vtx = .true. ! not important
      call occvtx4contr(1,occ_vtx,contr,op_info)
      call topo_contr(ifac,reo,vtx_reo,contr,occ_vtx,fix_vtx)
      call canon_contr(contr,reo,vtx_reo)
      deallocate(occ_vtx,vtx_reo,fix_vtx)

      if (ntest.ge.100) then
        write(luout,*) 'on output:'
        call prt_contr2(luout,contr,op_info)
      end if

      return
      end
