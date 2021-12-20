*----------------------------------------------------------------------*
      subroutine transpose_contr(contr,op_info,multi)
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
      logical, intent(in), optional ::
     &     multi

      integer ::
     &     nvtx, narc, nxarc, ivtx, jvtx, ihelp, isuper, idx,
     &     idxst, idxnd, njoined, njoined_res, iarc, ifac
      logical ::
     &     skip, reo, use_multi, scal

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
      logical, external ::
     &     occ_is_diag_blk

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'transpose_contr')
        write(lulog,*) 'on input:'
        call prt_contr2(lulog,contr,op_info)
      end if

      if (present(multi)) then
        use_multi = multi
      else
        use_multi = .false. 
      endif

      nvtx = contr%nvtx
      vertex => contr%vertex
      svertex => contr%svertex
      narc = contr%narc
      arc => contr%arc
      nxarc = contr%nxarc
      xarc => contr%xarc
      
      ! A new logical 'scal' has been used of dagger is not necessary
      ! while transposing the formula for operator 
      ! corresponding to CI coefficients
      if(use_multi) then 
        scal = op_info%op_arr(contr%idx_res)%op%n_occ_cls.gt.1
      else 
        scal = .true.
      end if

      ! set transposition info for result
      if (scal) contr%dagger = .not.contr%dagger

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
c dbg
c        print *,'ivtx',ivtx
c dbg
        ! handle super-vertex only once
        isuper = svertex(ivtx)          
        skip = .false.
        do jvtx = 1, ivtx-1 ! previously visited?
          skip = skip.or.(isuper.eq.svertex(jvtx))
        end do
        if (skip) cycle

        if (vertex(ivtx)%dagger) then
          ! if the vertex was transpose, we know that iblk
          ! refers to the un-transpose operator, so we merely
          ! change the dagger label
          vertex(ivtx)%dagger = .false.
          do jvtx = ivtx+1, nvtx
            if (svertex(jvtx).ne.isuper) cycle
            vertex(jvtx)%dagger = .false.
          end do
        else
          ! find out whether the transposed block exists:
          op => op_info%op_arr(vertex(ivtx)%idx_op)%op
          njoined = op%njoined
          idxnd = vertex(ivtx)%iblk_op
          idxst = idxnd-njoined+1
c          idxst = vertex(ivtx)%iblk_op
c          idxnd = idxst+njoined-1
          occ => op%ihpvca_occ(1:ngastp,1:2,idxst:idxnd)
c dbg
c          print *,'idxst, idxnd: ',idxst,idxnd
c          call wrt_occ_n(6,occ,op%njoined)
c dbg
          ! ... but only if we have off-diagonal blocks
          ! *and* if we know that the operator is hermitian
          if (.not.occ_is_diag_blk(occ,njoined).and.
     &         abs(op%hermitian).eq.1)
     &         then
            idx = iblk_occ(occ,.true.,op,
     &                     op%blk_version((idxst-1)/njoined+1))
            ! we have to change the overall sign, if the
            ! operator is anti-hermitian:
            if (op%hermitian.eq.-1)
     &           contr%fac = -contr%fac
          else
            idx = -1
          end if
c dbg
c          print *,'idx = ',idx,occ_is_diag_blk(occ,njoined),op%hermitian
c dbg
          if (idx.gt.0) then
            idx = (idx-1)*njoined+1
            do jvtx = ivtx, nvtx
              if (svertex(jvtx).ne.isuper) cycle
              vertex(jvtx)%iblk_op = idx
              idx = idx+1
            end do
          else
            do jvtx = ivtx, nvtx
              if (svertex(jvtx).ne.isuper) cycle
              vertex(jvtx)%dagger = .true.
            end do
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
c        call prt_contr2(lulog,contr,op_info)
c dbg      
      allocate(occ_vtx(ngastp,2,nvtx),vtx_reo(nvtx),fix_vtx(nvtx))
      fix_vtx = .true. ! not important
      call occvtx4contr(1,occ_vtx,contr,op_info)
      call topo_contr(ifac,reo,vtx_reo,contr,occ_vtx,fix_vtx)
      call canon_contr(contr,reo,vtx_reo)
      deallocate(occ_vtx,vtx_reo,fix_vtx)

      if (ntest.ge.100) then
        write(lulog,*) 'on output:'
        call prt_contr2(lulog,contr,op_info)
      end if

      return
      end
