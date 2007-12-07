*----------------------------------------------------------------------*
      subroutine expand_term(fl_expand,nterms,
     &                       njoined,f_term,fpl_intm,op_info)
*----------------------------------------------------------------------*
*     expand O1.O2...Int...On (on f_term as formula list) to
*            O1.O2...I1aI1b...On + O1.O2...I2...On + ...
*     where the definition Int = I1aI1b... + I2... + ...
*     is given as formula pointer list
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_array.h'
      include 'def_formula_item_list.h'

      integer, parameter ::
     &     ntest = 00
      
      type(formula_item), target, intent(out) ::
     &     fl_expand
      integer, intent(out) ::
     &     nterms
      type(formula_item), target, intent(in) ::
     &     f_term
      type(formula_item_list), target, intent(in) ::
     &     fpl_intm
      type(operator_info) ::
     &     op_info

      type(contraction) ::
     &     proto
      integer ::
     &     nvtx, narc, narc0, ivtx, jvtx, kvtx, iarc,
     &     iop_intm, iblk_intm, iblk, iadd, njoined
      type(contraction), pointer ::
     &     term, intm
      type(formula_item_list), pointer ::
     &     fpl_intm_pnt
      type(formula_item), pointer ::
     &     fl_expand_pnt
      type(operator), pointer ::
     &     op_intm
      integer, pointer ::
     &     ivtx_term_reo(:), ivtx_intm_reo(:),
     &     occ_vtx(:,:,:), svmap(:), vtxmap(:),
     &     ipos_vtx(:)
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     vtx_in_contr, ifac, idxlist

      if (ntest.ge.100) then
        write(luout,*) '========================='
        write(luout,*) ' Here speaks expand_term'
        write(luout,*) '========================='
      end if

      iop_intm  = fpl_intm%item%contr%idx_res
      iblk_intm = fpl_intm%item%contr%iblk_res
      op_intm => op_info%op_arr(iop_intm)%op
      njoined = op_intm%njoined

      allocate(ipos_vtx(njoined))

      term => f_term%contr
      call get_vtx_in_contr(ipos_vtx,iop_intm,njoined,term)

      iblk = term%vertex(ipos_vtx(1))%iblk_op
      if (njoined.gt.1)
     &     iblk = (iblk-1)/njoined + 1

      if (iblk.ne.iblk_intm)
     &     call quit(1,'expand_term','inconsistency')

      call init_contr(proto)
      
      nterms = 0
      fpl_intm_pnt => fpl_intm
      fl_expand_pnt => fl_expand
      ! loop over intermediate blocks
      do

        intm => fpl_intm_pnt%item%contr

        ! estimate new number of vertices and arcs:
        nvtx = term%nvtx-njoined + intm%nvtx
        narc = max(term%narc,(term%nvtx-njoined)*(term%nvtx-njoined-1))
     &       + intm%narc

        ! make a map: which vertex goes where ...

        ! get supter-vertex map for intermediate
        allocate(svmap(intm%nvtx))
        if (njoined.eq.1) then
          svmap(1:intm%nvtx) = 1
        else
          allocate(occ_vtx(ngastp,2,intm%nvtx+njoined))
          call occvtx4contr(0,occ_vtx,intm,op_info)
          call svmap4contr(svmap,intm,occ_vtx,njoined)
          deallocate(occ_vtx)
        end if

        allocate(vtxmap(nvtx))

        call joinmap4contr(vtxmap,term,
     &                     -1,ipos_vtx,
     &                     svmap,intm%nvtx,njoined)

        deallocate(svmap)

        ! assemble proto contraction
        call resize_contr(proto,nvtx,narc,0,0)

        if (term%nvtx.gt.0) allocate(ivtx_term_reo(term%nvtx))
        if (intm%nvtx.gt.0)  allocate(ivtx_intm_reo(intm%nvtx))
        if (term%nvtx.gt.0) ivtx_term_reo(1:term%nvtx) = 0
        if (intm%nvtx.gt.0)  ivtx_intm_reo(1:intm%nvtx) = 0

        ! copy general info
        proto%idx_res = term%idx_res
        proto%iblk_res = term%iblk_res
        proto%fac = term%fac*intm%fac
        proto%nvtx = nvtx
        proto%nfac = 0

        ! set vertices and reordering arrays
        do ivtx = 1, nvtx
          if (vtxmap(ivtx).gt.0) then
            proto%vertex(ivtx) = term%vertex(vtxmap(ivtx))
            proto%svertex(ivtx) = term%svertex(vtxmap(ivtx))
            ivtx_term_reo(vtxmap(ivtx)) = ivtx
          else
            proto%vertex(ivtx) = intm%vertex(-vtxmap(ivtx))
            proto%svertex(ivtx) =
     &           term%nsupvtx + intm%svertex(-vtxmap(ivtx))
            ivtx_intm_reo(-vtxmap(ivtx)) = ivtx
          end if
        end do

        ! set up correct super-vertex info
        call update_svtx4contr(proto)

        ! copy arcs
        narc = 0
        do iarc = 1, term%narc
          if (idxlist(term%arc(iarc)%link(1),ipos_vtx,njoined,1).gt.0
     &    .or.idxlist(term%arc(iarc)%link(2),ipos_vtx,njoined,1).gt.0)
     &         cycle
          narc = narc+1
          proto%arc(narc)%link(1)=ivtx_term_reo(term%arc(iarc)%link(1))
          proto%arc(narc)%link(2)=ivtx_term_reo(term%arc(iarc)%link(2))
          proto%arc(narc)%occ_cnt = term%arc(iarc)%occ_cnt
        end do
        narc0 = narc
        do iarc = 1, intm%narc
          narc = narc+1
          proto%arc(narc)%link(1)=ivtx_intm_reo(intm%arc(iarc)%link(1))
          proto%arc(narc)%link(2)=ivtx_intm_reo(intm%arc(iarc)%link(2))
          proto%arc(narc)%occ_cnt = intm%arc(iarc)%occ_cnt
        end do
        ! all new connections must be between term and intm; 
        ! we have to add 0-contractions to avoid intra-term connections:
        do jvtx = 1, nvtx
          if (vtxmap(jvtx).lt.0) cycle
          kloop: do kvtx = jvtx+1, nvtx
            if (vtxmap(kvtx).lt.0) cycle kloop
            do iarc = 1, narc0
              if (proto%arc(iarc)%link(1).eq.jvtx.and.
     &            proto%arc(iarc)%link(2).eq.kvtx) cycle kloop
            end do
            narc = narc+1
            proto%arc(narc)%link(1) = jvtx
            proto%arc(narc)%link(2) = kvtx
            proto%arc(narc)%occ_cnt = 0
          end do kloop
        end do

        proto%narc = narc

        if (ntest.ge.100) then
          write(luout,*) 'generated proto-contraction:'
          call prt_contr2(luout,proto,op_info)
        end if
        
        ! a bit of bureaucracy ...
        allocate(occ_vtx(ngastp,2,proto%nvtx+1),fix_vtx(proto%nvtx))
        fix_vtx = .true.     ! "fix" all vertices -> ieqvfac will be 1
        call occvtx4contr(0,occ_vtx,proto,op_info)

        ! ... and go! get all possible connections
        call gen_contr2(fl_expand_pnt,proto,fix_vtx,occ_vtx,op_info)
        do
          if (fl_expand_pnt%command.eq.command_end_of_formula) exit
          nterms = nterms+1
          fl_expand_pnt => fl_expand_pnt%next
        end do
        if (ntest.ge.100) then
          write(luout,*) 'currently ',nterms,' terms'
        end if
        deallocate(occ_vtx,fix_vtx)

        deallocate(vtxmap)
        if (term%nvtx.gt.0) deallocate(ivtx_term_reo)
        if (intm%nvtx.gt.0)  deallocate(ivtx_intm_reo)

        if (.not.associated(fpl_intm_pnt%next)) exit
        fpl_intm_pnt => fpl_intm_pnt%next
      end do

      call dealloc_contr(proto)

      deallocate(ipos_vtx)

      return
      end
