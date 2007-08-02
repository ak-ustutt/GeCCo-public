*----------------------------------------------------------------------*
      subroutine expand_term(fl_expand,nterms,f_term,fpl_intm,op_info)
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
     &     iop_intm, iblk_intm, iadd
      type(contraction), pointer ::
     &     term, intm
      type(formula_item_list), pointer ::
     &     fpl_intm_pnt
      type(formula_item), pointer ::
     &     fl_expand_pnt
      integer, pointer ::
     &     ivtx_reo(:), occ_vtx(:,:,:)
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     vtx_in_contr, ifac

      if (ntest.ge.100) then
        write(luout,*) '========================='
        write(luout,*) ' Here speaks expand_term'
        write(luout,*) '========================='
      end if

      iop_intm  = fpl_intm%item%contr%idx_res
      iblk_intm = fpl_intm%item%contr%iblk_res

      term => f_term%contr
      ivtx = vtx_in_contr(iop_intm,term)

      if (term%vertex(ivtx)%iblk_op.ne.iblk_intm)
     &     call quit(1,'expand_term','inconsistency')

      call init_contr(proto)
      
      nterms = 0
      fpl_intm_pnt => fpl_intm
      fl_expand_pnt => fl_expand
      ! loop over intermediate blocks
      do

        intm => fpl_intm_pnt%item%contr

        ! assemble proto contraction
        nvtx = term%nvtx-1 + intm%nvtx
        narc = max(term%narc,ifac(term%nvtx-1)) + intm%narc
        call resize_contr(proto,nvtx,narc,0)

        ! copy general info
        proto%idx_res = term%idx_res
        proto%iblk_res = term%iblk_res
        proto%fac = term%fac*intm%fac
        proto%nvtx = nvtx
        proto%nfac = 0
        ! copy vertices
        do jvtx = 1, ivtx-1
          proto%vertex(jvtx) = term%vertex(jvtx)
        end do
        kvtx = 0
        do jvtx = ivtx, ivtx-1+intm%nvtx
          kvtx = kvtx+1
          proto%vertex(jvtx) = intm%vertex(kvtx)
        end do
        kvtx = ivtx
        do jvtx = ivtx+intm%nvtx, nvtx
          kvtx = kvtx+1
          proto%vertex(jvtx) = term%vertex(kvtx)
        end do
        ! copy arcs
        narc = 0
        do iarc = 1, term%narc
          if (term%arc(iarc)%link(1).eq.ivtx .or.
     &        term%arc(iarc)%link(2).eq.ivtx) cycle
          narc = narc+1
          iadd = 0
          if (term%arc(iarc)%link(1).gt.ivtx) iadd = intm%nvtx-1
          proto%arc(narc)%link(1) = term%arc(iarc)%link(1)+iadd
          iadd = 0
          if (term%arc(iarc)%link(2).gt.ivtx) iadd = intm%nvtx-1
          proto%arc(narc)%link(2) = term%arc(iarc)%link(2)+iadd
          proto%arc(narc)%occ_cnt = term%arc(iarc)%occ_cnt
        end do
        narc0 = narc
        do iarc = 1, intm%narc
          narc = narc+1
          proto%arc(narc)%link(1) = intm%arc(iarc)%link(1)+ivtx-1
          proto%arc(narc)%link(2) = intm%arc(iarc)%link(2)+ivtx-1
          proto%arc(narc)%occ_cnt = intm%arc(iarc)%occ_cnt
        end do
        ! all new connections must be between term and intm; 
        ! we have to add 0-contractions to avoid intra-term connections:
        do jvtx = 1, ivtx-1
          kloop: do kvtx = ivtx+intm%nvtx, nvtx
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

        if (.not.associated(fpl_intm_pnt%next)) exit
        fpl_intm_pnt => fpl_intm_pnt%next
      end do

      call dealloc_contr(proto)

      return
      end
