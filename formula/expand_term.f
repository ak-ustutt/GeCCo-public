*----------------------------------------------------------------------*
      subroutine expand_term(fl_expand,nterms,
     &                       njoined,f_term,fpl_intm,force,op_info)
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
      include 'ifc_operators.h'

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
      logical, intent(in)::
     &     force
      integer, intent(in) ::
     &     njoined

      type(contraction) ::
     &     proto
      logical ::
     &     adj_intm, ok
      integer ::
     &     nvtx, narc, narc0, ivtx, jvtx, kvtx, iarc,
     &     iop_intm, iblk_intm, iblk, iadd, nj_tgt, ieqvfac,
     &     idx, ioff
      integer ::
     &     occ_temp(ngastp,2)
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
     &     ipos_vtx(:), ol_map(:)
      logical, pointer ::
     &     fix_vtx(:), fix_vtx_intm(:), found(:)
      integer, pointer ::
     &     svertex(:), neqv(:), idx_eqv(:,:)
      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:)

      integer, external ::
     &     vtx_in_contr, ifac, idxlist, ielsqsum, get_eqvfac

      if (force) call quit(1,'expand_term','obsolete "force" feature?')

      if (ntest.ge.100) then
        write(luout,*) '========================='
        write(luout,*) ' Here speaks expand_term'
        write(luout,*) '========================='
      end if

      iop_intm  = fpl_intm%item%contr%idx_res
      iblk_intm = fpl_intm%item%contr%iblk_res
      adj_intm  = fpl_intm%item%contr%dagger
      op_intm => op_info%op_arr(iop_intm)%op
c      njoined = op_intm%njoined

      allocate(ipos_vtx(njoined))

      term => f_term%contr
      call get_vtx_in_contr(ipos_vtx,iop_intm,adj_intm,njoined,term)

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
c        narc = max(term%narc,(term%nvtx-njoined)*(term%nvtx-njoined-1))
c     &       + intm%narc
        narc = max(term%narc,(nvtx+2)*(nvtx+1)) 

        ! make a map: which vertex goes where ...

        ! get supter-vertex map for intermediate
        allocate(svmap(intm%nvtx))
        if (njoined.eq.1) then
          svmap(1:intm%nvtx) = 1
        else
          call svmap4contr2(svmap,intm,ok)
          if (.not.ok) call quit(1,'expand_term',
     &             'cannot handle non-unique svmap yet')
        end if

        allocate(vtxmap(nvtx))

        call joinmap4contr(vtxmap,term,
     &                     -1,ipos_vtx,
     &                     svmap,intm%nvtx,njoined)

        deallocate(svmap)

        ! assemble proto contraction
        call resize_contr(proto,nvtx+2,narc,0,0)

        if (term%nvtx.gt.0) allocate(ivtx_term_reo(term%nvtx))
        if (intm%nvtx.gt.0) allocate(ivtx_intm_reo(intm%nvtx))
        if (term%nvtx.gt.0) ivtx_term_reo(1:term%nvtx) = 0
        if (intm%nvtx.gt.0) ivtx_intm_reo(1:intm%nvtx) = 0

        ! copy general info
        proto%idx_res = term%idx_res
        proto%iblk_res = term%iblk_res
        proto%fac = term%fac*intm%fac
        proto%nvtx = nvtx+2
        proto%nfac = 0

        ! result vertex: first and last vertex
        nj_tgt = op_info%op_arr(proto%idx_res)%op%njoined
        ioff = (proto%iblk_res-1)*nj_tgt 
        if (nj_tgt.ne.1)
     &      call quit(1,'expand_term','nj>1 for res. not debugged yet')
        proto%vertex(1)%idx_op = proto%idx_res
        proto%vertex(1)%iblk_op = proto%iblk_res ! dummy
        proto%vertex(1)%dagger  =.false.
        proto%svertex(1) = 1
        proto%vertex(nvtx+2)%idx_op = proto%idx_res
        proto%vertex(nvtx+2)%iblk_op = proto%iblk_res+1 ! dummy
        proto%vertex(nvtx+2)%dagger  =.false.
        proto%svertex(nvtx+2) = 1

        ! set vertices and reordering arrays
        do ivtx = 1, nvtx
          if (vtxmap(ivtx).gt.0) then
            proto%vertex(ivtx+1) = term%vertex(vtxmap(ivtx))
            proto%svertex(ivtx+1) = term%svertex(vtxmap(ivtx))+1
            ivtx_term_reo(vtxmap(ivtx)) = ivtx
          else
            proto%vertex(ivtx+1) = intm%vertex(-vtxmap(ivtx))
            proto%svertex(ivtx+1) =
     &           term%nsupvtx + intm%svertex(-vtxmap(ivtx))+1
            ivtx_intm_reo(-vtxmap(ivtx)) = ivtx
          end if
        end do

        ! set up correct super-vertex info
        call update_svtx4contr(proto)

        ! copy arcs
        narc = 0
        do iarc = 1, term%narc
          if (idxlist(term%arc(iarc)%link(1),ipos_vtx,njoined,1).gt.0)
     &         then
            narc = narc+1
            proto%arc(narc)%link(1)=0
            proto%arc(narc)%link(2)=
     &           ivtx_term_reo(term%arc(iarc)%link(2))+1
            proto%arc(narc)%occ_cnt = term%arc(iarc)%occ_cnt
          else if
     &       (idxlist(term%arc(iarc)%link(2),ipos_vtx,njoined,1).gt.0)
     &           then
            narc = narc+1
            proto%arc(narc)%link(1)=
     &           ivtx_term_reo(term%arc(iarc)%link(1))+1
            proto%arc(narc)%link(2)=0
            proto%arc(narc)%occ_cnt = term%arc(iarc)%occ_cnt
          else
c          if (idxlist(term%arc(iarc)%link(1),ipos_vtx,njoined,1).gt.0
c     &      .or.idxlist(term%arc(iarc)%link(2),ipos_vtx,njoined,1).gt.0)
c     &      cycle  
          narc = narc+1
          proto%arc(narc)%link(1)=ivtx_term_reo(term%arc(iarc)%link(1))
     &                            +1
          proto%arc(narc)%link(2)=ivtx_term_reo(term%arc(iarc)%link(2))
     &                            +1
          proto%arc(narc)%occ_cnt = term%arc(iarc)%occ_cnt
          end if
        end do
        narc0 = narc
        do iarc = 1, intm%narc
          narc = narc+1
          proto%arc(narc)%link(1)=ivtx_intm_reo(intm%arc(iarc)%link(1))
     &                            +1
          proto%arc(narc)%link(2)=ivtx_intm_reo(intm%arc(iarc)%link(2))
     &                            +1
          proto%arc(narc)%occ_cnt = intm%arc(iarc)%occ_cnt
        end do
        ! all new connections must be between term and intm; 
        ! we have to add 0-contractions to avoid intra-term connections:
        do jvtx = 1, nvtx
          if (vtxmap(jvtx).lt.0) cycle
          kloop: do kvtx = jvtx+1, nvtx
            if (vtxmap(kvtx).lt.0) cycle kloop
            do iarc = 1, narc0
              if (proto%arc(iarc)%link(1).eq.jvtx+1.and.
     &            proto%arc(iarc)%link(2).eq.kvtx+1) cycle kloop
            end do
            narc = narc+1
            proto%arc(narc)%link(1) = jvtx+1
            proto%arc(narc)%link(2) = kvtx+1
            proto%arc(narc)%occ_cnt = 0
          end do kloop
        end do

        ! external arcs from term (not from intermediate, though
        ! this can be a problem if the intermediate is an nj<>0 op.)
        ! and [0] contractions to avoid unwanted external arcs
        occ_temp = 0
        do idx = 1, nj_tgt
          occ_temp = occ_temp 
     &             + op_info%op_arr(proto%idx_res)%op%ihpvca_occ(
     &                     1:ngastp,1:2,ioff+idx)
        end do
        if (iocc_nonzero(occ_temp)) then
          allocate(found(nvtx))
          found = .false.
          do idx = 1, term%nxarc
            ! ignore open lines of intermediate for now
            ! might be problematic for njoined.gt.1 !!!
            if (idxlist(term%xarc(idx)%link(1),ipos_vtx,njoined,1)
     &          .gt.0) cycle
            ivtx = ivtx_term_reo(term%xarc(idx)%link(1))
            if (found(ivtx)) call quit(1,'expand_term',
     &                  'more than one external arc for one vertex?')
            found(ivtx) = .true.
            ! connect excitation part to first vertex ...
            occ_temp = iocc_xdn(1,term%xarc(idx)%occ_cnt)
            narc = narc + 1
            proto%arc(narc)%link(1) = 1
            proto%arc(narc)%link(2) = 1 + ivtx
            proto%arc(narc)%occ_cnt = iocc_dagger(occ_temp)
            ! ... and deexcitation part to last vertex
            occ_temp = iocc_xdn(2,term%xarc(idx)%occ_cnt)
            narc = narc + 1
            proto%arc(narc)%link(1) = 1 + ivtx
            proto%arc(narc)%link(2) = nvtx + 2
            proto%arc(narc)%occ_cnt = occ_temp
          end do
          do jvtx = 1, term%nvtx
            ivtx = ivtx_term_reo(jvtx)
            if (ivtx.ne.0.and..not.found(ivtx)) then
              narc = narc + 1
              proto%arc(narc)%link(1) = 1
              proto%arc(narc)%link(2) = 1 + ivtx
              proto%arc(narc)%occ_cnt = 0
              narc = narc + 1
              proto%arc(narc)%link(1) = 1 + ivtx
              proto%arc(narc)%link(2) = nvtx + 2
              proto%arc(narc)%occ_cnt = 0
            end if
          end do
          deallocate(found)
        end if

        proto%narc = narc

        if (ntest.ge.100) then
          write(luout,*) 'generated proto-contraction:'
          call prt_contr2(luout,proto,op_info)
        end if
        
        ! a bit of bureaucracy ...
        allocate(occ_vtx(ngastp,2,nvtx+2),fix_vtx(nvtx+2),
     &           ol_map(nvtx+2))
        ol_map(1) = 1
        ol_map(nvtx+2) = -1
        ol_map(2:nvtx+1) = 0
        fix_vtx = .true.     ! "fix" all vertices -> ieqvfac will be 1

        ! get info about equivalent operators of the intermediate,
        ! "unfix" 2nd,3rd,... vertex of multivertex operators
        ! and apply resulting correction factor to proto contraction
        ! --> generates less (identical) terms
        allocate(vtx(intm%nvtx),topo(intm%nvtx,intm%nvtx),
     &           xlines(intm%nvtx,njoined),svertex(intm%nvtx),
     &           neqv(intm%nvtx),idx_eqv(intm%nvtx,intm%nvtx),
     &           fix_vtx_intm(intm%nvtx))
        call pack_contr(svertex,vtx,topo,xlines,intm,njoined)
        do ivtx = 1, intm%nvtx
          fix_vtx_intm(ivtx) = 
     &          idxlist(svertex(ivtx),svertex,intm%nvtx,1).eq.ivtx
          fix_vtx(ivtx_intm_reo(ivtx)+1) = fix_vtx_intm(ivtx)
        end do
        call set_eqv_map(neqv,idx_eqv,
     &                  vtx,svertex,topo,xlines,intm%nvtx,njoined)
        ieqvfac = get_eqvfac(neqv,fix_vtx_intm,intm%nvtx)
        proto%fac = proto%fac*dble(ieqvfac)
        deallocate(vtx,topo,xlines,svertex,neqv,idx_eqv,fix_vtx_intm)

        call occvtx4contr(1,occ_vtx,proto,op_info)

        ! daggered excitation and deexcitation parts of result operator
        occ_vtx(1:ngastp,1:2,1) = iocc_dagger(iocc_xdn(1,
     &       op_info%op_arr(proto%idx_res)%op%ihpvca_occ(
     &       1:ngastp,1:2,proto%iblk_res)))
        occ_vtx(1:ngastp,1:2,nvtx+2) = iocc_dagger(iocc_xdn(2,
     &       op_info%op_arr(proto%idx_res)%op%ihpvca_occ(
     &       1:ngastp,1:2,proto%iblk_res)))

        if(force)then
          occ_temp(1:ngastp,1)=0
          occ_temp(1:ngastp,2)=occ_vtx(1:ngastp,2,2)
          if(ielsqsum(occ_temp,ngastp*2).ne.0)then
            narc=narc+1
            proto%arc(narc)%link(1)=2
            proto%arc(narc)%link(2)=nvtx+1
            proto%arc(narc)%occ_cnt=occ_temp
            proto%narc=narc
            call prt_contr2(luout,proto,op_info)
          endif
        endif

        ! ... and go! get all possible connections
        call gen_contr4(.false.,fl_expand_pnt,proto,
     &       fix_vtx,occ_vtx,ol_map,op_info)
        do
          if (fl_expand_pnt%command.eq.command_end_of_formula) exit
          nterms = nterms+1
          fl_expand_pnt => fl_expand_pnt%next
        end do
        if (ntest.ge.100) then
          write(luout,*) 'currently ',nterms,' terms'
        end if
        deallocate(occ_vtx,fix_vtx,ol_map)

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
