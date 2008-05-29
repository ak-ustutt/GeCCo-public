*----------------------------------------------------------------------*
      subroutine join_contr2(contr_abc,contr_ac,contr_b,
     &     idxop_abc,iblk_abc,op_info)
*----------------------------------------------------------------------*
*     join contractions ac and b into one with result occupation abc
*
*     intermediate solution that can in principle generate
*     all possible cases if new connection is not unique but
*     currently restricts to the same behaviour as join_contr()
*     
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 000

      type(contraction), intent(out) ::
     &     contr_abc

      type(contraction), intent(in) ::
     &     contr_ac, contr_b
      integer, intent(in) ::
     &     idxop_abc, iblk_abc
      type(operator_info), intent(in) ::
     &     op_info
 
      logical ::
     &     reo
      integer ::
     &     nvtx_abc, nvtx_ac, nvtx_a, nvtx_b, nvtx_c,
     &     narc_abc, narc_abc0, narc_ac, narc_b, 
     &     idx, ivtx_abc, iarc, ivtx, jvtx, jvtx_last,
     &     nproto_ac, nproto_b, idxsuper,
     &     nsuper, njoined, isuper, njoined_abc
      type(formula_item) ::
     &     wrap
      type(operator), pointer ::
     &     opres
     
      integer, pointer ::
     &     ivtx_ac_reo(:), ivtx_b_reo(:),
     &     occ_vtx(:,:,:), svmap(:), ivtx_old(:)
      logical, pointer ::
     &     fix_vtx(:)

      type(cntr_arc), pointer ::
     &     arc(:)

      integer, external ::
     &     ifndmax, imltlist, idxlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'This is join_contr 2')
        write(luout,*) 'joining: AC, B'
        call prt_contr2(luout,contr_ac,op_info)
        call prt_contr2(luout,contr_b,op_info)
      end if

      nvtx_ac = contr_ac%nvtx
      narc_ac = contr_ac%narc

      ! count proto-vertices in AC:
      nproto_ac = 0
      do ivtx = 1, nvtx_ac
        if (contr_ac%vertex(ivtx)%idx_op.eq.0)
     &       nproto_ac = nproto_ac+1
      end do

      nvtx_b = contr_b%nvtx
      narc_b = contr_b%narc

      ! count proto-vertices in B (there should be none):
      nproto_b = 0
      do ivtx = 1, nvtx_b
        if (contr_b%vertex(ivtx)%idx_op.eq.0)
     &       nproto_b = nproto_b+1
      end do

      if (nvtx_b.le.0.or.nproto_b.ne.0)
     &     call quit(1,'join_contr2',
     &     'inserted contraction fragment must be non-empty and '//
     &     'must not contain proto-vertices!')

      ! get the super-vertex map for B
      allocate(svmap(nvtx_b))
      opres => op_info%op_arr(contr_b%idx_res)%op
      njoined = opres%njoined
      if (njoined.eq.1) then
        svmap(1:nvtx_b) = 1
      else
        allocate(occ_vtx(ngastp,2,nvtx_b+njoined))
        call occvtx4contr(0,occ_vtx,contr_b,op_info)
        call svmap4contr(svmap,contr_b,occ_vtx,njoined)
        deallocate(occ_vtx)
      end if
      ! largest index = number of super vertices
      nsuper = ifndmax(svmap,1,nvtx_b,1)
c dbg
c      print *,'nvtx_b: ',nvtx_b
c      print *,'svmap:  ',svmap(1:nvtx_b)
c      print *,'njoined: ',njoined
c      print *,'nsuper, nproto_ac: ',nsuper,nproto_ac
c dbg

      if (nsuper.ne.nproto_ac) then
        write(luout,*) 'joining: AC, B'
        call prt_contr2(luout,contr_ac,op_info)
        call prt_contr2(luout,contr_b,op_info)
        call quit(1,'join_contr2','incompatible contractions !')
      end if

      nvtx_abc = nvtx_ac-nproto_ac+nvtx_b

      ! generate a map: which vertex goes where
      allocate(ivtx_old(nvtx_abc))

      call joinmap4contr(ivtx_old,contr_ac,
     &                   0,-1,
     &                   svmap,nvtx_b,njoined)

      if (ntest.ge.1000) then
        write(luout,'(3x,a,10i5)') 'ivtx_old: ',ivtx_old(1:nvtx_abc)
      end if

      ! make some assumptions about the number of arcs in the 
      ! proto contraction
      narc_abc = min(
     &     narc_ac+narc_b+nvtx_ac*(nvtx_ac-1)+nvtx_b*(nvtx_b-1),
     &     nvtx_abc*(nvtx_abc-1))

      if (ntest.ge.1000) then
        write(luout,*) 'nvtx_ac, nvtx_b, nvtx_abc: ',
     &       nvtx_ac, nvtx_b, nvtx_abc
        write(luout,*) 'narc_ac, narc_b, narc_abc: ',
     &       narc_ac, narc_b, narc_abc
      end if

      call resize_contr(contr_abc,nvtx_abc,narc_abc,0,0)

      if (nvtx_ac.gt.0) allocate(ivtx_ac_reo(nvtx_ac))
      if (nvtx_b.gt.0)  allocate(ivtx_b_reo(nvtx_b))
      if (nvtx_ac.gt.0) ivtx_ac_reo(1:nvtx_ac) = 0
      if (nvtx_b.gt.0)  ivtx_b_reo(1:nvtx_b) = 0

      ! set prefactor
      contr_abc%fac = contr_ac%fac*contr_b%fac
      ! set result
      contr_abc%idx_res = idxop_abc
      contr_abc%iblk_res = iblk_abc
      contr_abc%dagger = contr_ac%dagger
      njoined_abc = op_info%op_arr(idxop_abc)%op%njoined

      ! set vertices and reordering arrays
      do ivtx = 1, nvtx_abc
        if (ivtx_old(ivtx).gt.0) then
          contr_abc%vertex(ivtx) = contr_ac%vertex(ivtx_old(ivtx))
          contr_abc%svertex(ivtx) = contr_ac%svertex(ivtx_old(ivtx))
          ivtx_ac_reo(ivtx_old(ivtx)) = ivtx
        else
          contr_abc%vertex(ivtx) = contr_b%vertex(-ivtx_old(ivtx))
          contr_abc%svertex(ivtx) =
     &         contr_ac%nsupvtx + contr_b%svertex(-ivtx_old(ivtx))
          ivtx_b_reo(-ivtx_old(ivtx)) = ivtx
        end if
      end do

      ! set up correct super-vertex info
      contr_abc%nvtx = nvtx_abc
      call update_svtx4contr(contr_abc)

      narc_abc = 0
      ! add all arcs from A and C, except the external ones
      do idx = 1, narc_ac
        if (contr_ac%arc(idx)%link(1).le.0.or.
     &      contr_ac%arc(idx)%link(2).le.0) cycle
        narc_abc = narc_abc + 1
        contr_abc%arc(narc_abc)%link(1) =
     &       ivtx_ac_reo(contr_ac%arc(idx)%link(1))
        contr_abc%arc(narc_abc)%link(2) =
     &       ivtx_ac_reo(contr_ac%arc(idx)%link(2))
        contr_abc%arc(narc_abc)%occ_cnt =
     &       contr_ac%arc(idx)%occ_cnt
      end do
      ! add all arcs from B
      do idx = 1, narc_b
        narc_abc = narc_abc + 1
        contr_abc%arc(narc_abc)%link(1) =
     &       ivtx_b_reo(contr_b%arc(idx)%link(1))
        contr_abc%arc(narc_abc)%link(2) =
     &       ivtx_b_reo(contr_b%arc(idx)%link(2))
        contr_abc%arc(narc_abc)%occ_cnt =
     &       contr_b%arc(idx)%occ_cnt
      end do
      ! add the external arcs from A and C
      ! which make gen_contr consider only special connections
      do idx = 1, narc_ac
        if (contr_ac%arc(idx)%link(1).le.0.or.
     &      contr_ac%arc(idx)%link(2).le.0) then
          narc_abc = narc_abc + 1
          if (contr_ac%arc(idx)%link(1).eq.0) then
            contr_abc%arc(narc_abc)%link(1) = 0
          else if (contr_ac%arc(idx)%link(1).lt.0) then
            ! fix for unique re-substitutions:
            idxsuper = -contr_ac%arc(idx)%link(1)
            if (imltlist(idxsuper,svmap,nvtx_b,1).eq.1) then
              contr_abc%arc(narc_abc)%link(1) =
     &             ivtx_b_reo(idxlist(idxsuper,svmap,nvtx_b,1))
            else
              contr_abc%arc(narc_abc)%link(1) = 0
            end if
          else
            contr_abc%arc(narc_abc)%link(1) = 
     &           ivtx_ac_reo(contr_ac%arc(idx)%link(1))
          end if
          if (contr_ac%arc(idx)%link(2).eq.0) then
            contr_abc%arc(narc_abc)%link(2) = 0
          else if (contr_ac%arc(idx)%link(2).lt.0) then
            ! fix for unique re-substitutions:
            idxsuper = -contr_ac%arc(idx)%link(2)
            if (imltlist(idxsuper,svmap,nvtx_b,1).eq.1) then
              contr_abc%arc(narc_abc)%link(2) =
     &             ivtx_b_reo(idxlist(idxsuper,svmap,nvtx_b,1))
            else
              contr_abc%arc(narc_abc)%link(2) = 0
            end if
          else
            contr_abc%arc(narc_abc)%link(2) =
     &           ivtx_ac_reo(contr_ac%arc(idx)%link(2))
          end if
          contr_abc%arc(narc_abc)%occ_cnt =
     &         contr_ac%arc(idx)%occ_cnt
        end if
      end do
      ! add [0] connections for any other intra-AC connection
      ! as these must not be generated
      narc_abc0 = narc_abc
      do ivtx = 1, nvtx_abc
        if (ivtx_old(ivtx).lt.0) cycle
        jloop: do jvtx = ivtx+1, nvtx_abc
          if (ivtx_old(jvtx).lt.0) cycle jloop
          do iarc = 1, narc_abc0
            if (contr_abc%arc(iarc)%link(1).eq.ivtx.and.
     &          contr_abc%arc(iarc)%link(2).eq.jvtx) cycle jloop
          end do
          narc_abc = narc_abc+1
          contr_abc%arc(narc_abc)%link(1) = ivtx
          contr_abc%arc(narc_abc)%link(2) = jvtx
          contr_abc%arc(narc_abc)%occ_cnt = 0
        end do jloop
      end do
 
      contr_abc%nvtx = nvtx_abc
      ! preliminary setting
      contr_abc%narc = narc_abc
      contr_abc%nfac = 0

      if (ntest.ge.1000) then
        write(luout,*) 'generated proto-contraction:'
        call prt_contr2(luout,contr_abc,op_info)
      end if

      ! make a "wrap" formula list for gen_contr output:
      call init_formula(wrap)

      ! set fix_vtx and occ_vtx arrays
      allocate(fix_vtx(nvtx_abc),occ_vtx(ngastp,2,nvtx_abc+njoined_abc))
      fix_vtx = .true. ! "fix" all vertices -> ieqvfac will be 1
      call occvtx4contr(0,occ_vtx,contr_abc,op_info)

      ! generate all possible contractions (hopefully only 1)
c      call gen_contr2(wrap,contr_abc,fix_vtx,occ_vtx,op_info)
      call gen_contr3(wrap,contr_abc,fix_vtx,
     &                occ_vtx,njoined_abc,op_info)

      ! none at all?
      if (wrap%command.eq.command_end_of_formula) then
        write(luout,*) 'proto-contraction:'
        call prt_contr2(luout,contr_abc,op_info)
        call quit(1,'join_contr2',
     &       'no possible connection found')
      end if

      ! more than one?
      if (wrap%next%command.ne.command_end_of_formula) then
        write(luout,*) 'proto-contraction:'
        call prt_contr2(luout,contr_abc,op_info)
        write(luout,*) 'generated terms:'
        call print_form_list(luout,wrap,op_info)
        call quit(1,'join_contr2',
     &       'version 2beta allows only unique recombinations!')
      end if

      ! copy the "wrapped" final contraction to output place
      call copy_contr(wrap%contr,contr_abc)

      if (nvtx_ac.gt.0) deallocate(ivtx_ac_reo)
      if (nvtx_b.gt.0)  deallocate(ivtx_b_reo)
      deallocate(fix_vtx,occ_vtx,ivtx_old)
      call dealloc_formula_list(wrap)

      deallocate(svmap)

      if (ntest.ge.100) then
        write(luout,*) 'generated contraction:'
        call prt_contr2(luout,contr_abc,op_info)
      end if
      
      return
      end
