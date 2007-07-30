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
c      include 'def_operator.h'
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
     &     idx, idx_abc, iarc, ivtx, jvtx
c     &     nconnect_a, nconnect_c,
c     &     nconnect_b_dx, nconnect_b_ex,
c     &     ieqvfac
      type(formula_item) ::
     &     wrap
c      integer ::
c     &     iocc(ngastp,2), jocc(ngastp,2),
c     &     iocc_a_dx(ngastp,2), iocc_b_dx(ngastp,2),
c     &     iocc_b_ex(ngastp,2), iocc_c_ex(ngastp,2)
     
      integer, pointer ::
     &     ivtx_a(:), ivtx_b_dx(:), ivtx_b_ex(:), ivtx_c(:),
     &     ivtx_ac_reo(:), ivtx_b_reo(:), ivtx_reo(:),
     &     occ_vtx(:,:,:)
      logical, pointer ::
     &     fix_vtx(:)

      type(cntr_arc), pointer ::
     &     arc(:)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'This is join_contr')
        write(luout,*) 'joining: AC, B'
        call prt_contr2(luout,contr_ac,op_info)
        call prt_contr2(luout,contr_b,op_info)
      end if

      nvtx_ac = contr_ac%nvtx
      narc_ac = contr_ac%narc
      ! count all vertices up to first "0" entry
      nvtx_a = 0
      do idx = 1, nvtx_ac
        if (contr_ac%vertex(idx)%idx_op.ne.0) then
          nvtx_a = nvtx_a+1
c          vtxlist(idx) = idx
        else
          exit
        end if
      end do

      nvtx_c = nvtx_ac-nvtx_a-1

      if (contr_b%nvtx.le.0)
     &     call quit(1,'join_contr',
     &     'inserted contraction fragment must be non-empty')

      nvtx_b = contr_b%nvtx
      narc_b = contr_b%narc

      nvtx_abc = nvtx_a+nvtx_b+nvtx_c
      narc_abc = narc_ac+narc_b+nvtx_a*(nvtx_b-1)+nvtx_c*(nvtx_b-1)
     &           +nvtx_a*nvtx_c

      if (ntest.ge.1000) then
        write(luout,*) 'nvtx_a, nvtx_b, nvtx_c: ',nvtx_a,nvtx_b,nvtx_c
        write(luout,*) 'narc_ac, narc_b, narc_abc: ',
     &       narc_ac, narc_b, narc_abc
      end if

      call resize_contr(contr_abc,nvtx_abc,narc_abc,0)

      if (nvtx_ac.gt.0) allocate(ivtx_ac_reo(nvtx_ac))
      if (nvtx_b.gt.0)  allocate(ivtx_b_reo(nvtx_b))
      if (nvtx_ac.gt.0) ivtx_ac_reo(1:nvtx_ac) = 0
      if (nvtx_b.gt.0)  ivtx_b_reo(1:nvtx_b) = 0

      ! set prefactor
      contr_abc%fac = contr_ac%fac*contr_b%fac
      ! set result
      contr_abc%idx_res = idxop_abc
      contr_abc%iblk_res = iblk_abc

      ! set vertices and reordering arrays
      ! collect vertices from A
      do idx = 1, nvtx_a
        contr_abc%vertex(idx) = contr_ac%vertex(idx)
        ivtx_ac_reo(idx) = idx
      end do
      ! collect vertices from B
      idx_abc = nvtx_a
      do idx = 1, nvtx_b
        idx_abc = idx_abc+1
        contr_abc%vertex(idx_abc) = contr_b%vertex(idx)
        ivtx_b_reo(idx) = idx_abc        
      end do
      ! collect vertices from C
      do idx = nvtx_a+2,nvtx_ac
        idx_abc = idx_abc+1
        contr_abc%vertex(idx_abc) = contr_ac%vertex(idx)
        ivtx_ac_reo(idx) = idx_abc
      end do
c dbg
c      print *,'ivtx_b_reo: ',ivtx_b_reo
c      print *,'ivtx_ac_reo: ',ivtx_ac_reo
c dbg

      narc_abc = 0
      ! add all arcs from A and C, except the external ones
      do idx = 1, narc_ac
        if (contr_ac%arc(idx)%link(1).eq.0.or.
     &      contr_ac%arc(idx)%link(2).eq.0) cycle
        narc_abc = narc_abc + 1
        contr_abc%arc(narc_abc)%link(1) =
     &       ivtx_ac_reo(contr_ac%arc(idx)%link(1))
        contr_abc%arc(narc_abc)%link(2) =
     &       ivtx_ac_reo(contr_ac%arc(idx)%link(2))
c dbg
c        print *,' ori: ',contr_ac%arc(idx)%link(1),
c     &                   contr_ac%arc(idx)%link(2)
c        print *,' reo: ',ivtx_ac_reo(contr_ac%arc(idx)%link(1)),
c     &                   ivtx_ac_reo(contr_ac%arc(idx)%link(2))
c dbg
        contr_abc%arc(narc_abc)%occ_cnt =
     &       contr_ac%arc(idx)%occ_cnt
      end do
c dbg
c      contr_abc%narc = narc_abc
c        write(luout,*) 'generated proto-contraction (-2):'
c        call prt_contr2(luout,contr_abc,op_info)
c dbg
      ! add all arcs from B
c      idx_abc = narc_abc !narc_ac
      do idx = 1, narc_b
c        idx_abc = idx_abc+1
        narc_abc = narc_abc + 1
        contr_abc%arc(narc_abc)%link(1) =
     &       ivtx_b_reo(contr_b%arc(idx)%link(1))
        contr_abc%arc(narc_abc)%link(2) =
     &       ivtx_b_reo(contr_b%arc(idx)%link(2))
        contr_abc%arc(narc_abc)%occ_cnt =
     &       contr_b%arc(idx)%occ_cnt
      end do
c dbg
c      contr_abc%narc = narc_abc
c        write(luout,*) 'generated proto-contraction (-1):'
c        call prt_contr2(luout,contr_abc,op_info)
c dbg
      ! add the external arcs from A and C
      ! which make gen_contr consider only special connections
      do idx = 1, narc_ac
        if (contr_ac%arc(idx)%link(1).eq.0.or.
     &      contr_ac%arc(idx)%link(2).eq.0) then
          narc_abc = narc_abc + 1
          if (contr_ac%arc(idx)%link(1).eq.0) then
            contr_abc%arc(narc_abc)%link(1) = 0
          else
            contr_abc%arc(narc_abc)%link(1) = 
     &           ivtx_ac_reo(contr_ac%arc(idx)%link(1))
          end if
          if (contr_ac%arc(idx)%link(2).eq.0) then
            contr_abc%arc(narc_abc)%link(2) = 0
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
c dbg
c      contr_abc%narc = narc_abc
c        write(luout,*) 'generated proto-contraction (0):'
c        call prt_contr2(luout,contr_abc,op_info)
c dbg
      narc_abc0 = narc_abc
      do ivtx = 1, nvtx_a
        jloop: do jvtx = nvtx_a+nvtx_b+1, nvtx_abc
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
      allocate(fix_vtx(nvtx_abc),occ_vtx(ngastp,2,nvtx_abc+1))
      fix_vtx = .true. ! "fix" all vertices -> ieqvfac will be 1
      call occvtx4contr(occ_vtx,contr_abc,op_info)

      ! generate all possible contractions (hopefully only 1)
      call gen_contr2(wrap,contr_abc,fix_vtx,occ_vtx,op_info)

      ! none at all?
      if (wrap%command.eq.command_end_of_formula) then
        write(luout,*) 'proto-contraction:'
        call prt_contr2(luout,contr_abc,op_info)
        call quit(1,'join_contr',
     &       'no possible connection found')
      end if

      ! more than one?
      if (wrap%next%command.ne.command_end_of_formula) then
        write(luout,*) 'proto-contraction:'
        call prt_contr2(luout,contr_abc,op_info)
        write(luout,*) 'generated terms:'
        call print_form_list(luout,wrap,op_info)
        call quit(1,'join_contr',
     &       'version 2beta allows only unique recombinations!')
      end if

      ! copy the "wrapped" final contraction to output place
      call copy_contr(wrap%contr,contr_abc)

      deallocate(fix_vtx,occ_vtx)
      call dealloc_formula_list(wrap)

      if (ntest.ge.100) then
        write(luout,*) 'generated contraction:'
        call prt_contr2(luout,contr_abc,op_info)
      end if
      
      return
      end
