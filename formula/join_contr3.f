*----------------------------------------------------------------------*
      subroutine join_contr3(contr_abc,contr_ac,contr_b,
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
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 0000

      type(contraction), intent(out) ::
     &     contr_abc

      type(contraction), intent(in) ::
     &     contr_ac, contr_b
      integer, intent(in) ::
     &     idxop_abc, iblk_abc
      type(operator_info), intent(in) ::
     &     op_info
 
      logical ::
     &     reo, ok
      integer ::
     &     nvtx_abc, nvtx_ac, nvtx_a, nvtx_b, nvtx_c,
     &     narc_abc, narc_abc0, narc_ac, narc_b, nxarc_b,
     &     idx, ivtx_abc, iarc, ivtx, jvtx, jvtx_last,
     &     nproto_ac, nproto_b, ivtx1, ivtx2, jvtx1, jvtx2,
     &     nsuper, njoined, isuper, njoined_abc, ncnt_b, icnt_b,
     &     icnt_ac, ncnt_ac
      type(formula_item) ::
     &     wrap
      type(operator), pointer ::
     &     opres
     
      integer, pointer ::
     &     ivtx_ac_reo(:), ivtx_b_reo(:),
     &     occ_vtx(:,:,:), svmap(:), ivtx_old(:),
     &     connect_map(:,:),
     &     occ_ex(:,:,:), occ_dx(:,:,:), occ_ac(:,:,:), idx_proto(:),
     &     list_b(:), list_ac(:)
      logical, pointer ::
     &     fix_vtx(:)

      type(cntr_arc), pointer ::
     &     arc(:)

      integer, external ::
     &     ifndmax, imltlist, idxlist
      call quit(1,'join_contr3','call to obsolete routine')

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'This is join_contr')
        write(lulog,*) 'joining: AC, B'
        call prt_contr2(lulog,contr_ac,op_info)
        call prt_contr2(lulog,contr_b,op_info)
      end if

      nvtx_ac = contr_ac%nvtx
      narc_ac = contr_ac%narc

      allocate(idx_proto(nvtx_ac))
      ! count proto-vertices in AC:
      nproto_ac = 0
      do ivtx = 1, nvtx_ac
        if (contr_ac%vertex(ivtx)%idx_op.eq.0) then
          nproto_ac = nproto_ac+1
          idx_proto(nproto_ac) = ivtx
        end if
      end do

      nvtx_b = contr_b%nvtx
      narc_b = contr_b%narc
      nxarc_b = contr_b%nxarc

      ! count proto-vertices in B (there should be none):
      nproto_b = 0
      do ivtx = 1, nvtx_b
        if (contr_b%vertex(ivtx)%idx_op.eq.0)
     &       nproto_b = nproto_b+1
      end do

      if (nvtx_b.le.0.or.nproto_b.ne.0)
     &     call quit(1,'join_contr',
     &     'inserted contraction fragment must be non-empty and '//
     &     'must not contain proto-vertices!')

      ! get the super-vertex map for B
      allocate(svmap(nvtx_b))
      opres => op_info%op_arr(contr_b%idx_res)%op
      njoined = opres%njoined
      call svmap4contr2(svmap,contr_b,ok)
      if (.not.ok) call quit(1,'join_contr3',
     &        'not prepared for non-unique svmap!')

      ! largest index = number of super vertices
      nsuper = ifndmax(svmap,1,nvtx_b,1)

      if (nsuper.ne.nproto_ac)
     &     call quit(1,'join_contr','incompatible contractions!')
      
      nvtx_abc = nvtx_ac-nproto_ac+nvtx_b

      ! generate a map: which vertex goes where
      allocate(ivtx_old(nvtx_abc))

      call joinmap4contr(ivtx_old,contr_ac,nvtx_abc,
     &                   0,-1,
     &                   svmap,nvtx_b,njoined)

      if (ntest.ge.1000) then
        write(lulog,'(3x,a,10i5)') 'ivtx_old: ',ivtx_old(1:nvtx_abc)
        write(lulog,'(3x,a,10i5)') 'idx_proto:',idx_proto(1:nproto_ac)
      end if

      ! make some assumptions about the number of arcs in the 
      ! proto contraction
      narc_abc = min(
     &     narc_ac+narc_b+nvtx_ac*(nvtx_ac-1)+nvtx_b*(nvtx_b-1),
     &     nvtx_abc*(nvtx_abc-1))

      if (ntest.ge.1000) then
        write(lulog,*) 'nvtx_ac, nvtx_b, nvtx_abc: ',
     &       nvtx_ac, nvtx_b, nvtx_abc
        write(lulog,*) 'narc_ac, narc_b, narc_abc: ',
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

      if (ntest.ge.1000) then
        write(lulog,*) 'ivtx_ac_reo:',ivtx_ac_reo(1:nvtx_ac)
        write(lulog,*) 'ivtx_b_reo: ',ivtx_b_reo(1:nvtx_b)
      end if

      ! set up correct super-vertex info
      contr_abc%nvtx = nvtx_abc
      call update_svtx4contr(contr_abc)

      narc_abc = 0
      ! add all arcs from AC which are not connected to proto vertices
      do idx = 1, narc_ac
        ivtx1 = ivtx_ac_reo(contr_ac%arc(idx)%link(1))
        ivtx2 = ivtx_ac_reo(contr_ac%arc(idx)%link(2))
        if (ivtx1.eq.0.or.ivtx2.eq.0) cycle
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

      ! set remaining vertices
      allocate(connect_map(nvtx_b,nvtx_ac))
      connect_map = 0
      ! loop over all AC vertices (excluding proto vertices)
      do ivtx = 1, nvtx_ac
        ivtx1 = ivtx_ac_reo(ivtx)
        if (ivtx1.le.0) cycle
        ! look for connections to proto vertices
        do iarc = 1, narc_ac
          ivtx1 = contr_ac%arc(iarc)%link(1)
          ivtx2 = contr_ac%arc(iarc)%link(2)
          if (ivtx1.ne.ivtx.and.ivtx2.ne.ivtx) cycle
          jvtx1 = ivtx_ac_reo(ivtx1)
          jvtx2 = ivtx_ac_reo(ivtx2)
          if (jvtx1.ne.0.and.jvtx2.ne.0) cycle
          if (jvtx1.eq.0) then
            isuper = idxlist(ivtx1,idx_proto,nsuper,1)
            ivtx1 = ivtx2
            jvtx1 = jvtx2
c            b_before_ac = .true.
          else
            isuper = idxlist(ivtx2,idx_proto,nsuper,1)
c            b_before_ac = .false.
          end if
          ! look for appropriate vertex on B
          do ivtx2 = 1, nvtx_b
            if (svmap(ivtx2).eq.isuper) then
              connect_map(ivtx2,ivtx1) = isuper
            end if
          end do
            

        end do

      end do

      if (ntest.ge.1000) then
        write(lulog,*) 'connect_map:'
        call iwrtma(connect_map,nvtx_b,nvtx_ac,nvtx_b,nvtx_ac)
      end if

      allocate(occ_ex(ngastp,2,nvtx_b),
     &         occ_dx(ngastp,2,nvtx_b),list_b(nvtx_b),
     &         occ_ac(ngastp,2,nvtx_ac),list_ac(nvtx_ac))

      ! loop over super vertices
      do isuper = 1, nsuper
        ! how many vertices in B do contribute?
        ncnt_b = imltlist(isuper,svmap,nvtx_b,1)        
        ! get external lines
        occ_ex = 0
        occ_dx = 0
        list_b = 0
        icnt_b = 0
        do idx = 1, nxarc_b
          if (contr_b%xarc(idx)%link(2).ne.isuper) cycle
          icnt_b = icnt_b+1
          list_b(icnt_b) = idx
          occ_ex(1:ngastp,1:2,icnt_b) =
     &         occ_ex(1:ngastp,1:2,icnt_b) +
     &         iocc_xdn(1,contr_b%xarc(idx)%occ_cnt)
          occ_dx(1:ngastp,1:2,icnt_b) =
     &         occ_dx(1:ngastp,1:2,icnt_b) +
     &         iocc_xdn(2,contr_b%xarc(idx)%occ_cnt)
        end do
        if (ntest.ge.1000) then
          write(lulog,*) 'list_b: ',list_b(1:ncnt_b)
          write(lulog,*) 'EX:'
          call wrt_occ_n(lulog,occ_ex,ncnt_b)
          write(lulog,*) 'DX:'
          call wrt_occ_n(lulog,occ_dx,ncnt_b)
        end if
        ! loop over vertices inserted before present B nodes
        icnt_ac = 0
        do ivtx = 1, idx_proto(isuper)
          ! look for matching arcs
          ivtx1 = ivtx
          ivtx2 = idx_proto(isuper)
c dbg
          write(6,*) 'looking for ',ivtx1,ivtx2,' arcs'
c dbg
          do iarc = 1, narc_ac
            if (contr_ac%arc(iarc)%link(1).ne.ivtx1 .or.
     &          contr_ac%arc(iarc)%link(2).ne.ivtx2) cycle
c dbg
            write(6,*) 'found arc # ',iarc
c            wrt_occ_n(lulog,contr_ac%arc(iarc)%occ_cnt,1)
c dbg            
            icnt_ac = icnt_ac+1
            list_ac(icnt_ac) = ivtx1
            occ_ac(1:ngastp,1:2,icnt_ac) = contr_ac%arc(iarc)%occ_cnt
          end do
        end do
        ncnt_ac = icnt_ac
        if (ntest.ge.1000) then
          write(lulog,*) 'BEFORE'
          write(lulog,*) 'list_ac: ',list_ac(1:ncnt_ac)
          write(lulog,*) 'OCCs:'
          call wrt_occ_n(lulog,occ_ac,ncnt_ac)
        end if
        write(6,*) 'do something'

        ! loop over vertices inserted after present B nodes
        icnt_ac = 0
        do ivtx = idx_proto(isuper), nvtx_ac
          ! look for matching arcs
          ivtx1 = idx_proto(isuper)
          ivtx2 = ivtx
c dbg
          write(6,*) 'looking for ',ivtx1,ivtx2,' arcs'
c dbg
          do iarc = 1, narc_ac
            if (contr_ac%arc(iarc)%link(1).ne.ivtx1 .or.
     &          contr_ac%arc(iarc)%link(2).ne.ivtx2) cycle
c dbg
            write(6,*) 'found arc # ',iarc
c            wrt_occ_n(lulog,contr_ac%arc(iarc)%occ_cnt,1)
c dbg            
            icnt_ac = icnt_ac+1
            list_ac(icnt_ac) = ivtx2
            occ_ac(1:ngastp,1:2,icnt_ac) = contr_ac%arc(iarc)%occ_cnt
          end do
        end do
        ncnt_ac = icnt_ac
        if (ntest.ge.1000) then
          write(lulog,*) 'AFTER'
          write(lulog,*) 'list_ac: ',list_ac(1:ncnt_ac)
          write(lulog,*) 'OCCs:'
          call wrt_occ_n(lulog,occ_ac,ncnt_ac)
        end if
        write(6,*) 'do something'

c        ! how many vertices in AC do contribute?
c        ncnt_ac = imltlist(isuper,
      end do
c      ! check connectivity type
c      do ivtx = 1, nvtx_ac
c        nconn = sum(connect_map(1:nvtx_b,ivtx))
c        if (nconn.eq.0) cycle
c        do jvtx = 1, nvtx_b
c          if (connect_map(jvtx,ivtx).eq.0) cycle
c          nconn2 = sum(connect_map(jvtx,1:nvtx_ac))
c          ! case one: single B vtx connected to one or more AC vtx's
c          if (nconn2.eq.1) then
c            ! treat this case only once
c            if (sum(connect_map(jvtx,1:ivtx-1)).gt.0) cycle
c            problem = .false. ! null problemo
c            isuper = svmap(jvtx)
c            
c
c          ! case two: more than one B vtx's connected to single AC vtx
c          else if (nconn.eq.1) then
c
c            
c          ! case three: no unique recombination possible
c          else
c            problem = .true.
c          end if
c          if (problem) then
c            write(lulog,*) 'AC'
c            call prt_contr2(lulog,contr_ac,op_info)
c            write(lulog,*) 'B'
c            call prt_contr2(lulog,contr_b,op_info)
c            call quit(1,'join_contr3','no unique recombination')
c          end if
c        end do
c      end do

      stop 'testing'

      if (nvtx_ac.gt.0) deallocate(ivtx_ac_reo)
      if (nvtx_b.gt.0)  deallocate(ivtx_b_reo)
      deallocate(svmap)

      if (ntest.ge.100) then
        write(lulog,*) 'generated contraction:'
        call prt_contr2(lulog,contr_abc,op_info)
      end if
      
      return
      end
