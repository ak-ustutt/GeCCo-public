*----------------------------------------------------------------------*
      subroutine join_contr2a(mode,fl_abc,nterms,contr_abc,contr_ac,
     &     contr_b,idxop_abc,iblk_abc,op_info)
*----------------------------------------------------------------------*
*     join contractions ac and b into one with result occupation abc
*
*     extended intermediate solution that can generate
*     all possible cases if new connection is not unique 
*
*     mode=0: single output contraction is returned on contr_abc
*     mode=1: output formula is returned on fl_abc, #of terms on nterms
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
     &     ntest = 000

      type(formula_item), intent(out), target ::
     &     fl_abc
      type(contraction), intent(out) ::
     &     contr_abc
      integer, intent(out) ::
     &     nterms

      type(contraction), intent(in) ::
     &     contr_ac, contr_b
      integer, intent(in) ::
     &     idxop_abc, iblk_abc, mode
      type(operator_info), intent(in) ::
     &     op_info
      type(formula_item), pointer ::
     &     fl_abc_pnt
 
      logical ::
     &     reo, unique, done
      integer ::
     &     nvtx_abc, nvtx_ac, nvtx_a, nvtx_b, nvtx_c,
     &     narc_abc, narc_abc0, narc_ac, narc_b, nxarc_ac,
     &     idx, jdx, ivtx_abc, iarc, ivtx, jvtx, jvtx_last,
     &     nproto_ac, nproto_b, idxsuper,
     &     nsuper, nsuper_non0, njoined, isuper, njoined_abc,
     &     occ_x(ngastp,2), occ_over(ngastp,2), icnt, ioff
      integer(8) ::
     &     base, overlap
      type(operator), pointer ::
     &     opres
     
      integer, pointer ::
     &     ivtx_ac_reo(:), ivtx_b_reo(:),
     &     occ_vtx(:,:,:), svmap(:), ivtx_old(:), svtx(:), ol_map(:)
      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:)
      logical, pointer ::
     &     fix_vtx(:), found(:)

      type(cntr_arc), pointer ::
     &     arc(:)

      integer, external ::
     &     ifndmax, imltlist, idxlist, int8_expand
      integer(8), external ::
     &     int8_pack

      if (mode.ne.0.and.mode.ne.1) call quit(1,'join_contr2a',
     &        'mode must be 0 or 1')

      base = pack_base

      if (ntest.ge.50) then
        call write_title(luout,wst_dbg_subr,'This is join_contr 2a')
        write(luout,*) 'joining: AC, B'
        call prt_contr2(luout,contr_ac,op_info)
        call prt_contr2(luout,contr_b,op_info)
      end if

      nvtx_ac = contr_ac%nvtx
      narc_ac = contr_ac%narc
      nxarc_ac = contr_ac%nxarc

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
     &     call quit(1,'join_contr2a',
     &     'inserted contraction fragment must be non-empty and '//
     &     'must not contain proto-vertices!')

      ! get the super-vertex map for B
      allocate(svmap(nvtx_b))
      opres => op_info%op_arr(contr_b%idx_res)%op
      njoined = opres%njoined
      if (njoined.eq.1) then
        svmap(1:nvtx_b) = 1
        unique = .true.
      else if (mode.eq.0) then
        if (njoined.ne.nvtx_b)
     &   call quit(1,'join_contr2a',
     &   'expected B to be single operator intermediate')
        ! CAUTION: not clear if this works in the general case
        svmap = 0
        ioff = (contr_b%iblk_res-1)*njoined
        do idx = 1, njoined
          if (iocc_nonzero(op_info%op_arr(
     &        contr_b%idx_res)%op%ihpvca_occ(1:ngastp,1:2,ioff+idx)))
     &        svmap(idx) = idx
        end do
        unique = .true.
c        ! need to call old svmap4contr because svmap4contr2 cannot
c        ! deal with single operator intermediates where open lines
c        ! are not explicitly given
c        allocate(occ_vtx(ngastp,2,nvtx_b+njoined))
c        call occvtx4contr(0,occ_vtx,contr_b,op_info)
c        call svmap4contr(svmap,contr_b,occ_vtx,njoined)
c        deallocate(occ_vtx)
      else
c        allocate(occ_vtx(ngastp,2,nvtx_b+njoined))
c        call occvtx4contr(0,occ_vtx,contr_b,op_info)
c dbg
c        print *,'opres: ',trim(opres%name)
c        print *,'njoined = ',njoined
c        print *,'call in join_contr2a'
c dbg
        call svmap4contr2(svmap,contr_b,unique)

        ! quick fix: middle zero vertex would not be accounted for
        if (njoined.eq.3.and.svmap(3).eq.3.and.svmap(2).eq.0)
     &          unique = .false.

        if (.not.unique) call pseudo_svmap(svmap,contr_b,njoined)
c        deallocate(occ_vtx)
      end if
      ! largest index = number of super vertices (at least 1)
c      nsuper = ifndmax(svmap,1,nvtx_b,1)
      nsuper = max(1,ifndmax(svmap,1,nvtx_b,1))

      ! check for zero occupations in B
      nsuper_non0 = nsuper
c      if (nsuper.ne.nproto_ac) then
cc dbg
c        print *,'unique= ',unique
c        print *,'svmap=',svmap
cc dbg
c        nsuper_non0 = 0
c        ioff = (contr_b%iblk_res-1)*njoined
c        do idx = 1, njoined
c          if (iocc_nonzero(op_info%op_arr(
c     &        contr_b%idx_res)%op%ihpvca_occ(1:ngastp,1:2,ioff+idx)))
c     &    then
c            nsuper_non0 = nsuper_non0+1
c          else if (.not.unique) then
c            do jdx = idx, nvtx_b
c              svmap(jdx) = svmap(jdx)-1
c            end do
c          end if
c        end do
c      end if

      if (nsuper.ne.nproto_ac.and.nsuper_non0.gt.nproto_ac) then
        write(luout,*) 'join_contr2a: joining: AC, B'
        call prt_contr2(luout,contr_ac,op_info)
        call prt_contr2(luout,contr_b,op_info)
        write(luout,*) 'nsuper, nproto_ac: ',nsuper, nproto_ac
        write(luout,*) 'nsuper_non0: ',nsuper_non0
        write(luout,*) 'svmap: ',svmap
        write(luout,*) 'mode:  ',mode
        call quit(1,'join_contr2a','incompatible contractions !')
      end if

      if (ntest.ge.50) write(luout,*) 'nsuper, nsuper_non0: ',
     &                                 nsuper, nsuper_non0
      if (ntest.ge.50) write(luout,*) 'svmap: ',
     &                                 svmap

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
     &     narc_ac+narc_b+nvtx_ac*(nvtx_ac-1)+nvtx_b*(nvtx_b-1)
     &     +2*nvtx_ac,(nvtx_abc+2)*(nvtx_abc+1))

      if (ntest.ge.1000) then
        write(luout,*) 'nvtx_ac, nvtx_b, nvtx_abc: ',
     &       nvtx_ac, nvtx_b, nvtx_abc
        write(luout,*) 'narc_ac, narc_b, narc_abc: ',
     &       narc_ac, narc_b, narc_abc
      end if

      call resize_contr(contr_abc,nvtx_abc+2,narc_abc,0,0)

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

      ! result vertex: first and last vertex
      njoined_abc = op_info%op_arr(idxop_abc)%op%njoined
      ioff = (iblk_abc-1)*njoined_abc
      contr_abc%vertex(1)%idx_op = idxop_abc
      contr_abc%vertex(1)%iblk_op = ioff+1 ! dummy
      contr_abc%vertex(1)%dagger = .false.
      contr_abc%svertex(1) = 1
      contr_abc%vertex(nvtx_abc+2)%idx_op = idxop_abc
      contr_abc%vertex(nvtx_abc+2)%iblk_op = ioff+2 ! dummy
      contr_abc%vertex(nvtx_abc+2)%dagger = .false.
      contr_abc%svertex(nvtx_abc+2) = 1

      ! set vertices and reordering arrays
      do ivtx = 1, nvtx_abc
        if (ivtx_old(ivtx).gt.0) then
          contr_abc%vertex(ivtx+1) =
     &             contr_ac%vertex(ivtx_old(ivtx))
          contr_abc%svertex(ivtx+1) =
     &             contr_ac%svertex(ivtx_old(ivtx))+1
          ivtx_ac_reo(ivtx_old(ivtx)) = ivtx
        else
          contr_abc%vertex(ivtx+1) =
     &             contr_b%vertex(-ivtx_old(ivtx))
          contr_abc%svertex(ivtx+1) =
     &         contr_ac%nsupvtx + contr_b%svertex(-ivtx_old(ivtx))+1
          ivtx_b_reo(-ivtx_old(ivtx)) = ivtx
        end if
      end do

      ! set up correct super-vertex info
      contr_abc%nvtx = nvtx_abc+2
      call update_svtx4contr(contr_abc)

      narc_abc = 0
      ! add all arcs from A and C, except the external ones
      do idx = 1, narc_ac
        if (contr_ac%arc(idx)%link(1).le.0.or.
     &      contr_ac%arc(idx)%link(2).le.0) cycle
        narc_abc = narc_abc + 1
        contr_abc%arc(narc_abc)%link(1) = 1 +
     &       ivtx_ac_reo(contr_ac%arc(idx)%link(1))
        contr_abc%arc(narc_abc)%link(2) = 1 +
     &       ivtx_ac_reo(contr_ac%arc(idx)%link(2))
        contr_abc%arc(narc_abc)%occ_cnt =
     &       contr_ac%arc(idx)%occ_cnt
      end do
      ! add all arcs from B
      do idx = 1, narc_b
        narc_abc = narc_abc + 1
        contr_abc%arc(narc_abc)%link(1) = 1 +
     &       ivtx_b_reo(contr_b%arc(idx)%link(1))
        contr_abc%arc(narc_abc)%link(2) = 1 +
     &       ivtx_b_reo(contr_b%arc(idx)%link(2))
        contr_abc%arc(narc_abc)%occ_cnt =
     &       contr_b%arc(idx)%occ_cnt
      end do

      ! external arcs from A and C
      ! and [0] contractions to avoid unwanted external arcs
      occ_x = 0
      do idx = 1, njoined_abc
        occ_x = occ_x + op_info%op_arr(idxop_abc)%op%ihpvca_occ(
     &                   1:ngastp,1:2,ioff+idx)
      end do
      if (iocc_nonzero(occ_x)) then
        allocate(found(nvtx_abc))
        found = .false.
        do idx = 1, nxarc_ac
          ivtx = ivtx_ac_reo(contr_ac%xarc(idx)%link(1))
          if (found(ivtx)) call quit(1,'join_contr2a',
     &                'more than one external arc for one vertex?')
          found(ivtx) = .true.
          ! connect excitation part to first vertex ...
          occ_x = iocc_xdn(1,contr_ac%xarc(idx)%occ_cnt)
          narc_abc = narc_abc + 1
          contr_abc%arc(narc_abc)%link(1) = 1
          contr_abc%arc(narc_abc)%link(2) = 1 + ivtx
          contr_abc%arc(narc_abc)%occ_cnt = iocc_dagger(occ_x)
          ! ... and deexcitation part to last vertex
          occ_x = iocc_xdn(2,contr_ac%xarc(idx)%occ_cnt)
          narc_abc = narc_abc + 1
          contr_abc%arc(narc_abc)%link(1) = 1 + ivtx
          contr_abc%arc(narc_abc)%link(2) = nvtx_abc + 2
          contr_abc%arc(narc_abc)%occ_cnt = occ_x
        end do
        do jvtx = 1, nvtx_ac
          ivtx = ivtx_ac_reo(jvtx)
          if (ivtx.ne.0.and..not.found(ivtx)) then
            narc_abc = narc_abc + 1
            contr_abc%arc(narc_abc)%link(1) = 1
            contr_abc%arc(narc_abc)%link(2) = 1 + ivtx
            contr_abc%arc(narc_abc)%occ_cnt = 0
            narc_abc = narc_abc + 1
            contr_abc%arc(narc_abc)%link(1) = 1 + ivtx
            contr_abc%arc(narc_abc)%link(2) = nvtx_abc + 2
            contr_abc%arc(narc_abc)%occ_cnt = 0
          end if
        end do
        deallocate(found)
      end if

      if (.not.unique) then
        allocate(svtx(nvtx_b),vtx(nvtx_b),topo(nvtx_b,nvtx_b),
     &           xlines(nvtx_b,njoined))
        call pack_contr(svtx,vtx,topo,xlines,contr_b,njoined)
        if (ntest.ge.100) then
          write(luout,*) 'no unique svmap! Using xlines instead:'
          call prt_contr_p(luout,svtx,vtx,topo,xlines,nvtx_b,njoined)
        end if
      end if
      ! add the external arcs from A and C
      ! which make gen_contr consider only special connections
      do idx = 1, narc_ac
        if (contr_ac%arc(idx)%link(1).le.0.or.
     &      contr_ac%arc(idx)%link(2).le.0) then
          narc_abc = narc_abc + 1
          done = .false.
          if (contr_ac%arc(idx)%link(1).eq.0) then
            contr_abc%arc(narc_abc)%link(1) = 0
          else if (contr_ac%arc(idx)%link(1).lt.0.and.unique) then
            ! fix for unique re-substitutions:
            idxsuper = -contr_ac%arc(idx)%link(1)
            if (imltlist(idxsuper,svmap,nvtx_b,1).eq.1) then
              contr_abc%arc(narc_abc)%link(1) = 1 +
     &             ivtx_b_reo(idxlist(idxsuper,svmap,nvtx_b,1))
            else
              contr_abc%arc(narc_abc)%link(1) = 0
            end if
          else if (contr_ac%arc(idx)%link(1).lt.0) then
            idxsuper = -contr_ac%arc(idx)%link(1)
            do ivtx = 1, nvtx_b
              if (xlines(ivtx,idxsuper).eq.0) cycle
              occ_x = 0
              icnt = int8_expand(xlines(ivtx,idxsuper),base,occ_x)
              occ_over = iocc_overlap(occ_x,.false.,
     &                     contr_ac%arc(idx)%occ_cnt,.false.)
              overlap = int8_pack(occ_over,ngastp*2,base)
              if (overlap.ne.0) then
                contr_abc%arc(narc_abc)%link(1) = 1 +
     &                     ivtx_b_reo(ivtx)
                contr_abc%arc(narc_abc)%link(2) = 1 +
     &               ivtx_ac_reo(contr_ac%arc(idx)%link(2))
                contr_abc%arc(narc_abc)%occ_cnt = occ_over
                xlines(ivtx,idxsuper) = xlines(ivtx,idxsuper) - overlap
                narc_abc = narc_abc + 1
                done = .true.
              end if
            end do
          else
            contr_abc%arc(narc_abc)%link(1) = 1 +
     &           ivtx_ac_reo(contr_ac%arc(idx)%link(1))
          end if
          if (contr_ac%arc(idx)%link(2).eq.0) then
            contr_abc%arc(narc_abc)%link(2) = 0
          else if (contr_ac%arc(idx)%link(2).lt.0.and.unique) then
            ! fix for unique re-substitutions:
            idxsuper = -contr_ac%arc(idx)%link(2)
            if (imltlist(idxsuper,svmap,nvtx_b,1).eq.1) then
              contr_abc%arc(narc_abc)%link(2) = 1 +
     &             ivtx_b_reo(idxlist(idxsuper,svmap,nvtx_b,1))
            else
              contr_abc%arc(narc_abc)%link(2) = 0
            end if
          else if (contr_ac%arc(idx)%link(2).lt.0) then
            idxsuper = -contr_ac%arc(idx)%link(2)
            do ivtx = 1, nvtx_b
              if (xlines(ivtx,idxsuper).eq.0) cycle
              occ_x = 0
              icnt = int8_expand(xlines(ivtx,idxsuper),base,occ_x)
              occ_over = iocc_overlap(occ_x,.false.,
     &                     contr_ac%arc(idx)%occ_cnt,.true.)
              overlap = int8_pack(occ_over,ngastp*2,base)
              if (overlap.ne.0) then
                occ_over = iocc_overlap(occ_x,.true.,
     &                     contr_ac%arc(idx)%occ_cnt,.false.)
                contr_abc%arc(narc_abc)%link(1) = 1 +
     &               ivtx_ac_reo(contr_ac%arc(idx)%link(1))
                contr_abc%arc(narc_abc)%link(2) = 1 +
     &                     ivtx_b_reo(ivtx)
                contr_abc%arc(narc_abc)%occ_cnt = occ_over
                xlines(ivtx,idxsuper) = xlines(ivtx,idxsuper) - overlap
                narc_abc = narc_abc + 1
                done = .true.
              end if
            end do
          else
            contr_abc%arc(narc_abc)%link(2) = 1 +
     &           ivtx_ac_reo(contr_ac%arc(idx)%link(2))
          end if
          if (unique.and..not.done) then
            contr_abc%arc(narc_abc)%occ_cnt =
     &         contr_ac%arc(idx)%occ_cnt
          else
            narc_abc = narc_abc - 1
          end if
        end if
      end do

      if (.not.unique) then
        if (.not.all(xlines.eq.0)) then
          ! it seems, something did not work ...
          write(luout,*) 'what I generated so far seems inconsistent:'
          call prt_contr2(luout,contr_abc,op_info)
          call quit(1,'join_contr2a','trap!')
        end if
        deallocate(svtx,vtx,topo,xlines)
      end if

      ! add [0] connections for any other intra-AC connection
      ! as these must not be generated
      narc_abc0 = narc_abc
      do ivtx = 1, nvtx_abc
        if (ivtx_old(ivtx).lt.0) cycle
        jloop: do jvtx = ivtx+1, nvtx_abc
          if (ivtx_old(jvtx).lt.0) cycle jloop
          do iarc = 1, narc_abc0
            if (contr_abc%arc(iarc)%link(1).eq.ivtx+1.and.
     &          contr_abc%arc(iarc)%link(2).eq.jvtx+1)
     &         cycle jloop
          end do
          narc_abc = narc_abc+1
          contr_abc%arc(narc_abc)%link(1) = ivtx + 1
          contr_abc%arc(narc_abc)%link(2) = jvtx + 1
          contr_abc%arc(narc_abc)%occ_cnt = 0
        end do jloop
      end do

      ! preliminary setting
      contr_abc%narc = narc_abc
      contr_abc%nxarc = 0
      contr_abc%nfac = 0

      if (ntest.ge.1000) then
        write(luout,*) 'generated proto-contraction:'
        call prt_contr2(luout,contr_abc,op_info)
      end if

      ! set fix_vtx and occ_vtx arrays
      allocate(fix_vtx(nvtx_abc+2),
     &         occ_vtx(ngastp,2,nvtx_abc+2),
     &         ol_map(nvtx_abc+2))
      ol_map = 0
      fix_vtx = .true. ! "fix" all vertices -> ieqvfac will be 1
      call occvtx4contr(1,occ_vtx,contr_abc,op_info)

      ! daggered excitation and deexitation parts of result operator
      occ_vtx(1:ngastp,1:2,1) = 0
      occ_vtx(1:ngastp,1:2,nvtx_abc+2) = 0
      do idx = 1, njoined_abc
        occ_x = iocc_xdn(1,op_info%op_arr(idxop_abc)%op%ihpvca_occ(
     &                   1:ngastp,1:2,ioff+idx))
        if (iocc_nonzero(occ_x)) then
          if (ol_map(1).ne.0) then
            call quit(1,'join_contr2a','not yet adapted for this')
          else
            occ_vtx(1:ngastp,1:2,1) = iocc_dagger(occ_x)
            ol_map(1) = idx
          end if
        end if
        occ_x = iocc_xdn(2,op_info%op_arr(idxop_abc)%op%ihpvca_occ(
     &                   1:ngastp,1:2,ioff+idx))
        if (iocc_nonzero(occ_x)) then
          if (ol_map(nvtx_abc+2).ne.0) then
            call quit(1,'join_contr2a','not yet adapted for this')
          else
            occ_vtx(1:ngastp,1:2,nvtx_abc+2) = iocc_dagger(occ_x)
            ol_map(nvtx_abc+2) = -idx
          end if
        end if
      end do
      if (ol_map(1).eq.0) ol_map(1) = 1
      if (ol_map(nvtx_abc+2).eq.0) ol_map(nvtx_abc+2) = -1

      ! generate all possible contractions
      call gen_contr4(.false.,fl_abc,contr_abc,fix_vtx,
     &                occ_vtx,ol_map,op_info)

      ! sum up duplicate terms
      if (mode.eq.1) call sum_terms(fl_abc,op_info)

      ! count terms
      nterms = 0
      fl_abc_pnt => fl_abc
      do
        if (fl_abc_pnt%command.eq.command_end_of_formula) exit
        nterms = nterms+1
        fl_abc_pnt => fl_abc_pnt%next
      end do

      ! none at all? or more than one for mode=0?
      if (nterms.eq.0.or.mode.eq.0.and.nterms.gt.1) then
        write(luout,*) 'proto-contraction:'
        call prt_contr2(luout,contr_abc,op_info)
        if (nterms.eq.0) call quit(1,'join_contr2a',
     &       'no possible connection found')
        write(luout,*) 'generated terms:'
        call print_form_list(luout,fl_abc,op_info)
        call quit(1,'join_contr2a',
     &       'mode=0 allows only unique recombinations!')
      end if

      ! copy the "wrapped" final contraction to output place
      if (mode.eq.0) call copy_contr(fl_abc%contr,contr_abc)

      if (nvtx_ac.gt.0) deallocate(ivtx_ac_reo)
      if (nvtx_b.gt.0)  deallocate(ivtx_b_reo)
      deallocate(fix_vtx,occ_vtx,ivtx_old,ol_map)

      deallocate(svmap)

      if (ntest.ge.50) then
        write(luout,*) 'generated contraction(s):'
        call print_form_list(luout,fl_abc,op_info)
      end if
      
      return
      end
