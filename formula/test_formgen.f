      subroutine test_formgen(op_info)

      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'ifc_operators.h'

      type(operator_info) ::
     &     op_info

      type(contraction) ::
     &     proto_contr
      type(formula_item), pointer ::
     &     form_head, current

      logical ::
     &     init, dag, fix_vtx(20)
      integer, allocatable ::
     &     occ_vtx(:,:,:)
      integer ::
     &     nvtx, narc, ivtx, idx, iblk, npart, itest

      integer ::
     &     occ_conn(ngastp,2,10), occ_dist(ngastp,2),
     &     occ_min(ngastp,2), occ_max(ngastp,2)

      logical, external ::
     &     next_part_connection

c      write(luout,*) 'testing next_part_connection:'
c      occ_min = 0
c      occ_max = 0
c      occ_max(1,2) = 2
c      occ_max(2,1) = 1
c      occ_max(4,1) = 1
c      occ_dist = occ_max
c      write(luout,*) 'occ, min, max: '
c      call wrt_occ(luout,occ_dist)
c      call wrt_occ(luout,occ_min)
c      call wrt_occ(luout,occ_max)
c
c      idx = 0
c      npart = 3
c      init = .true.
c      do while(next_part_connection(occ_conn,
c     &           init,occ_dist,npart,
c     &           occ_min,occ_max))
c        init = .false.
c        idx = idx+1
c        write(luout,*) '# ',idx
c        write(luout,'("|",4i4,"|",4i4,"|",4i4,"|",4i4,"|")')
c     &       occ_conn(1:ngastp,1,1:npart)
c        write(luout,'("|",4i4,"|",4i4,"|",4i4,"|",4i4,"|")')
c     &       occ_conn(1:ngastp,2,1:npart)
c        if (idx==100) exit
c      end do
c      stop 'test test'

      allocate(form_head)
      form_head%prev => null()
      form_head%next => null()
      form_head%contr => null()
      form_head%interm => null()

      call init_contr(proto_contr)
      itest = 2
      
      select case(itest)
      case(1)
      nvtx=3
      narc=1
      call resize_contr(proto_contr,nvtx,narc,0)
      ! ------------------------------
      ! set up a prototype-contraction
      ! ------------------------------
      ! set target occupation, op index, block number
      proto_contr%idx_res = 0
      proto_contr%iblk_res = 0
      proto_contr%fac = 1d0
      proto_contr%nvtx = nvtx
      proto_contr%narc = narc
      ! set for each vertex: 
      !  operator index, block number 
      proto_contr%vertex(1)%idx_op = 3
      proto_contr%vertex(1)%iblk_op = 2
      proto_contr%vertex(2)%idx_op = 1
      proto_contr%vertex(2)%iblk_op = 6
      proto_contr%vertex(3)%idx_op = 2
      proto_contr%vertex(3)%iblk_op = 2
      ! set contraction rules:
      !  for each pair for which an arc is set: at least one connection
      !  if arc%occ_cnt is set: use exactly that connection
      !  put (at least) first element of arc%occ_cnt to -1 to signal
      !  that contraction is arbitrary
      !  a zero contraction would mean: no contraction at all
      proto_contr%arc(1)%link(1) = 2
      proto_contr%arc(1)%link(2) = 3
      proto_contr%arc(1)%occ_cnt = -1
      fix_vtx(1:6) = .false.
      case(2)
      nvtx=5
      narc=3
      call resize_contr(proto_contr,nvtx,narc,0)
      proto_contr%idx_res = 0
      proto_contr%iblk_res = 0
      proto_contr%fac = 1d0
      proto_contr%nvtx = nvtx
      proto_contr%narc = narc
      proto_contr%vertex(1)%idx_op = 3
      proto_contr%vertex(1)%iblk_op = 2
      proto_contr%vertex(2)%idx_op = 1
      proto_contr%vertex(2)%iblk_op = 8
      proto_contr%vertex(3)%idx_op = 2
      proto_contr%vertex(3)%iblk_op = 1
      proto_contr%vertex(4)%idx_op = 2
      proto_contr%vertex(4)%iblk_op = 1
      proto_contr%vertex(5)%idx_op = 2
      proto_contr%vertex(5)%iblk_op = 2
      proto_contr%arc(1)%link(1) = 2
      proto_contr%arc(1)%link(2) = 3
      proto_contr%arc(1)%occ_cnt = -1
      proto_contr%arc(2)%link(1) = 2
      proto_contr%arc(2)%link(2) = 4
      proto_contr%arc(2)%occ_cnt = -1
      proto_contr%arc(3)%link(1) = 2
      proto_contr%arc(3)%link(2) = 5
      proto_contr%arc(3)%occ_cnt = -1
      case(3)
      nvtx=5
      narc=2
      fix_vtx(1:6) = .false.
      call resize_contr(proto_contr,nvtx,narc,0)
      proto_contr%idx_res = 0
      proto_contr%iblk_res = 0
      proto_contr%fac = 1d0
      proto_contr%nvtx = nvtx
      proto_contr%narc = narc
      proto_contr%vertex(1)%idx_op = 3
      proto_contr%vertex(1)%iblk_op = 2
      proto_contr%vertex(2)%idx_op = 3
      proto_contr%vertex(2)%iblk_op = 2
      proto_contr%vertex(3)%idx_op = 1
      proto_contr%vertex(3)%iblk_op = 6
      proto_contr%vertex(4)%idx_op = 2
      proto_contr%vertex(4)%iblk_op = 2
      proto_contr%vertex(5)%idx_op = 2
      proto_contr%vertex(5)%iblk_op = 2
      proto_contr%arc(1)%link(1) = 3
      proto_contr%arc(1)%link(2) = 4
      proto_contr%arc(1)%occ_cnt = -1
      proto_contr%arc(2)%link(1) = 3
      proto_contr%arc(2)%link(2) = 5
      proto_contr%arc(2)%occ_cnt = -1
      fix_vtx(1:6) = .false.
      case(4)
      nvtx=6
      narc=4
      call resize_contr(proto_contr,nvtx,narc,0)
      proto_contr%idx_res = 0
      proto_contr%iblk_res = 0
      proto_contr%fac = 1d0
      proto_contr%nvtx = nvtx
      proto_contr%narc = narc
      proto_contr%vertex(1)%idx_op = 3
      proto_contr%vertex(1)%iblk_op = 2
      proto_contr%vertex(2)%idx_op = 1
      proto_contr%vertex(2)%iblk_op = 8
      proto_contr%vertex(3)%idx_op = 1
      proto_contr%vertex(3)%iblk_op = 6
      proto_contr%vertex(4)%idx_op = 2
      proto_contr%vertex(4)%iblk_op = 2
      proto_contr%vertex(5)%idx_op = 1
      proto_contr%vertex(5)%iblk_op = 6
      proto_contr%vertex(6)%idx_op = 2
      proto_contr%vertex(6)%iblk_op = 2
      proto_contr%arc(1)%link(1) = 2
      proto_contr%arc(1)%link(2) = 4
      proto_contr%arc(1)%occ_cnt = -1
      proto_contr%arc(2)%link(1) = 2
      proto_contr%arc(2)%link(2) = 6
      proto_contr%arc(2)%occ_cnt = -1
      proto_contr%arc(3)%link(1) = 3
      proto_contr%arc(3)%link(2) = 4
      proto_contr%arc(3)%occ_cnt(1:ngastp,1:2) = 0
      proto_contr%arc(3)%occ_cnt(1,1) = 2
      proto_contr%arc(4)%link(1) = 5
      proto_contr%arc(4)%link(2) = 6
      proto_contr%arc(4)%occ_cnt(1:ngastp,1:2) = 0
      proto_contr%arc(4)%occ_cnt(1,1) = 2
      fix_vtx(1:6) = .false.
      fix_vtx(3) = .true.
      fix_vtx(5) = .true.
      end select

      ! for convenience in gen_contr(): set occupations for target and
      ! each vertex
      allocate(occ_vtx(ngastp,2,nvtx+1))
      if (proto_contr%idx_res.eq.0) then
        occ_vtx(1:ngastp,1:2,1) = 0
      else
        idx = proto_contr%idx_res
        iblk = proto_contr%iblk_res
        dag = op_info%op_arr(idx)%op%dagger
        if (.not.dag) then
          occ_vtx(1:ngastp,1:2,ivtx+1) =
     &       op_info%op_arr(idx)%op%ihpvca_occ(1:ngastp,1:2,iblk)
        else
          occ_vtx(1:ngastp,1:2,ivtx+1) = iocc_dagger(
     &       op_info%op_arr(idx)%op%ihpvca_occ(1:ngastp,1:2,iblk))
        end if
      end if
      do ivtx = 1, nvtx
        idx = proto_contr%vertex(ivtx)%idx_op
        iblk = proto_contr%vertex(ivtx)%iblk_op
        dag = op_info%op_arr(idx)%op%dagger
        if (.not.dag) then
          occ_vtx(1:ngastp,1:2,ivtx+1) =
     &       op_info%op_arr(idx)%op%ihpvca_occ(1:ngastp,1:2,iblk)
        else
          occ_vtx(1:ngastp,1:2,ivtx+1) = iocc_dagger(
     &       op_info%op_arr(idx)%op%ihpvca_occ(1:ngastp,1:2,iblk))
        end if
      end do

      write(luout,*) 'occ_vtx:'
      call wrt_occ_n(luout,occ_vtx,nvtx+1)

      call gen_contr(form_head,proto_contr,fix_vtx,occ_vtx,op_info)

      write(luout,*) 'generated list:'
      call print_form_list(luout,form_head,op_info)

      deallocate(occ_vtx)

      call dealloc_contr(proto_contr)

      call dealloc_formula_list(form_head)
      deallocate(form_head)

      call quit(1,'test_formgen','test exit')

      return
      end
