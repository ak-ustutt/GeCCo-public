*----------------------------------------------------------------------*
      subroutine fact_cost2(possible,cost,iscale,
     &     contr,occ_vtx,irestr_vtx,info_vtx,iarc,
     &     op_info,str_info,ihpvgas,ngas,nsym)
*----------------------------------------------------------------------*
*     calculate the computational cost of a given contraction
*     for each contraction step, the number of flops is summed up
*     and the maximum size of the occurring intermediates is evaluated
*     returned on cost(1:3): cost(1)<-flops, 
*                            cost(2)<-memory(tot)
*                            cost(3)<-memory(block)
*     the maximum scaling of contraction and intermediates is returned
*     on iscale(ngastp,1:2): iscale(,1)<-contr, iscale(,2)<-mem
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'multd2h.h'
      
      integer, parameter ::
     &     ntest = 00
      
      integer, intent(in) ::
     &     ngas, ihpvgas(ngas), nsym, iarc
      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,contr%nvtx),
     &     irestr_vtx(2,ngas,2,2,contr%nvtx),
     &     info_vtx(2,contr%nvtx)
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      real(8), intent(out) ::
     &     cost(3)
      integer, intent(out) ::
     &     iscale(ngastp,2)
      logical, intent(out) ::
     &     possible

      integer ::
     &     nvtx, narc,
     &     np_int, nh_int, nx_int, np_cnt, nh_cnt, nx_cnt
      integer ::
     &     iocc_cnt(ngastp,2), iocc_ext(ngastp,2,2),
     &     irestr_res(2,ngas,2,2),
     &     iocc_op(ngastp,2,2), iocc_int(ngastp,2),
     &     irestr_int(2,ngas,2,2), igr_int(ngastp,2),
     &     irestr_op(2,ngas,2,2,2), 
     &     mst_op(2), mst_int,
     &     igamt_op(2), igamt_int,
     &     idar1(2), idar2(2)
      real(8) ::
     &     flops, xmemtot, xmemblk

      integer, external ::
     &     idxlist, int_expand, int_pack, maxxlvl_op

      if (ntest.gt.0) then
        call write_title(luout,wst_dbg_subr,'this is fact_cost')
      end if

      nvtx = contr%nvtx
      narc = contr%narc

      ! preliminary version (valid for standard CC only):
      ! restriction on result:
c dbg
c      print *,'call to set_restr_prel'
c dbg
      call set_restr_prel(irestr_res,contr,op_info,ihpvgas,ngas)
 
      call get_bc_info2(idar1,idar2,
     &     iocc_ext,iocc_cnt,
     &     iocc_op,iocc_int,
     &     irestr_op,irestr_int,
     &     mst_op,mst_int,
     &     igamt_op,igamt_int,
     &     contr,occ_vtx,irestr_vtx,info_vtx,iarc,
     &     irestr_res,ihpvgas,ngas)

      ! count particle, hole, (active) spaces involved:
      ! in intermediate
      nh_int = iocc_ext(ihole,1,1)+iocc_ext(ihole,2,1)+
     &         iocc_ext(ihole,1,2)+iocc_ext(ihole,2,2)
      np_int = iocc_ext(ipart,1,1)+iocc_ext(ipart,2,1)+
     &         iocc_ext(ipart,1,2)+iocc_ext(ipart,2,2)
      nx_int = iocc_ext(iextr,1,1)+iocc_ext(iextr,2,1)+
     &         iocc_ext(iextr,1,2)+iocc_ext(iextr,2,2)
      ! in contraction length
      nh_cnt = iocc_cnt(ihole,1)+iocc_cnt(ihole,2)
      np_cnt = iocc_cnt(ipart,1)+iocc_cnt(ipart,2)
      nx_cnt = iocc_cnt(iextr,1)+iocc_cnt(iextr,2)

      if (ntest.ge.100) then
        write(luout,*) 'op 1:'
        call wrt_occ(luout,iocc_op(1,1,1))
c dbg
c        call wrt_rstr(luout,irestr_op,ngas)
c dbg
        write(luout,*) 'op 2:'
        call wrt_occ(luout,iocc_op(1,1,2))
        write(luout,*) 'externals 1:'
        call wrt_occ(luout,iocc_ext(1,1,1))
        write(luout,*) 'externals 2:'
        call wrt_occ(luout,iocc_ext(1,1,2))
        write(luout,*) 'contraction:'
        call wrt_occ(luout,iocc_cnt)
        write(luout,*) 'intermediate/result:'
        call wrt_occ(luout,iocc_int)
      end if

      if (ntest.ge.50) then
        write(luout,'(x,a,i2,a,i2,a)')
     &         'Contraction scales as  H^{',nh_int+nh_cnt,
     &                                   '}P^{',np_int+np_cnt,'}'
        write(luout,'(x,a,i2,a,i2,a)')
     &             'Intermediate scales as H^{',nh_int,
     &                                   '}P^{',np_int,'}'
      end if

      ! set iscale
      ! maximum is taken over total scaling, so H^2 .gt. P^1 
      if ( nh_int+np_int+nx_int+nh_cnt+np_cnt+nx_cnt.gt.
     &           iscale(ihole,1)+iscale(ipart,1)+iscale(iextr,1) .or.
     &    (nh_int+np_int+nx_int+nh_cnt+np_cnt+nx_cnt.eq.
     &           iscale(ihole,1)+iscale(ipart,1)+iscale(iextr,1).and.
     &     np_int+np_cnt.gt.iscale(ipart,1) ) ) then
        iscale(ihole,1) = nh_int+nh_cnt
        iscale(ipart,1) = np_int+np_cnt
        iscale(iextr,1) = nx_int+nx_cnt
      end if
      if ( nh_int+np_int+nx_int.gt.
     &     iscale(ihole,2)+iscale(ipart,2)+iscale(iextr,2) .or.
     &    (nh_int+np_int+nx_int.eq.
     &     iscale(ihole,2)+iscale(ipart,2)+iscale(iextr,2).and.
     &     np_int.gt.iscale(ipart,2) ) ) then
        iscale(ihole,2) = nh_int
        iscale(ipart,2) = np_int
        iscale(iextr,2) = nx_int
      end if

      possible = .true.     
      ! check whether intermediate can be addressed by
      ! the available graphs (preliminary fix)
      call get_grph4occ(igr_int,iocc_int,irestr_int,
     &       str_info,ihpvgas,ngas,.false.)
      if (ntest.ge.100) then
        write(luout,*) 'igr_int:'
        call wrt_occ(luout,igr_int)
      end if
      ! if not: 
      if (igr_int(1,1).lt.0) then
        cost(1:3) = huge(cost(1))
        ! do not allow this factorization
        possible = .false.
      end if

      if (possible) then
        call dummy_contr(flops,xmemtot,xmemblk,
     &       iocc_op,iocc_ext,iocc_int,iocc_cnt,
     &       irestr_op,irestr_int,
     &       mst_op,mst_int,igamt_op,igamt_int,
     &       str_info,ngas,ihpvgas,nsym)

        cost(1) = cost(1)+flops
        cost(2) = max(cost(2),xmemtot)
        cost(3) = max(cost(3),xmemblk)

        if (ntest.ge.10)
     &       write(luout,'(x,a,3(g20.10,x))') '$$ cost: ',cost(1:3)
        if (ntest.ge.10)
     &       write(luout,'(x,a,2(2i4,2x))') '$$ scal: ',iscale(1:2,1:2)
      end if

      return
      end

