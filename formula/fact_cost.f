*----------------------------------------------------------------------*
      subroutine fact_cost(possible,cost,iscale,ifact,nfact,
     &     cost_prev,
     &     contr,op_info,str_info,ihpvgas,ngas,nsym)
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
     &     ngas, ihpvgas(ngas), nsym, nfact, ifact(4,nfact)
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      real(8), intent(in) ::
     &     cost_prev(3)
      real(8), intent(out) ::
     &     cost(3)
      integer, intent(out) ::
     &     iscale(ngastp,2)
      logical, intent(out) ::
     &     possible

      integer ::
     &     idx, jdx, iops, ivtx, nvtx, narc,
     &     idx_int, icnt, ninter, idxop, iblkop,
     &     np_int, nh_int, np_cnt, nh_cnt,
     &     maxex, ica, igas, igastp, idum1(2), idum2(2)
      integer ::
     &     iocc(ngastp,2,2), iocc_cnt(ngastp,2), iocc_ext(ngastp,2,2),
     &     interm(nfact), iocc_int(ngastp,2,nfact),
     &     irestr_int(2,ngas,2,2,nfact),
     &     irestr_res(2,ngas,2,2), igr_int(ngastp,2),
     &     ivtx_expand(contr%nvtx,2), nexpand(2),
     &     irestr(2,ngas,2,2,2), iscr(contr%nvtx),
     &     mstop(2), mstint(nfact),
     &     igamtop(2), igamtint(nfact)
      real(8) ::
     &     flops, xmemtot, xmemblk

      integer, external ::
     &     idxlist, int_expand, int_pack, maxxlvl_op

      if (ntest.gt.0) then
        write(luout,*) '==================='
        write(luout,*) ' this is fact_cost'
        write(luout,*) '==================='
      end if

      cost(1:3) = 0d0
      iscale(1:ngastp,1:2) = 0
      nvtx = contr%nvtx
      narc = contr%narc
      ninter = 0

      ! preliminary version (valid for standard CC only):
      ! restriction on result:
c dbg
c      print *,'call to set_restr_prel'
c dbg
      call set_restr_prel(irestr_res,contr,op_info,ihpvgas,ngas)
 
      ! let's assume the current factorization is allowed
      possible = .true.
     
      ! loop over operation sequence
      do idx = 1, nfact
        if (ntest.ge.50) then
          write(luout,*) 'binary operation #',idx
        end if

        ! set up info for binary contraction
        call get_bc_info(iocc,irestr,
     &       mstop,igamtop,idum1,idum2,
     &       iocc_ext,iocc_cnt,
     &       ifact(1,idx),contr,op_info,irestr_res,
     &       interm,ninter,iocc_int,irestr_int, mstint,igamtint,
     &       ihpvgas,ngas)

        ! count particle, hole, (active) spaces involved:
        ! in intermediate
        nh_int = iocc_ext(1,1,1)+iocc_ext(1,2,1)+
     &           iocc_ext(1,1,2)+iocc_ext(1,2,2)
        np_int = iocc_ext(2,1,1)+iocc_ext(2,2,1)+
     &           iocc_ext(2,1,2)+iocc_ext(2,2,2)
        ! in contraction length
        nh_cnt = iocc_cnt(1,1)+iocc_cnt(1,2)
        np_cnt = iocc_cnt(2,1)+iocc_cnt(2,2)

        if (ntest.ge.100) then
          write(luout,*) 'externals 1:'
          call wrt_occ(luout,iocc_ext(1,1,1))
          write(luout,*) 'externals 2:'
          call wrt_occ(luout,iocc_ext(1,1,2))
          write(luout,*) 'contraction:'
          call wrt_occ(luout,iocc_cnt)
          write(luout,*) 'intermediate/result:'
          call wrt_occ(luout,iocc_int(1,1,ninter))
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
        if (nh_int+np_int+nh_cnt+np_cnt.gt.iscale(1,1)+iscale(2,1)) then
          iscale(1,1) = nh_int+nh_cnt
          iscale(2,1) = np_int+np_cnt
        end if
        if (nh_int+np_int.gt.iscale(1,2)+iscale(2,2)) then
          iscale(1,2) = nh_int
          iscale(2,2) = np_int
        end if

        ! check whether intermediate can be addressed by
        ! the available graphs (preliminary fix)
        call get_grph4occ(igr_int,iocc_int(1,1,ninter),
     &       irestr_int(1,1,1,1,ninter),
     &       str_info,ihpvgas,ngas,.false.)
        ! if not: 
        if (igr_int(1,1).lt.0) then
          cost(1:3) = huge(cost(1))
          ! do not allow this factorization
          possible = .false.
          exit
        end if

        call dummy_contr(flops,xmemtot,xmemblk,
     &       iocc,iocc_ext,iocc_int(1,1,ninter),iocc_cnt,
     &       irestr,irestr_int(1,1,1,1,ninter),
     &       mstop,mstint,igamtop,igamtint,
     &       str_info,ngas,ihpvgas,nsym)

        cost(1) = cost(1)+flops
        cost(2) = max(cost(2),xmemtot)
        cost(3) = max(cost(3),xmemblk)

        ! if the cost is already higher than before: exit
        if (cost(1).gt.cost_prev(1)) exit

      end do

      if (ntest.ge.10)
     &     write(luout,'(x,a,3(g20.10,x))') '$$ cost: ',cost(1:3)
      if (ntest.ge.10)
     &     write(luout,'(x,a,2(2i4,2x))') '$$ scal: ',iscale(1:2,1:2)

      return
      end

