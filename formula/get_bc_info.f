*----------------------------------------------------------------------*
      subroutine get_bc_info(iocc,irestr,
     &     mstop, igamtop, idxop, iblkop,
     &     iocc_ext, iocc_cnt, 
     &     ifact,contr,op_info,irestr_res,
     &     interm, ninter, iocc_int, irestr_int, mstint, igamtint,
     &     ihpvgas,ngas)
*----------------------------------------------------------------------*
*     set up information for binary contraction
*     given: the present factorization info and the contraction struct
*     return: occ's, restr's, ms, IRREP, block numbers etc. of
*             involved operators
*             occ etc. of intermediate is on record ninter of iocc_int etc. 
*     and:  keep track of information about (temporary) intermediates
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ifact(4)
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     ngas,ihpvgas(ngas),irestr_res(2,ngas,2,2)
      integer, intent(out) ::
     &     iocc(ngastp,2,2), irestr(2,ngas,2,2,2),
     &     iocc_ext(ngastp,2,2), iocc_cnt(ngastp,2),
     &     mstop(2), igamtop(2), idxop(2), iblkop(2)
      integer, intent(inout) ::
     &     ninter, interm(*), iocc_int(ngastp,2,*),
     &     irestr_int(2,ngas,2,2,*), mstint(*), igamtint(*)

      integer ::
     &     iops, ivtx, nvtx, narc, idx_int, icnt, jdx
      integer ::
     &     ivtx_expand(contr%nvtx,2), nexpand(2), iscr(contr%nvtx)
      type(operator), pointer ::
     &     cur_op

      integer, external ::
     &     idxlist, int_expand, int_pack

      nvtx = contr%nvtx
      narc = contr%narc

      ! loop over involved operators (=vertices)
      do iops = 1, 2
        ivtx = ifact(iops) 
        ! primitive vertex or result of preceeding contraction?
        !  if ivtx .gt. contr%nvtx, look up on interm-table
        if (ivtx.gt.nvtx) then
          ! get index of intermediate
          idx_int = idxlist(ivtx,interm,ninter,1)
          if (idx_int.le.0)
     &         call quit(1,'get_bc_info','no intermediate found')
          ! get occ, restr:
          iocc(1:ngastp,1:2,iops) = iocc_int(1:ngastp,1:2,idx_int)
          irestr(1:2,1:ngas,1:2,1:2,iops) =
     &           irestr_int(1:2,1:ngas,1:2,1:2,idx_int)

          mstop(iops) = mstint(idx_int)
          igamtop(iops) = igamtint(idx_int)

          ! expand name (original vertices)
          nexpand(iops) = int_expand(ivtx,nvtx+1,ivtx_expand(1,iops))

          ! not a primary operator
          idxop(iops) = -idx_int
          iblkop(iops) = 1

        else
          ! get info on primitive vertex
          ivtx_expand(1,iops) = ivtx
          nexpand(iops) = 1
          ! need: (index), block of operator
          !     get occupation + restriction
          idxop(iops) = contr%vertex(ivtx)%idx_op
          iblkop(iops) = contr%vertex(ivtx)%iblk_op
          cur_op => op_info%op_arr(idxop(iops))%op
          iocc(1:ngastp,1:2,iops) =
     &          cur_op%ihpvca_occ(1:ngastp,1:2,iblkop(iops))
          irestr(1:2,1:ngas,1:2,1:2,iops) =
     &         cur_op%igasca_restr(1:2,1:ngas,1:2,1:2,iblkop(iops))
          mstop(iops) = cur_op%mst
          igamtop(iops) = cur_op%gamt
        end if
        if (ntest.ge.100) then
          write(luout,*) 'info on operand #',iops
          write(luout,*) ' involved vertices: ',
     &         ivtx_expand(1:nexpand(iops),iops)
          write(luout,*) 'occupation:'
          call wrt_occ(luout,iocc(1,1,iops))
          write(luout,*) 'restrictions:'
          call wrt_rstr(luout,irestr(1,1,1,1,iops),ngas)
        end if
      end do

      ! get contraction info
      icnt = ifact(4)
      iocc_cnt(1:ngastp,1:2) = 0
      do while (icnt.gt.0)
        jdx = mod(icnt,narc+1)
        icnt = icnt/(narc+1)
        ! test whether contraction was meant for op(1)*op(2) or op(2)*op(1)
        if (idxlist(contr%arc(jdx)%link(1),
     &              ivtx_expand(1,1),nexpand(1),1).gt.0) then
c dbg
          if (idxlist(contr%arc(jdx)%link(2),
     &         ivtx_expand(1,2),nexpand(2),1).le.0) then
            write(luout,*) 'arc%link: ',contr%arc(jdx)%link(1:2)
            write(luout,*) ' list 1: ',ivtx_expand(1:nexpand(1),1)
            write(luout,*) ' list 2: ',ivtx_expand(1:nexpand(2),2)
            write(luout,*) 'inconsistency!'
            call quit(1,'get_bc_info','inconsistency 1')
          end if
c dbg            
          ! add occupation information
          iocc_cnt = iocc_add(+1,iocc_cnt,.false.,+1,
     &           contr%arc(jdx)%occ_cnt,.false.)
        else
c dbg            
          if (idxlist(contr%arc(jdx)%link(2),
     &         ivtx_expand(1,1),nexpand(1),1).le.0 .or.
     &         idxlist(contr%arc(jdx)%link(1),
     &         ivtx_expand(1,2),nexpand(2),1).le.0) then
            write(luout,*) 'inconsistency (2)!'
            write(luout,*) 'arc%link: ',contr%arc(jdx)%link(1:2)
            write(luout,*) ' list 1: ',ivtx_expand(1:nexpand(1),1)
            write(luout,*) ' list 2: ',ivtx_expand(1:nexpand(2),2)
            call quit(1,'get_bc_info','inconsistency 2')
          end if
c dbg            
          ! add the complex conjugate:
          iocc_cnt = iocc_add(+1,iocc_cnt,.false.,+1,
     &                contr%arc(jdx)%occ_cnt,.true.)
        end if
      end do

      ! partition occupations into contraction and external indices
      iocc_ext(1:ngastp,1:2,1) = iocc_add(1,iocc(1,1,1),.false.,
     &                                -1,iocc_cnt,.false.)
      iocc_ext(1:ngastp,1:2,2) = iocc_add(1,iocc(1,1,2),.false.,
     &                                -1,iocc_cnt,.true.)

      ! set name of intermediate:
      iscr(1:nexpand(1)) = ivtx_expand(1:nexpand(1),1)
      iscr(nexpand(1)+1:nexpand(1)+nexpand(2)) =
     &     ivtx_expand(1:nexpand(2),2)
      call isort(iscr,nexpand(1)+nexpand(2),+1)
      ninter = ninter+1
      ! store occupation of intermediate
      iocc_int(1:ngastp,1:2,ninter) =
     &       iocc_add(1,iocc_ext(1,1,1),.false.,
     &                      1,iocc_ext(1,1,2),.false.)
      ! and its total Ms
      mstint(ninter) = mstop(1)+mstop(2)
      igamtint(ninter) = multd2h(igamtop(1),igamtop(2))

      ! restriction on intermediate
      ! preliminary version:
      ! defined by restriction on result
      call fit_restr(irestr_int(1,1,1,1,ninter),iocc_int(1,1,ninter),
     &     irestr_res,ihpvgas,ngas)
      if (ntest.ge.100) then
        write(luout,*) 'restriction on intermediate:'
        call wrt_rstr(luout,irestr_int,ngas)
      end if
      ! end preliminary

      ! store name of intermediate
      interm(ninter) = int_pack(iscr,nexpand(1)+nexpand(2),nvtx+1)

      return
      end
