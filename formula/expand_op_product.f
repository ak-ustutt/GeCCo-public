*----------------------------------------------------------------------*
      subroutine expand_op_product(form_list,idx_res,
     &                             fac,nops,idx_op,
     &                             iblk_min_in, iblk_max_in,
     &                             connect,nconnect,force,
     &                             op_info)
*----------------------------------------------------------------------*
*     given a list of operator indices and a result operator generate 
*     all contractions arising from
*
*        fac * <C_res|Op(1)Op(2) ....Op(n)|A_res>
*
*     idx_res: operator describing resulting "shape" of contraction
*     idx_op(nops): operators to be contracted
*     connect(2,nconnect): list of operator pairs that must be connected
*                          e.g. (1,3) if first and third operator should
*                          be connected.
*     force is set true if we want the connection of the first and last 
*     vertices in a non-standard way. If this is so then the first 
*     element of connect should be (1,nconnect).
*
*     andreas, june 2007
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'ifc_operators.h'
      include 'ifc_baserout.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(inout), target ::
     &     form_list
      real(8), intent(in) ::
     &     fac
      integer, intent(in) ::
     &     nops, nconnect,
     &     idx_res, idx_op(nops), iblk_min_in(nops), iblk_max_in(nops),
     &     connect(2,nconnect)
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in)::
     &     force

      type(formula_item), pointer ::
     &     form_pnt
      type(contraction) ::
     &     proto
      type(operator), pointer ::
     &     op_res, op
      logical ::
     &     fix_vtx(nops)
      logical ::
     &     ok
      integer ::
     &     nvtx, narc, iarc, iblk_res, iop, jop
      integer ::
     &     iblk_min(nops), iblk_max(nops), iblk_op(nops),
     &     occ_test(ngastp,2), occ_temp(ngastp,2)
      integer, pointer ::
     &     occ_vtx(:,:,:), neqv(:), idx_eqv(:,:), iop_typ(:)

      integer, external ::
     &     vtx_type

      if (ntest.ge.100) then
        write(luout,*) '============================='
        write(luout,*) ' Info from expand_op_product'
        write(luout,*) '============================='
        write(luout,*) ' idx_res = ',idx_res
        write(luout,*) ' fac, nops, nconnect: ',fac,nops,nconnect
        write(luout,*) ' idx_op: ',idx_op(1:nops)
        write(luout,*) ' connect:',connect(1:2,1:nconnect)
      end if

      form_pnt => form_list

      ! initialize proto-contraction
      call init_contr(proto)
      nvtx = nops
      narc = nconnect
      call resize_contr(proto,nvtx,narc,0,0)
      proto%nvtx = nvtx
      ! currently, we expand primitive operators only
      call set_primitive_contr(proto)
      proto%narc = narc
      proto%fac  = fac
      proto%idx_res = idx_res
      ! set connectivity info on proto-contr.
      do iarc = 1, nconnect
        proto%arc(iarc)%link(1) = min(connect(1,iarc),connect(2,iarc)) 
        proto%arc(iarc)%link(2) = max(connect(1,iarc),connect(2,iarc)) 
        proto%arc(iarc)%occ_cnt = -1
      end do

      ! aux-array needed for gen_contr: occupations of vertices
      allocate(occ_vtx(ngastp,2,nops+1))
      
      fix_vtx = .false. ! do not fix anything

      op_res => op_info%op_arr(idx_res)%op

      ! allocate and initialize equivalence table
      allocate(neqv(nops),idx_eqv(nops,nops),iop_typ(nops))
      neqv(1:nops) = 1
      do iop = 1, nops
        idx_eqv(1,iop) = iop
        idx_eqv(2:nops,iop) = 0
      end do

      do iop = 1, nops
        iop_typ(iop) = vtx_type(op_info%op_arr(idx_op(iop))%op)
      end do

      ! look for equivalent commuting operators and set info on proto-contr.
      do iop = 1, nops
        op => op_info%op_arr(idx_op(iop))%op
        if (iblk_max_in(1).lt.0) then
          iblk_max(iop) = op%n_occ_cls
        else if (iblk_max_in(iop).le.0) then
          iblk_max(iop) = op%n_occ_cls
        else
          iblk_max(iop) = iblk_max_in(iop)
        end if
        proto%vertex(iop)%idx_op = idx_op(iop)
        do jop = 1, iop-1
          if (neqv(jop).lt.0) cycle
          if (idx_op(iop).eq.idx_op(jop) .and.
     &        iop_typ(iop).eq.iop_typ(jop) .and.
     &       (iop_typ(iop).eq.vtxtyp_ph_ex .or.
     &        iop_typ(iop).eq.vtxtyp_ph_dx) ) then
            neqv(jop) = neqv(jop)+1
            neqv(iop) = -1
            idx_eqv(neqv(jop),jop) = iop
          end if
        end do
      end do
      if (iblk_min_in(1).le.0) then
        iblk_min(1:nops) = 1
      else
        iblk_min(1:nops) = iblk_min_in(1:nops)
      end if

      if (ntest.ge.100) then
        write(luout,*) 'iop_typ: ',iop_typ(1:nops)
        write(luout,*) 'iblk_min:',iblk_min(1:nops)
        write(luout,*) 'iblk_max:',iblk_max(1:nops)
        write(luout,*) 'neqv:    ',neqv(1:nops)
        write(luout,*) 'idx_eqv: '
        call iwrtma(idx_eqv,nops,nops,nops,nops)
      end if

      ! --------------------------
      !  loop over result blocks
      ! --------------------------
      do iblk_res = 1, op_res%n_occ_cls
        proto%iblk_res = iblk_res
        occ_vtx(1:ngastp,1:2,1) =
     &       op_res%ihpvca_occ(1:ngastp,1:2,iblk_res)
        if (op_res%dagger) then
          occ_vtx(1:ngastp,1:2,1) = iocc_dagger(occ_vtx(1:ngastp,1:2,1))
        end if

        ! first distribution of blocks
        !  the others will be generated by next_dist2() at
        !  the bottom of the loop
c        iblk_op(1:nops) = 1
        iblk_op(1:nops) = iblk_min(1:nops)
        do

          if (ntest.ge.100) then
            write(luout,*) 'current dist: ',iblk_op(1:nops)
          end if

          ! check that equivalent commuting operators 
          ! are in ascending order
          ok = .true.
          do iop = 1, nops
            if (neqv(iop).gt.1) then
              do jop = 1, neqv(iop)-1
                ok = ok.and.iblk_op(idx_eqv(jop  ,iop)).le.
     &                      iblk_op(idx_eqv(jop+1,iop))
              end do
              if (.not.ok) exit
            end if
          end do

          if (ntest.ge.100) then
            write(luout,*) 'check1: ',ok
          end if

          if (ok)  then

            ! check whether contraction is possible:
            ! [EX_res]-[DX_res]^dag = sum_i [EX_op(i)] - [DX_op(i)]^dag
            occ_test = -iocc_xdn(1,occ_vtx(1:ngastp,1:2,1))
     &         + iocc_dagger(iocc_xdn(2,occ_vtx(1:ngastp,1:2,1)))
            do iop = 1, nops
              op => op_info%op_arr(idx_op(iop))%op
              if (.not.op%dagger) then
                occ_vtx(1:ngastp,1:2,iop+1) =
     &               op%ihpvca_occ(1:ngastp,1:2,iblk_op(iop))
              else
                occ_vtx(1:ngastp,1:2,iop+1) = iocc_dagger(
     &               op%ihpvca_occ(1:ngastp,1:2,iblk_op(iop)) )
              end if
              occ_test = occ_test
     &           + iocc_xdn(1,occ_vtx(1:ngastp,1:2,iop+1))
     &           - iocc_dagger(iocc_xdn(2,occ_vtx(1:ngastp,1:2,iop+1)))
            end do
            ok = .not.iocc_nonzero(occ_test)

            if (ntest.ge.100) then
              write(luout,*) 'occ_vtx:'
              call wrt_occ_n(luout,occ_vtx,nops+1)
              write(luout,*) 'test:'
              call wrt_occ(luout,occ_test)
              write(luout,*) 'check2: ',ok
            end if

          end if

          if (ok) then
            ! set remaining proto-contraction info
            do iop = 1, nops
              proto%vertex(iop)%iblk_op = iblk_op(iop)
            end do

            if(force)then
              occ_temp(1:ngastp,1)=0
              occ_temp(1:ngastp,2)=occ_vtx(1:ngastp,2,2)
              if(ielsqsum(occ_temp,ngastp*2).ne.0)then
                proto%arc(1)%occ_cnt = occ_temp
c     &               form_pnt%contr%arc(1)%occ_cnt+occ_temp  
              endif  
            endif  

            ! generate contractions
            call gen_contr2(form_pnt,proto,fix_vtx,occ_vtx,op_info)

            ! advance pointer
            do while(form_pnt%command.ne.command_end_of_formula)
              form_pnt => form_pnt%next
            end do
          end if

          if (.not.next_dist2(iblk_op,nops,iblk_min,iblk_max,+1)) exit
        end do

      end do

      call dealloc_contr(proto)
      deallocate(neqv,idx_eqv,occ_vtx,iop_typ)

      return
      end

