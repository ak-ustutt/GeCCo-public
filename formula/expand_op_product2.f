*----------------------------------------------------------------------*
      subroutine expand_op_product2(form_list,idx_res,
     &                             fac,nvtx,nops,
     &                             idx_op_vtx,idx_sv_vtx,
     &                             iblk_min_in, iblk_max_in,
     &                             connect,nconnect,
     &                             avoid,navoid,
     &                             inproj,ninproj,
     &                             fix_in,op_info)
*----------------------------------------------------------------------*
*     given a list of operator indices and a result operator generate 
*     all contractions arising from
*
*        fac * <0|Op(1)Op(2) ....Op(n)|0>
*
*     idx_res: operator describing resulting "shape" of contraction
*     idx_op_vtx(nvtx): operators to be contracted, including the place
*         holders for the positions where open lines may occur
*       operators with more than one vertex must appear n-times in the
*       definition
*     idx_sv_vtx(nvtx): super-vertex information, needed for operators with
*       more than one vertex, and for open line definitions
*     examples: no open lines (closed diagram):
*              <0|Op(1)Op(2)Op(3)|0>  -> (/1,2,3/) -> idx_op_vtx
*                                        (/1,2,3/) -> idx_sv_vtx
*               open diagram, no open lines in between:
*              <0|XOp(1)Op(2)Op(3)X|0> -> (/X,1,2,3,X/) -> idx_op_vtx
*                                         (/1,2,3,4,1/) -> idx_sv_vtx
*               closed diagram, open lines in between:
*              <0|Op(1)XOp(2)Op(3)|0> -> (/1,X,X,2,3/) -> idx_op_vtx
*                                        (/1,2,2,3,4/) -> idx_sv_vtx
*               open diagram, open lines in between:
*              <0|XOp(1)XOp(2)Op(3)X|0> -> (/X,1,X,X,2,3,X/) -> idx_op_vtx
*                                          (/1,2,1,1,3,4,1/) -> idx_sv_vtx
*       where X is the number of the result operator
*
*     connect(2,nconnect): list of operator pairs that must be connected
*                          e.g. (1,3) if first and third operator should
*                          be connected.
*     inproj(4,ninproj): list of operator pairs that must be connected 
*                  (1st and 2nd index) but whose
*                  connection results from an rank n (3rd index) 
*                  inner projection (needed for R12)
*                  type of projector on 4th index:
*
*                   type     rank 1       rank 2        rank n
*                 ----------------------------------------------
*                    0        all          all           all
*                    1        |H><H|       Q(ans.1)      not impl.
*                    2        |P><P|       P(ans.1)
*                    3        |V><V|       Q(ans.3)
*                    4        |X><X|       P(ans.3)           
*
*     andreas, june 2007
*     new version with general result shapes, march 2008
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
      logical, intent(in) ::
     &     fix_in
      integer, intent(in) ::
     &     nvtx, nops, nconnect, navoid, ninproj,
     &     idx_res, idx_op_vtx(nvtx), idx_sv_vtx(nvtx),
     &     iblk_min_in(nops), iblk_max_in(nops),
     &     connect(2,nconnect), avoid(2,navoid), inproj(4,ninproj)
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     form_pnt
      type(contraction) ::
     &     proto
      type(operator), pointer ::
     &     op_res
      type(operator_array), pointer ::
     &     ops(:)
      logical ::
     &     ok
      integer ::
     &     narc, iarc, iblk_res, iop, jop, ivtx, nopen, idx, ioff, inc,
     &     num_res, nvtx_res, njoined_res, iblk_res_min, iblk_res_max,
     &     nterms, sv_open
      integer ::
     &     iblk_min(nops), iblk_max(nops), iblk_op(nops),
     &     occ_test(ngastp,2), occ_temp(ngastp,2)
      integer, pointer ::
     &     occ_vtx(:,:,:), neqv(:), idx_eqv(:,:), iop_typ(:), ol_map(:),
     &     idx_op(:), svertex(:), joined(:,:)
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     vtx_type

      if (ntest.ge.100) then
        write(lulog,*) '============================='
        write(lulog,*) ' Info from expand_op_product'
        write(lulog,*) '============================='
        write(lulog,*) ' idx_res = ',idx_res
        write(lulog,*) ' fac: ',fac
        write(lulog,*) ' nvtx, nops, nconnect: ',nvtx,nops,nconnect
        write(lulog,*) ' idx_op_vtx: ',idx_op_vtx(1:nvtx)
        write(lulog,*) ' idx_sv_vtx: ',idx_sv_vtx(1:nvtx)
        write(lulog,*) ' connect:',connect(1:2,1:nconnect)
        write(lulog,*) ' avoid:',avoid(1:2,1:navoid)
        write(lulog,*) ' inproj: ', inproj(1:2,1:ninproj)
      end if

      if (idx_res.lt.0) then
        write(lulog,*)
     &       'idx_res < 0 detected; if you meant to set up a formula'
        write(lulog,*)
     &       'for the adjoint operator: set up the formula for the'
        write(lulog,*)
     &       'actual operator and use transpose_formula()'
        call quit(1,'expand_op_product2','incorrect usage')
      end if

      form_pnt => form_list

      ! initialize proto-contraction
      call init_contr(proto)
c      nvtx = nops
      narc = nconnect+navoid+ninproj
      call resize_contr(proto,nvtx,narc,0,0)
      proto%nvtx = nvtx
c      ! currently, we expand primitive operators only
      proto%svertex = idx_sv_vtx
      call update_svtx4contr(proto)
      ! actual number of operators (=super-vertices), set in prev. routine
      if (nops.ne.proto%nsupvtx) then
        call quit(1,'expand_op_product2',
     &              'number of operators does not conform with '//
     &              'super-vertex information')
      end if
      proto%narc = narc
      proto%fac  = fac
      proto%idx_res = idx_res

      do ivtx = 1, nvtx
        ! set info
        proto%vertex(ivtx)%idx_op = abs(idx_op_vtx(ivtx))
        proto%vertex(ivtx)%dagger = idx_op_vtx(ivtx).lt.0
      end do

      ! set connectivity info on proto-contr.
      do iarc = 1, nconnect
        proto%arc(iarc)%link(1) = min(connect(1,iarc),connect(2,iarc)) 
        proto%arc(iarc)%link(2) = max(connect(1,iarc),connect(2,iarc)) 
        proto%arc(iarc)%occ_cnt = -1
      end do

      ioff = nconnect
      do iarc = 1, navoid
        proto%arc(ioff+iarc)%link(1) = min(avoid(1,iarc),avoid(2,iarc)) 
        proto%arc(ioff+iarc)%link(2) = max(avoid(1,iarc),avoid(2,iarc)) 
        proto%arc(ioff+iarc)%occ_cnt = 0
      end do

      allocate(ol_map(nvtx),idx_op(nops),fix_vtx(nvtx))
      fix_vtx = fix_in

      ! identify open-line vertices
      nopen = 0
      sv_open = 0
      do ivtx = 1, nvtx
        iop = idx_op_vtx(ivtx)
        ! open line vertex? (Should have same supervtx number)
        if (iop.eq.idx_res.and.sv_open.eq.0)
     &     sv_open = idx_sv_vtx(ivtx)
        if (iop.eq.idx_res.and.idx_sv_vtx(ivtx).eq.sv_open) then
          nopen = nopen+1
          idx = (nopen+1)/2
          if (mod(nopen+1,2).eq.1) idx = -idx
          ol_map(ivtx) = idx
        else
          ol_map(ivtx) = 0
        end if
      end do

      if (mod(nopen,2).ne.0)
     &     call quit(1,'expand_op_product2',
     &     'expected even number of open line vertices')

      ! assign numbers 1 .. nops to operators
      do ivtx = 1, nvtx        
        idx_op(proto%svertex(ivtx)) = idx_op_vtx(ivtx)
      end do

      ! number of result operator
      num_res = 0 ! no open lines?
      do ivtx = 1, nvtx
        if (ol_map(ivtx).ne.0) then
          num_res = proto%svertex(ivtx)
        end if
      end do

      ! we use svertex instead of input idx_sv
      svertex => proto%svertex
      joined => proto%joined

      ! set some pointers to operators
      op_res => op_info%op_arr(idx_res)%op
      allocate(ops(nops))
      do iop = 1, nops
c dbg*
c        print *,'iop, idx_op(iop): ',iop, idx_op(iop)
c dbg*
        ops(iop)%op => op_info%op_arr(abs(idx_op(iop)))%op
      end do

      if (num_res.eq.0) then
        nvtx_res = 0
      else
        nvtx_res = joined(0,num_res)
        if (nvtx_res.ne.2*op_res%njoined) then
          write(lulog,*) 'nvtx_res, njoined: ',nvtx_res,op_res%njoined
          write(lulog,*) 'nvtx_res must be 2 times njoined!'
          call quit(1,'expand_op_product2','inconsistency')
        end if
        njoined_res = nvtx_res/2
      end if

c dbg*
c      print *,'nvtx_res',nvtx_res
c      print *,'njoined_res',njoined_res
c      print *,'ol_map: ',ol_map
c      print *,'svertx: ',proto%svertex
c      print *,'joined: '
c      do iop = 1, nops
c        print '(x,i3," > ",i3," : ",10i3)',iop,proto%joined(0,iop),
c     &       proto%joined(1:proto%joined(0,iop),iop)
c      end do
c      print *,'idx_op: ',idx_op
c      print *,'num_res:',num_res
c dbg*
      ! aux-array needed for gen_contr: occupations of vertices
      allocate(occ_vtx(ngastp,2,nvtx))

      ! allocate and initialize equivalence table
      allocate(neqv(nops),idx_eqv(nops,nops),iop_typ(nops))
      neqv(1:nops) = 1
      do iop = 1, nops
        idx_eqv(1,iop) = iop
        idx_eqv(2:nops,iop) = 0
      end do

      do iop = 1, nops
c dbg*
c        print *,'iop = ',iop
c dbg*
        iop_typ(iop) = vtx_type(ops(iop)%op)
      end do

      ! look for equivalent commuting operators in order to avoid
      ! calls to gen_contr for unwanted combinations
      ! we need not be 100% perfect here, as gen_contr has the
      ! full machinery to decide whether a certain contribution
      ! is allowed or not, but we do not want to bother it too
      ! much with crappy combinations
      ! extended to non-commuting op.s if explicitly declared as
      ! exitation or deexcitation op.s
      do iop = 1, nops
        if (iop.eq.num_res) then
          iblk_max(iop) = 1 ! do not advance result blocks
                            ! this is done in a separate loop
          cycle
        end if
        if (iblk_max_in(1).lt.0) then
          iblk_max(iop) = ops(iop)%op%n_occ_cls
c dbg*
c          print *,' iop = ',iop,
c     &         ' "',trim(ops(iop)%op%name),'" >',iblk_max(iop)
c dbg*
        else if (iblk_max_in(iop).le.0) then
          iblk_max(iop) = ops(iop)%op%n_occ_cls
        else
          iblk_max(iop) = min(iblk_max_in(iop),ops(iop)%op%n_occ_cls)
        end if
      end do

      do iop = 1, nops
        do jop = 1, iop-1
          ! fix_in: don't assume that we need all permutations of commuting ops
          if (fix_in.or.neqv(jop).lt.0) cycle
          if (idx_op_vtx(iop).eq.idx_op_vtx(jop) .and.
     &        (iop_typ(iop).eq.iop_typ(jop) .and.
     &         (iop_typ(iop).eq.vtxtyp_ph_ex .or.
     &          iop_typ(iop).eq.vtxtyp_ph_dx) .or.
               ! also op.s with valence indices if explicitly declared as
     &         ops(iop)%op%species.eq.ops(jop)%op%species .and.
     &         (ops(iop)%op%species.eq.1 .or.     ! excitation or
     &          ops(iop)%op%species.eq.2))) then  ! deexcitation op.
            neqv(jop) = neqv(jop)+1
            neqv(iop) = -1
            idx_eqv(neqv(jop),jop) = iop
          end if
        end do
      end do
      if (iblk_min_in(1).lt.0) then
        iblk_min(1:nops) = 1
      else
        iblk_min(1:nops) = iblk_min_in(1:nops)
        ! do not advance result blocks:
        if (num_res.gt.0) iblk_min(num_res) = 1
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'iop_typ: ',iop_typ(1:nops)
        write(lulog,*) 'iblk_min:',iblk_min(1:nops)
        write(lulog,*) 'iblk_max:',iblk_max(1:nops)
        write(lulog,*) 'neqv:    ',neqv(1:nops)
        write(lulog,*) 'idx_eqv: '
        call iwrtma(idx_eqv,nops,nops,nops,nops)
      end if

      ! --------------------------
      !  loop over result blocks
      ! --------------------------
      iblk_res_min = 1
      iblk_res_max = 1
      if (num_res.gt.0) then
        if (iblk_min_in(1).gt.0 .and.
     &      iblk_min_in(num_res).gt.0)
     &       iblk_res_min = iblk_min_in(num_res)
        iblk_res_max = op_res%n_occ_cls
        if (iblk_max_in(1).gt.0 .and.
     &      iblk_max_in(num_res).gt.0)
     &       iblk_res_max = min(iblk_max_in(num_res),op_res%n_occ_cls)
      end if

      do iblk_res = iblk_res_min, iblk_res_max
        proto%iblk_res = iblk_res

        ioff = (iblk_res-1)*njoined_res
        do ivtx = 1, nvtx_res
          idx = (ivtx+1)/2
          if (mod(ivtx+1,2).eq.0) then
            occ_vtx(1:ngastp,1:2,joined(ivtx,num_res)) =
     &           iocc_dagger(iocc_xdn(1,
     &             op_res%ihpvca_occ(1:ngastp,1:2,ioff+idx)))
          else
            occ_vtx(1:ngastp,1:2,joined(ivtx,num_res)) =
     &           iocc_dagger(iocc_xdn(2,
     &             op_res%ihpvca_occ(1:ngastp,1:2,ioff+idx)))
          end if
          proto%vertex(joined(ivtx,num_res))%iblk_op = ioff+ivtx
c          proto%vertex(joined(ivtx,num_res))%iblk_op = iblk_res
        end do

        ! first distribution of blocks
        !  the others will be generated by next_dist2() at
        !  the bottom of the loop
        iblk_op(1:nops) = iblk_min(1:nops)
        do

          if (ntest.ge.100) then
            write(lulog,*) 'current dist: ',iblk_op(1:nops)
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
            write(lulog,*) 'check1: ',ok
          end if

          if (ok)  then

            ! set up occ_vtx
            do iop = 1, nops
              if (iop.eq.num_res) cycle
              if (idx_op(iop).gt.0) then
                ioff = (iblk_op(iop)-1)*ops(iop)%op%njoined
                do ivtx = 1, joined(0,iop)
                  occ_vtx(1:ngastp,1:2,joined(ivtx,iop)) =
     &              ops(iop)%op%ihpvca_occ(1:ngastp,1:2,ioff+ivtx)
                end do
              else
                ioff = (iblk_op(iop))*ops(iop)%op%njoined+1
                do ivtx = 1, joined(0,iop)
                  occ_vtx(1:ngastp,1:2,joined(ivtx,iop)) =
     &                 iocc_dagger(                  
     &              ops(iop)%op%ihpvca_occ(1:ngastp,1:2,ioff-ivtx))
                end do
              end if
c dbg*
c              print *,'iop, ioff, n: ',iop, ioff, joined(0,iop)
c              print *,'              ',joined(1:joined(0,iop),iop)
c              print *,'proto'
c              call prt_contr2(lulog,proto,op_info)
c              write(lulog,*) 'occ_vtx test:'
c              call wrt_occ_n(lulog,occ_vtx,nvtx)
c dbg*
              
            end do

            if (ninproj.gt.0) then
              ! Handle inner projection
              ! reset number of arcs
              proto%narc = nconnect+navoid
              call set_inproj(proto,occ_vtx,ok,inproj,ninproj)
            else
              ! check whether contraction is possible:
              ! sum_i [EX_op(i)] - [DX_op(i)]^dag = 0
              occ_test = 0
              do ivtx = 1, nvtx
                occ_test = occ_test
     &             + iocc_xdn(1,occ_vtx(1:ngastp,1:2,ivtx))
     &             - iocc_dagger(iocc_xdn(2,occ_vtx(1:ngastp,1:2,ivtx)))
              end do
              ok = iocc_zero(occ_test)
            end if


            if (ntest.ge.100) then
              write(lulog,*) 'occ_vtx:'
              call wrt_occ_n(lulog,occ_vtx,nvtx)
              write(lulog,*) 'test:'
              call wrt_occ(lulog,occ_test)
              write(lulog,*) 'check2: ',ok
            end if

          end if

          if (ok) then
            ! set remaining proto-contraction info
            do iop = 1, nops
              if (iop.eq.num_res) cycle
              ioff = (iblk_op(iop)-1)*ops(iop)%op%njoined
              if (idx_op(iop).gt.0) then
                do ivtx = 1, joined(0,iop)
                  proto%vertex(joined(ivtx,iop))%iblk_op = ioff+ivtx
C??                proto%vertex(joined(ivtx,iop))%iblk_op = iblk_op(iop)
                end do
              else
                ioff = ioff+ops(iop)%op%njoined+1
                do ivtx = 1, joined(0,iop)
                  proto%vertex(joined(ivtx,iop))%iblk_op = ioff-ivtx
                end do
              end if
            end do

            if (ntest.ge.1000) then
              write(lulog,*) 'passing proto contraction:'
              call prt_contr2(lulog,proto,op_info)
            end if

            ! generate contractions
            call gen_contr4(.false.,form_pnt,proto,
     &           fix_vtx,occ_vtx(1,1,1),ol_map,op_info)

            ! advance pointer
            nterms = 0
            do while(form_pnt%command.ne.command_end_of_formula)
              nterms = nterms+1
              form_pnt => form_pnt%next
            end do

            if (ntest.ge.100)
     &           write(lulog,*) nterms,' new terms ...'

          end if

          if (.not.next_dist2(iblk_op,nops,iblk_min,iblk_max,+1)) exit
        end do

      end do

      call dealloc_contr(proto)
      deallocate(neqv,idx_eqv,occ_vtx,iop_typ,ops,ol_map,idx_op,fix_vtx)

      return
      end

