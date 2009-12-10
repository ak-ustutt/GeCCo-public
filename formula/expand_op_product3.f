*----------------------------------------------------------------------*
      subroutine expand_op_product3(form_list,idx_res,
     &                             fac,nvtx,nops,
     &                             idx_op_vtx,idx_sv_vtx,
     &                             iblk_min_in, iblk_max_in,
     &                             connect,nconnect,
     &                             avoid,navoid,
     &                             descr,ndescr,
     &                             op_info)
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
*
*     descr: a descriptor string of the following form
*        <vtx1>,<vtx2>,<descr_C>,<descr_A>
*     
*         where <vtx1>,<vtx2> vertex numbers (both set: contraction,
*                                        only vtx1 set: restriction)
*
*         and <descr_C>,<descr_A> C/A descriptor strings, a short hand
*            form for occupations, best explaind by examples:
*
*              HP  :     equivalent to (1,1,0,0)
*              HHX :     equivalent to (2,0,0,1)
*              H[PX]:    equivalent to (1,1,0,0),(1,0,0,1)
*
*        so e.g. 2,4,[HPX],H requests a contraction of vertices 2 and 4
*                            over all orbitals in vtx1(C) (=vtx2(A))
*                            and over H in vtx1(A) (=vtx2(C))
*                2,-,,HH     constrains the allowed occupations of
*                            vertex 2 to those having two hole indices in A
*
*     andreas, june 2007
*     new version with general result shapes, march 2008
*     new version with inproj replaced by descriptors, nov 2009
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_occ_list.h'
      include 'ifc_operators.h'
      include 'ifc_baserout.h'

      integer, parameter ::
     &     ntest = 00
      character, parameter ::
     &     i_am = 'expand_op_product3'

      integer, parameter ::
     &     maxlist_raw = 50

      type(formula_item), intent(inout), target ::
     &     form_list
      real(8), intent(in) ::
     &     fac
      integer, intent(in) ::
     &     nvtx, nops, nconnect, navoid, 
     &     idx_res, idx_op_vtx(nvtx), idx_sv_vtx(nvtx),
     &     iblk_min_in(nops), iblk_max_in(nops),
     &     connect(2,nconnect), avoid(2,navoid), ndescr
      character(len=*) ::
     &     descr(ndescr)
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     form_pnt
      type(contraction) ::
     &     proto0, proto
      type(operator), pointer ::
     &     op_res
      type(operator_array), pointer ::
     &     ops(:)
      type(occ_list) ::
     &     contr_list, occrs_list
      logical ::
     &     ok, vtx_ok
      integer ::
     &     narc, iarc, iblk_res, iop, jop, ivtx, vtx1, vtx2, nopen,
     &     idx, jdx, ioff, inc, istart, iend,
     &     num_res, nvtx_res, njoined_res, iblk_res_min, iblk_res_max,
     &     nterms, n_occrs, n_cntpairs, nlist
      integer ::
     &     iblk_min(nops), iblk_max(nops), iblk_op(nops),
     &     occ_test(ngastp,2), occ_temp(ngastp,2)
      integer, pointer ::
     &     occ_vtx(:,:,:), neqv(:), idx_eqv(:,:), iop_typ(:), ol_map(:),
     &     idx_op(:), svertex(:), joined(:,:), occ_list_raw(:,:,:),
     &     vtx_occrs(:,:), vtx_cntrq(:,:), occrs(:,:,:), cntrq(:,:,:),
     &     cnt_comb(:), cnt_comb_mnmx(:,:)

      integer, external ::
     &     vtx_type

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,i_am)
        write(luout,*) ' idx_res = ',idx_res
        write(luout,*) ' fac: ',fac
        write(luout,*) ' nvtx, nops, nconnect: ',nvtx,nops,nconnect
        write(luout,*) ' idx_op_vtx: ',idx_op_vtx(1:nvtx)
        write(luout,*) ' idx_sv_vtx: ',idx_sv_vtx(1:nvtx)
        write(luout,*) ' connect:',connect(1:2,1:nconnect)
        write(luout,*) ' avoid:',avoid(1:2,1:navoid)
        write(luout,*) ' dscr:'
        do idx = 1, ndescr
          write(luout,*) ' : ', trim(descr(idx))
        end do
      end if

      if (idx_res.lt.0) then
        write(luout,*)
     &       'idx_res < 0 detected; if you meant to set up a formula'
        write(luout,*)
     &       'for the adjoint operator: set up the formula for the'
        write(luout,*)
     &       'actual operator and use transpose_formula()'
        call quit(1,i_am,'incorrect usage')
      end if

      form_pnt => form_list

      ! initialize proto-contraction
      call init_contr(proto0)
      call init_contr(proto)
c      nvtx = nops
      narc = nconnect+navoid
      call resize_contr(proto0,nvtx,narc,0,0)
      proto0%nvtx = nvtx
c      ! currently, we expand primitive operators only
      proto0%svertex = idx_sv_vtx
      call update_svtx4contr(proto0)
      ! actual number of operators (=super-vertices), set in prev. routine
      if (nops.ne.proto0%nsupvtx) then
        call quit(1,i_am,
     &              'number of operators does not conform with '//
     &              'super-vertex information')
      end if
      proto0%narc = narc
      proto0%fac  = fac
      proto0%idx_res = idx_res

      do ivtx = 1, nvtx
        ! set info
        proto0%vertex(ivtx)%idx_op = abs(idx_op_vtx(ivtx))
        proto0%vertex(ivtx)%dagger = idx_op_vtx(ivtx).lt.0
      end do

      ! set connectivity info on proto-contr.
      do iarc = 1, nconnect
        proto0%arc(iarc)%link(1) = min(connect(1,iarc),connect(2,iarc)) 
        proto0%arc(iarc)%link(2) = max(connect(1,iarc),connect(2,iarc)) 
        proto0%arc(iarc)%occ_cnt = -1
      end do

      ioff = nconnect
      do iarc = 1, navoid
        proto0%arc(ioff+iarc)%link(1) = min(avoid(1,iarc),avoid(2,iarc)) 
        proto0%arc(ioff+iarc)%link(2) = max(avoid(1,iarc),avoid(2,iarc)) 
        proto0%arc(ioff+iarc)%occ_cnt = 0
      end do

      if (ndescr.gt.0) then
        call init_occ_list(contr_list)
        call init_occ_list(occrs_list)
        allocate(occ_list_raw(ngastp,2,maxlist_raw))
      end if

      ! process descriptor info
      do idx = 1, ndescr
        ! decode and merge
        call process_descriptor(vtx1,vtx2,occ_list_raw,nlist,
     &                                    descr(idx),maxlist_raw)
        if (vtx2.gt.0) then
          call update_occ_list(contr_list,
     &                         vtx1,vtx2,occ_list_raw,nlist)
        else
          call update_occ_list(occrs_list,
     &                         vtx1,0,occ_list_raw,nlist)
        end if
      end do

      if (ndescr.gt.0) then
        n_occrs = occrs_list%n_vtx_inf
        vtx_occrs => occrs_list%vtx_inf
        occrs => occrs_list%occ

        n_cntpairs = contr_list%n_vtx_inf
        vtx_cntrq => contr_list%vtx_inf
        cntrq => contr_list%occ

        deallocate(occ_list_raw)
      end if

      if (n_cntpairs.gt.0) then
        allocate(cnt_comb(n_cntpairs),cnt_comb_mnmx(2,n_cntpairs))
        do idx = 1, n_cntpairs
          cnt_comb_mnmx(1,idx) = vtx_cntrq(3,idx)
          cnt_comb_mnmx(2,idx) = vtx_cntrq(4,idx)
        end do
      else
        allocate(cnt_comb(1),cnt_comb_mnmx(2,1))
        cnt_comb_mnmx(1:2,1) = 1
      end if

      allocate(ol_map(nvtx),idx_op(nops))

      ! identify open-line vertices
      nopen = 0
      do ivtx = 1, nvtx
        iop = idx_op_vtx(ivtx)
        ! open line vertex?
        if (iop.eq.idx_res) then
          nopen = nopen+1
          idx = (nopen+1)/2
          if (mod(nopen+1,2).eq.1) idx = -idx
          ol_map(ivtx) = idx
        else
          ol_map(ivtx) = 0
        end if
      end do

      if (mod(nopen,2).ne.0)
     &     call quit(1,i_am,
     &     'expected even number of open line vertices')

      ! assign numbers 1 .. nops to operators
      do ivtx = 1, nvtx        
        idx_op(proto0%svertex(ivtx)) = idx_op_vtx(ivtx)
      end do

      ! number of result operator
      num_res = 0 ! no open lines?
      do ivtx = 1, nvtx
        if (ol_map(ivtx).ne.0) then
          num_res = proto0%svertex(ivtx)
        end if
      end do

      ! we use svertex instead of input idx_sv
      svertex => proto0%svertex
      joined => proto0%joined

      ! set some pointers to operators
      op_res => op_info%op_arr(idx_res)%op
      allocate(ops(nops))
      do iop = 1, nops
        ops(iop)%op => op_info%op_arr(abs(idx_op(iop)))%op
      end do

      if (num_res.eq.0) then
        nvtx_res = 0
      else
        nvtx_res = joined(0,num_res)
        if (nvtx_res.ne.2*op_res%njoined) then
          write(luout,*) 'nvtx_res, njoined: ',nvtx_res,op_res%njoined
          call quit(1,i_am,'inconsistency')
        end if
        njoined_res = nvtx_res/2
      end if

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
        iop_typ(iop) = vtx_type(ops(iop)%op)
      end do

      ! look for equivalent commuting operators in order to avoid
      ! calls to gen_contr for unwanted combinations
      ! we need not be 100% perfect here, as gen_contr has the
      ! full machinery to decide whether a certain contribution
      ! is allowed or not, but we do not want to bother it too
      ! much with crappy combinations
      do iop = 1, nops
        if (iop.eq.num_res) then
          iblk_max(iop) = 1 ! do not advance result blocks
                            ! this is done in a separate loop
          cycle
        end if
        if (iblk_max_in(1).lt.0) then
          iblk_max(iop) = ops(iop)%op%n_occ_cls
        else if (iblk_max_in(iop).le.0) then
          iblk_max(iop) = ops(iop)%op%n_occ_cls
        else
          iblk_max(iop) = iblk_max_in(iop)
        end if
      end do

      do iop = 1, nops
        do jop = 1, iop-1
          if (neqv(jop).lt.0) cycle
          if (idx_op_vtx(iop).eq.idx_op_vtx(jop) .and.
     &        iop_typ(iop).eq.iop_typ(jop) .and.
     &       (iop_typ(iop).eq.vtxtyp_ph_ex .or.
     &        iop_typ(iop).eq.vtxtyp_ph_dx) ) then
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
      iblk_res_min = 1
      iblk_res_max = 1
      if (num_res.gt.0) then
        if (iblk_min_in(1).gt.0 .and.
     &      iblk_min_in(num_res).gt.0)
     &       iblk_res_min = iblk_min_in(num_res)
        iblk_res_max = op_res%n_occ_cls
        if (iblk_max_in(1).gt.0 .and.
     &      iblk_max_in(num_res).gt.0)
     &       iblk_res_max = iblk_max_in(num_res)
      end if

      do iblk_res = iblk_res_min, iblk_res_max
        proto0%iblk_res = iblk_res

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
          proto0%vertex(joined(ivtx,num_res))%iblk_op = ioff+idx
c          proto0%vertex(joined(ivtx,num_res))%iblk_op = ioff+ivtx
c          proto0%vertex(joined(ivtx,num_res))%iblk_op = iblk_res
        end do

        ! first distribution of blocks
        !  the others will be generated by next_dist2() at
        !  the bottom of the loop
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
              
            end do

          end if

          ! check that occupations conform with restrictions
          ! from descriptors
          ok = .true.
          do idx = 1, n_occrs
            ivtx = vtx_occrs(1,idx)
            istart = vtx_occrs(3,idx)
            iend = vtx_occrs(4,idx)
            vtx_ok = .false.
            do jdx = istart, iend
              vtx_ok = vtx_ok.or.(iocc_equal(occrs(:,:,jdx),.false.,
     &                                      occ_vtx(:,:,ivtx), .false.))
            end do
            ok = ok.and.vtx_ok
          end do

          if (ok) then
 
            ! set remaining proto-contraction info
            do iop = 1, nops
              if (iop.eq.num_res) cycle
              ioff = (iblk_op(iop)-1)*ops(iop)%op%njoined
              if (idx_op(iop).gt.0) then
                do ivtx = 1, joined(0,iop)
                  proto0%vertex(joined(ivtx,iop))%iblk_op = ioff+ivtx
                end do
              else
                ioff = ioff+ops(iop)%op%njoined+1
                do ivtx = 1, joined(0,iop) 
                  proto0%vertex(joined(ivtx,iop))%iblk_op = ioff-ivtx
                end do
              end if
            end do

           
            ! loop over all contraction descriptors
            ! (or a single turn, if no descriptors are given)
            if (n_cntpairs.gt.0)
     &           cnt_comb(1:n_cntpairs) = cnt_comb_mnmx(1,1:n_cntpairs)
            do

              ! get copy of present proto-contraction
              call copy_contr(proto0,proto)
              
              if (n_cntpairs.gt.0) then
                ! set proto contraction accordingly
                do idx = 1, n_cntpairs
                  vtx1 = vtx_cntrq(1,idx)
                  vtx2 = vtx_cntrq(2,idx)
                  call contr_add_arc(proto,vtx1,vtx2,
     &                               cntrq(:,:,cnt_comb(idx)))
                  if (ntest.ge.100) then
                    write(luout,*) 'trying ',vtx1,vtx2
                    call wrt_occ(luout,cntrq(:,:,cnt_comb(idx)))
                  end if
                end do
                if (ntest.ge.100) then
                  write(luout,*) 'current proto:'
                  call prt_contr2(luout,proto,op_info)
                end if

                ! check that the current request can be satisfied
                ok = .true.
                occ_temp = 0
                do ivtx = 1, proto%nvtx
                  call get_unconnected4vertex
     &                 (occ_test,ivtx,proto,op_info)
                  ok = ok.and.iocc_bound('>=',occ_test ,.false.,
     &                                        occ_temp,.false.)
                  if (ntest.ge.100) then
                    write(luout,*) 'vtx # ',ivtx,'  ok: ',ok
                    call wrt_occ(luout,occ_test)
                  end if
                end do

              else  ! works this way only, if no contractions are pre-def.d
              ! check whether contraction is possible:
              ! sum_i [EX_op(i)] - [DX_op(i)]^dag = 0
                occ_test = 0
                do ivtx = 1, nvtx
                  occ_test = occ_test
     &               + iocc_xdn(1,occ_vtx(1:ngastp,1:2,ivtx))
     &               - iocc_dagger(
     &                 iocc_xdn(2,occ_vtx(1:ngastp,1:2,ivtx)))
                end do
                ok = iocc_zero(occ_test)
              end if


              if (ntest.ge.100) then
                write(luout,*) 'occ_vtx:'
                call wrt_occ_n(luout,occ_vtx,nvtx)
                write(luout,*) 'test:'
                call wrt_occ(luout,occ_test)
                write(luout,*) 'check2: ',ok
              end if

              if (ok) then
c              ! set remaining proto-contraction info
c                do iop = 1, nops
c                  if (iop.eq.num_res) cycle
c                  ioff = (iblk_op(iop)-1)*ops(iop)%op%njoined
c                  if (.not.proto%vertex(iop)%dagger) then
c                    do ivtx = 1, joined(0,iop)
c                      proto%vertex(joined(ivtx,iop))%iblk_op = ioff+ivtx
c                    end do
c                  else
c                    ioff = ioff+ops(iop)%op%njoined+1
c                    do ivtx = 1, joined(0,iop)
c                      proto%vertex(joined(ivtx,iop))%iblk_op = ioff-ivtx
c                    end do
c                  end if
c                end do

                if (ntest.ge.1000) then
                  write(luout,*) 'passing proto contraction:'
                  call prt_contr2(luout,proto,op_info)
                end if

                ! generate contractions
                call gen_contr4(.false.,form_pnt,proto,
     &               occ_vtx(1,1,1),ol_map,op_info)

                ! advance pointer
                nterms = 0
                do while(form_pnt%command.ne.command_end_of_formula)
                  nterms = nterms+1
                  form_pnt => form_pnt%next
                end do

                if (ntest.ge.100)
     &               write(luout,*) nterms,' new terms ...'

              end if

              if (n_cntpairs.eq.0) exit
              if (.not.next_dist(cnt_comb,n_cntpairs,cnt_comb_mnmx,+1))
     &             exit

            end do ! loop over poss. contr.

          end if

          if (.not.next_dist2(iblk_op,nops,iblk_min,iblk_max,+1)) exit
        end do

      end do

      call dealloc_contr(proto0,proto)
      deallocate(neqv,idx_eqv,occ_vtx,iop_typ,ops,ol_map,idx_op)
      deallocate(cnt_comb,cnt_comb_mnmx)

      if (ndescr.gt.0) then
        call clean_occ_list(contr_list)
        call clean_occ_list(occrs_list)
      end if

      return
      end

