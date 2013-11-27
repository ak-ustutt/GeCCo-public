*----------------------------------------------------------------------*
      subroutine contr_insert
     &     (contr,op_info,nvtx,vtxinlist,idxins,iblkins,idxres)
*----------------------------------------------------------------------*
*     inserts block iblkins of operator idxins between the vertices
*     indicated by vtxinlist (as often as possible).
*
*     matthias, May 2010
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_contraction_list.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(inout) ::
     &     contr
      integer, intent(in) ::
     &     nvtx, idxins, iblkins, idxres
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     vtxinlist(nvtx)

      logical ::
     &     ok
      integer ::
     &     iins, jins, nvtx_new, nins, ivtx, jvtx,
     &     njoined_res, njoined_ins, ins_cur, ieqvfac, nins_old
      type(operator_array), pointer ::
     &     op_arr(:)

      type(operator), pointer ::
     &     op_ins, op_res
      integer(8) ::
     &     base, idxnew, occ_ins_x, occ_ins_d, overlap

      integer ::
     &     iocc_ins_x(ngastp,2), iocc_ins_d(ngastp,2)

      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:),
     &     vtx_new(:), topo_new(:,:), xlines_new(:,:),
     &     occ_sum_x(:), occ_sum_d(:), occ_tmp(:)
      integer, pointer ::
     &     svertex(:), svertex_new(:),
     &     occ_ins(:,:,:), insert(:), ins_at_vtx(:),
     &     neqv(:), idx_eqv(:,:)
      logical, pointer ::
     &     vil_new(:) !vtxinlist_new
      
      integer(8), external ::
     &     pack_vtx, int8_pack, occ_overlap_p
      integer, external ::
     &     ifac
      logical, external ::
     &     occ_bound_p

      base = pack_base

      if (ntest.ge.100) then
        write(lulog,*) '====================='
        write(lulog,*) ' contr_insert at work'
        write(lulog,*) '====================='
        write(lulog,*) ' Contraction on input:'
        call prt_contr2(lulog,contr,op_info)
        write(lulog,*) 'vtxinlist = ',vtxinlist
        write(lulog,*) 'idxins = ',idxins
        write(lulog,*) 'iblkins= ',iblkins
      end if

      op_arr => op_info%op_arr

      op_ins => op_arr(idxins)%op
      op_res => op_arr(idxres)%op

      njoined_ins = op_ins%njoined
      njoined_res = op_res%njoined

      if (njoined_ins.ne.1.or.njoined_res.ne.1)
     &       call quit(1,'contr_insert',
     &       'works only for njoined=1 operators so far')

      occ_ins => op_ins%ihpvca_occ
      iocc_ins_x = iocc_xdn(1,occ_ins(1:ngastp,1:2,iblkins))
      iocc_ins_d = iocc_xdn(2,occ_ins(1:ngastp,1:2,iblkins))
      occ_ins_x = int8_pack(iocc_ins_x,ngastp*2,base)
      occ_ins_d = int8_pack(iocc_ins_d,ngastp*2,base)

      allocate(vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,njoined_res))
      allocate(svertex(nvtx))
      allocate(occ_sum_x(nvtx),occ_sum_d(nvtx),insert(nvtx),
     &         ins_at_vtx(nvtx))

      call pack_contr(svertex,vtx,topo,xlines,contr,njoined_res)

      if (any(xlines(1:nvtx,1:njoined_res).ne.0))
     &        call quit(1,'contr_insert',
     &       'only for scalar target operator!')

      if (ntest.ge.100) then
        write(lulog,*) 'contraction in topo form:'
        call prt_contr_p(lulog,svertex,vtx,topo,xlines,nvtx,njoined_res)
      end if

      occ_sum_x(1:nvtx) = 0
      occ_sum_d(1:nvtx) = 0

      ! scalar?
      if (occ_ins_x.eq.0.and.occ_ins_d.eq.0) then
        ! insert scalar only if it is our only hope...
        if (.not.all(vtxinlist(1:nvtx))) then
          deallocate(vtx,topo,xlines,svertex,occ_sum_x,occ_sum_d,
     &               insert,ins_at_vtx)
          return
        end if
        ! ...but insert not more than one
        nins = 1
        ins_at_vtx(1) = 1
        insert(1) = 1
      else
      allocate(occ_tmp(nvtx))

      ! sum excitation (up) and deexcitation (down) arcs
      do jvtx = 1, nvtx
        if (vtxinlist(jvtx)) then
          do ivtx = 1, jvtx-1
            if (vtxinlist(ivtx)) occ_sum_d(jvtx) =
     &                           occ_sum_d(jvtx) + topo(ivtx,jvtx)
          end do
          do ivtx = jvtx+1, nvtx
            if (vtxinlist(ivtx)) occ_sum_x(jvtx) =
     &                           occ_sum_x(jvtx) + topo(ivtx,jvtx)
          end do
        end if
      end do
c dbg
c      print *,'occ_ins_x/d:',occ_ins_x,occ_ins_d
c      print *,'occ_sum_x:',occ_sum_x
c      print *,'occ_sum_d:',occ_sum_d
c dbgend

      nins = 0
      insert(1:nvtx) = 0
      ins_at_vtx(1:nvtx) = 0
      ! (first check excitation parts)
      do ivtx = nvtx, 1, -1
       if (.not.vtxinlist(ivtx)) cycle
       ! try to insert as many times as possible
       nins_old = nins-1
       do while (nins_old.lt.nins)
        nins_old = nins
        ! op. to insert must be fully connected with one op. via exc. part ...
        if (occ_overlap_p(occ_sum_x(ivtx),occ_ins_x).eq.occ_ins_x) then
          ! ... and with (several) others via deexcitation part ...
          occ_tmp(1:nvtx) = 0
          do jvtx = ivtx+1, nvtx
            if (vtxinlist(jvtx))
     &         occ_tmp(jvtx) = occ_overlap_p(topo(ivtx,jvtx),occ_ins_d)
c dbg
c            ! also fully connected with one op via deexcitation part
c            if (occ_tmp(jvtx).ne.0) exit
c dbgend
          end do
          if (sum(occ_tmp(1:nvtx)).eq.occ_ins_d) then
            ok = .true.
            do jvtx = ivtx+1, nvtx
             ok = ok.and.occ_bound_p('<=',occ_tmp(jvtx),occ_sum_d(jvtx))
            end do
            if (ok) then
              nins = nins + 1
              insert(nins) = ivtx
              ins_at_vtx(nins) = ivtx
              occ_sum_x(ivtx) = occ_sum_x(ivtx) - occ_ins_x
              do jvtx = ivtx+1, nvtx
                occ_sum_d(jvtx) = occ_sum_d(jvtx) - occ_tmp(jvtx)
              end do
            end if
          end if
        end if
       end do
      end do
      ! (now check deexcitation parts)
      do ivtx = 1, nvtx
       if (.not.vtxinlist(ivtx)) cycle
       ! try to insert as many times as possible (CAUTION: not debugged)
       nins_old = nins-1
       do while (nins_old.lt.nins)
        nins_old = nins
        ! ... or with one op. via deexc. part ...
        if (occ_overlap_p(occ_sum_d(ivtx),occ_ins_d).eq.occ_ins_d) then
          ! ... and with (several) others via excitation part.
          occ_tmp(1:nvtx) = 0
          do jvtx = 1, ivtx-1
            if (vtxinlist(jvtx))
     &         occ_tmp(jvtx) = occ_overlap_p(topo(ivtx,jvtx),occ_ins_x)
c dbg
c            ! also fully connected with one op via excitation part
c            if (occ_tmp(jvtx).ne.0) exit
c dbgend
          end do
          if (sum(occ_tmp(1:nvtx)).eq.occ_ins_x) then
            ok = .true.
            do jvtx = 1, ivtx-1
             ok = ok.and.occ_bound_p('<=',occ_tmp(jvtx),occ_sum_x(jvtx))
            end do
            if (ok) then
              nins = nins + 1
              insert(nins) = ivtx-1
              ins_at_vtx(nins) = ivtx
              occ_sum_d(ivtx) = occ_sum_d(ivtx) - occ_ins_d
              do jvtx = 1, ivtx-1
                occ_sum_x(jvtx) = occ_sum_x(jvtx) - occ_tmp(jvtx)
              end do
            end if
          end if
        end if
       end do
      end do
      if (nins.gt.nvtx) call quit(1,'contr_insert',
     &        'increase dimension of insert, ins_at_vtx')
      deallocate(occ_tmp)
      end if
c dbg
c      print *,'nins: ',nins
c      print *,'insert:     ',insert(1:nins)
c      print *,'ins_at_vtx: ',ins_at_vtx(1:nins)
c      print *,'occ_sum_x:',occ_sum_x
c      print *,'occ_sum_d:',occ_sum_d
c dbgend

      nvtx_new = nvtx + nins
      if (nins.gt.0) then
        allocate(svertex_new(nvtx_new),vtx_new(nvtx_new),
     &     topo_new(nvtx_new,nvtx_new),xlines_new(nvtx_new,njoined_res),
     &     vil_new(nvtx_new))

        svertex_new = 0
        vtx_new = 0
        topo_new = 0
        xlines_new = 0
        vil_new = .false.
        svertex_new(1:nvtx) = svertex(1:nvtx)
        vtx_new(1:nvtx) = vtx(1:nvtx)
        topo_new(1:nvtx,1:nvtx) = topo(1:nvtx,1:nvtx)
        vil_new(1:nvtx) = vtxinlist(1:nvtx)
c dbg
c          write(lulog,*) 'modified contraction (1):'
c          call prt_contr_p(lulog,svertex_new,vtx_new,topo_new,
c     &         xlines_new,nvtx_new,njoined_res)
c dbgend
      end if

      ! do insertions
      do iins = 1, nins

        ins_cur = insert(iins)
        do ivtx = nvtx+iins, ins_cur+2, -1
          svertex_new(ivtx) = svertex_new(ivtx-1)
          vtx_new(ivtx) = vtx_new(ivtx-1)
          vil_new(ivtx) = vil_new(ivtx-1)
          topo_new(ivtx,1:nvtx_new) = topo_new(ivtx-1,1:nvtx_new)
        end do
        vil_new(ins_cur+1) = .false.
        topo_new(ins_cur+1,1:nvtx_new) = 0
        do jvtx = nvtx+iins, ins_cur+2, -1
          topo_new(1:nvtx_new,jvtx) = topo_new(1:nvtx_new,jvtx-1)
        end do
        topo_new(1:nvtx_new,ins_cur+1) = 0
        do jins = iins, nins
          if (insert(jins).ge.ins_cur) insert(jins) = insert(jins) + 1
          if (ins_at_vtx(jins).gt.ins_cur)
     &           ins_at_vtx(jins) = ins_at_vtx(jins) + 1
        end do
c dbg
c        write(lulog,*) 'modified contraction for iins=',iins
c        call prt_contr_p(lulog,svertex_new,vtx_new,topo_new,
c     &       xlines_new,nvtx_new,njoined_res)
c        print *,'insert:     ',insert(1:nins)
c        print *,'ins_at_vtx: ',ins_at_vtx(1:nins)
c dbgend

        svertex_new(insert(iins)) = nvtx+iins
        idxnew = pack_vtx(idxins,iblkins,.false.)
        vtx_new(insert(iins)) = idxnew

        ! shift arcs from ins_at_vtx() to insert()
        if (insert(iins).gt.ins_at_vtx(iins)) then
          ! deexcitation part
          do ivtx = insert(iins)+1, nvtx+iins
            if (vil_new(ivtx)) then
              overlap = 
     &         occ_overlap_p(topo_new(ins_at_vtx(iins),ivtx),occ_ins_d)
              topo_new(insert(iins),ivtx) = overlap
              topo_new(ins_at_vtx(iins),insert(iins)) =
     &         topo_new(ins_at_vtx(iins),insert(iins)) + overlap
              topo_new(ins_at_vtx(iins),ivtx) =
     &         topo_new(ins_at_vtx(iins),ivtx) - overlap
c dbg
c              ! first shift must be only one
c              if (overlap.ne.0) exit
c dbgend
            end if
          end do
          ! excitation part
          do ivtx = insert(iins)+1, nvtx+iins
            if (vil_new(ivtx)) then
              overlap =      
     &         occ_overlap_p(topo_new(ivtx,ins_at_vtx(iins)),occ_ins_x)
              topo_new(ivtx,insert(iins)) = overlap
              topo_new(insert(iins),ins_at_vtx(iins)) =
     &         topo_new(insert(iins),ins_at_vtx(iins)) + overlap
              topo_new(ivtx,ins_at_vtx(iins)) =
     &         topo_new(ivtx,ins_at_vtx(iins)) - overlap
c dbg
c              ! first shift must be only one
c              if (overlap.ne.0) exit
c dbgend
            end if
          end do
        else if (insert(iins).lt.ins_at_vtx(iins)) then
          ! excitation part
          do ivtx = insert(iins)-1, 1, -1
            if (vil_new(ivtx)) then
              overlap =      
     &         occ_overlap_p(topo_new(ins_at_vtx(iins),ivtx),occ_ins_x)
              topo_new(insert(iins),ivtx) = overlap
              topo_new(ins_at_vtx(iins),insert(iins)) =
     &         topo_new(ins_at_vtx(iins),insert(iins)) + overlap
              topo_new(ins_at_vtx(iins),ivtx) =
     &         topo_new(ins_at_vtx(iins),ivtx) - overlap
c dbg
c              ! first shift must be only one
c              if (overlap.ne.0) exit
c dbgend
            end if
          end do
          ! deexcitation part
          do ivtx = insert(iins)-1, 1, -1
            if (vil_new(ivtx)) then
              overlap =
     &         occ_overlap_p(topo_new(ivtx,ins_at_vtx(iins)),occ_ins_d)
              topo_new(ivtx,insert(iins)) = overlap
              topo_new(insert(iins),ins_at_vtx(iins)) =
     &         topo_new(insert(iins),ins_at_vtx(iins)) + overlap
              topo_new(ivtx,ins_at_vtx(iins)) =
     &         topo_new(ivtx,ins_at_vtx(iins)) - overlap
c dbg
c              ! first shift must be only one
c              if (overlap.ne.0) exit
c dbgend
            end if
          end do
        else
          call quit(1,'contr_insert','this should not happen')
        end if

      end do

      if (nins.gt.0) then

        if (ntest.ge.100) then
          write(lulog,*) 'modified contraction:'
          call prt_contr_p(lulog,svertex_new,vtx_new,topo_new,
     &         xlines_new,nvtx_new,njoined_res)
        end if

        ! for insertion of more than one identical vertices,
        ! we have to divide by a permutation factor
        allocate(neqv(nvtx_new),idx_eqv(nvtx_new,nvtx_new))
        call set_eqv_map(neqv,idx_eqv,vtx_new,svertex_new,
     &                   topo_new,xlines_new,nvtx_new,njoined_res)
        ieqvfac = 1
        do iins = 1, nins
          if (neqv(insert(iins)).lt.0) cycle
          ieqvfac = ieqvfac*ifac(neqv(insert(iins)))
        end do
        deallocate(neqv,idx_eqv)
        if (ntest.ge.100)
     &       write(lulog,*) 'Divide contraction factor by ',ieqvfac
        contr%fac = contr%fac/dble(ieqvfac)

        ! set result
        call unpack_contr(contr,
     &                svertex_new,vtx_new,topo_new,xlines_new,
     &                nvtx_new,njoined_res)
        call update_svtx4contr(contr)

      end if

      deallocate(vtx,topo,xlines,svertex)
      if (nins.gt.0)
     &     deallocate(svertex_new,vtx_new,topo_new,xlines_new,
     &                vil_new)

      if (ntest.ge.100) then
        write(lulog,*) 'Generated term:'
        call prt_contr2(lulog,contr,op_info)
      end if

      return
      end
