*----------------------------------------------------------------------*
      subroutine contr_deriv3
     &     (conder,nderiv,contr,op_info,idxder,idxmlt,idxres)
*----------------------------------------------------------------------*
*     get derivative of contraction on contr
*     idxder is the index of the operator with respect to which the
*     derivative has to be taken
*     idxmlt is the index of the operator which the derivative is 
*     multiplied with (0 if only the derivative is taken)
*     idxres is the index of the resulting operator (0, if scalar)
*     
*     new version taking advantage of topo and xline matrices
*
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
      logical, parameter ::
     &     strict = .false.

      type(contraction_list), intent(out), target ::
     &     conder
      integer, intent(out) ::
     &     nderiv

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     idxder, idxmlt, idxres
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     allowed
      integer ::
     &     nvtx, idx_0, ideriv, svtx_last, ij, ivtx, jvtx,
     &     nder_actually, nvtx_new,
     &     idx_mlt_poss, iblkmlt, iblkres,
     &     njoined_0, njoined_res, njoined_der, njoined_mlt,
     &     ipcr_0, ipcr_res, ipcr_der, ipcr_mlt
      type(contraction_list), pointer ::
     &     cur_conder
      type(operator_array), pointer ::
     &     op_arr(:)

      type(operator), pointer ::
     &     op_der, op_mlt, op_res, op_res0
      integer(8) ::
     &     idxnew

      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:),
     &     vtx_new(:), topo_new(:,:), xlines_new(:,:),
     &     occ_res_p(:), occ_mlt_p(:)
      integer, pointer ::
     &     svertex(:), svertex_new(:),
     &     ivtxder(:,:), norder(:), iblkder(:), occ_mlt(:,:,:),
     &     occ_res(:,:,:), neqv(:), idx_eqv(:,:)
      
      integer, external ::
     &     iblk_occ, iblk_corresp, rank_occ
      integer(8), external ::
     &     pack_vtx
      

      if (ntest.ge.100) then
        write(luout,*) '====================='
        write(luout,*) ' contr_deriv3 at work'
        write(luout,*) '====================='
        write(luout,*) ' Contraction on input:'
        call prt_contr2(luout,contr,op_info)
        write(luout,*) 'idxder = ',idxder
        write(luout,*) 'idxmlt = ',idxmlt
        write(luout,*) 'idxres = ',idxres
      end if

      op_arr => op_info%op_arr

      idx_0 = contr%idx_res
      op_der => op_arr(abs(idxder))%op
      op_res => op_arr(idxres)%op
      op_res0 => op_arr(idx_0)%op
      if (idxmlt.gt.0) op_mlt => op_arr(idxmlt)%op

      njoined_0 = op_res0%njoined
      njoined_res = op_res%njoined
      njoined_der = op_der%njoined
      if (idxmlt.gt.0) njoined_mlt = op_mlt%njoined

      ipcr_0   = rank_occ('C-A',
     &     op_res0%ihpvca_occ,njoined_0)
      ipcr_der   = rank_occ('C-A',
     &     op_der%ihpvca_occ(1,1,1),njoined_der)
      if (idxder.lt.0) ipcr_der = -ipcr_der ! allow der. wrt. daggered op.
      ipcr_res   = rank_occ('C-A',
     &     op_res%ihpvca_occ(1,1,1),njoined_res)
      ipcr_mlt = 0
      if (idxmlt.gt.0) then
        occ_mlt => op_mlt%ihpvca_occ
        ipcr_mlt   = rank_occ('C-A',
     &       occ_mlt,njoined_mlt)
      end if

      if (ntest.ge.100) then
        write(luout,'(1x,a10,2i4)')
     &       trim(op_res0%name),njoined_0,ipcr_0
        write(luout,'(1x,a10,2i4)')
     &       trim(op_der%name),njoined_der,ipcr_der
        if (idxmlt.gt.0)
     &    write(luout,'(1x,a10,2i4)')
     &       trim(op_mlt%name),njoined_mlt,ipcr_mlt
        write(luout,'(1x,a10,2i4)')
     &       trim(op_res%name),njoined_res,ipcr_res
      end if

      if (idxmlt.gt.0.and.njoined_der.ne.njoined_mlt)
     &     call quit(1,'contr_deriv3',
     &     'mismatch of derivative and multiplication '//
     &     'operator (njoined)')
      if (idxmlt.eq.0.and.njoined_0+njoined_der.lt.njoined_res)
     &     call quit(1,'contr_deriv3',
     &     'mismatch of derivative and initial/final result '//
     &     'operator (njoined)')

      if (ipcr_0-ipcr_der+ipcr_mlt.ne.ipcr_res) then
        write(luout,*) ipcr_0,' - ',ipcr_der,' + ',ipcr_mlt,
     &       ' != ',ipcr_res
        call quit(1,'contr_deriv3',
     &     'particle creation ranks do not match')
      end if

      nvtx = contr%nvtx
      allocate(vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,njoined_0))
      allocate(svertex(nvtx),neqv(nvtx),idx_eqv(nvtx,nvtx))

      call pack_contr(svertex,vtx,topo,xlines,contr,njoined_0)

      if (ntest.ge.100) then
        write(luout,*) 'contraction in topo form:'
        call prt_contr_p(luout,svertex,vtx,topo,xlines,nvtx,njoined_0)
      end if

      call set_eqv_map(neqv,idx_eqv,
     &                 vtx,svertex,topo,xlines,nvtx,njoined_0)

      ! count number of appearances of operator in contraction
      ! we rely on the fact that the super vertices are numbered
      ! in increasing sequence, such that each new super vertex 
      ! has a larger index than the previous one
      nderiv = 0
      svtx_last = -1
      do ivtx = 1, nvtx
        if (neqv(ivtx).lt.0) cycle
        if (contr%vertex(ivtx)%idx_op.eq.abs(idxder).and.
     &     (contr%vertex(ivtx)%dagger.eqv.(idxder.lt.0)).and.
     &      svertex(ivtx).gt.svtx_last) then
          nderiv = nderiv+1
          svtx_last = svertex(ivtx)
        end if
      end do

      ! order for each block
      if (nderiv.gt.0)
     &     allocate(norder(nderiv),iblkder(nderiv),
     &     ivtxder(njoined_der,nderiv))
      ideriv = 0
      svtx_last = -1
      do ivtx = 1, nvtx
        if (neqv(ivtx).lt.0) cycle
        if (contr%vertex(ivtx)%idx_op.eq.abs(idxder).and.
     &     (contr%vertex(ivtx)%dagger.eqv.(idxder.lt.0)).and.
     &      svertex(ivtx).gt.svtx_last) then
          ideriv = ideriv+1
          iblkder(ideriv) = contr%vertex(ivtx)%iblk_op
          norder(ideriv) = neqv(ivtx)
          ! collect all vertices (if njoined_der > 1)
          ivtxder(1,ideriv) = ivtx
          if (njoined_der.gt.1) then
            ij = 1
            svtx_last = svertex(ivtx)
            do jvtx = ivtx+1, nvtx
              if (svertex(jvtx).eq.svtx_last) then
                ij = ij+1
                ivtxder(ij,ideriv) = jvtx
              end if
            end do
          end if
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'initially: nderiv = ',nderiv
        if (nderiv.gt.0) then
          write(luout,*) 'ivtxder: ',ivtxder
          write(luout,*) 'iblkder: ',iblkder
          write(luout,*) 'norder:  ',norder
        end if
      end if

      nder_actually = 0

      nvtx_new = nvtx
      if (idxmlt.le.0) nvtx_new = nvtx - njoined_der
      if (nderiv.gt.0)
     &     allocate(svertex_new(nvtx_new),vtx_new(nvtx_new),
     &     topo_new(nvtx_new,nvtx_new),xlines_new(nvtx_new,njoined_res),
     &     occ_res_p(max(1,njoined_res)),
     &     occ_res(ngastp,2,max(1,njoined_res)) )
      if (nderiv.gt.0.and.idxmlt.gt.0)
     &     allocate(occ_mlt_p(njoined_mlt))

      cur_conder => conder
      nullify(cur_conder%prev)
      nullify(cur_conder%next)
      ! assemble derivative contractions:
      do ideriv = 1, nderiv

        idx_mlt_poss = 1
        do while(idx_mlt_poss.gt.0)

          if (idxmlt.gt.0) then
            ! have to find the corresponding block
            iblkmlt = iblk_corresp(idx_mlt_poss,
     &           iblkder(ideriv),op_der,op_mlt,idxder.lt.0)

            if (iblkmlt.eq.0) call quit(1,'contr_deriv3',
     &           'corresponding operator block not found')

            svertex_new= svertex
            vtx_new    = vtx
            topo_new   = topo
            xlines_new = xlines

            if (njoined_mlt.ne.1)
     &           call quit(1,'contr_deriv3',
     &           'it''s time to reconsider handling of iblk for nj>1')
            idxnew = pack_vtx(idxmlt,iblkmlt,.false.)

            ! replace by new operator
            call topo_rename_vtxs(svertex_new,vtx_new,
     &           nvtx+1,0,idxnew,1,
     &           ivtxder(1,ideriv),nvtx,njoined_mlt)

            ! if the new operator has not the full rank:
            ! check consequences (may only reduce rank of xlines)
            if (ipcr_der.ne.ipcr_mlt) then
              call pack_occ(occ_mlt_p,
     &             occ_mlt(1:,1:,(iblkmlt-1)*njoined_mlt+1),njoined_mlt)
              call topo_rank_change(allowed,xlines_new,topo_new,
     &             occ_mlt_p,
     &             ivtxder(1,ideriv),njoined_mlt,nvtx,njoined_res)
              if (.not.allowed) cycle
            end if

          else
            ! only one possibility anyway:
            idx_mlt_poss = -1   ! toggle loop exit
            ! remove operator and place contractions to that operator
            ! in xlines
            ! what if there are xlines already?
            call topo_remove_vtxs(
     &           svertex_new,vtx_new,topo_new,xlines_new,
     &           svertex,    vtx,    topo,    xlines,
     &           ivtxder(1,ideriv),njoined_der,
     &           nvtx,nvtx_new,njoined_0,njoined_res)

          end if

          if (ntest.ge.100) then
            write(luout,*) 'modified contraction:'
            call prt_contr_p(luout,svertex_new,vtx_new,topo_new,
     &           xlines_new,nvtx_new,njoined_res)
          end if

          ! determine block of result op
          call occ_res_from_xlines(occ_res_p,
     &         xlines_new,nvtx_new,njoined_res)
          call unpack_occ(occ_res,occ_res_p,njoined_res)

          iblkres = iblk_occ(occ_res,.false.,op_res,
     &                       op_der%blk_version((iblkder(ideriv)-1)/
     &                       njoined_der+1))

          if (ntest.ge.100) then
            write(luout,*) 'result occupation'
            write(luout,'(1x,8i8.8)') occ_res_p
            call wrt_occ_n(luout,occ_res,njoined_res)
            write(luout,*) '-> iblkres = ',iblkres
          end if


          if (iblkres.le.0) then
            if (strict)
     &           call quit(1,'contr_deriv3','block not found on op_der')
            cycle
          end if

          nder_actually = nder_actually+1
          ! convert to contr
          if (nder_actually.gt.1) then
            allocate(cur_conder%next)
            cur_conder%next%prev => cur_conder
            cur_conder => cur_conder%next
            nullify(cur_conder%next)
          end if
          allocate(cur_conder%contr)
          call init_contr(cur_conder%contr)

          ! idx and block of result
          cur_conder%contr%idx_res = idxres
          cur_conder%contr%iblk_res = iblkres
          cur_conder%contr%dagger   = contr%dagger 

          ! new pre-factor
          cur_conder%contr%fac = contr%fac*dble(norder(ideriv))

          ! set result
          call unpack_contr(cur_conder%contr,
     &                  svertex_new,vtx_new,topo_new,xlines_new,
     &                  nvtx_new,njoined_res)

        end do ! further possibilities

      end do ! ider

      deallocate(vtx,topo,xlines,svertex,neqv,idx_eqv)
      if (nderiv.gt.0)
     &     deallocate(norder,iblkder,ivtxder)
      if (nderiv.gt.0)
     &     deallocate(svertex_new,vtx_new,topo_new,xlines_new,
     &     occ_res_p,occ_res)
      if (nderiv.gt.0.and.idxmlt.gt.0)
     &     deallocate(occ_mlt_p)

      nderiv = nder_actually

      if (ntest.ge.100) then
        write(luout,*) 'Generated derivative terms: ',nderiv
        cur_conder => conder
        do ideriv = 1, nderiv
          write(luout,*) 'term #',ideriv
          call prt_contr2(luout,cur_conder%contr,op_info)
          call check_xarcs(cur_conder%contr,op_info)
          if (ideriv.lt.nderiv) cur_conder => cur_conder%next
        end do
      end if

      return
      end
