*----------------------------------------------------------------------*
      subroutine contr_deriv2
     &     (conder,nder,contr,op_info,idxder,idxmlt,idxres)
*----------------------------------------------------------------------*
*     get derivative of contraction on contr
*     idxder is the index of the operator with respect to which the
*     derivative has to be taken
*     idxmlt is the index of the operator which the derivative is 
*     multiplied with (0 if only the derivative is taken)
*     idxres is the index of the resulting operator (0, if scalar)
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

      type(contraction_list), intent(out), target ::
     &     conder
      integer, intent(out) ::
     &     nder

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     idxder, idxmlt, idxres
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     nvtx, ntup, itup, ieqvfac, njoined_ori, njoined,
     &     iblk_last, idx, ider, ivtx, iarc, jvtx, jarc, il
      integer, pointer ::
     &     iocc(:,:,:)
      type(contraction_list), pointer ::
     &     cur_conder
      type(operator_array), pointer ::
     &     op_arr(:)

      integer, pointer ::
     &     nord(:), iblk(:), ivtxder(:), occ_vtx(:,:,:),
     &     topomap(:,:),eqv_map(:),
     &     neqv(:),idx_eqv(:,:),neqv_tup(:),idx_tup(:)
      
      integer, external ::
     &     iblk_occ

      if (ntest.ge.100) then
        write(luout,*) '====================='
        write(luout,*) ' contr_deriv at work'
        write(luout,*) '====================='
        write(luout,*) ' Contraction on input:'
        call prt_contr2(luout,contr,op_info)
        write(luout,*) 'idxder = ',idxder
        write(luout,*) 'idxmlt = ',idxmlt
        write(luout,*) 'idxres = ',idxres
      end if

      op_arr => op_info%op_arr

      nvtx = contr%nvtx
      njoined_ori = op_arr(contr%idx_res)%op%njoined
      njoined = op_arr(idxres)%op%njoined
      allocate(topomap(nvtx,nvtx),eqv_map(nvtx),
     &         neqv(nvtx),idx_eqv(nvtx,nvtx),
     &         neqv_tup(nvtx),idx_tup(nvtx),
     &         occ_vtx(ngastp,2,nvtx+max(njoined_ori,njoined)),
     &         iocc(ngastp,2,njoined))

      call occvtx4contr(1,occ_vtx,contr,op_info)

      call topomap4contr(2,topomap,eqv_map,neqv,idx_eqv,contr,occ_vtx)

      call eqvfac4contr(ieqvfac,neqv_tup,idx_tup,ntup,
     &                  topomap,eqv_map,neqv,idx_eqv,contr)
      nder = 0
      do itup = 1, ntup
        idx = idx_tup(itup)
        if (contr%vertex(idx)%idx_op.eq.idxder) then
          nder = nder+1
        end if
      end do

      ! order for each block
      if (nder.gt.0) allocate(nord(nder),iblk(nder),ivtxder(nder))
      ider = 0
      do itup = 1, ntup
        idx = idx_tup(itup)
        if (contr%vertex(idx)%idx_op.eq.idxder) then
          if (contr%joined(0,contr%svertex(idx)).gt.1)
     &         call quit(1,'contr_deriv',
     &         'cannot take derivative wrt. non-primitive vertices')
          ider = ider+1
          iblk(ider) = contr%vertex(idx)%iblk_op
          nord(ider) = neqv_tup(itup)
          ivtxder(ider) = idx
        end if
      end do

      deallocate(topomap)
      deallocate(eqv_map)
      deallocate(neqv)
      deallocate(idx_eqv)
      deallocate(neqv_tup)
      deallocate(idx_tup)

      cur_conder => conder
      nullify(cur_conder%prev)
      nullify(cur_conder%next)
      ! assemble derivative contractions:
      do ider = 1, nder
        if (ider.gt.1) then
          allocate(cur_conder%next)
          cur_conder%next%prev => cur_conder
          cur_conder => cur_conder%next
          nullify(cur_conder%next)
        end if
        allocate(cur_conder%contr)
        call init_contr(cur_conder%contr)
        call resize_contr(cur_conder%contr,contr%nvtx,contr%narc,0)
c        allocate(cur_conder%contr%vertex(contr%nvtx))
c        allocate(cur_conder%contr%arc(contr%narc))
c        ! needed for deallocation routine:
c        cur_conder%contr%mxvtx=contr%nvtx
c        cur_conder%contr%mxarc=contr%narc
c        cur_conder%contr%mxfac=0

        ! resulting operator (block: see below)
        cur_conder%contr%idx_res = idxres
c        if (idxmlt.ne.0) then
c          cur_conder%contr%iblk_res = contr%iblk_res
c        else
c          if (contr%idx_res.ne.0) then
c            iocc = iocc_add(1,op_arr(contr%idx_res)%op%
c     &                        ihpvca_occ(1,1,contr%iblk_res),.false.,
c     &                  1,op_arr(idxder)%op%ihpvca_occ(1,1,iblk(ider)),
c     &                  .not.op_arr(idxder)%op%dagger)
c          else
c            if (.not.op_arr(idxder)%op%dagger) then
c              iocc = iocc_dagger(op_arr(idxder)%op%
c     &             ihpvca_occ(1,1,iblk(ider)))
c            else
c              iocc(1:ngastp,1:2) =
c     &             op_arr(idxder)%op%ihpvca_occ(1:ngastp,1:2,iblk(ider))
c            end if
c          end if
c          idx = iblk_occ(iocc,.false.,op_arr(idxres)%op)
c          if (idx.le.0) then
c            write(luout,*) 'occupation not found for operator ',idxres
c            call wrt_occ(luout,iocc)
c            call quit(0,'contr_deriv','occupation not found') 
c          end if
c          cur_conder%contr%iblk_res = idx
c        end if

        ! new pre-factor
        cur_conder%contr%fac = contr%fac*dble(nord(ider))
        
        ! copy and update vertex information
c        ivtxder = 0 
        jvtx = 0
        do ivtx = 1, contr%nvtx
          ! vertex of derivative target ?
          if (ivtx.eq.ivtxder(ider).and.
     &        contr%vertex(ivtx)%idx_op.eq.idxder.and.
     &        contr%vertex(ivtx)%iblk_op.eq.iblk(ider)) then
            ! remember that vertex and skip
            if (idxmlt.ne.0) then
              jvtx = jvtx+1
              ! copy super-vertex info ...
              cur_conder%contr%svertex(jvtx) = contr%svertex(ivtx)
            end if
          else
            jvtx = jvtx + 1 ! counter for copy target
            ! just copy information
            cur_conder%contr%vertex(jvtx)%idx_op =
     &                                        contr%vertex(ivtx)%idx_op
            cur_conder%contr%vertex(jvtx)%iblk_op =
     &                                        contr%vertex(ivtx)%iblk_op
            cur_conder%contr%svertex(jvtx) = contr%svertex(ivtx)
          end if
        end do
        ! add info on new operator to contract derivative with
        ! (if any)
        if (idxmlt.ne.0) then
          cur_conder%contr%vertex(ivtxder(ider))%idx_op = idxmlt
          ! get occupation of original operator at vertex
          ! and search for the corresponding block of new operator
          ! (usually it will be the same)
          idx = iblk_occ(op_arr(contr%vertex(ivtxder(ider))%idx_op)%op%
     &         ihpvca_occ(1,1,contr%vertex(ivtxder(ider))%iblk_op),
     &         .false.,
     &         op_arr(idxder)%op)
          cur_conder%contr%vertex(ivtxder(ider))%iblk_op = idx
          cur_conder%contr%nvtx = contr%nvtx
          cur_conder%contr%nsupvtx = contr%nsupvtx
        else
          cur_conder%contr%nvtx = contr%nvtx-1
          cur_conder%contr%nsupvtx = contr%nsupvtx-1
        end if

        ! rectify super-vertex info
        call update_svtx4contr(cur_conder%contr)

        ! copy and update arc information
        jarc = 0
        arc_loop: do iarc = 1, contr%narc
          ! skip, if only derivative is taken and 
          ! vertex ivtxder is involved:
          if ((contr%arc(iarc)%link(1).eq.ivtxder(ider).or.
     &         contr%arc(iarc)%link(2).eq.ivtxder(ider)).and.
     &       idxmlt.eq.0)
     &         cycle arc_loop
          jarc = jarc+1
          ! change vertex numbers
          do il = 1, 2
            if (contr%arc(iarc)%link(il).lt.ivtxder(ider).or.
     &           idxmlt.gt.0) then
              ! if vertex was replaced only or if
              ! vertices below ivtxder: unchanged
              cur_conder%contr%arc(jarc)%link(il) =
     &             contr%arc(iarc)%link(il)
            else
              ! all vertices above ivtxder: shift by -1
              cur_conder%contr%arc(jarc)%link(il) =
     &             contr%arc(iarc)%link(il)-1
            end if
          end do ! il
          cur_conder%contr%arc(jarc)%occ_cnt(1:ngastp,1:2) =
     &         contr%arc(iarc)%occ_cnt(1:ngastp,1:2)
        end do arc_loop 
        cur_conder%contr%narc = jarc

        ! still, no factorization info
        cur_conder%contr%nfac = 0

        ! get occupation of target ...
        call occvtx4contr(1,occ_vtx,cur_conder%contr,op_info)
        call occ_contr(iocc,cur_conder%contr,occ_vtx,njoined)

        ! ... and the corresponding block in the target operator
        cur_conder%contr%iblk_res = iblk_occ(iocc,.false.,
     &                                          op_arr(idxres)%op)

        if (cur_conder%contr%iblk_res.le.0) then
          write(luout,*) 'resulting occupation:'
          call wrt_occ_n(luout,iocc,njoined)
          call quit(1,'contr_deriv',
     &         'result occupation not defined for operator '//
     &         trim(op_arr(idxres)%op%name))
        end if

      end do ! ider

      if (nder.gt.0) deallocate(nord,iblk)

      deallocate(iocc,occ_vtx)

      if (ntest.ge.100) then
        write(luout,*) 'Generated derivative terms: ',nder
        cur_conder => conder
        do ider = 1, nder
          write(luout,*) 'term #',ider
          call prt_contr2(luout,cur_conder%contr,op_info)
          if (ider.lt.nder) cur_conder => cur_conder%next
        end do
      end if

      return
      end
