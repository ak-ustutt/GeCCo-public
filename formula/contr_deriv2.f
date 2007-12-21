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
     &     ntest = 100
      logical, parameter ::
     &     strict = .false.

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
     &     iblk_last, idx, ider, ivtx, iarc, jvtx, jarc, il, ierr,
     &     nder_actually, idx_0, ipcr_0, ipcr_res, ipcr_mlt, ipcr_der
      integer, pointer ::
     &     iocc(:,:,:),
     &     iocc2(:,:,:)
      type(contraction_list), pointer ::
     &     cur_conder
      type(operator_array), pointer ::
     &     op_arr(:)

      integer, pointer ::
     &     nord(:), iblk(:), ivtxder(:), occ_vtx(:,:,:),
     &     topomap(:,:),eqv_map(:),
     &     neqv(:),idx_eqv(:,:),neqv_tup(:),idx_tup(:)
      
      integer, external ::
     &     iblk_occ, ielsum

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

      ! check for particle creation rank
      ! we assume that operators have a defined rank, so we
      ! check the first occupation only
      ! o  lala :::> njoined is needed as well
      print *,'njoined ignored!'
      idx_0    = contr%idx_res
      ipcr_0   = ielsum(op_arr(idx_0 )%op%ihpvca_occ(1,1,1),ngastp)-
     &           ielsum(op_arr(idx_0 )%op%ihpvca_occ(1,2,1),ngastp)
      ipcr_der = ielsum(op_arr(idxder)%op%ihpvca_occ(1,1,1),ngastp)-
     &           ielsum(op_arr(idxder)%op%ihpvca_occ(1,2,1),ngastp)
      ipcr_res = ielsum(op_arr(idxres)%op%ihpvca_occ(1,1,1),ngastp)-
     &           ielsum(op_arr(idxres)%op%ihpvca_occ(1,2,1),ngastp)
      ipcr_mlt = 0
      if (idxmlt.gt.0) then
        ipcr_mlt = ielsum(op_arr(idxmlt)%op%ihpvca_occ(1,1,1),ngastp)-
     &             ielsum(op_arr(idxmlt)%op%ihpvca_occ(1,2,1),ngastp)
      end if
      if (ipcr_0-ipcr_der+ipcr_mlt.ne.ipcr_res) then
        write(luout,*) ipcr_0,' - ',ipcr_der,' + ',ipcr_mlt,
     &       ' != ',ipcr_res
        call quit(1,'contr_deriv2',
     &     'particle creation ranks do not match')
      end if
      
      nvtx = contr%nvtx
      njoined_ori = op_arr(contr%idx_res)%op%njoined
      njoined = op_arr(idxres)%op%njoined
      allocate(topomap(nvtx,nvtx),eqv_map(nvtx),
     &         neqv(nvtx),idx_eqv(nvtx,nvtx),
     &         neqv_tup(nvtx),idx_tup(nvtx),
     &         occ_vtx(ngastp,2,nvtx+max(njoined_ori,njoined)),
     &         iocc(ngastp,2,njoined),
     &         iocc2(ngastp,2,njoined))

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

      ierr = 0
      nder_actually = 0

      cur_conder => conder
      nullify(cur_conder%prev)
      nullify(cur_conder%next)
      ! assemble derivative contractions:
      do ider = 1, nder
        if (ider.gt.1.and.ierr.eq.0) then
          allocate(cur_conder%next)
          cur_conder%next%prev => cur_conder
          cur_conder => cur_conder%next
          nullify(cur_conder%next)
        end if
        allocate(cur_conder%contr)
        call init_contr(cur_conder%contr)
        call resize_contr(cur_conder%contr,contr%nvtx,contr%narc,0)

        ! resulting operator (block: see below)
        cur_conder%contr%idx_res = idxres

        ! new pre-factor
        cur_conder%contr%fac = contr%fac*dble(nord(ider))
        
        ! copy and update vertex information
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
          iocc(1:ngastp,1:2,1) =
     &         op_arr(contr%vertex(ivtxder(ider))%idx_op)%op%
     &         ihpvca_occ(1:ngastp,1:2,
     &                    contr%vertex(ivtxder(ider))%iblk_op)
          ! exception for rank change (QUICK FIX):
c dbg-QUICK FIX:::
          if (ipcr_mlt.eq.-1.and.ipcr_der.eq.0) then
            iocc(2,1,1) = iocc(2,1,1)-1
            ! print a nasty mark to remind us of this bad line:
            print *,'QUICK FIX active'
c            call wrt_occ_n(luout,iocc,1)
          else if (ipcr_mlt.ne.0.or.ipcr_der.ne.0) then
            call quit(1,'contr_deriv',
     &           'QUICK FIX is not prepared for this')
          end if
          idx = iblk_occ(iocc,.false.,op_arr(idxmlt)%op)
          if (idx.lt.0)
     &         call quit(1,'contr_dervi2','idx<0??')
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
        call occ_contr(iocc,ierr,cur_conder%contr,occ_vtx,njoined)

c dbg-QUICK FIX2:::
        if (ierr.eq.0.and.ipcr_mlt.eq.-1.and.ipcr_der.eq.0) then
          print *,'QUICK FIX2 active:'
          call get_unconnected4vertex(iocc2,ivtxder(ider),
     &         cur_conder%contr,op_info)
c          call wrt_occ_n(luout,iocc2,1)
          if (.not.iocc_bound('>=',iocc2,.false.,
     &                             (/0,0,0,0,0,0,0,0/),.false.)) ierr=66
          if (ierr.ne.0) print *,'skipping a contribution'
        end if

        if (ierr.ne.0) then
          if (strict) then
            write(luout,*) 'derivative gives njoined = ',-ierr-1
            write(luout,*) 'expected:                  ',njoined
            call quit(1,'contr_deriv2','incompatible result operator')
          end if
        else        
          ! ... and the corresponding block in the target operator
          cur_conder%contr%iblk_res = iblk_occ(iocc,.false.,
     &                                          op_arr(idxres)%op)
c dbg
          print *,'idx, occ: ',cur_conder%contr%iblk_res
          call wrt_occ_n(luout,iocc,njoined)
c dbg          

          if (cur_conder%contr%iblk_res.le.0) then
            if (strict) then
              write(luout,*) 'resulting occupation:'
              call wrt_occ_n(luout,iocc,njoined)
              call quit(1,'contr_deriv',
     &           'result occupation not defined for operator '//
     &           trim(op_arr(idxres)%op%name))
            end if
            ierr = +1
          end if
        end if
        
        if (ierr.ne.0) then
          call dealloc_contr(cur_conder%contr)
          deallocate(cur_conder%contr)
        else
          nder_actually = nder_actually + 1
        end if

      end do ! ider

      if (nder.gt.0) deallocate(nord,iblk)

      deallocate(iocc2,iocc,occ_vtx)

      nder = nder_actually

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
