*----------------------------------------------------------------------*
      subroutine contr_deriv2
     &     (conder,nder,contr,op_arr,idxder,idxmlt,idxres)
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
      include 'def_operator.h'
      include 'def_operator_array.h'
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
      type(operator_array), intent(in) ::
     &     op_arr(*)


      integer ::
     &     iblk_last, idx, ider, ivtx, iarc, ivtxder, jvtx, jarc, il
      integer ::
     &     iocc(ngastp,2)
      type(contraction_list), pointer ::
     &     cur_conder

      integer, allocatable ::
     &     nord(:), iblk(:)
      
      integer, external ::
     &     iblk_occ

      if (ntest.ge.100) then
        write(luout,*) '====================='
        write(luout,*) ' contr_deriv at work'
        write(luout,*) '====================='
        write(luout,*) ' Contraction on input:'
        call prt_contr2(luout,contr,op_arr)
        write(luout,*) 'idxder = ',idxder
        write(luout,*) 'idxmlt = ',idxmlt
        write(luout,*) 'idxres = ',idxres
      end if

      ! count appearence of operator in contraction
      ! and appearence of different blocks of that op
      ! we assume that only derivatives of standard-ordered
      ! contractions are taken, so identical blocks come 
      ! subsequently
      iblk_last = -1
      nder = 0
      do idx = 1, contr%nvtx
        if (contr%vertex(idx)%idx_op.eq.idxder) then
          if (contr%vertex(idx)%iblk_op.ne.iblk_last) nder = nder+1
          iblk_last = contr%vertex(idx)%iblk_op
        end if
      end do

      ! order for each block
      allocate(nord(nder),iblk(nder))
      iblk_last = -1
      nord(1:nder) = 0
      ider = 0
      do idx = 1, contr%nvtx
        if (contr%vertex(idx)%idx_op.eq.idxder) then
          if (contr%vertex(idx)%iblk_op.ne.iblk_last) then
            ider = ider+1
            iblk(ider) = contr%vertex(idx)%iblk_op
          end if
          iblk_last = contr%vertex(idx)%iblk_op
          nord(ider) = nord(ider)+1
        end if
      end do

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
        allocate(cur_conder%contr%vertex(contr%nvtx))
        allocate(cur_conder%contr%arc(contr%narc))
        ! needed for deallocation routine:
        cur_conder%contr%mxvtx=contr%nvtx
        cur_conder%contr%mxarc=contr%narc
        cur_conder%contr%mxfac=0

        ! resulting operator and block:
        cur_conder%contr%idx_res = idxres
        if (idxmlt.ne.0) then
          cur_conder%contr%iblk_res = contr%iblk_res
        else
          if (contr%idx_res.ne.0) then
            iocc = iocc_add(1,op_arr(contr%idx_res)%op%
     &                        ihpvca_occ(1,1,contr%iblk_res),.false.,
     &                  1,op_arr(idxder)%op%ihpvca_occ(1,1,iblk(ider)),
     &                  .not.op_arr(idxder)%op%dagger)
          else
            if (.not.op_arr(idxder)%op%dagger) then
              iocc = iocc_dagger(op_arr(idxder)%op%
     &             ihpvca_occ(1,1,iblk(ider)))
            else
              iocc(1:ngastp,1:2) =
     &             op_arr(idxder)%op%ihpvca_occ(1:ngastp,1:2,iblk(ider))
            end if
          end if
          idx = iblk_occ(iocc,.false.,op_arr(idxres)%op)
          if (idx.le.0) then
            write(luout,*) 'occupation not found for operator ',idxres
            call wrt_occ(luout,iocc)
            call quit(0,'contr_deriv','occupation not found') 
          end if
          cur_conder%contr%iblk_res = idx
        end if

        ! new pre-factor
        cur_conder%contr%fac = contr%fac/dble(nord(ider))
        
        ! copy and update vertex information
        ivtxder = 0 
        jvtx = 0
        do ivtx = 1, contr%nvtx
          ! find first vertex of derivative target
          if (ivtxder.eq.0.and.
     &        contr%vertex(ivtx)%idx_op.eq.idxder.and.
     &        contr%vertex(ivtx)%iblk_op.eq.iblk(ider)) then
            ! remember that vertex and skip
            ivtxder = ivtx
            if (idxmlt.ne.0) jvtx = jvtx+1
          else
            jvtx = jvtx + 1 ! counter for copy target
            ! just copy information
            cur_conder%contr%vertex(jvtx)%idx_op =
     &                                        contr%vertex(ivtx)%idx_op
            cur_conder%contr%vertex(jvtx)%iblk_op =
     &                                        contr%vertex(ivtx)%iblk_op
          end if
        end do
        ! add info on new operator to contract derivative with
        ! (if any)
        if (idxmlt.ne.0) then
          cur_conder%contr%vertex(ivtxder)%idx_op = idxmlt
          ! get occupation of original operator at vertex
          ! and search for the corresponding block of new operator
          ! (usually it will be the same)
          idx = iblk_occ(op_arr(contr%vertex(ivtxder)%idx_op)%op%
     &         ihpvca_occ(1,1,contr%vertex(ivtxder)%iblk_op),.false.,
     &         op_arr(idxder)%op)
          cur_conder%contr%vertex(ivtxder)%iblk_op = idx
          cur_conder%contr%nvtx = contr%nvtx
        else
          cur_conder%contr%nvtx = contr%nvtx-1
        end if

        ! copy and update arc information
        jarc = 0
        arc_loop: do iarc = 1, contr%narc
          ! skip, if only derivative is taken and 
          ! vertex ivtxder is involved:
          if ((contr%arc(iarc)%link(1).eq.ivtxder.or.
     &         contr%arc(iarc)%link(2).eq.ivtxder).and.idxmlt.eq.0)
     &         cycle arc_loop
          jarc = jarc+1
          ! change vertex numbers
          do il = 1, 2
            if (contr%arc(iarc)%link(il).lt.ivtxder.or.idxmlt.gt.0) then
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

      end do ! ider

      deallocate(nord,iblk)

      if (ntest.ge.100) then
        write(luout,*) 'Generated derivative terms: ',nder
        cur_conder => conder
        do ider = 1, nder
          write(luout,*) 'term #',ider
          call prt_contr2(luout,cur_conder%contr,op_arr)
          if (ider.lt.nder) cur_conder => cur_conder%next
        end do
      end if

      return
      end
