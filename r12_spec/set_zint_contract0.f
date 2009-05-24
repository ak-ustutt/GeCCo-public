      subroutine set_zint_contract0(flist,ansatz,
     &     idx_opsin,nopsin,max_ext_in_J,max_ext_in_K,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     Routine used to form the various contractions necessary to form
*     the CABS approximation to the Z-intermediate.
*     Here for contracted Z
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_orbinf.h'

      type(formula_item), target, intent(inout) ::
     &     flist
      integer, intent(in) ::
     &     ansatz, nopsin, idx_opsin(nopsin)
      integer, intent(in) ::
     &     max_ext_in_J, max_ext_in_K
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     idx, idx_shape, nops ,nvtx, idx_temp, ndef, next,
     &     idx_rg, idx_rr, idx_gr
      real(8) ::
     &     fac
      character ::
     &     name*(form_maxlen_label*2)

      integer, pointer ::
     &     idxarr(:), svtxarr(:), occ_def(:,:,:)
      character(4), parameter ::
     &     op_temp = 'TEMP'

      type(operator), pointer ::
     &     op
      type(formula_item), pointer ::
     &     form_pnt

      integer, external ::
     &     idx_oplist2


      if(ntest.ge.100)then
        write(luout,*)'Z-Intermediate Contraction'
        write(luout,*)'Constituent operators: '
        do idx = 1, nopsin
          write(luout,*)trim(op_info%op_arr(idx_opsin(idx))%op%name)
        enddo
      endif

      ! Get indices of input operators.
      idx_shape = idx_opsin(1)

      ! Point to the formula and move to the end of the list.
      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

c dbg
c      goto 100 ! Exchange parts only
c dbg

      ! Add the G^{p'q}_{km}.FF_{p'l}^{ij} terms.
      call expand_op_product2(form_pnt,idx_shape,
     &     0.5d0,6,3,
c     &     (/idx_shape,idx_opsin(3),idx_opsin(3),idx_shape,idx_shape,
c     &     idx_opsin(4),idx_shape,idx_shape,idx_opsin(4),idx_shape/),
c     &     (/        1,           2,           2,        1,        1,
c     &                3,        1,        1,           3,        1/),
     &     (/idx_shape,idx_opsin(3),idx_opsin(3),
     &     idx_opsin(4),idx_opsin(4),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     (/2,5/),1,
     &     (/2,4,1,0/),1,
     &     op_info)

      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the FF_{p'q}^{km}.G^{p'l}_{ij} terms.
      call expand_op_product2(form_pnt,idx_shape,
     &     0.5d0,6,3,
c     &     (/idx_shape,idx_opsin(3),idx_opsin(3),idx_shape,idx_shape,
c     &     idx_opsin(4),idx_shape,idx_shape,idx_opsin(4),idx_shape/),
c     &     (/        1,           2,           2,        1,        1,
c     &                3,        1,        1,           3,        1/),
     &     (/idx_shape,idx_opsin(4),idx_opsin(4),
     &     idx_opsin(3),idx_opsin(3),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     (/2,4/),1,
     &     (/3,4,1,0/),1,
     &     op_info)

      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      !  - R.P.G.R
      do idx_rg = 1, 4
        if (idx_rg.eq.3) cycle
        do idx_rr = 1, 4
          if (idx_rr.eq.3) cycle
          if ((idx_rg.eq.4.and.idx_rr.ne.1).or.
     &        (idx_rr.eq.4.and.idx_rg.ne.1)) cycle
          do idx_gr = 1, 4
            if (idx_gr.eq.3) cycle

            next = 0
            if (idx_rg.eq.4) next = next+1
            if (idx_rr.eq.4) next = next+1
            if (idx_gr.eq.4) next = next+1

            if (next.gt.max_ext_in_J) cycle

            call expand_op_product2(form_pnt,idx_shape,
     &           -1d0,6,4,
     &           (/idx_shape,-idx_opsin(2),idx_opsin(3),
     &           idx_opsin(3),idx_opsin(2),idx_shape/),
     &           (/        1,            2,             3,
     &           3,          4,        1/),
     &           -1,-1,
     &           0,0,
     &           0,0,
     &           (/2,3,1,idx_rg,2,5,1,idx_rr,3,5,1,idx_gr/),3,
     &           op_info)

            do while(associated(form_pnt%next))
              form_pnt => form_pnt%next
            enddo

          end do
        end do
      end do

      !  - R.G.R.R
      do idx_rg = 1, 4
        if (idx_rg.eq.3) cycle
        do idx_rr = 1, 4
          if (idx_rr.eq.3) cycle
          do idx_gr = 1, 4
            if (idx_gr.eq.3) cycle
            if ((idx_gr.eq.4.and.idx_rr.ne.1).or.
     &          (idx_rr.eq.4.and.idx_gr.ne.1)) cycle

            next = 0
            if (idx_rg.eq.4) next = next+1
            if (idx_rr.eq.4) next = next+1
            if (idx_gr.eq.4) next = next+1

            if (next.gt.max_ext_in_J) cycle

            call expand_op_product2(form_pnt,idx_shape,
     &           -1d0,6,4,
     &           (/idx_shape,-idx_opsin(2),idx_opsin(3),
     &           idx_opsin(3),idx_opsin(2),idx_shape/),
     &           (/        1,            2,             3,
     &           3,          4,        1/),
     &           -1,-1,
     &           0,0,
     &           0,0,
     &           (/2,3,1,idx_rg,2,5,1,idx_rr,3,5,1,idx_gr/),3,
     &           op_info)

            do while(associated(form_pnt%next))
              form_pnt => form_pnt%next
            enddo

          end do
        end do
      end do


      !  + R.P.G.P.R
      do idx_rg = 1, 4
        if (idx_rg.eq.3) cycle
        do idx_rr = 1, 4
          if (idx_rr.eq.3) cycle
          if ((idx_rg.eq.4.and.idx_rr.ne.1).or.
     &        (idx_rr.eq.4.and.idx_rg.ne.1)) cycle
          do idx_gr = 1, 4
            if (idx_gr.eq.3) cycle
            if ((idx_gr.eq.4.and.idx_rr.ne.1).or.
     &          (idx_rr.eq.4.and.idx_gr.ne.1)) cycle

            next = 0
            if (idx_rg.eq.4) next = next+1
            if (idx_rr.eq.4) next = next+1
            if (idx_gr.eq.4) next = next+1

            if (next.gt.max_ext_in_J) cycle

            call expand_op_product2(form_pnt,idx_shape,
     &           1d0,6,4,
     &           (/idx_shape,-idx_opsin(2),idx_opsin(3),
     &           idx_opsin(3),idx_opsin(2),idx_shape/),
     &           (/        1,            2,             3,
     &           3,          4,        1/),
     &           -1,-1,
     &           0,0,
     &           0,0,
     &           (/2,3,1,idx_rg,2,5,1,idx_rr,3,5,1,idx_gr/),3,
     &           op_info)

            do while(associated(form_pnt%next))
              form_pnt => form_pnt%next
            enddo

          end do
        end do
      end do

      ! exchange terms
      !  - R.Q.K.Q.R
      do idx_rg = 2, 4
        if (idx_rg.eq.3) cycle        
        do idx_rr = 2, 4
          if (idx_rr.eq.3) cycle
          if (idx_rg.eq.2.and.idx_rr.eq.2) cycle
          do idx_gr = 2, 4
            if (idx_gr.eq.3) cycle
            if (idx_gr.eq.2.and.idx_rr.eq.2) cycle

            next = 0
            if (idx_rg.eq.4) next = next+1
            if (idx_rr.eq.4) next = next+1
            if (idx_gr.eq.4) next = next+1

            if (next.gt.max_ext_in_K) cycle

            call expand_op_product2(form_pnt,idx_shape,
     &           1d0,6,4,
     &           (/idx_shape,-idx_opsin(2),idx_opsin(3),
     &           idx_opsin(3),idx_opsin(2),idx_shape/),
     &           (/        1,            2,             3,
     &           3,          4,        1/),
     &           -1,-1,
     &           0,0,
     &           0,0,
     &           (/2,4,1,idx_rg,2,5,1,idx_rr,3,5,1,idx_gr/),3,
     &           op_info)

            do while(associated(form_pnt%next))
              form_pnt => form_pnt%next
            enddo

          end do
        end do
      end do
      
      if (ntest.ge.100) then
        write(luout,*)'formula before summing: Z-Int.'
        call print_form_list(luout,flist,op_info)
      end if

      call sum_terms(flist,op_info)

      if (ntest.ge.100) then
        write(luout,*)'Final formula: Z-Int.'
        call print_form_list(luout,flist,op_info)
      end if
c dbg
c      stop 'test ex'
c dbg

      return
      end
