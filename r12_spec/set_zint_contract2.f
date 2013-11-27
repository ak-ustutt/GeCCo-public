      subroutine set_zint_contract2(flist,ansatz,
     &     idx_opsin,nopsin,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     Routine used to form the various contractions necessary to form
*     the CABS approximation to the Z-intermediate.
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

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
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info


      integer ::
     &     idx, idx_shape, nops ,nvtx, idx_temp, ndef, idx_prj,
     &     idx_prj2, idx_prj3
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
        write(lulog,*)'Z-Intermediate Contraction'
        write(lulog,*)'Constituent operators: '
        do idx = 1, nopsin
          write(lulog,*)trim(op_info%op_arr(idx_opsin(idx))%op%name)
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
c dbg
c      call warn('set_zint_contract2','HEAVY DEBUGGING!')
c      goto 50
c dbg

      ! Add the G^{p'q}_{km}.FF_{p'l}^{ij} terms.
      idx_prj = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     0.5d0,10,3,
     &     (/idx_shape,idx_opsin(3),idx_opsin(3),idx_shape,idx_shape,
     &     idx_opsin(4),idx_shape,idx_shape,idx_opsin(4),idx_shape/),
     &     (/        1,           2,           2,        1,        1,
     &                3,        1,        1,           3,        1/),
     &     -1,-1,
     &     (/3,7/),1,
     &     0,0,
     &     (/3,5,1,idx_prj,2,6,1,idx_prj,4,6,1,idx_prj/),3,
     &     .false.,op_info)

      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      idx_prj = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     0.5d0,10,3,
     &     (/idx_shape,idx_opsin(3),idx_opsin(3),idx_shape,idx_shape,
     &     idx_opsin(4),idx_shape,idx_shape,idx_opsin(4),idx_shape/),
     &     (/        1,           2,           2,        1,        1,
     &                3,        1,        1,           3,        1/),
     &     -1,-1,
     &     (/3,7,2,4/),2,
     &     0,0,
     &     (/2,6,1,idx_prj,4,6,1,idx_prj/),2,
     &     .false.,op_info)

      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

c dbg
c 50   call warn('HEAVY','DEBUGGING !')
c dbg

      do idx = 1, 2             ! Full
c      do idx = 1, 1 ! SA
        idx_prj  = 1
        idx_prj2 = 2*idx
        call expand_op_product2(form_pnt,idx_shape,
     &       0.5d0,10,3,
     &       (/idx_shape,idx_opsin(3),idx_opsin(3),idx_shape,idx_shape,
     &       idx_opsin(4),idx_shape,idx_shape,idx_opsin(4),idx_shape/),
     &       (/        1,           2,           2,        1,        1,
     &                  3,        1,        1,           3,        1/),
     &       -1,-1,
     &       (/3,7/),1,
     &       0,0,
     &       (/3,5,1,idx_prj,2,6,1,idx_prj2,4,6,1,idx_prj/),3,
     &       .false.,op_info)

        form_pnt => flist
        do while(associated(form_pnt%next))
          form_pnt => form_pnt%next
        enddo

        call expand_op_product2(form_pnt,idx_shape,
     &       0.5d0,10,3,
     &       (/idx_shape,idx_opsin(3),idx_opsin(3),idx_shape,idx_shape,
     &       idx_opsin(4),idx_shape,idx_shape,idx_opsin(4),idx_shape/),
     &       (/        1,           2,           2,        1,        1,
     &                  3,        1,        1,           3,        1/),
     &       -1,-1,
     &       (/3,7,2,4/),2,
     &       0,0,
     &       (/2,6,1,idx_prj2,4,6,1,idx_prj/),2,
     &       .false.,op_info)

        form_pnt => flist
        do while(associated(form_pnt%next))
          form_pnt => form_pnt%next
        enddo
      enddo
c dbg
c 50   call warn('HEAVY','DEBUGGING !')
c dbg

      do idx = 1,3 ! Full
c      do idx = 1,2 ! SA
        idx_prj = 1
        idx_prj2 = 2*idx-min(idx,2)
        call expand_op_product2(form_pnt,idx_shape,
     &      0.5d0,10,3,
     &      (/idx_shape,-idx_opsin(4),idx_shape,idx_shape,idx_opsin(3),
     &      -idx_opsin(4),idx_opsin(3),idx_shape,idx_shape,idx_shape/),
     &      (/        1,            2,        1,        1,           3,
     &                 2,           3,        1,        1,        1/),
     &      -1,-1,
     &      (/5,8/),1,
     &      0,0,
     &      (/6,7,1,idx_prj2,6,9,1,idx_prj,7,9,1,idx_prj/),3,
     &      .false.,op_info)

        ! Point to the formula and move to the end of the list.
        form_pnt => flist
        do while(associated(form_pnt%next))
          form_pnt => form_pnt%next
        enddo
      enddo

c dbg
c        ! terminates after first loop
c        call warn('HEAVY','DEBUGGING !')
c        goto 200
c dbg

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

c dbg
c      print *,'skipping all R.G.R type terms'
c      goto 200 ! skip the rest
c dbg


c dbg
c 50   call warn('HEAVY','DEBUG')
c dbg
      ! Add the F_{lm}^{pq}.G_{pm}^{nk}.R_{nq}^{ij}.
      idx_prj  = 2
      idx_prj2 = 1
      ! #1
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      idx_prj  = 2
      idx_prj2 = 1
      ! #2
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj,6,9,1,idx_prj2/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #3
      idx_prj  = 2
      idx_prj2 = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj2,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #4
      idx_prj  = 2
      idx_prj2 = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #5
      idx_prj  = 2
      idx_prj2 = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj2,6,9,1,idx_prj2/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #6
      idx_prj  = 2
      idx_prj2 = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj,6,9,1,idx_prj2/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #7
      idx_prj  = 2
      idx_prj2 = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj2,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #8
      idx_prj  = 2
      idx_prj2 = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj2,6,9,1,idx_prj2/),3,
     &     .false.,op_info)

c dbg
c      goto 200 ! SA
c dbg

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q}.G_{p"m}^{rp}.F_{rq}^{ij}.
      idx_prj  = 4
      idx_prj2 = 2
      idx_prj3 = 1
      ! #9
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj2,6,9,1,idx_prj2/),3,
     &     .false.,op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #10
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj2,6,9,1,idx_prj3/),3,
     &     .false.,op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #11
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj3,6,9,1,idx_prj2/),3,
     &     .false.,op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #12
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj3,6,9,1,idx_prj3/),3,
     &     .false.,op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q}.G_{p"m}^{rp}.R_{rq}^{ij}.
      idx_prj  = 4
      idx_prj2 = 2
      idx_prj3 = 1
      ! #13
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj2,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #14
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj3,2,9,1,idx_prj2,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #15
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj3,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! #16
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj3,2,9,1,idx_prj3,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q"}.G_{p"m}^{np}.R_{nq"}^{ij}.
      idx_prj = 4
      idx_prj2= 1
      ! #17
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj,6,9,1,idx_prj2/),3,
     &     .false.,op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      idx_prj = 4
      idx_prj2= 1
      ! #18
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      ! Extra terms from my derivation.
      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"o}.G_{p"m}^{r"p}.R_{r"o}^{ij}.
      idx_prj = 4
      idx_prj2= 1
      ! #19
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
c     &     (/2,6,1,idx_prj,2,9,1,idx_prj,5,9,1,idx_prj2/),3,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj2,6,9,1,idx_prj/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{lm}^{oq"}.G_{on}^{hk}.F_{hq"}^{ij}.
      idx_prj = 4
      idx_prj2= 1
      ! #20
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj,6,9,1,idx_prj2/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{lm}^{oq"}.G_{on}^{bk}.F_{bq"}^{ij}.
      idx_prj3 = 2
      ! #21
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj,6,9,1,idx_prj3/),3,
     &     .false.,op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{lm}^{aq"}.G_{an}^{ok}.R_{oq"}^{ij}.
      ! #22
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj3,2,9,1,idx_prj,6,9,1,idx_prj2/),3,
     &     .false.,op_info)

c dbg
c      call warn('HEAVY','DEBUG')
c      goto 200
c dbg
c dbg
c      goto 200 ! Only Coulomb terms needed.
c dbg


      ! EXCHANGE TYPE TERMS

      ! Point to the formula and move to the end of the list.
 100  do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{kl}^{p"q"}.G_{p"m}^{r"s}.F_{r"q"}^{ij}.
 110  idx_prj = 4
      call expand_op_product2(form_pnt,idx_shape,
     &     1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj,5,9,1,idx_prj/),3,
     &     .false.,op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"b}.G_{p"m}^{pr"}.R_{r"b}^{ij}.
 120  idx_prj = 4
      idx_prj2= 2
      call expand_op_product2(form_pnt,idx_shape,
     &     1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj2,5,9,1,idx_prj/),3,
     &     .false.,op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{aq"}.G_{am}^{pb}.F_{bq"}^{ij}.
 130  idx_prj = 4
      idx_prj2= 2
      call expand_op_product2(form_pnt,idx_shape,
     &     1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj,5,9,1,idx_prj2/),3,
     &     .false.,op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q"}.G_{p"m}^{pa}.R_{aq"}^{ij}.
      idx_prj = 4
      idx_prj2= 2
 140  call expand_op_product2(form_pnt,idx_shape,
     &     1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj,2,9,1,idx_prj,5,9,1,idx_prj2/),3,
     &     .false.,op_info)

c dbg
      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{aq"}.G_{am}^{pr"}.R_{r"q"}^{ij}.
      idx_prj = 4
      idx_prj2= 2
 150  call expand_op_product2(form_pnt,idx_shape,
     &     1d0,10,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_opsin(3),idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &                  3,        1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,6,1,idx_prj2,2,9,1,idx_prj,5,9,1,idx_prj/),3,
     &     .false.,op_info)
c dbg

 200  if(ntest.ge.100)then
        write(lulog,*)'Final formula: Z-Int.'
        call print_form_list(lulog,flist,op_info)
      endif

c dbg
c      stop
c dbg

      return
      end
