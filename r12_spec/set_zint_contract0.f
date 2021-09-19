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
      integer, intent(in) ::
     &     max_ext_in_J, max_ext_in_K
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     idx, idx_shape, nops ,nvtx, idx_temp, ndef, next,
     &     idx_rg, idx_rr, idx_gr, icase, nterms,
     &     prj(8)
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
        write(lulog,*) 'max EXTR in J: ',max_ext_in_J
        write(lulog,*) 'max EXTR in K: ',max_ext_in_K
      endif

      ! Get indices of input operators.
      idx_shape = idx_opsin(1)

      ! Point to the formula and move to the end of the list.
      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the G^{p'q}_{km}.FF_{p'l}^{ij} terms.
      ! due to the formal shifting of commuting operators, the
      ! terms are somewhat tricky and thus "hand taylored" for
      ! each kind of Z

      ! special for Z0:
      call set_ffg_for_z0()
        
      ! settings for Z1
      call set_ffg_for_z1()

      ! settings for Z2
      call set_ffg_for_z2()
      call set_ffg_for_z2_2()

      ! remove the projected terms
      !  - R.P.G.R - R.G.P.R + R.P.G.P.R
c dbg
c      if (ntest.ge.500) call warn('zint0','debug!!!')
c      if (ntest.ge.500) goto 999
c dbg
      call set_fjf_for_zn()

      ! and finally: the exchange terms
      !  - R.Q.K.Q.R
      if (max_ext_in_K.gt.0)
     &     call set_fkf_for_zn()
c dbg
c 999  print *,'and jumped'
c dbg

      if (ntest.ge.100) then
        write(lulog,*)'formula before summing: Z-Int.'
        call print_form_list(lulog,flist,op_info)
      end if

      ! we have set the fjf terms in by straight-forward
      ! expansion of the projector terms; some of them cancel
      ! so we have to ....
      call sum_terms(flist,nterms,op_info)

      if (ntest.ge.100) then
        write(lulog,*)'Final formula: Z-Int.'
        call print_form_list(lulog,flist,op_info)
c dbg
c        if (ntest.ge.500) stop 'testing'
c dbg
      end if

      return

      ! some embedded subroutines follow
      contains

      subroutine set_ffg_for_z0()

      implicit none

        call expand_op_product3(form_pnt,idx_shape,
     &    -0.5d0,6,3,
     &     (/idx_shape,idx_opsin(3),idx_opsin(3),
     &     idx_opsin(4),idx_opsin(4),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,       ','6,,,       ',
     &       '3,,H,H     ','2,,H,[HPX],',
     &       '4,,H[HPX], ','5,,,HH     ',
     &       '2,4,,[HPX] ','3,4,,H     '/),8,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo
        call expand_op_product3(form_pnt,idx_shape,
     &     0.5d0,6,3,
     &     (/idx_shape,idx_opsin(3),idx_opsin(3),
     &     idx_opsin(4),idx_opsin(4),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
c     &     (/'2,,H,H','3,,H,[HPX],','4,,H[HPX],','5,,,HH',
c     &       '3,4,,[HPX]','2,3,,H'/),6,
c  the above gives problems with trace_op, so we use instead:
     &     (/'1,,,       ','6,,,       ',
     &       '3,,H,H     ','2,,H,[HPX],',
     &       '4,,H[HPX], ','5,,,HH     ',
     &       '2,4,,[HPX] ','2,3,H,     '/),8,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo

        ! term group #3
        call expand_op_product3(form_pnt,idx_shape,
     &    -0.5d0,6,3,
     &     (/idx_shape,-idx_opsin(4),-idx_opsin(4),
     &    -idx_opsin(3),-idx_opsin(3),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,       ','6,,,       ',
c     &       '2,,HH,','3,,,H[HPX]','5,,H,H','4,,[HPX],H',
c     &       '3,4,,[HPX]','3,5,,H'/),8,
     &       '2,,HH,     ','3,,,H[HPX] ',
     &       '4,,H,H     ','5,,[HPX],H ',
     &       '3,5,,[HPX] ','3,4,,H     '/),8,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo
        ! term group #4
        call expand_op_product3(form_pnt,idx_shape,
     &     0.5d0,6,3,
     &     (/idx_shape,-idx_opsin(4),-idx_opsin(4),
     &    -idx_opsin(3),-idx_opsin(3),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,      ','6,,,      ',
     &       '2,,HH,    ','3,,,H[HPX]',
     &       '4,,H,H    ','5,,[HPX],H',
     &       '3,5,,[HPX]','2,3,H,    '/),8,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo

      end subroutine

      subroutine set_ffg_for_z1()

      implicit none

c dbg
c      print *,'call for Z #1'
c dbg
      ! G.FF terms
      call expand_op_product3(form_pnt,idx_shape,
     &     0.5d0,6,3,
     &     (/idx_shape,idx_opsin(3),idx_opsin(3),
     &     idx_opsin(4),idx_opsin(4),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,       ','6,,P,H     ',
c     &       '3,,H,P','2,,H,[HPX],','4,,H[HPX],','5,,,H[HP]',
     &       '3,,H,P     ','2,,H,[HPX],',
     &       '4,,H[HPX], ','5,,,HH     ',
     &       '2,4,,[HPX] '/),7,
     &     .false.,op_info)
      idx = 0
      do while(associated(form_pnt%next))
        idx = idx+1
        form_pnt => form_pnt%next
      enddo
      
      ! FF.G terms
      call expand_op_product3(form_pnt,idx_shape,
     &     0.5d0,6,3,
     &     (/idx_shape,-idx_opsin(4),-idx_opsin(4),
     &    -idx_opsin(3),-idx_opsin(3),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,       ','6,,P,H     ',
     &       '2,,H[HP],  ','3,,,H[HPX] ',
     &       '4,,H,P     ','5,,[HPX],H ',
     &       '3,5,,[HPX] '/),7,
     &     .false.,op_info)
      idx = 0
      do while(associated(form_pnt%next))
        idx = idx+1
        form_pnt => form_pnt%next
      enddo
      
      ! one FF.G terms is missing: here, due to the
      ! shifting of indices a normal ordered contraction
      ! becomes anti-normal order - so we have to enforce
      ! it (and reverse the sign!)
      call expand_op_product3(form_pnt,idx_shape,
     &     -0.5d0,6,3,
     &     (/idx_shape,-idx_opsin(4),-idx_opsin(4),
     &    -idx_opsin(3),-idx_opsin(3),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,       ','6,,P,H     ',
     &       '2,,H[HP],  ','3,,,H[HPX] ',
     &       '4,,H,P     ','5,,[HPX],H ',
     &       '3,5,,[HPX] ','3,4,,H     '/),8,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo

      end subroutine

      subroutine set_ffg_for_z2()

      implicit none

c dbg
c      print *,'call for Z #2'
c dbg
      ! G.FF terms
      call expand_op_product3(form_pnt,idx_shape,
     &     0.5d0,6,3,
     &     (/idx_shape,idx_opsin(3),idx_opsin(3),
     &     idx_opsin(4),idx_opsin(4),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,       ','6,,PP,HH   ',
     &       '3,,H,P     ','2,,H,[HPX],',
     &       '4,,H[HPX], ','5,,,H[HP]  ',
     &       '2,4,,[HPX] '/),7,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo
        
        ! FF.G terms
        call expand_op_product3(form_pnt,idx_shape,
     &     0.5d0,6,3,
     &     (/idx_shape,-idx_opsin(4),-idx_opsin(4),  ! with dagger??
     &     idx_opsin(3),idx_opsin(3),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,         ','6,,PP,HH     ',
     &       '2,,HH,       ','3,,,[HP][HPX]',
     &       '5,,H,P       ','4,,[HPX],[HP]',
     &       '3,4,,[HPX]   '/),7,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo
        ! enforce the self-contraction terms
        call expand_op_product3(form_pnt,idx_shape,
     &    -0.5d0,6,3, ! "-": anti-normal-order hole contraction
     &     (/idx_shape,-idx_opsin(4),-idx_opsin(4),
     &     idx_opsin(3),idx_opsin(3),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,         ','6,,PP,HH     ',
     &       '2,,HH,       ','3,,,[HP][HPX]',
     &       '5,,H,P       ','4,,[HPX],[HP]',
     &       '3,4,,[HPX]   ','4,5,,H       '/),8,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo
        
        call expand_op_product3(form_pnt,idx_shape,
     &     -0.5d0,6,3,
     &     (/idx_shape,idx_opsin(4),idx_opsin(4),
     &     idx_opsin(3),idx_opsin(3),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,       ','6,,PP,HH   ',
     &       '2,,HH,     ','3,,,H[HPX] ',
     &       '5,,H,P     ','4,,[HPX],P ',
     &       '3,4,,[HPX] ','3,5,,H     '/),8,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo

      end subroutine

      subroutine set_ffg_for_z2_2()

      implicit none

c dbg
c      print *,'call for Z #2-2'
c dbg
      ! G.FF terms
      call expand_op_product3(form_pnt,idx_shape,
     &     0.5d0,6,3,
     &     (/idx_shape,idx_opsin(3),idx_opsin(3),
     &     idx_opsin(4),idx_opsin(4),idx_shape/),
     &     (/        1,           2,           2,
     &                3,        3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/'1,,,P      ','6,,PP,H    ',
     &       '2,4,,[HPX] ','2,5,H,     '/),4,
     &     .false.,op_info)
        idx = 0
        do while(associated(form_pnt%next))
          idx = idx+1
          form_pnt => form_pnt%next
        enddo
        
        ! FF.G terms
!        call expand_op_product3(form_pnt,idx_shape,
!     &     0.5d0,6,3,
!     &     (/idx_shape,-idx_opsin(4),-idx_opsin(4),  ! with dagger??
!     &     idx_opsin(3),idx_opsin(3),idx_shape/),
!     &     (/        1,           2,           2,
!     &                3,        3,        1/),
!     &     -1,-1,
!     &     0,0,
!     &     0,0,
!     &     (/'1,,,P','6,,PP,H',
!     &       '2,,HH,','3,,,[HP][HPX]','5,,H,P','4,,[HPX],[HP]',
!     &       '3,4,,[HPX]'/),7,
!     &     .false.,op_info)
!        idx = 0
!        do while(associated(form_pnt%next))
!          idx = idx+1
!          form_pnt => form_pnt%next
!        enddo
!        ! enforce the self-contraction terms
!        call expand_op_product3(form_pnt,idx_shape,
!     &    -0.5d0,6,3, ! "-": anti-normal-order hole contraction
!     &     (/idx_shape,-idx_opsin(4),-idx_opsin(4),
!     &     idx_opsin(3),idx_opsin(3),idx_shape/),
!     &     (/        1,           2,           2,
!     &                3,        3,        1/),
!     &     -1,-1,
!     &     0,0,
!     &     0,0,
!     &     (/'1,,,P','6,,PP,H',
!     &       '2,,HH,','3,,,[HP][HPX]','5,,H,P','4,,[HPX],[HP]',
!     &       '3,4,,[HPX]','4,5,,H'/),8,
!     &     .false.,op_info)
!        idx = 0
!        do while(associated(form_pnt%next))
!          idx = idx+1
!          form_pnt => form_pnt%next
!        enddo
!        
!        call expand_op_product3(form_pnt,idx_shape,
!     &     -0.5d0,6,3,
!     &     (/idx_shape,idx_opsin(4),idx_opsin(4),
!     &     idx_opsin(3),idx_opsin(3),idx_shape/),
!     &     (/        1,           2,           2,
!     &                3,        3,        1/),
!     &     -1,-1,
!     &     0,0,
!     &     0,0,
!     &     (/'1,,,P','6,,PP,H',
!     &       '2,,HH,','3,,,H[HPX]','5,,H,P','4,,[HPX],P',
!     &       '3,4,,[HPX]','3,5,,H'/),8,
!     &     .false.,op_info)
!        idx = 0
!        do while(associated(form_pnt%next))
!          idx = idx+1
!          form_pnt => form_pnt%next
!        enddo
!
      end subroutine

      subroutine set_fjf_for_zn()

      implicit none

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
     &           (/3,4/),1,
c     &           (/2,4,1,idx_rg,2,5,1,idx_rr,4,5,1,idx_gr/),3,
     &           (/2,3,1,idx_rg,2,5,1,idx_rr,3,5,1,idx_gr/),3,
     &           .false.,op_info)

            idx = 0
            do while(associated(form_pnt%next))
              idx = idx+1
              form_pnt => form_pnt%next
            enddo
c dbg
            print *,' R.P.G.R',idx_rg,idx_rr,idx_gr,': ',idx,' terms'
c dbg            

          end do
        end do
      end do

      !  - R.G.P.R
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
     &           (/3,4/),1,
     &           (/2,3,1,idx_rg,2,5,1,idx_rr,3,5,1,idx_gr/),3,
c     &           (/2,4,1,idx_rg,2,5,1,idx_rr,4,5,1,idx_gr/),3,
     &           .false.,op_info)

            idx = 0
            do while(associated(form_pnt%next))
              idx = idx+1
              form_pnt => form_pnt%next
            enddo
c dbg
            print *,' R.G.P.R',idx_rg,idx_rr,idx_gr,': ',idx,' terms'
c dbg            

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
     &           (/3,4/),1,
     &           (/2,3,1,idx_rg,2,5,1,idx_rr,3,5,1,idx_gr/),3,
c     &           (/2,4,1,idx_rg,2,5,1,idx_rr,4,5,1,idx_gr/),3,
     &           .false.,op_info)

            idx = 0
            do while(associated(form_pnt%next))
              idx = idx+1
              form_pnt => form_pnt%next
            enddo
c dbg
            print *,' R.P.G.P.R',idx_rg,idx_rr,idx_gr,': ',idx,' terms'
c dbg            

          end do
        end do
      end do

      end subroutine

      subroutine set_fkf_for_zn()

      implicit none

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
     &           1d0,6,4,  ! must be a "+" due to phase of G_Z
     &           (/idx_shape,-idx_opsin(2),idx_opsin(3),
     &           idx_opsin(3),idx_opsin(2),idx_shape/),
     &           (/        1,            2,             3,
     &           3,          4,        1/),
     &           -1,-1,
     &           0,0,
     &           (/3,4/),1,
c     &           (/2,4,1,idx_rg,2,5,1,idx_rr,3,5,1,idx_gr/),3,
     &           (/2,3,1,idx_rg,2,5,1,idx_rr,4,5,1,idx_gr/),3,
     &           .false.,op_info)

            do while(associated(form_pnt%next))
              form_pnt => form_pnt%next
            enddo

          end do
        end do
      end do

      end subroutine

      end
