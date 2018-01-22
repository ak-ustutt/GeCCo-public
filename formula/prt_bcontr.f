*----------------------------------------------------------------------*
      subroutine prt_bcontr(lulog,bcontr)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lulog
      type(binary_contr) ::
     &     bcontr

      integer ::
     &     ij

      write(lulog,*) bcontr%fact
      if (.not.bcontr%tra_op1) write(lulog,'(x,a,i4)')
     &     trim(bcontr%label_op1),bcontr%iblk_op1
      if (     bcontr%tra_op1) write(lulog,'(x,a,i4)')
     &     trim(bcontr%label_op1)//'^+',bcontr%iblk_op1

      do ij = 1, bcontr%nj_op1
        call wrt_occ_rstr(lulog,ij,
     &       bcontr%occ_op1(1:,1:,ij),
     &       bcontr%rst_op1(1:,1:,1:,1:,1:,ij),bcontr%ngas,bcontr%nspin)
      end do

      if (bcontr%n_operands.eq.2) then
        if (.not.bcontr%tra_op2) write(lulog,'(x,a,i4)')
     &       trim(bcontr%label_op2),bcontr%iblk_op2
        if (     bcontr%tra_op2) write(lulog,'(x,a,i4)')
     &       trim(bcontr%label_op2)//'^+',bcontr%iblk_op2
        do ij = 1, bcontr%nj_op2
          call wrt_occ_rstr(lulog,ij,
     &         bcontr%occ_op2(1:,1:,ij),
     &         bcontr%rst_op2(1:,1:,1:,1:,1:,ij),
     &                                         bcontr%ngas,bcontr%nspin)
        end do
      end if

      if (bcontr%n_cnt.gt.0) then
        write(lulog,*) 'contracted by'
c        call wrt_occ_n(lulog,bcontr%occ_cnt,bcontr%n_cnt)
        do ij = 1, bcontr%n_cnt
          write(lulog,*) "Hello", bcontr%occ_cnt(1:,1:,ij)
          call wrt_occ_rstr(lulog,ij,
     &         bcontr%occ_cnt(1:,1:,ij),
     &         bcontr%rst_cnt(1:,1:,1:,1:,1:,ij),
     &         bcontr%ngas,bcontr%nspin)
        end do
        write(lulog,*) 'external indices of first operator'
        do ij = 1, bcontr%nj_op1
          call wrt_occ_rstr(lulog,ij,
     &         bcontr%occ_ex1(1:,1:,ij),
     &         bcontr%rst_ex1(1:,1:,1:,1:,1:,ij),
     &         bcontr%ngas,bcontr%nspin)
        end do
        if (bcontr%n_operands.eq.2) then
          write(lulog,*) 'external indices of second operator'
          do ij = 1, bcontr%nj_op2
            call wrt_occ_rstr(lulog,ij,
     &         bcontr%occ_ex2(1:,1:,ij),
     &         bcontr%rst_ex2(1:,1:,1:,1:,1:,ij),
     &         bcontr%ngas,bcontr%nspin)
          end do
        end if
        write(lulog,'(2x,a,6i3,/12x,6i3)') 'merge_op1 ',bcontr%merge_op1
        if (bcontr%n_operands.eq.2) 
     &     write(lulog,'(2x,a,6i3,/12x,6i3)') 
     &                                     'merge_op2 ',bcontr%merge_op2
        write(lulog,'(2x,a,6i3,/12x,6i3)')
     &                               'merge_op1op2 ',bcontr%merge_op1op2
        write(lulog,'(2x,a,6i3,/12x,6i3)') 
     &                               'merge_op2op1 ',bcontr%merge_op2op1
      end if

      if (.not.bcontr%tra_res) write(lulog,'(x,"==> ",a,i4)')
     &     trim(bcontr%label_res),bcontr%iblk_res
      if (     bcontr%tra_res) write(lulog,'(x,"==> ",a,i4)')
     &     trim(bcontr%label_res)//'^+',bcontr%iblk_res

      do ij = 1, bcontr%nj_res
        call wrt_occ_rstr(lulog,ij,
     &       bcontr%occ_res(1:,1:,ij),
     &       bcontr%rst_res(1:,1:,1:,1:,1:,ij),bcontr%ngas,bcontr%nspin)
      end do      

      return
      end
