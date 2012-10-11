      integer function cmp_bcontr_spc(bc1,bc2)
      ! compare binary contractions; special version for factor out

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(binary_contr), intent(in) ::
     &    bc1, bc2

      logical, external ::
     &    list_cmp

      cmp_bcontr_spc = 0
      if (bc1%n_cnt.ne.bc2%n_cnt) return

      if (bc1%iblk_res.ne.bc2%iblk_res) return
      if (trim(bc1%label_res).ne.trim(bc2%label_res)) return

      if (.not.list_cmp(bc1%occ_cnt,bc2%occ_cnt,bc1%n_cnt*ngastp*2))
     &     return
      if (.not.list_cmp(bc1%rst_cnt,bc2%rst_cnt,
     &                  bc1%n_cnt*bc1%ngas*bc1%nspin*8))
     &     return

      if (.not.list_cmp(bc1%occ_ex1,bc2%occ_ex1,bc1%nj_op1*ngastp*2))
     &     return
      if (.not.list_cmp(bc1%rst_ex1,bc2%rst_ex1,
     &                  bc1%nj_op1*bc1%ngas*bc1%nspin*8))
     &     return

      if (.not.list_cmp(bc1%occ_ex2,bc2%occ_ex2,bc1%nj_op2*ngastp*2))
     &     return
      if (.not.list_cmp(bc1%rst_ex2,bc2%rst_ex2,
     &                  bc1%nj_op2*bc1%ngas*bc1%nspin*8))
     &     return

      if (.not.list_cmp(bc1%merge_op1,
     &                               bc2%merge_op1,size(bc1%merge_op1)))
     &     return
      if (.not.list_cmp(bc1%merge_op2,
     &                               bc2%merge_op2,size(bc1%merge_op2)))
     &     return
      if (.not.list_cmp(bc1%merge_op1op2,
     &                            bc2%merge_op1op2,size(bc1%merge_op2)))
     &     return
      if (.not.list_cmp(bc1%merge_op2op1,
     &                            bc2%merge_op2op1,size(bc1%merge_op2)))
     &     return

      ! note that short term intermediates *never* are the
      ! same, so we exclude them from the comparison
      if (bc1%iblk_op1.eq.bc2%iblk_op1.and.
     &    (bc1%tra_op1.eqv.bc2%tra_op1)) then
        cmp_bcontr_spc = 1
        if (trim(bc1%label_op1).eq.trim(bc2%label_op1)
     &      .and.bc1%label_op1(1:5).ne.'_STIN') return
      end if

      if (bc1%iblk_op2.eq.bc2%iblk_op2.and.
     &    (bc1%tra_op2.eqv.bc2%tra_op2)) then
        cmp_bcontr_spc = 2
        if (trim(bc1%label_op2).eq.trim(bc2%label_op2)
     &      .and.bc1%label_op2(1:5).ne.'_STIN') return
      end if

      cmp_bcontr_spc = 0

      return
      end
