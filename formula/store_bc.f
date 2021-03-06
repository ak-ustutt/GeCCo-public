*----------------------------------------------------------------------*
      subroutine store_bc(fl_item,
     &       fact,fact_itf,
     &       label12,label1,label2,
     &       iblk12,iblk1,iblk2,
     &       tra12,tra1,tra2,
     &       nj12,nj1,nj2,
     &       iocc12,iocc1,iocc2,
     &       irst12,irst1,irst2,
     &       ex1,ex2,cnt,
     &       irstex1,irstex2,irstcnt,nj_cnt,
     &       merge1,merge2,
     &       merge12,merge21,
     &       itf_index_info,
     &       orb_info)
*----------------------------------------------------------------------*
*     store info on unary or binary contraction on fl_item
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(inout) ::
     &     fl_item
      type(orbinf), intent(in) ::
     &     orb_info
      real(8) ::
     &     fact, fact_itf  ! corrected factor for ITF post-proc.
      character(len=*), intent(in) ::
     &     label12,label1,label2
      integer, intent(in) ::
     &     iblk12,iblk1,iblk2,nj_cnt,
     &     nj12, nj1, nj2
      logical, intent(in) ::
     &     tra12,tra1,tra2
      integer, intent(in) ::
     &     ex1(ngastp,2,nj1),
     &     ex2(ngastp,2,nj2),
     &     cnt(ngastp,2,nj_cnt),
     &     iocc12(ngastp,2,nj12),
     &     iocc1(ngastp,2,nj1),
     &     iocc2(ngastp,2,nj2),
     &     irst12(2,orb_info%ngas,2,2,orb_info%nspin,nj12),
     &     irst1(2,orb_info%ngas,2,2,orb_info%nspin,nj1),
     &     irst2(2,orb_info%ngas,2,2,orb_info%nspin,nj2),
     &     irstex1(2,orb_info%ngas,2,2,orb_info%nspin,nj1),
     &     irstex2(2,orb_info%ngas,2,2,orb_info%nspin,nj2),
     &     irstcnt(2,orb_info%ngas,2,2,orb_info%nspin,nj_cnt),
     &     merge1(*), merge2(*), merge12(*), merge21(*),
     &     itf_index_info(*)

      integer ::
     &     n_operands, lenmap, ngas, nspin
      
      integer, external ::
     &     len_merge_map
      
      if (ntest.ge.100) write(lulog,*) 'storing binary contraction'

      ngas = orb_info%ngas
      nspin = orb_info%nspin

      n_operands = 1
      if (label2(1:3).ne.'---') n_operands = 2

      fl_item%bcontr%fact = fact
      fl_item%bcontr%fact_itf = fact_itf
      fl_item%bcontr%n_operands = n_operands
      fl_item%bcontr%n_cnt     = nj_cnt
      fl_item%bcontr%ngas      = ngas
      fl_item%bcontr%nspin     = nspin
      fl_item%bcontr%label_res = trim(label12)
      fl_item%bcontr%label_op1 = trim(label1)
      fl_item%bcontr%label_op2 = trim(label2)
      fl_item%bcontr%iblk_res = iblk12
      fl_item%bcontr%iblk_op1 = iblk1
      fl_item%bcontr%iblk_op2 = iblk2
      fl_item%bcontr%tra_res = tra12
      fl_item%bcontr%tra_op1 = tra1
      fl_item%bcontr%tra_op2 = tra2
      fl_item%bcontr%nj_res = nj12
      fl_item%bcontr%nj_op1 = nj1
      fl_item%bcontr%nj_op2 = nj2

      allocate(fl_item%bcontr%occ_res(ngastp,2,nj12))
      fl_item%bcontr%occ_res = iocc12
      allocate(fl_item%bcontr%rst_res(2,ngas,2,2,nspin,nj12))
      fl_item%bcontr%rst_res = irst12

      allocate(fl_item%bcontr%occ_op1(ngastp,2,nj1))
      fl_item%bcontr%occ_op1 = iocc1
      allocate(fl_item%bcontr%rst_op1(2,ngas,2,2,nspin,nj1))
      fl_item%bcontr%rst_op1 = irst1

      if (n_operands.eq.2) then
        allocate(fl_item%bcontr%occ_op2(ngastp,2,nj2))
        fl_item%bcontr%occ_op2 = iocc2
        allocate(fl_item%bcontr%rst_op2(2,ngas,2,2,nspin,nj2))
        fl_item%bcontr%rst_op2 = irst2
      end if

      if (nj_cnt.gt.0) then
        allocate(fl_item%bcontr%occ_cnt(ngastp,2,nj_cnt))
        fl_item%bcontr%occ_cnt = cnt
        allocate(fl_item%bcontr%rst_cnt(2,ngas,2,2,nspin,nj_cnt))
        fl_item%bcontr%rst_cnt = irstcnt
        allocate(fl_item%bcontr%occ_ex1(ngastp,2,nj1))
        fl_item%bcontr%occ_ex1 = ex1
        allocate(fl_item%bcontr%rst_ex1(2,ngas,2,2,nspin,nj1))
        fl_item%bcontr%rst_ex1 = irstex1

        lenmap = len_merge_map(merge1,nj1)
        allocate(fl_item%bcontr%merge_op1(lenmap))
        fl_item%bcontr%merge_op1 = merge1(1:lenmap)

        if (n_operands.eq.2) then
          allocate(fl_item%bcontr%occ_ex2(ngastp,2,nj2))
          fl_item%bcontr%occ_ex2 = ex2
          allocate(fl_item%bcontr%rst_ex2(2,ngas,2,2,nspin,nj2))
          fl_item%bcontr%rst_ex2 = irstex2
          lenmap = len_merge_map(merge2,nj2)
          allocate(fl_item%bcontr%merge_op2(lenmap))
          fl_item%bcontr%merge_op2 = merge2(1:lenmap)
        end if

        lenmap = len_merge_map(merge12,nj12)
        allocate(fl_item%bcontr%merge_op1op2(lenmap))
        fl_item%bcontr%merge_op1op2 = merge12(1:lenmap)

        lenmap = len_merge_map(merge21,nj12)
        allocate(fl_item%bcontr%merge_op2op1(lenmap))
        fl_item%bcontr%merge_op2op1 = merge21(1:lenmap)
      end if

      ! extra for itf
      lenmap = 3+sum(itf_index_info(1:3))
      allocate(fl_item%bcontr%itf_index_info(lenmap))
      fl_item%bcontr%itf_index_info(1:lenmap) = itf_index_info(1:lenmap)

      return
      end
