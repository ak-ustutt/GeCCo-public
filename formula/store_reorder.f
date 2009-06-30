      subroutine store_reorder(fl_item,
     &     label_out,label_in,
     &     iblk_out,iblk_in,
     &     sign_reo,occ_op0,
     &     from_to,occ_shift,nreo,
     &     occ_opout,rst_opout,nj_out,
     &     occ_opin, rst_opin, nj_in,
     &     merge1,merge1inv,merge2,merge2inv,
     &     orb_info)

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
      character(len=*), intent(in) ::
     &     label_out, label_in
      integer, intent(in) ::
     &     iblk_out, iblk_in,
     &     nreo, nj_out, nj_in
      integer, intent(in) ::
     &     sign_reo,
     &     from_to(2,nreo),
     &     occ_shift(ngastp,2,nreo),
     &     occ_opin(ngastp,2,nj_in),
     &     rst_opin(2,orb_info%ngas,2,2,orb_info%nspin,nj_in),
     &     occ_opout(ngastp,2,nj_out),
     &     rst_opout(2,orb_info%ngas,2,2,orb_info%nspin,nj_out),
     &     occ_op0(ngastp,2,nj_in),
     &     merge1(*),merge1inv(*),
     &     merge2(*),merge2inv(*)

      integer ::
     &     ngas, nspin, lenmap

      integer, external ::
     &     len_merge_map

      if (ntest.ge.100)
     &     write(luout,*) 'storing reordering'
      
      ngas = orb_info%ngas
      nspin = orb_info%nspin

      fl_item%reo%label_in  = label_in
      fl_item%reo%label_out = label_out
      fl_item%reo%nreo    = nreo
      fl_item%reo%ngas      = ngas
      fl_item%reo%nspin     = nspin
      fl_item%reo%nj_in   = nj_in
      fl_item%reo%nj_out  = nj_out
      fl_item%reo%iblk_in   = iblk_in
      fl_item%reo%iblk_out  = iblk_out
      fl_item%reo%sign      = sign_reo
      allocate(fl_item%reo%from_to(2,nreo))
      fl_item%reo%from_to = from_to
      allocate(fl_item%reo%occ_opin(ngastp,2,nj_in))
      fl_item%reo%occ_opin = occ_opin
      allocate(fl_item%reo%rst_opin(2,ngas,2,2,nspin,nj_in))
      fl_item%reo%rst_opin = rst_opin
      allocate(fl_item%reo%occ_opout(ngastp,2,nj_out))
      fl_item%reo%occ_opout = occ_opout
      allocate(fl_item%reo%rst_opout(2,ngas,2,2,nspin,nj_out))
      fl_item%reo%rst_opout = rst_opout
      allocate(fl_item%reo%occ_shift(ngastp,2,nreo))
      fl_item%reo%occ_shift = occ_shift
      allocate(fl_item%reo%occ_op0(ngastp,2,nj_in))
      fl_item%reo%occ_op0 = occ_op0

      lenmap = len_merge_map(merge1,nj_in)
      allocate(fl_item%reo%merge_stp1(lenmap))
      fl_item%reo%merge_stp1 = merge1(1:lenmap)

      lenmap = len_merge_map(merge1inv,nj_in)
      allocate(fl_item%reo%merge_stp1inv(lenmap))
      fl_item%reo%merge_stp1inv = merge1inv(1:lenmap)

      lenmap = len_merge_map(merge2,nj_out)
      allocate(fl_item%reo%merge_stp2(lenmap))
      fl_item%reo%merge_stp2 = merge2(1:lenmap)

      lenmap = len_merge_map(merge2inv,nj_out)
      allocate(fl_item%reo%merge_stp2inv(lenmap))
      fl_item%reo%merge_stp2inv = merge2inv(1:lenmap)

      return
      end
