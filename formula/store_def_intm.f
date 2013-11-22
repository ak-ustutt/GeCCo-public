*----------------------------------------------------------------------*
      subroutine store_def_intm(fl_item,
     &     label,occ,rst,nj,nblk,
     &     parent1,parent2,tra_new,tra1,tra2,
     &     orb_info)
*----------------------------------------------------------------------*
*     write info for creating new intermediate on formula list item
*     (fl_item)
*     the operator labels parent1 and parent2 are needed for 
*     obtaining the list-attributes when evaluating the formula
*     (e.g. the IRREP of the new intermediate is the direct product
*      of the IRREPs of parent1 and parent2; parent2 may be empty)
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(inout) ::
     &     fl_item
      character(len=*), intent(in) ::
     &     label, parent1, parent2
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     nj, nblk, 
     &     occ(ngastp,2,nj*nblk),
     &     rst(2,orb_info%ngas,2,2,orb_info%nspin,nj*nblk)
      logical, intent(in) ::
     &     tra_new, tra1, tra2

      if (ntest.ge.100) write(lulog,*) 'storing def of intermediate'

      fl_item%interm%n_occ_cls = nblk
      fl_item%interm%njoined = nj
      fl_item%interm%name = trim(label)
      call init_operator(fl_item%interm,orb_info)

      fl_item%interm%ihpvca_occ = occ
      fl_item%interm%igasca_restr = rst
      
      call occ2caocc(fl_item%interm%ica_occ,occ,nj,nblk)

      fl_item%interm%formal_blk = .false.

      fl_item%parent1 = trim(parent1)
      fl_item%parent2 = trim(parent2)

      fl_item%tra = tra_new
      fl_item%tra1 = tra1
      fl_item%tra2 = tra2

      return
      end
