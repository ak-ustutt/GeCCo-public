      subroutine fl_item_rename_intm(fl_item,label_new,set_incore,
     &                               orb_info)

      ! wrapper for renaming an operator on a fl_item:
      ! we need to do some detour via clone_operator to
      ! maintain a correct memory book-keeping ...

      implicit none

      include 'opdim.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'


      character(len=19) ::
     &     i_am = 'fl_item_rename_intm'

      type(formula_item), intent(inout) ::
     &     fl_item
      character(len=*), intent(in) ::
     &     label_new
      logical ::
     &     set_incore
      type(orbinf) ::
     &     orb_info


      type(operator), pointer ::
     &     intm_new


      if (.not.associated(fl_item%interm))
     &    call quit(1,i_am,'No intermediate allocated here ...')

      allocate(intm_new)
      intm_new%name = label_new
      call clone_operator(intm_new,fl_item%interm,.false.,orb_info)
      call dealloc_operator(fl_item%interm)
      deallocate(fl_item%interm)
      fl_item%interm => intm_new
      if (set_incore) fl_item%incore = 1

      end

