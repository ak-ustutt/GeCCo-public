*----------------------------------------------------------------------*
      subroutine factor_out(form_head,label_intm,
     &     form_info,op_info)
*----------------------------------------------------------------------*
*     given a (un-factorized) formula list, prepend the rules to
*     calculate the intermediate (label "label_intm") and replace the 
*     appropriate terms in the formula list with that intermediate
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_formula_item.h'
      include 'par_formnames_gen.h'

      type(formula_item), target ::
     &     form_head
      character(*), intent(in) ::
     &     label_intm
      type(formula_info) ::
     &     form_info
      type(operator_info) ::
     &     op_info

      integer ::
     &     idx_intm, nrpl

      integer, external ::
     &     idx_formlist

      type(formula_item), pointer ::
     &     form_link, form_ptr, fintm_head, fintm_tail

      idx_intm = idx_formlist(label_intm,form_info)

      if (idx_intm.le.0)
     &     call quit(1,'factor_out',
     &     'label was not found on list: '//trim(label_intm))

      allocate(form_link)

      ! we currently cannot pass back pointers to complex types, so
      ! copy (!) current form_head to new place in memory ...
      form_link = form_head

      ! ... and use form_head as head for Intm list
      fintm_head => form_head
      ! read Intm definition
      call read_form_list(form_info%form_arr(idx_intm)%form%fhand,
     &     fintm_head)

      form_ptr => form_head
      ! advance form_ptr to end of list
      do while(associated(form_ptr%next))
        form_ptr => form_ptr%next
      end do
      fintm_tail => form_ptr

      ! factor out the sub-expressions
      call factor_out_subexpr2(form_link,fintm_head,nrpl,op_info)

      ! remove [END] from Intm list (if any)
      if (fintm_tail%command.eq.command_end_of_formula) then
        fintm_tail => fintm_tail%prev
        deallocate(fintm_tail%next)
      end if
      ! link the lists
      fintm_tail%next => form_link
      form_link%prev => fintm_tail

      return
      end
