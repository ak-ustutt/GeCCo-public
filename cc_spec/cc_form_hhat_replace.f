*----------------------------------------------------------------------*
      subroutine cc_form_hhat_replace(form_head,
     &     form_info,op_info)!,str_info,orb_info)
*----------------------------------------------------------------------*
*     given a (un-factorized) formula list, prepend the rules to
*     calculate Hhat and replace the appropriate terms in the formula
*     list with the Hhat intermediate
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
      type(formula_info) ::
     &     form_info
      type(operator_info) ::
     &     op_info

      integer ::
     &     idx_form_hhat

      integer, external ::
     &     idx_formlist

      type(formula_item), pointer ::
     &     form_link, form_ptr, fhhat_head, fhhat_tail

c      idx_op_hhat = idx_oplist2(op_hhat,op_info)
      idx_form_hhat = idx_formlist(label_cchhat,form_info)
c dbg
c      print *,'idx_form_hhat = ',idx_form_hhat
c dbg
      if (idx_form_hhat.le.0)
     &     call quit(1,'cc_form_hhat_replace',
     &     'Hhat was not found on list')

      allocate(form_link)

      ! we currently cannot pass back pointers to complex types, so
      ! copy (!) current form_head to new place in memory ...
      form_link = form_head

      ! ... and use form_head as head for Hhat list
      fhhat_head => form_head
      ! read Hhat definition
      call read_form_list(form_info%form_arr(idx_form_hhat)%form%fhand,
     &     fhhat_head)

      form_ptr => form_head
      ! advance form_ptr to end of list
      do while(associated(form_ptr%next))
        form_ptr => form_ptr%next
      end do
      fhhat_tail => form_ptr

      ! factor out the sub-expressions
      call factor_out_subexpr(form_link,fhhat_head,op_info)

      ! remove [END] from Hhat list (if any)
      if (fhhat_tail%command.eq.command_end_of_formula) then
        fhhat_tail => fhhat_tail%prev
        deallocate(fhhat_tail%next)
      end if
      ! link the lists
      fhhat_tail%next => form_link
      form_link%prev => fhhat_tail

      return
      end
