*----------------------------------------------------------------------*
      subroutine form_select_terms(f_output,f_input,
     &                      label_res,
     &                      ninclude,label_incl,iblk_incl,
     &                      ninclude_or,label_incl_or,iblk_incl_or,
     &                      nexclude,label_excl,iblk_excl,
     &                      op_info)
*----------------------------------------------------------------------*
*     keep all terms that
*     - have operator and block as requested by "include"
*     - do not have operator and block as requested by "exclude"
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ninclude, ninclude_or, nexclude,
     &     iblk_incl(ninclude), iblk_incl_or(ninclude_or),
     &     iblk_excl(nexclude)
      character(len=*) ::
     &     label_res,
     &     label_incl(ninclude), label_incl_or(ninclude_or),
     &     label_excl(nexclude)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     on_list, same, delete
      integer ::
     &     idx, idxop_res
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     fl_pnt, fl_pnt_next
      type(contraction), pointer ::
     &     contr

      integer, pointer ::
     &     idxop_incl(:), idxop_incl_or(:), idxop_excl(:)

      integer, external ::
     &     idx_oplist2, imltlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &       'here speaks form_select_terms')
        write(luout,*) ' f_input  = ',trim(f_input%label)
        write(luout,*) ' f_output = ',trim(f_output%label)
      end if

      allocate(idxop_incl(ninclude),idxop_incl_or(ninclude_or),
     &         idxop_excl(nexclude))

      idxop_res = idx_oplist2(label_res,op_info)
      if (idxop_res.le.0)
     &       call quit(1,'form_select_terms',
     &       'label not on list: '//trim(label_res))
      do idx = 1, ninclude
        idxop_incl(idx) = idx_oplist2(label_incl(idx),op_info)
        if (idxop_incl(idx).le.0)
     &       call quit(1,'form_select_terms',
     &       'label not on list: '//trim(label_incl(idx)))
      end do
      do idx = 1, ninclude_or
        idxop_incl_or(idx) = idx_oplist2(label_incl_or(idx),op_info)
        if (idxop_incl_or(idx).le.0)
     &       call quit(1,'form_select_terms',
     &       'label not on list: '//trim(label_incl_or(idx)))
      end do
      do idx = 1, nexclude
        idxop_excl(idx) = idx_oplist2(label_excl(idx),op_info)
        if (idxop_excl(idx).le.0)
     &       call quit(1,'form_select_terms',
     &       'label not on list: '//trim(label_excl(idx)))
      end do

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist)

      call select_terms(flist,
     &     idxop_res,
     &     ninclude,idxop_incl,iblk_incl,
     &     ninclude_or,idxop_incl_or,iblk_incl_or,
     &     nexclude,idxop_excl,iblk_excl,
     &     op_info)

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = 'XXX'
      call write_form_list(f_output%fhand,flist,'XXX')

      call dealloc_formula_list(flist)

      deallocate(idxop_incl,idxop_excl)

      return
      end
