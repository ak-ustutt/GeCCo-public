*----------------------------------------------------------------------*
      subroutine form_r12exc_split(f_output,f_input,
     &                      label_res,
     &                      mode,
     &                      labels,nlabels,
     &                      op_info)
*----------------------------------------------------------------------*
*     split expression for <L|Hbar|R> into contributions
*     mode:
*       1 - <L1|Hhat|R1>          (no R12)
*       2 - <L1|HBar-Hhat|R1>   
*       3 - <L1|HBar|R2>
*       4 - <L2|HBar|R1>
*       5 - <L2|HBar|R2>
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
     &     mode, nlabels
      character(len=*) ::
     &     label_res,
     &     labels(nlabels)
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
     &     idxop(:)

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'here speaks form_select_terms')
        write(lulog,*) ' f_input  = ',trim(f_input%label)
        write(lulog,*) ' f_output = ',trim(f_output%label)
      end if

      allocate(idxop(nlabels))

      idxop_res = idx_oplist2(label_res,op_info)
      if (idxop_res.le.0)
     &       call quit(1,'form_r12exc_split',
     &       'label not on list: '//trim(label_res))
      do idx = 1, nlabels
        if (trim(labels(idx)).eq.'-') then
          idxop(idx) = -1
          cycle
        end if
        idxop(idx) = idx_oplist2(labels(idx),op_info)
        if (idxop(idx).le.0)
     &       call quit(1,'form_r12exc_split',
     &       'label not on list: '//trim(labels(idx)))
      end do

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      call r12exc_split_terms(flist,
     &     mode,idxop_res,
     &     idxop(1),idxop(2),idxop(3),idxop(4),
     &     idxop(5),idxop(6),idxop(7),
     &     op_info)

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = 'XXX'
      call write_form_list(f_output%fhand,flist,'XXX')

      call dealloc_formula_list(flist)

      deallocate(idxop)

      return
      end
