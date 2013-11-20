*----------------------------------------------------------------------*
      subroutine form_del_terms(f_output,f_input,
c     &                      title,
     &                      nterms,idxterms,mode,
     &                      op_info)
*----------------------------------------------------------------------*
*     remove all terms that
*     mode=-1: are listed in idxterms
*     mode=+1: are not listed in idxterms
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
     &     mode, nterms, idxterms(nterms)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     on_list, same, delete
      integer ::
     &     idx
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     fl_pnt, fl_pnt_next
      type(contraction), pointer ::
     &     contr

      integer, external ::
     &     idx_oplist2, imltlist

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'here speaks form_del_terms')
        write(lulog,*) ' f_input  = ',trim(f_input%label)
        write(lulog,*) ' f_output = ',trim(f_output%label)
      end if

      same = trim(f_input%label).eq.trim(f_output%label)


      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      fl_pnt => flist

      if (.not.associated(fl_pnt))
     &     call quit(1,'form_del_terms',
     &     'empty formula list? something is buggy')

      idx = 0
      fl_loop: do

        fl_pnt_next => fl_pnt%next
        if (fl_pnt%command.eq.command_end_of_formula)
     &       exit fl_loop

        if (fl_pnt%command.eq.command_add_contribution) then

          idx = idx+1
c dbg          
c          print  *,'present term ',idx
c          call prt_contr2(lulog,fl_pnt%contr,op_info)        
c dbg

          on_list = imltlist(idx,idxterms,nterms,1).gt.0
          delete = mode.lt.0.and.on_list .or. mode.gt.0.and..not.on_list
c dbg
c          print *,'delete = ',delete
c dbg

          if (delete) then
            if (iprlvl.ge.10) then
              write(lulog,*) 'removing term # ',idx
              call prt_contr2(lulog,fl_pnt%contr,op_info)
            end if
            ! deallocate contents and re-link the list
            call delete_fl_node(fl_pnt)
            ! deallocate the node itself
            deallocate(fl_pnt)
          end if
        end if

        fl_pnt => fl_pnt_next
        if (.not.associated(fl_pnt))
     &       call quit(1,'form_del_terms',
     &       'unexpected end of formula list')

      end do fl_loop

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = 'XXX'
      call write_form_list(f_output%fhand,flist,'XXX')

      call dealloc_formula_list(flist)

      return
      end
