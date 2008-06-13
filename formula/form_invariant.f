*----------------------------------------------------------------------*
      subroutine form_invariant(f_output,f_input,
     &                      title,label_opres,
     &                      ncmpnd,label_op,
     &                      op_info)
*----------------------------------------------------------------------*
*     collect those contractions on ffoutput that do not depend on
*     any operator on list label_op(1:ncmpnd)
*
*     andreas, april 2007 (new version: december 2007)
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
c      include 'def_contraction_list.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ncmpnd
      character(*), intent(in) ::
     &     title,
     &     label_opres, label_op(ncmpnd)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     ok, same, transpose
      integer ::
     &     nterms, idum, idxinp, idx, icmpnd, idxres, idxop(ncmpnd), len
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     fl_pnt, fl_pnt_next
      type(contraction), pointer ::
     &     contr

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'here speaks form_indep')
        write(luout,*) ' f_input  = ',trim(f_input%label)
        write(luout,*) ' f_output = ',trim(f_output%label)
        write(luout,*) ' op_res  = ',trim(label_opres)
        do icmpnd = 1, ncmpnd
          write(luout,*) 'compound # ',icmpnd
          write(luout,*) ' op  = ',trim(label_op(icmpnd))
        end do
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! get indices
      idxres = idx_oplist2(label_opres,op_info)
      if (idxres.lt.0)
     &     call quit(1,'form_invariant',
     &     'required operators are not yet defined? '//
     &       trim(label_opres))

      do icmpnd = 1, ncmpnd

        ! look for transposition label
        len = len_trim(label_op(icmpnd))
        transpose = (label_op(icmpnd)(len-1:len).eq.'^+') 
        if (transpose) len = len-2

        idxop(icmpnd) = idx_oplist2(label_op(icmpnd)(1:len),op_info)

        if (idxop(icmpnd).lt.0)
     &     call quit(1,'form_invariant',
     &     'required operators are not yet defined? '//
     &       label_op(icmpnd)(1:len))

        if (transpose) idxop(icmpnd) = -idxop(icmpnd)

      end do
c dbg
c      print *,idxop
c      print *,idxres
c dbg

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist)

      fl_pnt => flist

      if (.not.associated(fl_pnt))
     &     call quit(1,'form_invariant',
     &     'empty formula list? something is buggy')

      nterms = 0
      fl_loop: do

        fl_pnt_next => fl_pnt%next
        if (fl_pnt%command.eq.command_end_of_formula)
     &       exit fl_loop

        fl_pnt%target = idxres

        if (fl_pnt%command.eq.command_add_contribution) then

          ok = .true.
          contr => fl_pnt%contr
          cmp_loop: do idx = 1, contr%nvtx
            do icmpnd = 1, ncmpnd
              if (contr%vertex(idx)%idx_op.eq.abs(idxop(icmpnd)).and.
     &           (contr%vertex(idx)%dagger.eqv.(idxop(icmpnd).lt.0)))
     &        then
                ok = .false.
                exit cmp_loop
              end if
            end do
          end do cmp_loop

          if (ok) then
            nterms = nterms+1
            contr%idx_res = idxres ! not completely OK
            if (ntest.ge.100) then
              call prt_contr2(luout,contr,op_info)
            end if
          else
            ! deallocate contents and re-link the list
            call delete_fl_node(fl_pnt)
            ! deallocate the node itself
            deallocate(fl_pnt)
          end if
        else 
        end if

        fl_pnt => fl_pnt_next
        if (.not.associated(fl_pnt))
     &       call quit(1,'form_invariant',
     &       'unexpected end of formula list')

      end do fl_loop

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

c dbg
c      write(luout,*)'TeX list'
c      call tex_form_list(luout,flist,op_info)
c dbg

      call dealloc_formula_list(flist)

      if (ntest.ge.10) then
        write(luout,*) 'generated terms: ',nterms
      end if
      
      return
      end
