*----------------------------------------------------------------------*
      subroutine form_op_replace_drv(f_output,f_input,
     &                      title,
     &                      nreplace,label_replace,
     &                      op_info)
*----------------------------------------------------------------------*
*
*     driver for renaming operators in input formula
*
*     label_replace contains 2*nreplace labels, the first of each
*     pair is the operator to be replaced and the second is the
*     operator to be inserted instead
*
*     f_input and f_output may be identical
*
*     andreas, april 2008
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
c      include 'def_contraction_list.h'
      include 'def_formula_item.h'
      include 'mdef_formula_info.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     nreplace
      character(*), intent(in) ::
     &     title,
     &     label_replace(nreplace*2)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     same, transpose1, transpose2
      integer ::
     &     irepl, irepl1, irepl2, idx1, idx2, len1, len2
      character ::
     &     name*(form_maxlen_label*2)

      type(filinf), pointer ::
     &     ffintm
      type(formula_item) ::
     &     flist, fl_intm

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'form_op_rename reports')
        write(luout,*) ' f_input  = ',trim(f_input%label)
        write(luout,*) ' f_output = ',trim(f_output%label)
        do irepl = 1, nreplace
          write(luout,*) irepl,trim(label_replace(2*(irepl-1)+1)),
     &                         trim(label_replace(2*(irepl-1)+2))
        end do
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist)

      ! loop over intermediates
      do irepl = 1, nreplace

        irepl1 = (irepl-1)*2+1
        irepl2 = (irepl-1)*2+2

        ! look for transposition label of first operator
        len1 = len_trim(label_replace(irepl1))
        transpose1 = (label_replace(irepl1)(len1-1:len1).eq.'^+') 
        if (transpose1) len1 = len1-2         
        ! get index
        idx1 = idx_oplist2(label_replace(irepl1)(1:len1),op_info)
        if (idx1.le.0)
     &       call quit(1,'form_op_replace_drv',
     &    'operator label 1 not on list: '//trim(label_replace(irepl1)))
        
        ! look for transposition label of second operator
        len2 = len_trim(label_replace(irepl2))
c dbg fix by mh
c original line        transpose2 = (label_replace(irepl2)(len2-1:len2).eq.'^+') 
        if (len2.gt.1) then
          transpose2 = (label_replace(irepl2)(len2-1:len2).eq.'^+')
        else
          transpose2 = .false.
        end if
c dbg end fix
        if (transpose2) len2 = len2-2
        ! get index
        idx2 = idx_oplist2(label_replace(irepl2)(1:len2),op_info)
        if (idx2.le.0)
     &       call quit(1,'form_op_replace_drv',
     &    'operator label 2 not on list: '//trim(label_replace(irepl2)))
        
        if (ntest.ge.100)
     &       write(luout,*)
     &       'now replacing: ',
     &       label_replace(irepl1)(1:len1),
     &       ' transpose: ',transpose1,' by ',
     &       label_replace(irepl2)(1:len2),
     &       ' transpose: ',transpose2

        call op_replace(idx1,transpose1,
     &                  idx2,transpose2,
     &                  .false.,.true.,flist,op_info)
c     &                  .false.,.false.,flist,op_info)

      end do

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

      if (ntest.ge.100) then
        call write_title(luout,wst_around_double,'Modified formula:')
        call print_form_list(luout,flist,op_info)
      end if

      call dealloc_formula_list(flist)

      return
      end
