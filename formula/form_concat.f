*----------------------------------------------------------------------*
      subroutine form_concat(f_output,
     &                      title,
     &                      nform,label_form,nfac,fac_form,
     &                      op_info,form_info)
*----------------------------------------------------------------------*
*
*     driver for concatenating several formulae into one
*
*     the first nfac forms are scaled by fac_form (otherwise 1.0 is assumed)
*
*     f_output should not be identical to any of the input forms
*
*     andreas, sept 2021
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
      include 'ifc_formula.h' 

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     nform, nfac
      real(8), intent(in) ::
     &     fac_form(nfac)
      character(*), intent(in) ::
     &     title,
     &     label_form(nform)
      type(formula), intent(inout) ::
     &     f_output
      type(operator_info), intent(in) ::
     &     op_info
      type(formula_info), intent(in) ::
     &     form_info

      logical ::
     &     same, transpose
      integer ::
     &     ii, idx, len, nrpl, nspl
      character ::
     &     name*(form_maxlen_label*2)
      real(8) ::
     &     factor
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      type(filinf), pointer ::
     &     ffform
      type(formula_item), pointer ::
     &     flist, fl_append(:)  ! have to be pointers as they get part of the list
                                ! (otherwise dealloc_formula_list will fail)

      integer, external ::
     &     idx_formlist, form_count

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &     'form_concat at work')
        write(lulog,*) ' f_output = ',trim(f_output%label)
        do ii = 1, nform
          write(lulog,*) ii,trim(label_form(ii))
        end do
        do ii = 1, nfac
          write(lulog,*) ii,fac_form(ii)
        end do
      end if

      call atim_csw(cpu0,sys0,wall0)

      allocate(flist,fl_append(nform))
      ! read in input formula
      call init_formula(flist)

      ! loop over forms
      do ii = 1, nform

        ! look for transposition label
        len = len_trim(label_form(ii))
        transpose = (label_form(ii)(len-1:len).eq.'^+') 
        if (transpose) len = len-2
         
        ! get index
        idx = idx_formlist(label_form(ii)(1:len),form_info)

        if (idx.le.0)
     &       call quit(1,'form_factor_out',
     &       'formula label not on list: '//trim(label_form(ii)))
        
        ! get factor
        if (ii.le.nfac) then
          factor = fac_form(ii)
        else
          factor = 1d0
        end if

        if (ntest.ge.100)
     &       write(lulog,*)
     &       'now processing: ',trim(label_form(ii)),
     &       ' transpose: ',transpose, ' factor: ',factor


        ffform => form_info%form_arr(idx)%form%fhand
        if (.not.associated(ffform))
     &       call quit(1,'form_factor_out',
     &       'formula file does not exist for '//
     &       trim(form_info%form_arr(idx)%form%label))

        call init_formula(fl_append(ii))
        call read_form_list(ffform,fl_append(ii),.true.)

        if (transpose)
     &       call transpose_formula(fl_append(ii),op_info)

        call append_flist(flist,fl_append(ii),factor,op_info)

        call atim_csw(cpu,sys,wall)
        

        ! no need to do this: all items are either deleted in append_flist or appended to prev. list
        ! call dealloc_formula_list(fl_append)

      end do

      ! write result
      write(name,'(a,".fml")') trim(f_output%label)
      call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

      if (iprlvl.ge.2)
     &     call prtim(lulog,'appending',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      ! this should release the entire (contatenated) list
      call dealloc_formula_list(flist)
      deallocate(flist)


      write(lulog,'(" end")')


      return
      end
