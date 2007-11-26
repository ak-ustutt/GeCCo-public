*----------------------------------------------------------------------*
      subroutine form_opt(ffform_opt,
     &     nfcat,idxform,
     &     form_info,op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     given a list of formulae, concatenate them into one formula
*     file, find optimal factorization and intermediates
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_formula_item.h'
      include 'ifc_input.h'

      integer, intent(in) ::
     &     nfcat, idxform(nfcat)
      type(filinf), intent(inout) ::
     &     ffform_opt
      type(formula_info), intent(inout) ::
     &     form_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      
      type(formula_item), pointer ::
     &     form_head, form_tail, form_ptr
      type(filinf), pointer ::
     &     cur_ffile

      integer ::
     &     icat, iprint, lentitle, isim

      character ::
     &     title*(form_maxlen_comment)

      logical, external ::
     &     rd_formula

      iprint = max(ntest,iprlvl)

      call write_title(luout,wst_section,'Formula optimization')

      ! initialize list
      allocate(form_head)
      form_ptr => form_head
      form_ptr%command = command_end_of_formula
      nullify(form_ptr%next)
      nullify(form_ptr%contr)
      nullify(form_ptr%interm)

      lentitle = 0
      ! ----------------------
      ! read in formula files:
      ! ----------------------
      if (iprint.gt.0)
     &     write(luout,'(x,a)') 'Reading in:'

      do icat = 1, nfcat
c dbg
c        print *,'idxform(icat) = ',idxform(icat)
c        print *,'>',trim(form_info%form_arr(idxform(icat))%form%label)
c        print *,'>',
c     &       trim(form_info%form_arr(idxform(icat))%form%fhand%name)
c dbg
      if (iprint.gt.0)
     &     write(luout,'(2x,"--",x,a)')
     &       trim(form_info%form_arr(idxform(icat))%form%label)
        cur_ffile => form_info%form_arr(idxform(icat))%form%fhand
        if (lentitle.lt.form_maxlen_comment) then
          if (icat.eq.1) then
            title = form_info%form_arr(idxform(icat))%form%label
          else
            title = trim(title)//'/'//
     &           form_info%form_arr(idxform(icat))%form%label
          end if
        end if

        call read_form_list(cur_ffile,form_ptr)
c dbg
        print *,'raw formula'
        call print_form_list(luout,form_ptr,op_info)
c dbg

        ! advance form_ptr to end of list
        ! (not possible via call list due to ifort problems)
        do while(associated(form_ptr%next))
          form_ptr => form_ptr%next
        end do
      end do
      form_tail => form_ptr

      title = trim(title)//' -- optimized'

      if (is_keyword_set('method.CC')) then
        call get_argument_value('calculate.routes','simtraf',ival=isim)
        if (iprlvl.gt.0)
     &       write(luout,'(2x,a)')
     &       'I will factor out the e^-T1 H e^T1 intermediate ...'
        if (isim.gt.0)
     &       call cc_form_hhat_replace(form_head,
     &                                 form_info,op_info)
      end if

      ! ----------------------------------------
      ! find optimal factorization for each term
      ! ----------------------------------------
c      form_ptr => form_head
      if (iprint.gt.0)
     &     write(luout,'(2x,a)')
     &       'Now looking for the optimal factorization of terms ...'
      call factorize(form_head,op_info,str_info,orb_info)

      if (iprint.ge.10) then
        call write_title(luout,wst_around_double,'Optimized formula:')
        call print_form_list(luout,form_head,op_info)
      end if

      call write_form_list(ffform_opt,form_head,title)

      call dealloc_formula_list(form_head)
      deallocate(form_head)
      
      return
      end
