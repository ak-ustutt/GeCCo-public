*----------------------------------------------------------------------*
      subroutine form_opt(f_opt,
     &     nfcat,flabels,ninterm,finlabels,
     &     form_info,op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     given a list of formulae, concatenate them into one formula
*     file, find optimal factorization and intermediates
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'routes.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_formula_item.h'
      include 'ifc_input.h'

c      logical, parameter ::
cc     &     new = .true.
c     &     new = .false.
      integer, intent(in) ::
     &     nfcat, ninterm
      character(*), intent(in) ::
     &     flabels(nfcat), finlabels(ninterm)
      type(formula), intent(inout) ::
     &     f_opt
      type(formula_info), intent(inout) ::
     &     form_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      
      type(formula_item), pointer ::
     &     fl_head, fl_tail, fl_ptr, fl_new, fl_opt
      type(filinf), pointer ::
     &     cur_ffile

      integer ::
     &     icat, iint, iprint, lentitle, idxform

      character ::
     &     title*(form_maxlen_comment), name*(form_maxlen_label*2)

      integer, external ::
     &     idx_formlist

      iprint = max(ntest,iprlvl)

      call write_title(luout,wst_section,'Formula optimization')

c dbg
      call print_op_info(luout,'op',op_info)
      call print_op_info(luout,'mel',op_info)
c dbg
      ! initialize list
      allocate(fl_head)
      fl_ptr => fl_head
      call init_formula(fl_ptr)
c      fl_ptr%command = command_end_of_formula
c      nullify(fl_ptr%next)
c      nullify(fl_ptr%contr)
c      nullify(fl_ptr%interm)

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
        idxform = idx_formlist(flabels(icat),form_info)

        if (idxform.le.0)
     &       call quit(1,'form_opt','formula label not on list: '//
     &       trim(flabels(icat)))

        if (iprint.gt.0)
     &     write(luout,'(2x,"--",x,a)')
     &       trim(form_info%form_arr(idxform)%form%label)

        cur_ffile => form_info%form_arr(idxform)%form%fhand
        if (.not.associated(cur_ffile))
     &       call quit(1,'form_opt',
     &       'formula file does not exist for '//
     &       trim(form_info%form_arr(idxform)%form%label))

        if (lentitle.lt.form_maxlen_comment) then
          if (icat.eq.1) then
            title = form_info%form_arr(idxform)%form%label
          else
            title = trim(title)//'/'//
     &           form_info%form_arr(idxform)%form%label
          end if
        end if

        call read_form_list(cur_ffile,fl_ptr,icat.eq.1)
c dbg
c        print *,'raw formula'
c        call print_form_list(luout,form_ptr,op_info)
c dbg

        ! advance form_ptr to end of list
        ! (not possible via call list due to ifort problems)
        do while(associated(fl_ptr%next))
          fl_ptr => fl_ptr%next
        end do
      end do
      fl_tail => fl_ptr

      title = trim(title)//' -- optimized'


      ! -----------------------------------------------------------
      ! round one: loop over suggested intermediates and
      ! factor them out; the intermediate definition will be
      ! prepended the present list
      ! in effect, the intermediates may successively be factored
      ! as well
      ! -----------------------------------------------------------
      do iint = 1, ninterm
        if (iprlvl.gt.0)
     &       write(luout,'(2x,a)')
     &       'I will factor out the intermediate: '//
     &       trim(finlabels(iint))//' ...'
        call factor_out(fl_head,finlabels(iint),
     &       form_info,op_info)
      end do
      ! ----------------------------------------
      ! round two:
      ! find optimal factorization for each term
      ! ----------------------------------------
      if (iprint.gt.0)
     &     write(luout,'(2x,a)')
     &       'Now looking for the optimal factorization of terms ...'
     
      if (irt_sched.eq.0) then ! old scheduler
        call factorize(fl_head,op_info,str_info,orb_info,f_opt%label)
        fl_opt => fl_head
      else ! new scheduler needs new factorization
        allocate(fl_new)
        call init_formula(fl_new)
c dbg
        print *,'testing new factorization'
c dbg        
        call factorize_new(fl_new,fl_head,
     &       op_info,str_info,orb_info,f_opt%label)

        fl_opt => fl_new

      end if

      ! ----------------------------------------
      ! round three:
      ! automatic intermediates (to come ...)
      ! ----------------------------------------

      if (iprint.ge.10) then
        call write_title(luout,wst_around_double,'Optimized formula:')
        call print_form_list(luout,fl_opt,op_info)
      end if
cmh      if (lustat.gt.0) call print_form_list(lustat,fl_opt,op_info)

      write(name,'(a,".fml")') trim(f_opt%label)
      call file_init(f_opt%fhand,name,ftyp_sq_unf,0)      
      f_opt%comment = trim(title)
      call write_form_list(f_opt%fhand,fl_opt,title)

      call dealloc_formula_list(fl_head)
      deallocate(fl_head)
      
      if (irt_sched.gt.0) then
        call dealloc_formula_list(fl_new)
        deallocate(fl_new)
      end if

      return
      end
