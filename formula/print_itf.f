*----------------------------------------------------------------------*
      subroutine print_itf(itflog,fl_head,itin,op_info,print_form,
     &                     formlog,tasks,itf_names,itf_targets)
*----------------------------------------------------------------------*
*     Print ITF lines to itflog
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &     itflog,                              ! Output file for ITF algo code
     &     formlog                              ! Output of GeCco formulae
      type(formula_item), intent(in), target ::
     &     fl_head                              ! Linked list of formulae
      type(operator_info), intent(in) ::
     &     op_info                              ! Operator info for printing formulae
      type(tensor_names), intent(in) ::
     &     itf_names                            ! contains renaming information
      type(code_targets), intent(in) ::
     &     itf_targets                          ! contains info for generating ITF code sections
      logical, intent(in) ::
     &     itin,                                ! Create ITIN lines or symmetrise residual at the end
     &     print_form,                          ! Print to optional formulae file
     &     tasks

      type(formula_item), pointer ::
     &     fl_item                              ! Current formula_item
      type(operator), pointer ::
     &     op
      logical ::
     &     last_was_reo                         ! last command was a pure [REORDER]
      integer ::
     &     counter(4), cnt      ! Counter array, 1 helper variable
      integer, parameter ::
     &     maxreo = 1
      integer ::
     &     nreo,
     &     occ_shift(ngastp,2,maxreo), from_to(2,maxreo)
      real(8) ::
     &     fact_reo
      character(len=maxlen_bc_label) ::
     &     label_op_ori, label_op_reo
      integer ::
     &   inter_itype(MAXINT,INDEX_LEN),         ! Store intermediate index-type (itype) info from previous line
     &   ii, idx_code, idx
      type(x_inter) ::
     &   x_dict(MAXX)
      type(inter_spin_cases) ::
     &   inter_spin_dict        ! Store intermediate spin cases from previous line


      ! Point to start of linked list of formulae
      fl_item => fl_head

      counter(1) = 1    ! contraction number
      counter(2) = 1    ! intermediate number
      counter(3) = 1    ! k4e
      counter(4) = 1    ! x intermediate
      inter_itype = 0

      do ii = 1, MAXX
         x_dict(ii)%label = ''
         x_dict(ii)%ops = 0
      end do

      do ii = 1, MAX_SPIN_CASES
         inter_spin_dict%names(ii) = ''
      end do
      inter_spin_dict%ncase = 0

      idx_code = 0

      last_was_reo = .false.

      ! Loop over formula_items, end of the list points to NULL
      do while (associated(fl_item%next))

         ! Optionally print the formula items to another output file
         if (print_form) then
          write(formlog,*) "FORMULA NUMBER: ", counter(1)
          cnt = counter(1)-1 ! beware of print_form_item2: it increments the counter
          call print_form_item2(formlog,'LONG',cnt,fl_item,
     &                          op_info)
         end if

         if (fl_item%command==command_add_intm .or.
     &       fl_item%command==command_cp_intm .or.
     &       fl_item%command==command_add_bc .or.
     &       fl_item%command==command_bc .or.
     &       fl_item%command==command_bc_reo .or.
     &       fl_item%command==command_add_bc_reo .or.
     &       fl_item%command==command_add_reo ) then

            if (fl_item%command==command_bc_reo .or.
     &          fl_item%command==command_add_bc_reo .or.
     &          fl_item%command==command_add_reo) then
              ! patch the result of the contraction with the reordered intermediate info
              ! itf_index_info is already set up correctly
              if (fl_item%bcontr%nj_res < fl_item%reo%nj_out)
     &             call quit(1,'print_itf','dimension trouble')
              fl_item%bcontr%nj_res = fl_item%reo%nj_out
              fl_item%bcontr%occ_res(:,:,1:fl_item%reo%nj_out) =
     &             fl_item%reo%occ_opout(:,:,1:fl_item%reo%nj_out)
            end if

            ! check if there is info about previous reordering on the list
            if (last_was_reo) then
              last_was_reo = .false.
              !write(itflog,'(1x,"patching with: (nreo = ",i4,")")') nreo
              !call wrt_occ_n(itflog,occ_shift,nreo)
              call patch_bcontr_for_reo(fl_item%bcontr,
     &             label_op_ori, label_op_reo, from_to, occ_shift,
     &             nreo)
              !write(itflog,'(1x,"patched!")')
            end if

            call command_to_itf(fl_item%bcontr,itin,itflog,
     &                          fl_item%command, inter_itype,
     &                          itf_names,
     &                          counter,tasks,x_dict,
     &                          inter_spin_dict)

            ! Count the number of terms
            counter(1) = counter(1) + 1

         else if (fl_item%command==command_add_contribution) then
            write(itflog,*) '[CONTR]',fl_item%target  ! this case should not appear in an opt. formula
c         else if (fl_item%command==command_add_bc_reo) then
c            write(itflog,*) '[CONTRACT][REORDER][ADD]',
c     &           fl_item%target
c            call prt_bcontr(itflog,fl_item%bcontr)
c            call prt_reorder(itflog,fl_item%reo)
c            call warn('print_itf',
c     &           'uncovered case appeared: [CONTRACT][REORDER][ADD] ')
c         else if (fl_item%command==command_add_reo) then
c            write(itflog,*) '[REORDER][ADD]',
c     &           fl_item%target
c            call prt_bcontr(itflog,fl_item%bcontr)
c            call prt_reorder(itflog,fl_item%reo)
c            call warn('print_itf',
c     &           'uncovered case appeared: [REORDER][ADD] ')
         else if (fl_item%command==command_reorder) then
!     save here info about reordered operator
            if (last_was_reo) then
              call quit(1,'print_itf',
     &             'did not expect double reordering')
            end if
            last_was_reo = .true.
            label_op_ori = fl_item%reo%label_in
            label_op_reo = fl_item%reo%label_out
            if (fl_item%reo%nj_in.ne.fl_item%reo%nj_out) then
              call quit(1,'print_itf',
     &             'mismatch of nj!')
            end if
            nreo = fl_item%reo%nreo+fl_item%reo%nreo_i0
            if (nreo.gt.maxreo) then
              call quit(1,'print_ift',
     &             'maxreo restricted to max. tested case. Extend!')
            end if
            from_to(1:2,1:nreo) = fl_item%reo%from_to(1:2,1:nreo)
            occ_shift(1:ngastp,1:2,1:nreo) =
     &           fl_item%reo%occ_shift(1:ngastp,1:2,1:nreo)
            !fact_reo = dble(fl_item%reo%sign) factor should not change as the reordering step does not change our
                                             ! previous interpretation of the diagram
c            write(itflog,*) '[REORDER]',
c     &           fl_item%target
c            call prt_reorder(itflog,fl_item%reo)
c            call warn('print_itf',
c     &           'uncovered case appeared: [REORDER] ')
         else if (fl_item%command==command_symmetrise) then
            write(itflog,*) '[SYMMETRISE]',fl_item%target
            call warn('print_itf',
     &           'uncovered case appeared: [SYMMETRISE] ')
         else if (fl_item%command==command_set_target_init ) then
            ! try to find target on target list
            idx = -1
            do ii = 1, itf_targets%ntargets
              if (itf_targets%idx_target(ii)==fl_item%target) then
                idx = ii
                exit
              end if
            end do
            if (idx > 0) then
              ! check if this starts a new code section
              if (itf_targets%idx_code(idx).ne.idx_code) then
                idx_code = itf_targets%idx_code(idx)
                write(itflog,'("CODE_BLOCK: ",a)')
     &               trim(itf_targets%code_name(idx_code))
              end if
            else
              op => op_info%op_arr(fl_item%target)%op
              call warn('print_itf','undeclared target '//trim(op%name))
            end if
            !write(itflog,*) '[INIT TARGET] ',trim(op%name)
         else if (fl_item%command==command_set_target_update ) then
           call quit(1,'print_itf','not prepared for switching targets')
           !write(itflog,*) '[SET TARGET]',fl_item%target
         else if (fl_item%command==command_end_of_formula .or.
     &            fl_item%command==command_set_target_init .or.
     &            fl_item%command==command_set_target_update .or.
     &            fl_item%command==command_new_intermediate .or.
     &            fl_item%command==command_del_intermediate ) then
            ! Do nothing
            ! write(itflog,*) '[END]'
            ! write(itflog,*) '[INIT TARGET]',fl_item%target
            ! write(itflog,*) '[SET TARGET]',fl_item%target
            ! write(itflog,*) '[NEW INTERMEDIATE]',fl_item%target
            ! write(itflog,*) '[DELETE INTERMEDIATE]',fl_item%target
            ! write(itflog,*) '[REORDER]',fl_item%target
         else
            write(itflog,*) 'unknown command ',fl_item%command,
     &           fl_item%target
         end if

         ! Check if at the end of the list, if not, point to the next item
         fl_item => fl_item%next

      end do

      return

      contains

      subroutine patch_bcontr_for_reo(bcontr,
     &             label_op_ori, label_op_reo, from_to, occ_shift,
     &             nreo)
      ! patch bcontr such that the info about the reordered operator is transformed back to that of the original one
      
      type(binary_contr), intent(inout) ::
     &     bcontr
      character(len=*), intent(in) ::
     &     label_op_ori, label_op_reo
      integer, intent(in) ::
     &     nreo,
     &     from_to(2,nreo),
     &     occ_shift(ngastp,2,nreo)
      integer ::
     &     ii

      ! check, which of the operators to patch
      if (trim(label_op_reo)==trim(bcontr%label_op1)) then
        bcontr%label_op1 = label_op_ori
        do ii = 1, nreo
          bcontr%occ_op1(:,:,from_to(1,ii)) =
     &         bcontr%occ_op1(:,:,from_to(1,ii)) + occ_shift(:,:,ii)
          bcontr%occ_op1(:,:,from_to(2,ii)) =
     &         bcontr%occ_op1(:,:,from_to(2,ii)) - occ_shift(:,:,ii)
          bcontr%occ_ex1(:,:,from_to(1,ii)) =
     &         bcontr%occ_ex1(:,:,from_to(1,ii)) + occ_shift(:,:,ii)
          bcontr%occ_ex1(:,:,from_to(2,ii)) =
     &         bcontr%occ_ex1(:,:,from_to(2,ii)) - occ_shift(:,:,ii)
        end do
      else if (trim(label_op_reo)==trim(bcontr%label_op2)) then
        bcontr%label_op2 = label_op_ori
        do ii = 1, nreo
          bcontr%occ_op2(:,:,from_to(1,ii)) =
     &         bcontr%occ_op2(:,:,from_to(1,ii)) + occ_shift(:,:,ii)
          bcontr%occ_op2(:,:,from_to(2,ii)) =
     &         bcontr%occ_op2(:,:,from_to(2,ii)) - occ_shift(:,:,ii)
          bcontr%occ_ex2(:,:,from_to(1,ii)) =
     &         bcontr%occ_ex2(:,:,from_to(1,ii)) + occ_shift(:,:,ii)
          bcontr%occ_ex2(:,:,from_to(2,ii)) =
     &         bcontr%occ_ex2(:,:,from_to(2,ii)) - occ_shift(:,:,ii)
        end do
      else
        write(lulog,*) 'Error: no operator matches with ',
     &       trim(label_op_reo)
        write(lulog,*) 'op1, op2: ',
     &       trim(bcontr%label_op1),
     &       trim(bcontr%label_op2)
        write(lulog,*) 'see also last output on bcontr.tmp'
        call quit(1,'print_itf','problem with [REORDER]')
      end if

      end subroutine
      
      end
