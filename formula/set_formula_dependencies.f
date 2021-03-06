*----------------------------------------------------------------------*
      subroutine set_formula_dependencies(depend,flist,op_info)
*----------------------------------------------------------------------*
*     scan formula list flist and find out:
*       which ME-lists are updated -> set indices on idxlist
*       on which ME-list do these results depend ->> depends_on_idxlist
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_dependency_info.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, parameter ::
     &     maxscr = 1024,
     &     ntest = 00

      type(formula_item), intent(in), target ::
     &     flist
      type(dependency_info), intent(inout) ::
     &     depend
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     sorted, error
      integer ::
     &     ipass, iptr, iptr_last, itarget, nlist, maxlist,
     &     nvtx, ivtx, idx_op
      integer ::
     &     idxscr(maxscr)

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     fl_ptr
      integer, pointer ::
     &     op2list(:), idxlist(:), depends_on_idxlist(:,:)

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'this is set_formula_dependencies')
      end if

      error = .false.
      op2list => op_info%op2list
! first pass: determine dimensions of the lists
      ! second pass: allocate lists
      do ipass = 1, 2
        if (ntest.ge.100) write(lulog,*) 'pass: ',ipass
        fl_ptr => flist
        iptr = 0
        iptr_last = 0
        itarget = 0
        maxlist=0
        sorted = .false.
        do
          select case(fl_ptr%command)
          case(command_end_of_formula,command_set_target_init)
            ! new target starts (or end of the story)
            ! did we have any contributions to previous target?
            if (iptr.gt.0) then
              ! remember length and extract the unique indices
              nlist = iptr
              if (.not.sorted) call unique_list(idxscr,nlist)
              ! remember max. dimension
              maxlist = max(maxlist,nlist)
              if (ipass.eq.2.and.itarget.gt.0) then
                ! store indices (converted to list indices)
                do iptr = 1, nlist
                  depends_on_idxlist(iptr,itarget) =
     &                 op2list(idxscr(iptr))
                  ! safety trap
                  if (depends_on_idxlist(iptr,itarget).le.0) then
                    error = .true.
                    write(lulog,*) 'no list for operator '//
     &                   trim(op_info%op_arr(idxscr(iptr))%op%name)
                  end if
                end do
              end if
            end if
            ! end of formula list? exit
            if (fl_ptr%command.eq.command_end_of_formula)
     &           exit
            ! else:
            ! increment target counter
            itarget = itarget+1
            ! and start from new
            iptr = 0
            iptr_last = 0
            sorted = .false.
            if (ipass.eq.2) then
! store targets (converted to list indices)
               if (ntest .ge.100)then
                  write (lulog,*) "found target",
     &                 op_info%op_arr(fl_ptr%target)%op%name
                  write (lulog,*) "found target list",
     &                 op_info%mel_arr(op2list(fl_ptr%target))
     &                 %mel%label
               end if 
              idxlist(itarget) = op2list(fl_ptr%target)
            end if
          case(command_add_contribution)
            nvtx = fl_ptr%contr%nvtx
            vertex => fl_ptr%contr%vertex
            sorted = sorted.and.nvtx.eq.0
            if (iptr+nvtx.gt.maxscr)
     &           call quit(1,'set_formula_dependencies',
     &           'increase maxscr')
            do ivtx = 1, nvtx
              iptr = iptr+1
              idxscr(iptr) = vertex(ivtx)%idx_op
            end do
          ! --> adaptieren fuer neue flists
          case(command_add_intm,command_cp_intm,
     &         command_bc,command_add_bc,
     &         command_bc_reo,command_add_bc_reo)
            nvtx = fl_ptr%bcontr%n_operands
            if (iptr+nvtx.gt.maxscr)
     &           call quit(1,'set_formula_dependencies',
     &           'increase maxscr')
            idx_op = idx_oplist2(fl_ptr%bcontr%label_op1,op_info)
            if (idx_op.gt.0) then
              sorted = .false.
              iptr = iptr+1
              idxscr(iptr) = idx_op
            end if
            if (nvtx.gt.1) then
              idx_op = idx_oplist2(fl_ptr%bcontr%label_op2,op_info)
              if (idx_op.gt.0) then
                sorted = .false.
                iptr = iptr+1
                idxscr(iptr) = idx_op                
              end if
            end if
          end select

          ! if we have a decent number of indices collected
          ! --> process then
          if (iptr-iptr_last.gt.128) then
            sorted = .true.
            nlist = iptr
            call unique_list(idxscr,nlist)
            iptr = nlist
            iptr_last = iptr
          end if

          if (.not.associated(fl_ptr%next)) exit
          fl_ptr => fl_ptr%next
        end do
        if (ipass.eq.1) then
          depend%ntargets = itarget
          depend%ndepend  = maxlist
          if (ntest.ge.100) write(lulog,*) 'ntargets, ndepend:',
     &         itarget,maxlist
          allocate(depend%idxlist(itarget),
     &             depend%depends_on_idxlist(maxlist,itarget))
          idxlist => depend%idxlist
          depends_on_idxlist => depend%depends_on_idxlist
          idxlist = 0
          depends_on_idxlist = 0
        end if

      end do
      if (ntest.ge.100.or.error) then
        write(lulog,*) 'targets: dependencies'
        idxlist => depend%idxlist
        depends_on_idxlist => depend%depends_on_idxlist
        do itarget = 1, depend%ntargets
          write(lulog,'(x,i3," :",20i4)') idxlist(itarget),
     &         depends_on_idxlist(1:depend%ndepend,itarget)
        end do
        if (error)
     &       call quit(1,'set_formula_dependencies','undefined list')
      end if

      return
      end
