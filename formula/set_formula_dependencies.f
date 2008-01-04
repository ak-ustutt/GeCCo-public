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
     &     sorted
      integer ::
     &     ipass, iptr, iptr_last, itarget, nlist, maxlist,
     &     nvtx, ivtx
      integer ::
     &     idxscr(maxscr)

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     fl_ptr
      integer, pointer ::
     &     op2list(:), idxlist(:), depends_on_idxlist(:,:)


      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &       'this is set_formula_dependencies')
      end if

      op2list => op_info%op2list 
      do ipass = 1, 2

        if (ntest.ge.100) write(luout,*) 'pass: ',ipass

        fl_ptr => flist
        iptr = 0
        iptr_last = 0
        itarget = 0
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
            ! if we have a decent number of indices collected
            ! --> process then
            if (iptr-iptr_last.gt.128) then
              sorted = .true.
              nlist = iptr
              call unique_list(idxscr,nlist)
              iptr = nlist
              iptr_last = iptr
            end if

          end select

          if (.not.associated(fl_ptr%next)) exit
          fl_ptr => fl_ptr%next

        end do

        if (ipass.eq.1) then
          depend%ntargets = itarget
          depend%ndepend  = maxlist
          if (ntest.ge.100) write(luout,*) 'ntargets, ndepend:',
     &         itarget,maxlist
          allocate(depend%idxlist(itarget),
     &             depend%depends_on_idxlist(maxlist,itarget))
          idxlist => depend%idxlist
          depends_on_idxlist => depend%depends_on_idxlist
          idxlist = 0
          depends_on_idxlist = 0
        end if

      end do

      if (ntest.ge.100) then
        write(luout,*) 'targets: dependencies'
        idxlist => depend%idxlist
        depends_on_idxlist => depend%depends_on_idxlist
        do itarget = 1, depend%ntargets
          write(luout,'(x,i3," :",20i4)') idxlist(itarget),
     &         depends_on_idxlist(1:depend%ndepend,itarget)
        end do
        
      end if

      return
      end
