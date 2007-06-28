*----------------------------------------------------------------------*
      subroutine print_action_list(act_list)
*----------------------------------------------------------------------*
*     print the elements on list
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_action.h'
      include 'def_action_list.h'

      type(action_list), intent(in), target ::
     &     act_list

      type(action_list), pointer ::
     &     current
      integer ::
     &     nel

      ! point to first element
      current => act_list
      ! loop over elements
      nel = 0
      do while (associated(current%act))
      
        nel = nel+1
        write(luout,*) 'Action #',nel

        select case (current%act%action_type)
          case (iaction_import)
            write(luout,*) 'Import'
          case (iaction_evaluate)
            write(luout,*) 'Evaluate'
          case (iaction_setup_prc)
            write(luout,*) 'Setup diagonal preconditioner'
          case (iaction_solve_leq)
            write(luout,*) 'Solve system of linear equations'
            write(luout,*) ' # sets: ',current%act%nop_opt
          case (iaction_solve_nleq)
            write(luout,*) 'Solve system of non-linear equations'
            write(luout,*) ' # sets: ',current%act%nop_opt
          case (iaction_solve_evp)
            write(luout,*) 'Solve eigenvalue problem'
            write(luout,*) ' # sets: ',current%act%nop_opt
          case (iaction_solve_gevp)
            write(luout,*) 'Solve general eigenvalue problem'
            write(luout,*) ' # sets: ',current%act%nop_opt
          case default
            write(luout,*) 'action = ',current%act%action_type
            call quit(0,'print_action_list','unknown action')
        end select

        if (current%act%nop_in.gt.0) then
          write(luout,'(x,a,6(i4,x),/,10x,6(i4,x))')
     &         'Input operators (def.):      ',
     &         current%act%idxopdef_in(1:current%act%nop_in)
          write(luout,'(x,a,6(i4,i4,x),/,10x,6(i4,i4,x))')
     &         'Input operators (file rec.): ',
     &         current%act%idxopfile_in(1:2,1:current%act%nop_in)
        end if
        if (current%act%nop_out.gt.0) then
          write(luout,'(x,a,6(i4,x),/,10x,6(i4,x))')
     &         'Output operators (def.):     ',
     &         current%act%idxopdef_out(1:current%act%nop_out)
          write(luout,'(x,a,6(i4,i4,x),/,10x,6(i4,i4,x))')
     &         'Output operators (file rec.):',
     &         current%act%idxopfile_out(1:2,1:current%act%nop_out)
        end if
        if (current%act%nform.gt.0) then
          write(luout,'(x,a,6i4)') 'Formula files: ',
     &         current%act%idx_formula(1:current%act%nform)
        end if

        if (associated(current%next)) then
          current => current%next
        else
          exit
        end if

      end do

      return
      end 
