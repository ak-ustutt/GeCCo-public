*----------------------------------------------------------------------*
      subroutine add_action(act_list,nactions,
     &     action_type,nop_in,nop_out,
     &     idxopdef_in,idxopdef_out,
     &     idxopfile_in,idxopfile_out,
     &     nform,idx_formula)
*----------------------------------------------------------------------*
*     an 'action' is given by the parameters on the input list
*     compare to elements on act_list, and if the particular action
*     is not yet on the list, add it as last element
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ifc_baserout.h'
      include 'def_filinf.h'
      include 'def_action.h'
      include 'def_action_list.h'

      type(action_list), intent(inout), target ::
     &     act_list
      integer, intent(inout) ::
     &     nactions

      integer, intent(in) ::
     &     action_type
      integer, intent(in) ::
     &     nop_in,nop_out,
     &     idxopdef_in(*),idxopdef_out(*),
     &     idxopfile_in(2,*),idxopfile_out(2,*),
     &     nform,
     &     idx_formula(*)

      type(action_list), pointer ::
     &     current
      logical ::
     &     unique, equal

      unique = .true.
      ! point to first element
      current => act_list
      ! loop over elements and compare 
      do while (associated(current%act))
      
        if (action_type.eq.current%act%action_type) then
          equal = .true.
          if (nop_in.gt.0) then
            equal = equal.and.current%act%nop_in.eq.nop_in
            equal = equal.and.
     &          list_cmp(current%act%idxopdef_in,idxopdef_in,nop_in)
            equal = equal.and.
     &          list_cmp(current%act%idxopfile_in,idxopfile_in,2*nop_in)
          else
            equal = equal.and.current%act%nop_in.eq.0
          end if
          if (nop_out.gt.0) then
            equal = equal.and.current%act%nop_out.eq.nop_out
            equal = equal.and.
     &          list_cmp(current%act%idxopdef_out,idxopdef_out,nop_out)
            equal = equal.and.
     &          list_cmp(current%act%idxopfile_out,
     &                                         idxopfile_out,2*nop_out)
          else
            equal = equal.and.current%act%nop_out.eq.0
          end if
          if (nform.gt.0) then
            equal = equal.and.current%act%nform.eq.nform
            equal = equal.and.
     &          list_cmp(current%act%idx_formula,idx_formula,nform)
          else
            equal = equal.and.current%act%nform.eq.0
          end if

        end if

        unique = .not.equal
        if (equal) exit

        if (associated(current%next)) then
          current => current%next
        else
          exit
        end if

      end do             

      ! if new element is unique:
      if (unique) then
        ! add as last element
        nactions = nactions+1
        if (associated(current%act)) then
          allocate(current%next)
          current%next%prev => current
          current => current%next
          nullify(current%next)
        end if

        allocate(current%act)
        if (nop_in.gt.0)
     &       allocate(current%act%idxopdef_in(nop_in),
     &                current%act%idxopfile_in(2,nop_in))
        if (nop_out.gt.0)
     &       allocate(current%act%idxopdef_out(nop_out),
     &                current%act%idxopfile_out(2,nop_out))
        if (nform.gt.0)
     &       allocate(current%act%idx_formula(nform))

        current%act%action_type = action_type
        
        current%act%nop_in = 0
        if (nop_in.gt.0) then
          current%act%nop_in = nop_in
          current%act%idxopdef_in(1:nop_in) =
     &                idxopdef_in(1:nop_in)
          current%act%idxopfile_in(1:2,1:nop_in) =
     &                idxopfile_in(1:2,1:nop_in)
        end if

        current%act%nop_out = 0
        if (nop_out.gt.0) then
          current%act%nop_out = nop_out
          current%act%idxopdef_out(1:nop_out) =
     &                idxopdef_out(1:nop_out)
          current%act%idxopfile_out(1:2,1:nop_out) =
     &                idxopfile_out(1:2,1:nop_out)
        end if

        current%act%nform = 0
        if (nform.gt.0) then
          current%act%nform = nform
          current%act%idx_formula(1:nform) = idx_formula(1:nform)
        end if

      end if

      return
      end 
*----------------------------------------------------------------------*
*     inactivated version with optional arguments:
*     does not work with present intel compiler (up to v9.1)
*----------------------------------------------------------------------*
      subroutine add_action_with_optional_arguments(act_list,nactions,
     &     action_type,nop_in,nop_out,
     &     idxopdef_in,idxopdef_out,
     &     idxopfile_in,idxopfile_out,
     &     nform,
     &     idx_formula)
*----------------------------------------------------------------------*
*     an 'action' is given by the parameters on the input list
*     compare to elements on act_list, and if the particular action
*     is not yet on the list, add it as last element
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ifc_baserout.h'
      include 'def_filinf.h'
      include 'def_action.h'
      include 'def_action_list.h'

      type(action_list), intent(inout), target ::
     &     act_list
      integer, intent(inout) ::
     &     nactions

      integer, intent(in) ::
     &     action_type
      integer, intent(in), optional ::
     &     nop_in,nop_out,
     &     idxopdef_in(*),idxopdef_out(*),
     &     idxopfile_in(2,*),idxopfile_out(2,*),
     &     nform,
     &     idx_formula(*)

      type(action_list), pointer ::
     &     current
      logical ::
     &     unique, equal

      ! consistency checks
      if (present(nop_in) .and.
     &    (.not.present(idxopdef_in) .or.
     &     .not.present(idxopfile_in))) then
        write(luout,*) 'in:',present(idxopdef_in),present(idxopfile_in) 
        call quit(0,'add_action','inconsisten call list (nop_in)')
      end if
      if (present(nop_out) .and.
     &    (.not.present(idxopdef_out) .or.
     &     .not.present(idxopfile_out))) then
        write(luout,*) 'out:',present(idxopdef_out),
     &                        present(idxopfile_out) 
        call quit(0,'add_action','inconsisten call list (nop_out)')
      end if
      if (present(nform) .and.
     &    (.not.present(idx_formula))) then
        call quit(0,'add_action','inconsisten call list (nform)')
      end if
      

      unique = .true.
      ! point to first element
      current => act_list
      ! loop over elements and compare 
      do while (associated(current%act))
      
        if (action_type.eq.current%act%action_type) then
          equal = .true.
          if (present(nop_in)) then
            equal = equal.and.current%act%nop_in.eq.nop_in
            equal = equal.and.
     &          list_cmp(current%act%idxopdef_in,idxopdef_in,nop_in)
            equal = equal.and.
     &          list_cmp(current%act%idxopfile_in,idxopfile_in,2*nop_in)
          else
            equal = equal.and.current%act%nop_in.eq.0
          end if
          if (present(nop_out)) then
            equal = equal.and.current%act%nop_out.eq.nop_out
            equal = equal.and.
     &          list_cmp(current%act%idxopdef_out,idxopdef_out,nop_out)
            equal = equal.and.
     &          list_cmp(current%act%idxopfile_out,
     &                                         idxopfile_out,2*nop_out)
          else
            equal = equal.and.current%act%nop_out.eq.0
          end if
          if (present(nform)) then
            equal = equal.and.current%act%nform.eq.nform
            equal = equal.and.
     &          list_cmp(current%act%idx_formula,idx_formula,nform)
          else
            equal = equal.and.current%act%nform.eq.0
          end if

        end if

        unique = .not.equal
        if (equal) exit

        if (associated(current%next)) then
          current => current%next
        else
          exit
        end if

      end do             

      ! if new element is unique:
      if (unique) then
        ! add as last element
        nactions = nactions+1
        if (associated(current%act)) then
          allocate(current%next)
          current%next%prev => current
          current => current%next
          nullify(current%next)
        end if

        allocate(current%act)
        if (present(nop_in))
     &       allocate(current%act%idxopdef_in(nop_in),
     &                current%act%idxopfile_in(2,nop_in))
        if (present(nop_out))
     &       allocate(current%act%idxopdef_out(nop_out),
     &                current%act%idxopfile_out(2,nop_out))
        if (present(nform))
     &       allocate(current%act%idx_formula(nform))

        current%act%action_type = action_type
        
        current%act%nop_in = 0
        if (present(nop_in)) then
          current%act%nop_in = nop_in
          current%act%idxopdef_in(1:nop_in) =
     &                idxopdef_in(1:nop_in)
          current%act%idxopfile_in(1:2,1:nop_in) =
     &                idxopfile_in(1:2,1:nop_in)
        end if

        current%act%nop_out = 0
        if (present(nop_out)) then
          current%act%nop_out = nop_out
          current%act%idxopdef_out(1:nop_out) =
     &                idxopdef_out(1:nop_out)
          current%act%idxopfile_out(1:2,1:nop_out) =
     &                idxopfile_out(1:2,1:nop_out)
        end if

        current%act%nform = 0
        if (present(nform)) then
          current%act%nform = nform
          current%act%idx_formula(1:nform) = idx_formula(1:nform)
        end if

      end if

      return
      end 
