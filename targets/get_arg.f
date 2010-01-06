*----------------------------------------------------------------------*
      subroutine get_arg(arg_label,rule,tgt_info,
     &     val_label,val_label_list,val_log,val_log_list,
     &     val_int,val_int_list,val_occ,val_restr,
     &     val_rl8,val_rl8_list,val_str,ndim)
*----------------------------------------------------------------------*
*     get argument value from rule structure
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_target_info.h'

      type(action), intent(in), target ::
     &     rule
      type(target_info), intent(in) ::
     &     tgt_info
      character(len=*), intent(in) ::
     &     arg_label
      character(len=*), intent(out), optional ::
     &     val_label
      character(len=*), intent(out), optional ::
     &     val_label_list(:)
      logical, intent(out), optional ::
     &     val_log
      logical, intent(out), optional ::
     &     val_log_list(:)
      integer, intent(out), optional ::
     &     val_int
      integer, intent(out), optional ::
     &     val_int_list(:)
      integer, intent(out), optional ::
     &     val_occ(:,:,:)
      integer, intent(out), optional ::
     &     val_restr(:,:,:,:,:,:) 
      real(8), intent(out), optional ::
     &     val_rl8
      real(8), intent(out), optional ::
     &     val_rl8_list(:)
      character(len=*), intent(out), optional ::
     &     val_str
      integer, intent(out), optional ::
     &     ndim 

      integer ::
     &     idx_arg, idx_cmd, arg_dim, arg_type, ntype,
     &     idx_str_batch, n_str_batch, idxst, idxnd
      type(action_arg), pointer ::
     &     arg
      type(action), pointer ::
     &     proto_rule

      integer, external ::
     &     idx_command_arg, idx_action
     

      idx_arg = idx_command_arg(arg_label,rule)

      if (idx_arg.le.0) then
        ! try to find default value
        idx_cmd = idx_action(rule%command,tgt_info)

        if (idx_cmd.gt.0) then
          proto_rule => tgt_info%act_array(idx_cmd)%act 
          idx_arg = idx_command_arg(arg_label,proto_rule)
        else
          idx_arg = -1
        end if

        if (idx_arg.le.0)
     &     call quit(1,'get_arg',
     &     'argument is undefined: '//
     &     trim(rule%command)//':'//trim(arg_label))

        ! default
        arg => proto_rule%arg(idx_arg)
        
      else

        ! user value
        arg => rule%arg(idx_arg)

      end if

      arg_type = arg%type
      arg_dim  = arg%arg_dim

      if (present(ndim)) ndim = arg_dim

      select case(arg_type)
      case(aatype_label)
        if (present(val_label_list)) then
          if (size(val_label_list).lt.arg_dim)
     &         call quit(1,'get_arg','dimension too small (1)')
          val_label_list(1:arg_dim) =
     &         arg%val_label(1:arg_dim)
        else if (present(val_label).and.arg_dim.eq.1) then
          val_label =
     &         arg%val_label(1)
        else
          call quit(1,'get_arg','no output variable present (1)')          
        end if
      case(aatype_log)
        if (present(val_log_list)) then
          if (size(val_log_list).lt.arg_dim)
     &         call quit(1,'get_arg','dimension too small (2)')
          val_log_list(1:arg_dim) =
     &         arg%val_log(1:arg_dim)
        else if (present(val_log).and.arg_dim.eq.1) then
          val_log =
     &         arg%val_log(1)
        else
          call quit(1,'get_arg','no output variable present (2)')          
        end if
      case(aatype_int)
        if (present(val_int_list)) then
          if (size(val_int_list).lt.arg_dim)
     &         call quit(1,'get_arg','dimension too small (1)')
          val_int_list(1:arg_dim) =
     &         arg%val_int(1:arg_dim)
        else if (present(val_int).and.arg_dim.eq.1) then
          val_int =
     &         arg%val_int(1)
        else
          call quit(1,'get_arg','no output variable present (1)')          
        end if
      case(aatype_occ)
        if (.not.present(val_occ))
     &       call quit(1,'get_arg','val_occ not present')
        if (size(val_occ).lt.ngastp*2*arg_dim)
     &       call quit(1,'get_arg','dimension too small')
        val_occ(1:ngastp,1:2,1:arg_dim) =
     &       arg%val_occ(1:ngastp,1:2,1:arg_dim)
      case(aatype_restr)
        if (.not.present(val_restr))
     &       call quit(1,'get_arg','val_restr not present')
        if (size(val_restr).lt.size(arg%val_restr))
     &       call quit(1,'get_arg','dimension too small')
        val_restr(:,:,:,:,:,1:arg_dim) =
     &       arg%val_restr(:,:,:,:,:,1:arg_dim)
      case(aatype_rl8)
        if (present(val_rl8_list)) then
          if (size(val_rl8_list).lt.arg_dim)
     &         call quit(1,'get_arg','dimension too small (1)')
          val_rl8_list(1:arg_dim) =
     &         arg%val_rl8(1:arg_dim)
        else if (present(val_rl8).and.arg_dim.eq.1) then
          val_rl8 =
     &         arg%val_rl8(1)
        else
          call quit(1,'get_arg','no output variable present (1)')          
        end if
      case(aatype_str)
        if (.not.present(val_str))
     &       call quit(1,'get_arg','val_str not present')

        val_str(1:len(val_str)) = ' '
        idxst = 1
        do idx_str_batch = 1, arg%n_str_batch
          idxnd = idxst-1 + len(arg%val_str(idx_str_batch))
          if (idxnd.gt.len(val_str))
     &         call quit(1,'get_arg','dimension too small (7)')
          val_str(idxst:idxnd) = arg%val_str(idx_str_batch)
          idxst = idxnd+1
        end do

      end select

      return
      end
