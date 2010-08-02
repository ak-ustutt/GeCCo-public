*----------------------------------------------------------------------*
      subroutine set_arg(name_target,command,arg_label,arg_dim,tgt_info,
     &     val_label,val_log,val_int,val_occ,val_restr,val_rl8,val_str)
*----------------------------------------------------------------------*
*     add a new argument to (latest) rule with label "command"
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_target_info.h'

      type(target_info), intent(inout), target ::
     &     tgt_info
      character(len=*), intent(in) ::
     &     name_target, command, arg_label
      integer, intent(in) ::
     &     arg_dim
      character(len=*), intent(in), optional ::
     &     val_label(arg_dim)
      logical, intent(in), optional ::
     &     val_log(arg_dim)
      integer, intent(in), optional ::
     &     val_int(arg_dim)
      integer, intent(in), optional ::
     &     val_occ(ngastp,2,arg_dim)
      integer, intent(in), optional ::
     &     val_restr(:,:,:,:,:,:) 
      real(8), intent(in), optional ::
     &     val_rl8(arg_dim)
      character(len=*), intent(in), optional ::
     &     val_str

      integer ::
     &     idx_tgt, idx_cmd, idx_arg,
     &     arg_type, ntype,
     &     idx_str_batch, n_str_batch, idxst, idxnd,
     &     ndim(6)
      logical ::
     &     proto_mode
      type(target), pointer ::
     &     tgt
      type(action), pointer ::
     &     rule
      type(action_arg), pointer ::
     &     new_arg(:), last_arg

      integer, external ::
     &     idx_target, idx_target_command, idx_command_arg, idx_action
     
      proto_mode = .false.
      if (trim(name_target).eq.'_PROTO_') proto_mode = .true.

      if (.not.proto_mode) then
        ! look for target/rule in target list

        idx_tgt = idx_target(name_target,tgt_info)

        if (idx_tgt.le.0)
     &       call quit(1,'set_arg',
     &       'target not (yet) defined: '//trim(name_target))

        tgt => tgt_info%array(idx_tgt)%tgt
      
        ! get latest occurrence of rule with label "command"
        idx_cmd = idx_target_command(-1,command,tgt)

        if (idx_cmd.le.0)
     &     call quit(1,'set_arg',
     &     'rule not (yet) defined: '//trim(command))

        rule => tgt%rules(idx_cmd)

      else

        ! find prototype command
        idx_cmd = idx_action(command,tgt_info)

        rule => tgt_info%act_array(idx_cmd)%act

      end if

      idx_arg = idx_command_arg(arg_label,rule)

      if (idx_arg.gt.0)
     &     call quit(1,'set_arg',
     &     'argument is already defined: '//trim(arg_label))

      allocate(new_arg(rule%n_arguments+1))
      if (rule%n_arguments.gt.0) then
        new_arg(1:rule%n_arguments) =
     &       rule%arg(1:rule%n_arguments)
        deallocate(rule%arg)
      end if

      last_arg => new_arg(rule%n_arguments+1)
      last_arg%arg_label(1:len_target_name) = ' '
      last_arg%arg_label = trim(arg_label)
      last_arg%arg_dim = arg_dim
      last_arg%n_str_batch = 0

      ! determine type from given argument
      ntype = 0 ! should be 1 at the end of this section
      if (present(val_label)) then
        arg_type = aatype_label
        ntype = ntype+1
      end if
      if (present(val_log)) then
        arg_type = aatype_log
        ntype = ntype+1
      end if
      if (present(val_int)) then
        arg_type = aatype_int
        ntype = ntype+1
      end if
      if (present(val_occ)) then
        arg_type = aatype_occ
        ntype = ntype+1
      end if
      if (present(val_restr)) then
        arg_type = aatype_restr
        ntype = ntype+1
      end if
      if (present(val_rl8)) then
        arg_type = aatype_rl8
        ntype = ntype+1
      end if
      if (present(val_str)) then
        arg_type = aatype_str
        ntype = ntype+1
      end if
      
      ! error conditions
      if (ntype.eq.0)
     &     call quit(1,'set_arg','no optional argument given')
      if (ntype.gt.1)
     &     call quit(1,'set_arg',
     &                'more than one optional argument given')

      last_arg%type = arg_type

      select case(arg_type)
      case(aatype_label)
        allocate(last_arg%val_label(arg_dim))
        last_arg%val_label(1:arg_dim) =
     &           val_label(1:arg_dim)
      case(aatype_log)
        allocate(last_arg%val_log(arg_dim))
        last_arg%val_log(1:arg_dim) =
     &           val_log(1:arg_dim)
      case(aatype_int)
        allocate(last_arg%val_int(arg_dim))
        last_arg%val_int(1:arg_dim) =
     &           val_int(1:arg_dim)
      case(aatype_occ)
        allocate(last_arg%val_occ(ngastp,2,arg_dim))
        last_arg%val_occ(1:ngastp,1:2,1:arg_dim) =
     &           val_occ(1:ngastp,1:2,1:arg_dim)
      case(aatype_restr)
        ndim(1:6) = shape(val_restr)
        allocate(last_arg%val_restr(ndim(1),ndim(2),ndim(3),
     &                              ndim(4),ndim(5),arg_dim))
        last_arg%val_restr(:,:,:,:,:,1:arg_dim) =
     &           val_restr(:,:,:,:,:,1:arg_dim)
      case(aatype_rl8)
        allocate(last_arg%val_rl8(arg_dim))
        last_arg%val_rl8(1:arg_dim) =
     &           val_rl8(1:arg_dim)
      case(aatype_str)
        if (arg_dim.gt.1)
     &       call quit(1,'set_arg',
     &       'arg_dim can only be 1 for strings ('//
     &       trim(name_target)//','//trim(command)//
     &       ','//trim(arg_label)//')')
        n_str_batch = len_trim(val_str)/len_str_batch + 1
        allocate(last_arg%val_str(n_str_batch))
        idxst = 1
        do idx_str_batch = 1, n_str_batch
          idxnd = min(len_trim(val_str),idxst-1 + len_str_batch)
          last_arg%val_str(idx_str_batch)(1:len_str_batch) = ' '
          last_arg%val_str(idx_str_batch) = val_str(idxst:idxnd)
          idxst = idxnd+1
        end do
        last_arg%n_str_batch = n_str_batch
      end select

      rule%arg => new_arg
      rule%n_arguments = rule%n_arguments+1

      return
      end
