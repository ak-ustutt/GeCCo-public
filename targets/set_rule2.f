*----------------------------------------------------------------------*
      subroutine set_rule2(name_target,command,tgt_info)
*----------------------------------------------------------------------*
*     set entry "command" to target referenced by "name_target" 
*     (must exist on tgt_info)
*     initialize arguments arrays
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'

      type(target_info), intent(inout), target ::
     &     tgt_info
      character, intent(in) ::
     &     name_target*(*), command*(*)

      integer ::
     &     idx_tgt, idx, type, idx_cmd, iarg, nupd, nargs
      type(target), pointer ::
     &     tgt
      type(action), pointer ::
     &     new_rules(:), last_rule
      integer, external ::
     &     idx_target, idx_action

      idx_tgt = idx_target(name_target,tgt_info)

      if (idx_tgt.le.0)
     &     call quit(1,'set_rule',
     &     'target not (yet) defined: '//trim(name_target))
      if (len_trim(command).gt.len_command_name)
     &     call quit(1,'set_rule',
     &     'name of rule too long: '//trim(command))

      tgt => tgt_info%array(idx_tgt)%tgt

      ! check for existing command name
      idx_cmd = idx_action(command,tgt_info)
      ! error, if not available
      if (idx_cmd.le.0)
     &     call quit(1,'set_rule',
     &     'command not (yet) defined: '//trim(command))

      allocate(new_rules(tgt%n_rules+1))
      if (tgt%n_rules.gt.0) then
        new_rules(1:tgt%n_rules) =
     &       tgt%rules(1:tgt%n_rules)
        deallocate(tgt%rules)
      end if
      last_rule => new_rules(tgt%n_rules+1)
      last_rule%command(1:len_command_name) = ' '
      last_rule%command = trim(command)
      ! this routine works for the new behaviour
      last_rule%new = .true.
      ! get from prototype rule
      nargs = tgt_info%act_array(idx_cmd)%act%n_arguments
      type  = tgt_info%act_array(idx_cmd)%act%type
      nupd  = tgt_info%act_array(idx_cmd)%act%n_update

      last_rule%n_arguments = nargs
      if (nargs.eq.0) then
        last_rule%arg => null()
      else
        allocate(last_rule%arg(nargs))
        do iarg = 1, nargs
          last_rule%arg(iarg)%type =
     &         tgt_info%act_array(idx_cmd)%act%arg(iarg)%type
          last_rule%arg(iarg)%arg_label =
     &         tgt_info%act_array(idx_cmd)%act%arg(iarg)%arg_label
          last_rule%arg(iarg)%arg_dim = 0
          last_rule%arg(iarg)%val_label => null()
          last_rule%arg(iarg)%val_log => null() 
          last_rule%arg(iarg)%val_int => null() 
          last_rule%arg(iarg)%val_occ => null() 
          last_rule%arg(iarg)%val_restr => null()
          last_rule%arg(iarg)%val_rl8  => null()
          last_rule%arg(iarg)%val_str  => null()
        end do
      end if

      ! set type and n_update as given in rule prototype
      last_rule%type     = type
      last_rule%n_update = nupd

      tgt%rules => new_rules
      tgt%n_rules = tgt%n_rules+1

      return
      end
