*----------------------------------------------------------------------*
      subroutine set_rule(name_target,type,command,
     &                    labels,n_labels,n_update,
     &                    parameters,n_parameter_strings,
     &                    tgt_info)
*----------------------------------------------------------------------*
*     set entry "command" to target referenced by "name_target" 
*     (must exist on tgt_info)
*     n_labels is the number of target labels on "labels"
*     the first n_update of those will get update when the rule is
*       executed
*     n_parameter_strings is the number of strings "parameters"
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'

      type(target_info), intent(inout), target ::
     &     tgt_info
      integer, intent(in) ::
     &     n_parameter_strings, n_labels, n_update
      character, intent(in) ::
     &     name_target*(*), command*(*)
      character*(*) ::
     &     parameters(n_parameter_strings), labels(n_labels)

      integer ::
     &     idx_tgt, idx, type
      type(target), pointer ::
     &     tgt
      type(action), pointer ::
     &     new_rules(:), last_rule
      integer, external ::
     &     idx_target

      idx_tgt = idx_target(name_target,tgt_info)

      if (idx_tgt.le.0)
     &     call quit(1,'set_rule',
     &     'target not (yet) defined: '//trim(name_target))
      if (len_trim(command).gt.len_command_name)
     &     call quit(1,'set_rule',
     &     'name rule too long: '//trim(command))

      tgt => tgt_info%array(idx_tgt)%tgt

      allocate(new_rules(tgt%n_rules+1))
      if (tgt%n_rules.gt.0) then
        new_rules(1:tgt%n_rules) =
     &       tgt%rules(1:tgt%n_rules)
        deallocate(tgt%rules)
      end if
      last_rule => new_rules(tgt%n_rules+1)
      last_rule%command(1:len_command_name) = ' '
      last_rule%command = trim(command)
      ! this routine works for the old behaviour
      last_rule%new = .false.

      allocate(last_rule%labels(n_labels))
      allocate(last_rule%parameters(n_parameter_strings))

      do idx = 1, n_labels
        last_rule%labels(idx) = labels(idx)
      end do
      do idx = 1, n_parameter_strings
        last_rule%parameters(idx) = parameters(idx)
      end do
      last_rule%type     = type
      last_rule%n_update = n_update
      last_rule%n_labels = n_labels
      last_rule%n_parameter_strings = n_parameter_strings

      tgt%rules => new_rules
      tgt%n_rules = tgt%n_rules+1

      return
      end
