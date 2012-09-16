      integer function pt_handle_rules(tgt_info,cur_target,wlist,mode)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'def_word_list.h'
      include 'par_parse_targets.h'
      include 'ifc_targets.h'

      character(len=15), parameter ::
     &    i_am = 'pt_handle_rules'
      
      type(target_info), intent(inout) ::
     &    tgt_info
      type(word_list), intent(inout) ::
     &    wlist
      integer, intent(inout) ::
     &    mode
      character(len=*), intent(inout) ::
     &    cur_target

      integer ::
     &    error, arg_type, idx_arg, idx_cmd, ndim, idim, nargs
      character ::
     &    sep
      character(len=256) ::
     &    command, label, word

      type(action), pointer ::
     &    cmd_pnt
      type(action_arg), pointer ::
     &    arg_pnt
       
      integer, pointer ::
     &    list_int(:)
      real(8), pointer ::
     &    list_rl8(:)
      logical, pointer ::
     &    list_log(:), need_arg(:)
      character(len=len_target_name), pointer ::
     &    list_label(:)

      logical, external ::
     &    advance_word_list_entry
      integer, external ::
     &    idx_action, idx_command_arg

      error = no_error

      ! get current entry again
      call get_word_list_entry(command,sep,wlist)

      ! compare with all the defined rules
      call uppcas(command)
      idx_cmd = idx_action(command,tgt_info)

      if (idx_cmd.le.0) then
        error = error_unknown_comkey
      else if (mode.ne.scan_for_rules) then
        error = error_here_no_rule
      else if (sep.ne.'(') then
        error = error_exp_bob
      end if

      pt_handle_rules = error
      if (error.ne.no_error) return

      ! if these tests are passed, we can register the rule
      call set_rule2(cur_target,command,tgt_info)

      cmd_pnt => tgt_info%act_array(idx_cmd)%act

      ! make an array to identify required arguments
      nargs = cmd_pnt%n_arguments
      allocate(need_arg(max(1,nargs)))
      do idx_arg = 1, nargs
        arg_pnt => cmd_pnt%arg(idx_arg)
        need_arg(idx_arg) = arg_pnt%required
      end do

      ! now process the list of arguments
      if (.not.advance_word_list_entry(wlist,'d')) then
        call quit(1,i_am,'unexpected absence of d-branch in wlist')
      end if

      scan_arg: do
        call get_word_list_entry(label,sep,wlist)
        ! accept a line-break
        if (len_trim(label).eq.0.and.sep.eq.'E') then
          if (.not.advance_word_list_entry(wlist,' ')) then
            pt_handle_rules = error_unexp_eob
            return
          end if
          call get_word_list_entry(label,sep,wlist)
        end if

        call uppcas(label)
        ! check whether the argument label exists
        ! and get the type of expected argument
        idx_arg = idx_command_arg(label,cmd_pnt)
        if (idx_arg.lt.0) then
          pt_handle_rules = error_unknown_arglab
        end if
        arg_pnt => cmd_pnt%arg(idx_arg)
        arg_type = arg_pnt%type 
        ! remember that argument is set 
        need_arg(idx_arg) = .false.

        if (sep.eq.')'.and.trim(label).eq.'') exit scan_arg

        if (sep.ne.'=') then
          pt_handle_rules = error_exp_eq
          return
        end if

        ! get next entry
        if (.not.advance_word_list_entry(wlist,' ')) then
          error = error_exp_arg
        end if

        call get_word_list_entry(word,sep,wlist)
        ! an array upcoming?
        if (trim(word).eq.'') then
          if (sep.ne.'(') then
            pt_handle_rules = error_exp_bob
            return
          end if

          ! move down and count entries
          if (.not.advance_word_list_entry(wlist,'d')) then
           call quit(1,i_am,'unexpected error when descending in wlist')
          end if

          ndim = 0
          cnt_lp: do
            call get_word_list_entry(word,sep,wlist)
            if (sep.ne.','.and.sep.ne.')') then
              pt_handle_rules = error_exp_copa
              return
            end if
            ndim = ndim+1
            if (.not.advance_word_list_entry(wlist,' ')) exit cnt_lp
          end do cnt_lp

          ! return to upper level
          if (.not.advance_word_list_entry(wlist,'u')) then
            call quit(1,i_am,'unexpected error when ascending in wlist')
          end if
          ! and go back
          if (.not.advance_word_list_entry(wlist,'d')) then
           call quit(1,i_am,'unexpected error when descending in wlist')
          end if
          ! (we now are back where we started!)

        else ! single argument
          ndim = 1
        end if

        ! allocate space for arguments (if needed)
        select case(arg_type)
        case(aatype_label)
          allocate(list_label(ndim))
        case(aatype_log)
          allocate(list_log(ndim))
        case(aatype_int)
          allocate(list_int(ndim))
        case(aatype_rl8)
          allocate(list_rl8(ndim))
        case(aatype_str)
          if (ndim.gt.1) then
            pt_handle_rules = error_no_str_arr
            return
          end if
        case default
          call quit(1,i_am,
     &              'unsupported argument type for '//trim(label))
        end select

        ! read in argument lists
        do idim = 1, ndim
          call get_word_list_entry(word,sep,wlist)

          select case(arg_type)
          case(aatype_label)
            list_label(idim) = word
          case(aatype_log)
            read(word,*,err=101) list_log(idim)
          case(aatype_int)
            read(word,*,err=102) list_int(idim)
          case(aatype_rl8)
            read(word,*,err=103) list_rl8(idim)
          case(aatype_str)
            ! use word directly 
          end select
          
          if (idim.lt.ndim) then
            if (.not.advance_word_list_entry(wlist,' ')) 
     &        call quit(1,i_am,
     &                  'unexpected end of list when re-scanning') 
          end if

        end do

        ! register
        select case(arg_type)
        case(aatype_label)
          call set_arg(cur_target,command,label,ndim,tgt_info,
     &         val_label=list_label)
          deallocate(list_label)
        case(aatype_log)
          call set_arg(cur_target,command,label,ndim,tgt_info,
     &         val_log=list_log)
          deallocate(list_log)
        case(aatype_int)
          call set_arg(cur_target,command,label,ndim,tgt_info,
     &         val_int=list_int)
          deallocate(list_int)
        case(aatype_rl8)
          call set_arg(cur_target,command,label,ndim,tgt_info,
     &         val_rl8=list_rl8)
          deallocate(list_rl8)
        case(aatype_str)
          call set_arg(cur_target,command,label,1,tgt_info,
     &         val_str=word)
        end select

        ! re-position the pointer on wlist
        if (ndim.gt.1) then
          if (.not.advance_word_list_entry(wlist,'u')) then
            call quit(1,i_am,'unexpected error while ascending')
          end if
        end if

        if (sep.ne.','.and.sep.ne.')') then
          pt_handle_rules = error_exp_copa
          return
        end if
        if (.not.advance_word_list_entry(wlist,' ')) exit scan_arg
      end do scan_arg

      ! check whether all required targets were set
      nargs = cmd_pnt%n_arguments
      do idx_arg = 1, nargs
        if (need_arg(idx_arg)) then
          arg_pnt => cmd_pnt%arg(idx_arg)
          write(luout,*)'argument "',trim(arg_pnt%arg_label),'" not set'
          error = error_arg_not_set
        end if
      end do

      deallocate(need_arg)

      pt_handle_rules = error
      if (error.ne.no_error) return 

      if (.not.advance_word_list_entry(wlist,'u')) then
        call quit(1,i_am,'unexpected absence of u-branch in wlist')
      end if

      if (.not.advance_word_list_entry(wlist,' ')) then
        pt_handle_rules = error_unexp_eob
        return
      end if

      pt_handle_rules = no_error

      return

  101 pt_handle_rules = error_read_log
      return
  102 pt_handle_rules = error_read_int
      return
  103 pt_handle_rules = error_read_rl8
      return

      end function
