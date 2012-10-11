      subroutine parse_targets(tgt_info,wlist)
 
      implicit none
      
      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'def_word_list.h'
      include 'par_parse_targets.h'

      integer, parameter ::
     &     ntest = 00
      character(len=13), parameter ::
     &     i_am = 'parse_targets'

      type(target_info), intent(inout) ::
     &    tgt_info
      type(word_list) ::
     &    wlist

      character(len=256) ::
     &     word, cur_target, cur_rule
      character :: 
     &     sep
      integer ::
     &     mode, error, line, col


      logical, external ::
     &     advance_word_list_entry, 
     &     pt_handle_target, pt_handle_depend, pt_handle_rules

      if (ntest.gt.0) call write_title(luout,wst_dbg_subr,i_am)

      ! safty exit
      if (.not.associated(wlist%head)) 
     &      call quit(1,i_am,'buggy wlist')

      if (ntest.ge.100) then
        write(luout,*) 'word list:'
        call print_word_list(luout,wlist)
      end if

      ! set "current" pointer to head
      wlist%current => wlist%head

      mode = scan_for_targets
      error = no_error
      cur_target = '- '
      cur_rule = '- '
      main_loop: do

        ! read present entry
        call get_word_list_entry(word,sep,wlist)

        if (ntest.ge.100) then
          write(luout,*) 'word = ',trim(word),'   sep = ',sep
        end if

        if (mode.ne.scan_for_arguments) then
          select case(trim(word))
            ! nothing
            case('')
              ! only newline or close_block is an acceptable separator
              if (sep.ne.'E'.and.sep.ne.')') then
                error = error_unexp_sep
                exit main_loop
              end if
              if (sep.eq.')') then
                if (mode.ne.scan_for_rules) then
                  error = error_unexp_eob
                  exit main_loop
                end if
                ! current target definition has ended
                if (.not.advance_word_list_entry(wlist,'u')) then
                  call quit(1,i_am,'unexpected error for u-branch')
                end if
                mode = scan_for_targets
                cur_target = '- '
              end if
            ! a keyword
            !case('if')
            case('target')
              if (ntest.ge.100) write(luout,*) 'handle case: target'
              error = pt_handle_target(tgt_info,cur_target,wlist,mode)
              if (error.ne.no_error) exit main_loop
              cycle main_loop
            case('depend')
              if (ntest.ge.100) write(luout,*) 'handle case: depend'
              error = pt_handle_depend(tgt_info,cur_target,wlist,mode)
              if (error.ne.no_error) exit main_loop
              cycle main_loop
            case default
              ! check for rules 
              if (ntest.ge.100) write(luout,*) 'handle case: rules'
              error = pt_handle_rules(tgt_info,cur_target,wlist,mode)
              if (error.ne.no_error) exit main_loop
              cycle main_loop
          end select

        else

        end if

        if (.not.advance_word_list_entry(wlist,' ')) exit main_loop

      end do main_loop

      if (error.eq.no_error) return

      ! error handling:
      call get_word_list_entry(word,sep,wlist)
      call get_word_list_position(line,col,wlist)
      write(luout,'(/,x,a,i4,a,i3)') 
     &       'Input error near line ',line,' column ',col
      
      select case(error) 
      case(error_unexp_sep)
        write(luout,*) 'Unexpected separator "'//sep//'"'
      case(error_exp_bob)
        write(luout,*) 'I expected the beginning of a block: ('
      case(error_unexp_eob)
        write(luout,*) 'Unexpected end of block'
      case(error_here_no_tgt)
        write(luout,*) 'You cannot define a target in this context'
      case(error_here_no_dpd)
        write(luout,*) 'You cannot define a dependency in this context'
      case(error_req_nlsc)
        write(luout,*) 'I expected one of <newline> or ; '
      case(error_exp_copa)
        write(luout,*) 'I expected one of , or ) '
      case(error_label_too_long)
        write(luout,*) 'The label is too long'
      case(error_unknown_comkey)
        write(luout,*) 'Unknown command or keyword: "'//trim(word)//'"'
      case(error_here_no_rule)
        write(luout,*) 'You cannot define a rule in this context'
      case(error_unknown_arglab)
        write(luout,*) 'Unknown argument label: "'//trim(word)//'"'
      case(error_exp_eq)
        write(luout,*) 'Expected = sign'
      case(error_exp_arg)
        write(luout,*) 'Expected value for argument here'
      case(error_no_str_arr)
        write(luout,*) 'A string argument cannot be an array'
      case(error_read_log)
        write(luout,*) 'Error reading logical value from: "'
     &                  //trim(word)//'"'
      case(error_read_int)
        write(luout,*) 'Error reading integer value from: "'
     &                  //trim(word)//'"'
      case(error_read_rl8)
        write(luout,*) 'Error reading real value from: "'
     &                  //trim(word)//'"'
      case(error_arg_not_set)
        write(luout,*) 'One or more required arguments are missing!'
      case default
        write(luout,*) 'error code = ',error
      end select
       
      call quit(0,i_am,'Input error')
      end

