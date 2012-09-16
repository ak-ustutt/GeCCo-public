      integer function pt_handle_target(tgt_info,cur_target,wlist,mode)

      implicit none

      include 'mdef_target_info.h'
      include 'def_word_list.h'
      include 'par_parse_targets.h'

      character(len=16), parameter ::
     &    i_am = 'pt_handle_target'
      
      type(target_info), intent(inout) ::
     &    tgt_info
      type(word_list), intent(inout) ::
     &    wlist
      integer, intent(inout) ::
     &    mode
      character(len=*), intent(inout) ::
     &    cur_target

      integer ::
     &    error
      character ::
     &    sep
      character(len=256) ::
     &    label, word

      logical, external ::
     &    advance_word_list_entry

      error = no_error

      ! get separator after "target" keyword
      call get_word_list_entry(label,sep,wlist)
      if (mode.ne.scan_for_targets) then
        error = error_here_no_tgt
      ! only newline or white space my be separator
      else if (sep.ne.' '.and.sep.ne.'E') then
        error = error_unexp_sep
      end if
 
      if (error.ne.no_error) then
        pt_handle_target = error
        return
      end if

      ! get next entry (label expected here)
      if (.not.advance_word_list_entry(wlist,' ')) then
        pt_handle_target = error_unexp_eob
        return
      end if

      call get_word_list_entry(label,sep,wlist)

      ! check length
      if (len_trim(label).gt.len_target_name) then
        pt_handle_target = error_label_too_long
        return
      end if
     
      ! separator must be (
      if (sep.ne.'(') then      
        ! look into next line
        if (sep.eq.'E') then
          if (advance_word_list_entry(wlist,' ')) then
            call get_word_list_entry(word,sep,wlist)
            if (len_trim(word).ne.0.or.sep.ne.'(') then
              error = error_exp_bob
            end if
          else
            error = error_exp_bob
          end if
        else
          error = error_exp_bob
        end if
      end if

      if (error.ne.no_error) then
        pt_handle_target = error
        return
      end if

      ! now, that we made sure that there is a '(', we can branch down
      if (.not.advance_word_list_entry(wlist,'d')) then
        call quit(1,i_am,'unexpected absence of d-branch in wlist')
      end if

      call get_word_list_entry(word,sep,wlist)

      ! overread a line-break after the "("
      if (len_trim(word).eq.0.and.sep.eq.'E') then
        if (.not.advance_word_list_entry(wlist,' ')) then
          pt_handle_target = error_unexp_eob
          return
        end if
        call get_word_list_entry(word,sep,wlist)
      end if
      
      ! register target 
      if (trim(word).eq.'required') then
        call add_target2(trim(label),.true.,tgt_info)
        ! advance pointer once more
        if (sep.ne.'E'.and.sep.ne.';') then
          pt_handle_target = error_req_nlsc
          return
        end if
        if (.not.advance_word_list_entry(wlist,' ')) then
          pt_handle_target = error_unexp_eob
          return
        end if 
      else
        call add_target2(trim(label),.false.,tgt_info)
      end if

      cur_target = label
      mode = scan_for_rules

      end function
