      integer function pt_handle_depend(tgt_info,cur_target,wlist,mode)

      implicit none

      include 'mdef_target_info.h'
      include 'def_word_list.h'
      include 'par_parse_targets.h'

      character(len=16), parameter ::
     &    i_am = 'pt_handle_depend'
      
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
     &    label
      
      logical, external ::
     &    advance_word_list_entry

      error = no_error

      ! get separator after "depend" keyword
      call get_word_list_entry(label,sep,wlist)
      if (mode.ne.scan_for_rules) then
        error = error_here_no_dpd
      ! only white space or "(" may be separator
      !  "(" opens a block with labels
      else if (sep.ne.' '.and.sep.ne.'(') then
        error = error_unexp_sep
      end if
 
      if (error.ne.no_error) then
        pt_handle_depend = error
        return
      end if

      ! single depend
      if (sep.eq.' ') then
        ! get next entry (label expected here)
        if (.not.advance_word_list_entry(wlist,' ')) then
          pt_handle_depend = error_unexp_eob
          return
        end if
        ! get contents of new entry
        call get_word_list_entry(label,sep,wlist)
        ! register
        call set_dependency(cur_target,trim(label),tgt_info)
      else ! a list of depends
        if (.not.advance_word_list_entry(wlist,'d')) then
          call quit(1,i_am,'unexpected absence of d-branch in wlist')
        end if
        scan_dep: do
          call get_word_list_entry(label,sep,wlist)
          ! register
          call set_dependency(cur_target,trim(label),tgt_info)
          if (sep.ne.','.and.sep.ne.')') then
            pt_handle_depend = error_exp_copa
            return
          end if
          if (.not.advance_word_list_entry(wlist,' ')) exit scan_dep
        end do scan_dep
        if (.not.advance_word_list_entry(wlist,'u')) then
          call quit(1,i_am,'unexpected absence of u-branch in wlist')
        end if
  
      end if

      if (.not.advance_word_list_entry(wlist,' ')) then
        pt_handle_depend = error_unexp_eob
        return
      end if

      pt_handle_depend = no_error
      return

      end function
