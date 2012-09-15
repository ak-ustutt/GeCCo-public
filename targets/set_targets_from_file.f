      subroutine set_targets_from_file(tgt_info,input_name,
     &                              orb_info,env_type)

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'
      include 'def_word_list.h'

      character(len=13), parameter ::
     &     key_start = 'begin targets',
     &     key_end   = 'end targets  '

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      character(len=*), intent(in) ::
     &     input_name, env_type
      

      type(filinf) :: 
     &    input_file
      type(word_list) ::
     &    input_list, blocked_list

      logical, external ::
     &    file_exists

      call file_init(input_file,input_name,ftyp_sq_frm,0)

      if (.not.file_exists(input_file)) then
        call quit(0,'set_targets_from_file','file not found: "'//
     &              trim(input_name)//'"')     
      end if

      write(luout,*) 'reading file ',trim(input_name)

      call init_word_list(input_list)
      call lex_file(input_list,input_file,key_start,key_end)

      call init_word_list(blocked_list)
      call parse_blocks(blocked_list,input_list,luout)

      call clean_word_list(input_list)

      call parse_targets(tgt_info,blocked_list)
      !call print_word_list(luout,blocked_list)

      call clean_word_list(blocked_list)

      end

