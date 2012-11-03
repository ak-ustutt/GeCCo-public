      subroutine add_target3(string_arr,tgt_info)
*
*     most convenient version: use the parser to set a new target
*     the same syntax as for input files can be used, see documentation
*
*     andreas, october 2012
*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'
      include 'def_word_list.h'
      include 'par_special_chars.h'

      type(target_info), intent(inout) ::
     &     tgt_info
c      integer, intent(in) ::
c     &     nlines
      character(len=*), intent(in), dimension(:) ::
     &     string_arr
      
      integer ::
     &     ii
      type(word_list) ::
     &    input_list, blocked_list

c      if (.not.present(string_arr)) return

      call init_word_list(input_list)

      do ii = 1, size(string_arr)
        call lex_line(input_list,string_arr(ii),ii,sep,wspc,quot)
      end do

      call init_word_list(blocked_list)
      call parse_blocks(blocked_list,input_list,luout)
c dbg
c      call print_word_list(luout,blocked_list)
c dbg

      call clean_word_list(input_list)

      call parse_targets(tgt_info,blocked_list)

      call clean_word_list(blocked_list)

      end

