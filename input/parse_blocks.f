      subroutine parse_blocks(wlist_p,wlist,luerr)
      ! sort wlist into levels according to blocks indicated
      ! by the separators ( and )

      implicit none

      include 'def_word_list.h'

      character, parameter ::
     &      char_block_open  = '(',
     &      char_block_close = ')',
     &      char_comnt_1 = '#',
     &      char_comnt_2 = '!'

      type(word_list), intent(inout) ::
     &      wlist
      type(word_list), intent(inout) ::
     &      wlist_p
      integer, intent(in) ::
     &      luerr

      type(word_list_entry), pointer ::
     &      wl_pnt, wlp_pnt

      integer ::
     &      level
      logical ::
     &      branch 
      character ::
     &      sep_rem

      ! set to head
      wlist%current   => wlist%head
      wlist_p%current => wlist_p%head

      branch = .false.
      level = 0
      main_lp: do

        ! check whether this is a comment
        if (wlist%current%word(1:1).eq.char_comnt_1.or.
     &      wlist%current%word(1:1).eq.char_comnt_2) then

          ! skip until newline appears (marked by 'E')          
          comment_lp: do
            sep_rem = wlist%current%sep
            if (.not.associated(wlist%current%next)) exit main_lp
            wlist%current => wlist%current%next
            if (sep_rem.eq.'E') cycle main_lp
          end do comment_lp

        end if

! dbg 
!         print *,'make new entry',branch,level
! dbg
        call new_word_list_entry(wlist_p,branch)
        ! copy entry
! dbg 
!         print *,'and copy'
! dbg
        wlist_p%current%word = wlist%current%word
        wlist_p%current%sep  = wlist%current%sep 
        wlist_p%current%line = wlist%current%line
        wlist_p%current%col  = wlist%current%col 

! dbg
!        print *,'separator is: ',wlist%current%sep
! dbg
        ! open new block?
        if (wlist%current%sep.eq.char_block_open) then
          branch = .true.
          level = level+1
        else if (wlist%current%sep.eq.char_block_close) then
          branch = .false.
          level = level-1
          ! move parsed list pointer one level up
          wlist_p%current => wlist_p%current%up
        else
          branch = .false.
        end if

        if (level.lt.0) exit main_lp ! error!

        ! advance original linear list
        if (.not.associated(wlist%current%next)) exit main_lp
        wlist%current => wlist%current%next
 
      end do main_lp

      if (level.ne.0) then
        write(luerr,*) 'Error while parsing blocks:'
        if (level.lt.0)
     &    write(luerr,*) 'A ")" occurred without matching "("'
        if (level.gt.0)
     &    write(luerr,*) 'It seems that not all "(" are matched by ")"'
        write(luerr,'(x,a,i4,a,i3)') 'I noted this in line',
     &         wlist%current%line,' near column ',wlist%current%col
        call quit(0,'parse_blocks','Error while parsing blocks')
      end if

      end
