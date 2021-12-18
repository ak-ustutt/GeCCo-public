*----------------------------------------------------------------------------*
!>  Cut string into words
!>
!>  Cuts a string into smaller strings, separated by separators. 
!>  Also removes whitespaces and quotes.
!>  Note: if whitespacees are also separators they are only removed at the beginning.
!>  Otherwise the separator function takes precedence.
!>  @param wlist[out] wordlist: contains the separated words
!>  @param line[in] line (input string)
!>  @param lcnt[in] lcnt line number (for bookkeping in wlist)
!>  @param sep_list[in] list (string) of separators
!>  @param ws_list[in] list (string) of whitespace 
!>  @param qu_list[in] list (string) of quote markers
*----------------------------------------------------------------------------*
      subroutine lex_line(wlist,line,lcnt,sep_list,ws_list,qu_list)
*----------------------------------------------------------------------------*
      !
      ! cut line into words (stored on wlist), where the
      ! separators are on sep_list and the whitespace (to be
      ! cleaned away) is on ws_list
      !

      include 'def_word_list.h'

      type(word_list), intent(out) ::
     &     wlist      

      character(len=*), intent(in) ::
     &     line, sep_list, ws_list, qu_list
      integer, intent(in) ::
     &     lcnt

      logical ::
     &     ws_is_sep
      integer ::
     &     len_line, nsep, nws, nqu,
     &     idx, jdx, ipos, isep, istart, iend

      ! initialization
      len_line = len_trim(line)
      if (len_line.eq.0) return

      nsep = max(1,len_trim(sep_list))
      nws  = max(1,len_trim(ws_list))
      nqu  = max(1,len_trim(qu_list))

      ws_is_sep = .false.
      do idx = 1, nws
        do jdx = 1, nsep
          ws_is_sep = ws_is_sep.or.sep_list(jdx:jdx).eq.ws_list(idx:idx)
        end do
      end do
c dbg
c      print *,'ws_is_sep: ',ws_is_sep
c dbg

      ipos = 1
      lex_loop: do
c dbg
c        print *,'initial ipos = ',ipos
c        print *,' line: ',line(ipos:ipos)
c dbg
        if (ipos.gt.len_line) exit lex_loop

        ! skip leading white space
        do while (cmp_ch_list(line(ipos:ipos),ws_list,nws))
          ipos = ipos + 1
          if (ipos.gt.len_line) exit lex_loop
        end do

        istart = ipos

        ! check for quote
        if (cmp_ch_list(line(ipos:ipos),qu_list,nqu)) then
 
          ipos = ipos+1 
          istart = ipos ! ignore quote
          ! collect word until next quote is found
          do while (.not.cmp_ch_list(line(ipos:ipos),qu_list,nqu))
            ipos = ipos + 1
            if (ipos.gt.len_line) exit
          end do
          if (ipos.le.len_line) ipos = ipos+1
          iend = ipos-1
          isep = ipos
 
        else

          ! collect word until next separator is found
          do while (.not.cmp_ch_list(line(ipos:ipos),sep_list,nsep))
            ipos = ipos + 1
            if (ipos.gt.len_line) exit
          end do
          isep = ipos ! isep.gt.len_line -> EOL

        end if
c dbg
c        print *,'istart,isep,ipos: ',istart,isep,ipos
c dbg
        if (.not.ws_is_sep) then
c dbg
c          print *,'way 1'
c dbg
          ! remove white space 
          ! (only necessary if white space is not separator)
          iend = ipos-1
          do while (cmp_ch_list(line(iend:iend),ws_list,nws))
            iend = iend - 1
            if (iend.lt.istart) exit
          end do
        else if (nsep.gt.1) then          
c dbg
c          print *,'way 2'
c dbg
          ! find out whether another separator (instead of white space)
          ! exists (only if more than one separator active)
          iend = ipos-1
          ! skip further white space
          do while (cmp_ch_list(line(ipos:ipos),ws_list,nws))
            ipos = ipos + 1
            if (ipos.ge.len_line) exit
          end do
c dbg
c          print *,'after search ipos = ',ipos
c          print *,'after search isep = ',isep
c          print *,' line: ',line(ipos:ipos)
c dbg
          ! is the next character a separator?
          if (ipos.le.len_line) then
            if (cmp_ch_list(line(ipos:ipos),sep_list,nsep)) then
              isep = ipos
            else
              ipos = ipos - 1
            end if
          end if
c dbg
c          print *,'after check = ',isep
c dbg
        end if

! remove closing quote, if necessary:
        if (iend.gt.0)then
           if (cmp_ch_list(line(iend:iend),qu_list,nqu))
     &          iend = iend-1
           end if 
c dbg
c        print *,'iend = ',iend
c        print *,'call to new_word_list_entry'
c dbg
        call new_word_list_entry(wlist,.false.)
        ! put word to list
        if (iend.ge.istart)
     &       wlist%current%word = line(istart:iend)

        ! put separator to list
        if (ipos.le.len_line)
     &       wlist%current%sep = line(isep:isep)
        if (ipos.gt.len_line)
     &       wlist%current%sep = 'E'

        ! and also remember the position in the original file
        wlist%current%line = lcnt
        wlist%current%col  = istart

        if (isep.gt.len_line.or.ipos.gt.len_line) exit lex_loop

        ipos = ipos+1

      end do lex_loop

      ! if the last separator was not "E" add an empty record
      ! we need this later for a correct interpretation of comments
      if (wlist%current%sep.ne.'E') then
        call new_word_list_entry(wlist,.false.)
        wlist%current%sep = 'E'
        wlist%current%line = lcnt
        wlist%current%col  = istart
      end if

      ! point to head of list again
!      wlist%current => wlist%head
      
      return
      
      contains

      ! embedded function
      logical function cmp_ch_list(ch,chlist,nlist)

      integer ::
     &     nlist
      character(len=1) ::
     &     ch
      character(len=nlist) ::
     &     chlist
      
      integer ::
     &     ilist

      cmp_ch_list = .false.
      do ilist = 1, nlist
        cmp_ch_list = cmp_ch_list.or.ch.eq.chlist(ilist:ilist)
      end do

      end function

      end subroutine
