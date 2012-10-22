      subroutine lex_file(wlist,inp_file,keyword_start,keyword_end)

      implicit none

      include 'stdunit.h'
      include 'def_word_list.h'
      include 'def_filinf.h'
      include 'par_special_chars.h'


      type(word_list), intent(inout) ::
     &     wlist
      type(filinf), intent(inout) ::
     &     inp_file
      character(len=*), intent(in) ::
     &     keyword_start, keyword_end

      integer ::
     &     un, lcnt
      character(len=256) ::
     &     line

      call file_open(inp_file)

      un = inp_file%unit
      rewind un

      ! look for initial keyword to start (if requested):
      lcnt = 0
      if (trim(keyword_start).ne.'') then
        do
          read(un,'(a)',end=101) line
          lcnt = lcnt+1
          if (trim(line).eq.trim(keyword_start)) exit
        end do
      end if
          
      do 
         read(un,'(a)',end=99) line
         lcnt = lcnt+1
         if (trim(keyword_end).ne.''.and.
     &       trim(line).eq.trim(keyword_end)) goto 100
         call lex_line(wlist,line,lcnt,sep,wspc,quot)
      end do

 99   if (trim(keyword_end).ne.'') then
        write(luout,*) 'Error reading file '//trim(inp_file%name)
        write(luout,*) 'File ended before keyword '//trim(keyword_end)
     &                 //' was reached!' 
        call quit(0,'lex_file','File ended while scanning ...')
      end if

 100  call file_close_keep(inp_file)

      return

 101  write(luout,*) 'Error reading file '//trim(inp_file%name)
      write(luout,*) 'File ended before keyword '//trim(keyword_start)  
     &                 //' was found!'
        call quit(0,'lex_file','File ended while scanning ...')

      end
