*----------------------------------------------------------------------*
      subroutine clean_line(line,delimiter,n_delim)
*----------------------------------------------------------------------*
*     remove any superfluous spaces from line line
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     n_delim
      character, intent(in) ::
     &     delimiter(n_delim)
      character, intent(inout) ::
     &     line*(*)

      integer ::
     &     len,ipos,jpos,ispc,idx

      len = len_trim(line)

      ipos = 1
      jpos = 1
      ispc = 0
      do while(ipos.le.len)
        if (line(ipos:ipos).eq.' ') then
          ispc = ispc+1
          if (ispc.eq.1.and.ipos.lt.len.and.jpos.gt.1) then
            if (ipos.gt.jpos) line(jpos:jpos) = line(ipos:ipos)
            jpos = jpos+1
          end if
          ipos = ipos+1
        else
          ! a delimiter follows space?
          if (ispc.gt.0) then
            ispc = 0
            do idx = 1, n_delim
              if (delimiter(idx).eq.' ') cycle
              if (line(ipos:ipos).eq.delimiter(idx)) then
                ! remove space then
                ispc = 2
c                jpos = max(1,jpos-1)
                exit
              end if
            end do
          end if
          if (ipos.gt.jpos) line(jpos:jpos) = line(ipos:ipos)
          ipos = ipos+1
          jpos = jpos+1
        end if

      end do

      line(jpos:len) = ' '
      
      return
      end
