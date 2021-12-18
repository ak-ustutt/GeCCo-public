      integer function first_nonblank(str)

      implicit none

      character, intent(in) ::
     &     str*(*)

      integer ::
     &     ipos, len

      ipos = 1
      len = len_trim(str)

      do while(ipos.le.len)
        if (str(ipos:ipos).ne.' ') exit
        ipos = ipos+1
      end do

      if (ipos.gt.len) ipos = 0

      first_nonblank = ipos

      return
      end
