      logical function file_exists(ffhand)

      implicit none

      include 'def_filinf.h'

      type(filinf), intent(in) ::
     &     ffhand

      inquire(file=trim(ffhand%name),exist=file_exists)

      return
      end
