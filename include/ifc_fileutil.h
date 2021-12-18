      interface

      logical function file_exists(ffhand)
      implicit none
      include 'def_filinf.h'
      type(filinf), intent(in) ::
     &     ffhand
      end function

      end interface
