      subroutine hostname(host)

      implicit none

      character*(*), intent(out) ::
     &     host
      
      integer ::
     &     len

      call gethostname(host)
      len = index(host,char(0))
      host(len:len_trim(host)) = ' '

      return
      end
