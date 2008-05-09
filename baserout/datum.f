*----------------------------------------------------------------------*
      subroutine datum(str)
*----------------------------------------------------------------------*
*     properly formatted time_of_day output
*     adapted from Marco Kattanek's f90 routine for TURBOMOLE
*----------------------------------------------------------------------*
      implicit none

      character*(*), intent(out) ::
     &     str

      character(5) ::
     &     zone
      integer ::
     &     values(8)

      call date_and_time(str(1:8),str(9:18),zone,values)
      str = str(1:4)//'-'//str(5:6)//'-'//str(7:8)//' '//
     &     str(9:10)//':'//str(11:12)//':'//str(13:18)

      return
      end
