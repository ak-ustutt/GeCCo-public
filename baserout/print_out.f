!---------------------------------------------------------------------!
!> A subroutine to print an arbitrary string to luout and/or lulog
!!
!! @param string the string to be printed
!! @param out_unit string identifying the requested output unit.
!!+       "O": Print to luout only
!!+       "L": Print to lulog only
!!+       "B": Print to both lulog and luout

!----------------------------------------------------------------------!
      subroutine print_out(string, out_unit)
!----------------------------------------------------------------------!
      implicit none 

      include 'stdunit.h'

      character(len=*), parameter::
     &     i_am="print_out"
      
      character(len=*),intent(in)::
     &     string
      character(len=1),intent(in)::
     $     out_unit

      if (out_unit.EQ."O" .or. out_unit.EQ."B")
     &     write (luout,'(a)') trim(string)
      if (out_unit.EQ."L" .or. out_unit.EQ."B")
     &     write (lulog,'(a)') trim(string)
      if (out_unit.NE."O" .and. out_unit.NE."L" .and. out_unit.NE."B")
     &     call warn(i_am," unrecognized ouput unit: "//out_unit)
      
      end subroutine
