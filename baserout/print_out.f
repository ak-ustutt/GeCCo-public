!---------------------------------------------------------------------!
!> A subroutine to print an arbitrary string to luout and/or lulog
!!
!! @param string the string to be printed
!! @param out_unit string identifying the requested output unit.
!!+       "UOUT": Print to log and output
!!+       "ULOG": Print only to lulog
!----------------------------------------------------------------------!
      subroutine print_out(string, out_unit)
!----------------------------------------------------------------------!
      implicit none 

      include 'stdunit.h'

      character(len=*), parameter::
     &     i_am="print_out"
      
      character(len=*),intent(in)::
     &     string,out_unit

     
      write (lulog,'(a)') trim(string)

      if ( (lulog .ne. luout) .and.
     &     (trim(out_unit) .EQ. "UOUT"))
     & write (luout,'(a)') trim(string)


      if ((trim(out_unit) .NE. "UOUT") .and.
     &     (trim(out_unit) .ne. "ULOG")      )
     & call warn(i_am,"unrecognized ouput unit"//trim(out_unit))
      
      end subroutine
