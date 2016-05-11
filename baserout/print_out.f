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

      character(len=9), parameter::
     &     i_am="print_out"
      
      character(len=*),intent(in)::
     &     string,out_unit

     
      write (lulog,*) trim(string)

      if ( (lulog .ne. luout) .and.
     &     (trim(out_unit) .EQ. "UOUT"))
     & write (luout,*) trim(string)


      if ((trim(out_unit) .NE. "UOUT") .and.
     &(trim(out_unit) .EQ. "ULOG")      )
     & call warn(i_am,"unrecognized ouput unit"//out_unit)
      
      end subroutine
