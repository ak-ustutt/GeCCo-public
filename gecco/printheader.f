      subroutine printheader()

      implicit none

      include "stdunit.h"

      write(luout,'(6(8x,a,/))')
     &'*-----------------------------------------------------------*',
     &'|                                                           |',
     &'|                         G e C C o                         |',
     &'|          a string-based general contraction code          |',
     &'|                                                           |',
     &'*-----------------------------------------------------------*' 

      return
      end
