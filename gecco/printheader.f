      subroutine printheader()

      implicit none

      include "stdunit.h"

      write(lulog,'(13(8x,a,/))')
     &'*-----------------------------------------------------------*',
     &'|                                                           |',
     &'|                         G e C C o                         |',
     &'|          a string-based general contraction code          |',
     &'|                                                           |',
     &'|                                                           |',
     &'|   principal author:                                       |',
     &'|       andreas koehn (university of mainz, germany)        |',
     &'|   contributing authors:                                   |',
     &'|       jeppe olsen (university of aarhus, denmark)         |',
     &'|       matthias hanauer (university of mainz, germany)     |',
     &'|       gareth richings (university of mainz, germany)      |',
     &'|                                                           |',
     &'*-----------------------------------------------------------*' 

      return
      end
