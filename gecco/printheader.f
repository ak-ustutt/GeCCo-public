      subroutine printheader(lu)

      implicit none
      
      integer, intent(in) :: lu

      write(lu,'(16(8x,a,/))')
     &'*-----------------------------------------------------------*',
     &'|                                                           |',
     &'|                         G e C C o                         |',
     &'|          a string-based general contraction code          |',
     &'|                                                           |',
     &'|                                                           |',
     &'|   principal authors:                                      |',
     &'|       andreas koehn (university of mainz, germany)        |',
     &'|       matthias hanauer (university of mainz, germany)     |',
     &'|   contributing authors:                                   |',
     &'|       jeppe olsen (university of aarhus, denmark)         |',
     &'|       wenlan liu (university of mainz, germany)           |',
     &'|       gareth richings (university of mainz, germany)      |',
     &'|       pradipta k samanta (IACS kolkata, india)            |',
     &'|                                                           |',
     &'*-----------------------------------------------------------*' 

      return
      end
