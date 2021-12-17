      subroutine printheader(lu)

      implicit none
      
      integer, intent(in) :: lu

      write(lu,'(16(8x,a,/))')
     &'*-------------------------------------------------------------*',
     &'|                                                             |',
     &'|                          G e C C o                          |',
     &'|           a string-based general contraction code           |',
     &'|                                                             |',
     &'|                                                             |',
     &'|   principal author:                                         |',
     &'|       andreas koehn                                         |',
     &'|   with contributions from:                                  |',
     &'|       mathias hanauer, arne bargholz, yuri a aoto,          |',
     &'|       pradipta k samanta, joshua a black, alexander waigum, |',
     &'|       patrik zielinski, wenlan liu,                         |',
     &'|       gareth richings, jeppe olsen                          |',
     &'|                                                             |',
     &'*-------------------------------------------------------------*' 

      return
      end
