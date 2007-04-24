*----------------------------------------------------------------------*
      subroutine test_memman()
*----------------------------------------------------------------------*
*     a few dummy allocations and deallocations to test the memman
*     routines
*     we assume that the memory manager was initialized outside
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ifc_memman.h'

      integer ::
     &     ifree
      real(8), pointer ::
     &     xbuf1(:), xbuf2(:), xbuf3(:), xbuf4(:), xbuf5(:), xbuf6(:)
      integer, pointer ::
     &     ibuf1(:), ibuf2(:), ibuf3(:)

      write(luout,*) ' test_memman: '
      write(luout,*) ' memory map on entry'
      call mem_map(.true.)

      ifree = mem_setmark('section_1')
      ifree = mem_alloc_real(xbuf1,10000,'xbuf1')
      ifree = mem_alloc_real(xbuf2,30000,'xbuf2')
      ifree = mem_setmark('section_2')
      ifree = mem_setmark('section_3')
      ifree = mem_alloc_real(xbuf3,30000,'xbuf3')
      ifree = mem_setmark('section_4')
      ifree = mem_alloc_real(xbuf4,30000,'xbuf4')
      
      write(luout,*) 'memory map after first round'
      call mem_map(.true.)

      call mem_pushmark()
      ifree = mem_gotomark('section_2')
      ifree = mem_alloc_int(ibuf1,2000,'ibuf1')
      call mem_popmark()
      ifree = mem_alloc_int(ibuf2,2000,'ibuf2')

      write(luout,*) 'memory map after second round'
      call mem_map(.true.)

      ifree = mem_flushmark('section_3')
      ifree = mem_alloc_int(ibuf3,2000,'ibuf3')

      write(luout,*) 'memory map after second round'
      call mem_map(.true.)      

      ifree = mem_flushmark()
      ifree = mem_flushmark()
      ifree = mem_flushmark()

      write(luout,*) ' memory map on exit'
      call mem_map(.true.)

      call quit(1,'test_memman','test exit')

      return
      end
