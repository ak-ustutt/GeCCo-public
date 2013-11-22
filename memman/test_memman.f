*----------------------------------------------------------------------*
      subroutine test_memman()
*----------------------------------------------------------------------*
*     a few dummy allocations and deallocations to test the memman
*     routines
*     we assume that the memory manager was initialized outside
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'

      integer ::
     &     ifree, idx
      real(8), pointer ::
     &     xbuf1(:), xbuf2(:), xbuf3(:), xbuf4(:), xbuf5(:), xbuf6(:)
      integer, pointer ::
     &     ibuf1(:), ibuf2(:), ibuf3(:)

      type(filinf) ::
     &     ffbuf

      write(lulog,*) ' test_memman: '
      write(lulog,*) ' memory map on entry'
      call mem_map(.true.)

      ifree = mem_setmark('section_1')
      ifree = mem_alloc_real(xbuf1,10000,'xbuf1')
      ifree = mem_alloc_real(xbuf2,30000,'xbuf2')
      ifree = mem_setmark('section_2')
      ifree = mem_setmark('section_3')
      ifree = mem_alloc_real(xbuf3,30000,'xbuf3')
      ifree = mem_setmark('section_4')
      ifree = mem_alloc_real(xbuf4,30000,'xbuf4')
      
      write(lulog,*) 'memory map after first round'
      call mem_map(.true.)

      call mem_pushmark()
      ifree = mem_gotomark('section_2')
      ifree = mem_alloc_int(ibuf1,2000,'ibuf1')
      call mem_popmark()
      ifree = mem_alloc_int(ibuf2,2000,'ibuf2')

      write(lulog,*) 'memory map after second round'
      call mem_map(.true.)

      ifree = mem_flushmark('section_3')
      ifree = mem_alloc_int(ibuf3,2000,'ibuf3')

      write(lulog,*) 'memory map after second round'
      call mem_map(.true.)      

      ifree = mem_flushmark()
      ifree = mem_flushmark()
      ifree = mem_flushmark()

      write(lulog,*) ' memory map on exit'
      call mem_map(.true.)

      write(lulog,*) ' now testing virtual arrays'

      ifree = mem_alloc_int(ibuf1,1000,'ibuf1')
      ifree = mem_alloc_int(ibuf2,1000,'ibuf2')
      do idx = 1, 1000
        ibuf1(idx) = idx
      end do

      call file_init(ffbuf,'testbuffer',ftyp_da_unf,10)
      call file_open(ffbuf)
      call mem_init_vbuffer(ffbuf,'vbuffer1',50,10)

      write(lulog,*) 'put 1..35'
      call mem_iput(ffbuf,ibuf1,1,35)

      call mem_print_vbuffer(lulog,ffbuf)

      write(lulog,*) 'put 36..70'
      call mem_iput(ffbuf,ibuf1(36),36,70)

      call mem_print_vbuffer(lulog,ffbuf)

      write(lulog,*) 'and reading again ...'
      call mem_iget(ffbuf,ibuf2,1,70)
      write(lulog,*) 'read in:'
      write(lulog,'(10i6)') ibuf2(1:70)

      write(lulog,*) 'put 71..200'
      call mem_iput(ffbuf,ibuf1(71),71,200)

      write(lulog,*) ' memory map after that'
      call mem_map(.true.)

      call mem_print_vbuffer(lulog,ffbuf)

      write(lulog,*) 'and reading 101 ... 200'
      call mem_iget(ffbuf,ibuf2,101,200)
      write(lulog,*) 'read in:'
      write(lulog,'(10i6)') ibuf2(1:100)
      write(lulog,*) 'and reading 1 ... 150'
      call mem_iget(ffbuf,ibuf2,1,150)
      write(lulog,*) 'read in:'
      write(lulog,'(10i6)') ibuf2(1:150)
      write(lulog,*) 'and reading 76 ... 123'
      call mem_iget(ffbuf,ibuf2,76,123)
      write(lulog,*) 'read in:'
      write(lulog,'(10i6)') ibuf2(1:48)
      write(lulog,*) 'and reading 180 ... 200'
      call mem_iget(ffbuf,ibuf2,180,200)
      write(lulog,*) 'read in:'
      write(lulog,'(10i6)') ibuf2(1:21)

      call mem_clean_vbuffer('vbuffer1')
      call file_close_delete(ffbuf)

      call quit(1,'test_memman','test exit')

      return
      end
