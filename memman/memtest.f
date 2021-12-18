      subroutine memtest
c      use memman
      implicit none
      include 'stdunit.h'
      include 'ifc_memman.h'

      integer ::
     &     mem, mem_c

      integer, pointer ::
     &     iarr1(:), iarr2(:), iarr3(:)
      real(8), pointer ::
     &     xarr1(:), xarr2(:), xarr3(:)
      logical ::
     &     check

      check = .true.
      mem = 20000
      call mem_init(mem)

c      mem_c = memman_alloc(mtyp_int,300,'testarray i1',ipnt=iarr1)
      mem_c = mem_alloc_int(iarr1,300,'testarray i1')
      print *,'free = ',mem_c
c      mem_c = memman_alloc(mtyp_int,155,'testarray i2',ipnt=iarr2)
      mem_c = mem_alloc_int(iarr2,155,'testarray i2')
      print *,'free = ',mem_c
c      mem_c = memman_alloc(mtyp_rl8,4000,'testarray r2',xpnt=xarr2)
      mem_c = mem_alloc_real(xarr2,4000,'testarray r2')
      print *,'free = ',mem_c

      call mem_map(check)


      print *,'removing array i2'
c      mem_c = memman_dealloc('testarray i2')
      mem_c = mem_dealloc('testarray i2')
      print *,'free = ',mem_c
      call mem_map(check)


      print *,'adding section 2'
c      mem_c = memman_addsection('section 2')
      mem_c = mem_setmark('section 2')

      call mem_map(check)

      print *,'adding array r1'
c      mem_c = memman_alloc(mtyp_rl8,10000,'array r1',xpnt=xarr1)
      mem_c = mem_alloc_real(xarr1,10000,'array r1')
      print *,'free = ',mem_c

      call mem_map(check)
      print *,'adding array i2'
c      mem_c = memman_alloc(mtyp_int,300,'array i2',ipnt=iarr1)
      mem_c = mem_alloc_int(iarr1,300,'array i2')
      print *,'free = ',mem_c

      print *,'adding section 3'
c      mem_c = memman_addsection('section 3')
      mem_c = mem_setmark('section 3')
      print *,'free = ',mem_c

      call mem_map(check)

      print *,'removing section 2'
c      mem_c = memman_remsection('section 2')
      mem_c = mem_flushmark('section 2')
      print *,'free = ',mem_c

      call mem_map(check)


      call mem_clean()

      stop 'test success ??'
      end
