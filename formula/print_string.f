*----------------------------------------------------------------------*
*     a set of three subroutines to print a "string" as defined in
*     def_contraction.h
*
*     print_string: print [CA][HPVX] info
*     print_string_cnt: print contraction info
*     print_string_idx: print index info
*
*     arguments:
*     str - the string
*     nstr - length of string
*     nvtx - number of (primitive) vertices that it represents
*
*     andreas, Dec. 2020
*
*----------------------------------------------------------------------*
      subroutine print_string(str,nstr,nvtx)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      
      character, parameter ::
     &     ca_c(2)   = (/'C','A'/),
     &     hpvx_c(4) = (/'H','P','V','X'/) 
      
      integer :: nstr,nvtx
      type(string_element) :: str(nstr)

      integer :: ii, ivtx_last, iwatch

      iwatch = 0
      write(lulog,'(1x,"{")',advance='no')
      ivtx_last = 1
      do ii = 1, nstr
        do while (ivtx_last<str(ii)%vtx)
          write(lulog,'(" }{")',advance='no')
          ivtx_last = ivtx_last+1
          iwatch = iwatch+1
          if (iwatch.gt.100) exit
        end do
        if (.not.str(ii)%del)
     &    write(lulog,'(1x,a1,a1)',advance='no')
     &       ca_c(str(ii)%ca),hpvx_c(str(ii)%hpvx)
      end do
      do while (ivtx_last<nvtx)
        write(lulog,'(" }{")',advance='no')
        ivtx_last = ivtx_last+1
        iwatch = iwatch+1
        if (iwatch.gt.100) exit
      end do
      write(lulog,'(" }")')
      if (iwatch.gt.100)
     &     write(lulog,'(1x,"warning: input data seems corrupted")')
      
      end subroutine print_string
      subroutine print_string_cnt(str,nstr,nvtx)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      
      character, parameter ::
     &     cx_c(3)   = (/'/', 'c','x'/)
      
      integer :: nstr, nvtx
      type(string_element) :: str(nstr)

      integer :: ii, ivtx_last, idx, iwatch

      iwatch = 0
      write(lulog,'(1x,"{")',advance='no')
      ivtx_last = 1
      do ii = 1, nstr
        do while (ivtx_last<str(ii)%vtx)
          write(lulog,'(" }{")',advance='no')
          ivtx_last = ivtx_last+1
          iwatch = iwatch+1
          if (iwatch.gt.100) exit
        end do
        idx = 1
        if (str(ii)%cnt.gt.0) idx = 2
        if (str(ii)%ext) idx = 3
   
        if (.not.str(ii)%del)
     &    write(lulog,'(1x,a1,i1)',advance='no')
     &       cx_c(idx),str(ii)%cnt
      end do
      do while (ivtx_last<nvtx)
        write(lulog,'(" }{")',advance='no')
        ivtx_last = ivtx_last+1
        iwatch = iwatch+1
        if (iwatch.gt.100) exit
      end do
      write(lulog,'(" }")')
      if (iwatch.gt.100)
     &     write(lulog,'(1x,"warning: input data seems corrupted")')
      
      end subroutine print_string_cnt
      subroutine print_string_idx(str,nstr,nvtx)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      character, parameter ::
     &     hpvx_ic(4)   = (/'h', 'p', 'v' ,'x'/)
      
      integer :: nstr, nvtx
      type(string_element) :: str(nstr)

      integer :: ii, ivtx_last, iwatch

      iwatch = 0
      write(lulog,'(1x,"{")',advance='no')
      ivtx_last = 1
      do ii = 1, nstr
        do while (ivtx_last<str(ii)%vtx)
          write(lulog,'(" }{")',advance='no')
          ivtx_last = ivtx_last+1
          iwatch = iwatch+1
          if (iwatch.gt.100) exit
        end do

        if (.not.str(ii)%del)
     &       write(lulog,'(1x,a1,i1)',advance='no')
     &       hpvx_ic(str(ii)%hpvx),str(ii)%idx
      end do
      do while (ivtx_last<nvtx)
        write(lulog,'(" }{")',advance='no')
        ivtx_last = ivtx_last+1
        iwatch = iwatch+1
        if (iwatch.gt.100) exit
      end do
      write(lulog,'(" }")')
      if (iwatch.gt.100)
     &     write(lulog,'(1x,"warning: input data seems corrupted")')
      
      end subroutine print_string_idx
