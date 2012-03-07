*----------------------------------------------------------------------*
      integer function zirat()
*----------------------------------------------------------------------*
*
* Ratio between real and integer
*
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'

      integer(8) ::
     &     itest

      itest = Z'00FFFFFFFFFFFFFF'
      call iset_test(itest)
      if (itest.eq.0) then
        zirat = 1
      else if (itest.eq.Z'00000000FFFFFFFF'.or.
     &         itest.eq.Z'00FFFFFF00000000') then
        ! outcome depends on high-word/low-word ordering, but this
        ! is not important for us
        zirat = 2
      else
        write(luout,*) 'Silly outcome in ZIRAT, d''you run on a C64?',
     &       itest
        call quit(1,'zirat','silly result')
      end if

      return
      end
*----------------------------------------------------------------------*

      subroutine iset_test(itest)

      integer ::
     &     itest

      itest = 0

      return
      end
