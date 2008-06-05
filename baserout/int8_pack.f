*----------------------------------------------------------------------*
      integer(8) function int8_pack(int_exp,nexp,ibase)
*----------------------------------------------------------------------*
*     pack a list of small integers into one integer using base (ibase)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     nexp, int_exp(nexp), ibase
      integer(8) ::
     &     ib8, int, ibs, n, iexp

      int = 0
      ibs = 1
      ib8 = ibase
      do iexp = 1, nexp, 4
        n = int_exp(iexp)
        int = int + n*ibs
        ibs = ibs*ib8
        n = int_exp(iexp+1)
        int = int + n*ibs
        ibs = ibs*ib8
        n = int_exp(iexp+2)
        int = int + n*ibs
        ibs = ibs*ib8
        n = int_exp(iexp+3)
        int = int + n*ibs
        ibs = ibs*ib8
      end do
      do iexp = nexp - mod(nexp,4) + 1, nexp
        n = int_exp(iexp)
        int = int + n*ibs
        ibs = ibs*ib8
      end do

      int8_pack = int

      return

c 100  write(luout,*) 'intlist not suited for packing:'
c      write(luout,*) ' list = ',int_exp(1:nexp)
c      write(luout,*) ' base = ',ibase
c      call quit(1,'int8_pack','intlist not suited for packing')

      end
