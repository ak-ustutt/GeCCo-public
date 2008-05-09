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
      do iexp = 1, nexp
        if (ibs.le.0) goto 100
        n = int_exp(iexp)
        int = int + n*ibs
        ibs = ibs*ibase
      end do

      int8_pack = int

      return

 100  write(luout,*) 'intlist not suited for packing:'
      write(luout,*) ' list = ',int_exp(1:nexp)
      write(luout,*) ' base = ',ibase
      call quit(1,'int8_pack','intlist not suited for packing')

      end
