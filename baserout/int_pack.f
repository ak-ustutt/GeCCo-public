*----------------------------------------------------------------------*
      integer function int_pack(int_exp,nexp,ibase)
*----------------------------------------------------------------------*
*     pack a list of small integers into one integer using base (ibase)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     nexp, int_exp(nexp), ibase
      integer ::
     &     int, ibs, iexp

      int = 0
      ibs = 1
      do iexp = 1, nexp
        int = int + int_exp(iexp)*ibs
        ibs = ibs*ibase
        if (ibs.le.0) goto 100
      end do

      int_pack = int

      return

 100  write(luout,*) 'intlist not suited for packing:'
      write(luout,*) ' list = ',int_exp(1:nexp)
      write(luout,*) ' base = ',ibase
      call quit(1,'packint','intlist not suited for packing')

      end
