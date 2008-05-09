
*----------------------------------------------------------------------*
      integer function int8_expand(int,ibase,int_exp)
*----------------------------------------------------------------------*
*     expand a list of small integers stored in int (base ibase) 
*     into array int_exp(:)
*----------------------------------------------------------------------*
      implicit none

      integer(8), intent(in) ::
     &     int, ibase
      integer, intent(out) ::
     &     int_exp(*)

      integer(8) ::
     &     intwk
      integer ::
     &     icnt

      intwk = int
      icnt = 0
      do while (intwk.gt.0)
        icnt = icnt+1
        int_exp(icnt) = mod(intwk,ibase)
        intwk = intwk/ibase
      end do

      int8_expand = icnt

      return
      end
