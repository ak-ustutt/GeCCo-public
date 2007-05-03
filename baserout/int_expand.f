
*----------------------------------------------------------------------*
      integer function int_expand(int,ibase,int_exp)
*----------------------------------------------------------------------*
*     expand a list of small integers stored in int (base ibase) 
*     into array int_exp(:)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     int, ibase
      integer, intent(out) ::
     &     int_exp(*)

      integer ::
     &     intwk, icnt

      intwk = int
      icnt = 0
      do while (intwk.gt.0)
        icnt = icnt+1
        int_exp(icnt) = mod(intwk,ibase)
        intwk = intwk/ibase
      end do

      int_expand = icnt

      return
      end
