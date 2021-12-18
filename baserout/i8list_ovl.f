      subroutine i8list_ovl(ovl,list1,list2,nel,inc)

      implicit none
      integer(8), intent(out) ::
     &     ovl(nel)
      integer, intent(in) ::
     &     nel, inc
      integer(8), intent(in) ::
     &     list1(nel), list2(nel)

      integer ::
     &     idx, iel

      idx = 1
      do iel = 1, nel
        ovl(iel) = min(list1(idx),list2(idx))
        idx = idx+inc
      end do

      return
      end
