      subroutine list_ovl(ovl,list1,list2,nel,inc)

      implicit none
      integer, intent(out) ::
     &     ovl(nel)
      integer, intent(in) ::
     &     nel, inc, list1(nel), list2(nel)

      integer ::
     &     idx, iel

      idx = 1
      do iel = 1, nel
        ovl(iel) = min(list1(idx),list2(idx))
        idx = idx+inc
      end do

      return
      end
