      logical function list_cmp(list1,list2,nel)

      implicit none
      integer, intent(in) ::
     &     nel, list1(nel), list2(nel)

      integer ::
     &     iel

      list_cmp = .true.

      do iel = 1, nel
        list_cmp = list_cmp.and.list1(iel).eq.list2(iel)
      end do

      return
      end
