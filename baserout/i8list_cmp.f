      integer function i8list_cmp(list1,list2,nel)

      integer, intent(in) ::
     &     nel
      integer(8), intent(in) ::
     &     list1(nel), list2(nel)

      i8list_cmp = 0
      do iel = 1, nel
        if (list1(iel).gt.list2(iel)) i8list_cmp = +1
        if (list1(iel).lt.list2(iel)) i8list_cmp = -1
        if (i8list_cmp.ne.0) exit
      end do

      return
      end
