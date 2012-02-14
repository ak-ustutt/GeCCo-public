      integer function i8common_entry(list1,list2,n1,n2)

c     returns position of first entry of list1 which can be found on list2

      implicit none
      integer, intent(in) ::
     &     n1, n2
      integer(8), intent(in) ::
     &     list1(n1), list2(n2)

      integer ::
     &     idx, jdx

      i8common_entry = 0

      main_loop: do idx = 1, n1
        do jdx = 1, n2
          if (list1(idx).eq.list2(jdx)) then
            i8common_entry = idx
            exit main_loop
          end if
        end do
      end do main_loop

      return
      end
