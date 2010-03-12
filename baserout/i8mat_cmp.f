      integer function i8mat_cmp(list1,list2,nlin,ncol)
      ! not identical to i8list_cmp since here matrix is run
      ! through linewise (not columnwise)

      integer, intent(in) ::
     &     nlin, ncol
      integer(8), intent(in) ::
     &     list1(nlin,ncol), list2(nlin,ncol)

      i8mat_cmp = 0
      line_loop: do ilin = 1, nlin
        do icol = 1, ncol
          if (list1(ilin,icol).gt.list2(ilin,icol)) i8mat_cmp = +1
          if (list1(ilin,icol).lt.list2(ilin,icol)) i8mat_cmp = -1
          if (i8mat_cmp.ne.0) exit line_loop
        end do
      end do line_loop

      return
      end
