      logical function cmpxarr(xval,arr,nel)

      implicit none

      integer, intent(in) ::
     &     nel
      real(8), intent(in) ::
     &     arr(nel), xval

      integer ::
     &     idx

      cmpxarr = .true.
      do idx = 1, nel
        cmpxarr = cmpxarr.and.arr(idx).eq.xval
      end do
      return
      end
      logical function cmpiarr(ival,iarr,nel)

      implicit none

      integer, intent(in) ::
     &     nel,
     &     iarr(nel), ival

      integer ::
     &     idx

      cmpiarr = .true.
      do idx = 1, nel
        cmpiarr = cmpiarr.and.iarr(idx).eq.ival
      end do
      return
      end

