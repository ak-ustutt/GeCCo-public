*----------------------------------------------------------------------*
      logical function redundant_comb(comb,newcomb,redun,items,imax)
*----------------------------------------------------------------------*
*     checks if combination is redundant by comparing with redun
*     returns non-redundant combination on newcomb
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     imax, items, comb(items), redun(imax)

      integer, intent(out) ::
     &     newcomb(items)

      integer ::
     &     ii, counter(imax), temp(imax)

      counter = 0
      do ii = 1,items
        counter(comb(ii)) = counter(comb(ii)) + 1
      end do

      temp = redun
      redundant_comb = .false.

      ! check if redundant combination
      do ii = 1,items
        redundant_comb = (redundant_comb .or.
     &                    counter(redun(comb(ii))).eq.0)
      end do

      ! redefine redun pattern to allow for determination of newcomb
      do ii = 1,imax
        temp(ii) = temp(temp(ii))
      end do

      ! determine non-redundant combination
      do ii = 1,items
        newcomb(ii) = temp(comb(ii))
      end do

      return
      end
