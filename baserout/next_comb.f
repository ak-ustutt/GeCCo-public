*----------------------------------------------------------------------*
      logical function next_comb(comb,items,imax)
*----------------------------------------------------------------------*
*     produces next possible combination with ascending order
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     imax, items

      integer, intent(inout) ::
     &     comb(items)

      integer ::
     &     digit, ii, counter(items), jj, temp(items+1), old

      logical ::
     &     add, correct, skip

      next_comb = .false.

      correct = .true.
      do ii = 1,items-1
        if (comb(ii).ge.comb(ii+1)) correct = .false.
      end do

      if (correct) then
        temp(1:items) = comb(1:items)
        temp(items+1) = imax+1
        digit = items
        skip = .false.
        do ii = items,1,-1
          if ((temp(ii)+1.eq.temp(ii+1)).and.(.not.skip)) then
            digit = digit -1
          else
            skip = .true.
          end if
        end do
        if (skip) then
          temp(digit) = temp(digit) + 1
          do ii = digit+1,items
            temp(ii) = temp(ii-1) + 1
          end do
          next_comb = .true.
          comb(1:items) = temp(1:items)
        end if
      end if

      return
      end
