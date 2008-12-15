*----------------------------------------------------------------------*
      logical function next_comb(comb,items,imax,nreg)
*----------------------------------------------------------------------*
*     produces next possible combination with ascending order
*     nreg=1 : only one register (simple combination)
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     imax(nreg), items, nreg

      integer, intent(inout) ::
     &     comb(items)

      integer ::
     &     digit, ii, counter(items), jj, temp(items+1), old, ireg, tot

      logical ::
     &     add, correct, skip

      next_comb = .false.
      correct = .true.

      ! return false for zero dimension
      if (items.eq.0.or.maxval(imax).eq.0) then
        return
      end if

      ! if combination is zero vector, find initial combination
      if (all(comb.eq.0)) then
        correct = .false. ! correct but don't need further treatment
        do ireg = 1,nreg
          if (imax(ireg).ge.items) then
            do ii = 1,items
              comb(ii) = maxval(imax)*(ireg-1)+ii
            end do
            next_comb = .true.
            exit
          end if
        end do
      end if

      ! each item should be larger than the previous
      do ii = 1,items-1
        if (comb(ii).ge.comb(ii+1)) correct = .false.
      end do
      ! all items should be of the same register. determine register
      do ii = 1,items
        if (ii.eq.1) then
          ireg = (comb(ii)-1)/maxval(imax)+1
        else if (ireg.ne.(comb(ii)-1)/maxval(imax)+1) then
          correct = .false.
        end if
      end do

      if (correct) then
        ! check out next possible combination in current register
        temp(1:items) = comb(1:items)
        temp(items+1) = (ireg-1)*maxval(imax)+imax(ireg)+1
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
          ! set next combination
          temp(digit) = temp(digit) + 1
          do ii = digit+1,items
            temp(ii) = temp(ii-1) + 1
          end do
          next_comb = .true.
          comb(1:items) = temp(1:items)
        else
          ! scan size of following registers and take next possible one
          do while (.not.next_comb.and.ireg.lt.nreg)
            ireg = ireg + 1
            if (imax(ireg).ge.items) then
              do ii = 1,items
                comb(ii) = (ireg-1)*maxval(imax) + ii
              end do
              next_comb = .true.
            end if
          end do
        end if
      end if

      return
      end
