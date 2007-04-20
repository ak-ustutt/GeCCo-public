*----------------------------------------------------------------------*
      integer function ipermfac(ilist,len,ordered)
*----------------------------------------------------------------------*
*     obtain number of equivalent permutations of list ilist
*     i.e. for each n-tuple of equivalent elements on ilist, we
*     get a factor of n!
*     if ordered==.true., the list is in ordered sequence and
*     equivalent entries come subsequently
*     ... intended for small lists ...
*----------------------------------------------------------------------*
      
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     ordered
      integer, intent(in) ::
     &     len, ilist(len)

      logical ::
     &     succ
      integer ::
     &     idx, ihlp, iscr(len), ilast, nn
      integer, external ::
     &     ifac

      ! copy to scratch (do not change input)
      iscr(1:len) = ilist(1:len)
      ! bubble-sort (len is small, remember?) into ordered list
      ! OK, I am ashamed -- I will use insertion sort next time ....
      if (.not.ordered) then
        succ = .false.
        do while (.not.succ)
          succ = .true.
          do idx = 1, len-1
            if (iscr(idx).gt.iscr(idx+1)) then
              succ = .false.
              ihlp = iscr(idx)
              iscr(idx) = iscr(idx+1)
              iscr(idx+1) = ihlp
            end if
          end do
        end do
      end if

      ilast = iscr(1)
      nn = 1
      ipermfac = 1
      do idx = 2, len
        if (iscr(idx).eq.ilast) then
          nn = nn+1
        else
          ilast = iscr(idx)
          if (nn.gt.1) ipermfac = ipermfac*ifac(nn)
          nn = 1
        end if
      end do
      if (nn.gt.1) ipermfac = ipermfac*ifac(nn)

      if (ntest.ge.100) then
        write(luout,*) '-------------------'
        write(luout,*) ' ipermfac speaking'
        write(luout,*) '-------------------'
        write(luout,*) ' list:      ',ilist(1:len)
        if (.not.ordered) 
     &       write(luout,*) ' list(reo): ',iscr(1:len)
        write(luout,*) ' ipermfac = ',ipermfac
      end if 

      return
      end
