*----------------------------------------------------------------------*
!>  takes a sorts the lists of values in increasing value order
!!
!!  xlist is actually nopt lists of elements of the me_lists. Each sorted in ascending order.
!!  idxlist contains the corresponding indices of the elements on the lists.
!!  
!!  this subroutine merges the lists into xlist_all and idxlist_all where xlist_all is the values in ascending order
!!  idxlist_all are the indices of the elements where every element is indexed by a pair ( index of list (iopt),index on list )
!!   if choice is unequal zero only values of the respective operator are considered.
*----------------------------------------------------------------------*
      subroutine merge_min_lists(xlist_all,idxlist_all,ntrials_all,
     &     xlist,idxlist,nopt,maxtrials,ntrials,choice)

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nopt, maxtrials, ntrials(nopt), idxlist(maxtrials,nopt),
     &     choice
      real(8), intent(in) ::
     &     xlist(maxtrials,nopt)

      integer, intent(out) ::
     &     idxlist_all(2,maxtrials), ntrials_all
      real(8), intent(out) ::
     &     xlist_all(maxtrials)
      
      integer ::
     &     itry, idx(nopt), iopt, jopt
      real(8) ::
     &     xlow

      if (choice.eq.0) then

        idx(1:nopt) = 0 
  
        ntrials_all = min(maxtrials,sum(ntrials(1:nopt)))
  
        do itry = 1, ntrials_all
          
          ! look for next-lowest element in all lists 
          xlow = huge(xlow)
          jopt = -1
          do iopt = 1, nopt
            if (idx(iopt).lt.ntrials(iopt)) then
              if (xlist(idx(iopt)+1,iopt).lt.xlow) then
                xlow = xlist(idx(iopt)+1,iopt)
                jopt = iopt
              end if
            end if
          end do
  
          if (jopt.le.0)
     &         call quit(1,'merge_min_lists','something''s wrong !')
  
          ! increment this index
          idx(jopt) = idx(jopt)+1
          
          ! store information on new list
          xlist_all(itry) = xlist(idx(jopt),jopt)
          idxlist_all(1:2,itry) = (/jopt,idxlist(idx(jopt),jopt)/)
  
        end do
       elseif (choice.le.nopt) then
         ntrials_all = min(maxtrials,ntrials(choice))

         do itry = 1, ntrials_all

           xlist_all(itry) = xlist(itry,choice)
           idxlist_all(1:2,itry) = (/choice,idxlist(itry,choice)/)

         enddo
       else
         call quit(1,'merge_min_lists','wrong choice')
       endif

      if (ntest.ge.100) then
        write(lulog,*) 'new array:'
        do itry = 1, ntrials_all
          write(lulog,'(1x,f14.8,2x,i3,i10)')
     &         xlist_all(itry), idxlist_all(1:2,itry)
        end do
        write(lulog,*) 'end of array'

      end if

      end
