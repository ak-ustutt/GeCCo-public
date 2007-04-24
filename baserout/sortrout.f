*----------------------------------------------------------------------*
*     some sort routines
*----------------------------------------------------------------------*
      subroutine isort(ilist,nel,mode)
*----------------------------------------------------------------------*
*     wrapper
*     mode > 0 -- ascending order
*     mode < 0 -- descending order
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'

      integer, parameter ::
     &     nswitch = 10

      integer, intent(in) ::
     &     mode, nel
      integer, intent(inout) ::
     &     ilist(nel)
 

      if (mode.gt.0) then
        if (nel.lt.nswitch) then
          call isort_is_asc(ilist,nel)
        else
          call isort_qs_asc(ilist,nel)
        end if
      else if (mode.lt.0) then
        if (nel.lt.nswitch) then
          call isort_is_dsc(ilist,nel)
        else
          ! not checked yet
          call quit(0,'isort','check descending quicksort!')
c          call isort_qs_dsc(ilist,nel)
        end if
      else
        call quit(0,'isort','mode == 0')
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine isort_is_asc(ilist,nel)
*----------------------------------------------------------------------*
*     insertion sort (ascending order)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     ilist(nel)

      integer ::
     &     iel, jel, ihlp

      do iel = 2, nel
        ihlp = ilist(iel)
        jel = iel-1
        do while (jel.gt.0.and.ihlp.lt.ilist(jel))
          ilist(jel+1) = ilist(jel)
          jel = jel-1
        end do
        ilist(jel+1) = ihlp
      end do

      return
      end
*----------------------------------------------------------------------*
      subroutine isort_is_dsc(ilist,nel)
*----------------------------------------------------------------------*
*     insertion sort (descending order)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     ilist(nel)

      integer ::
     &     iel, jel, ihlp

      do iel = 2, nel
        ihlp = ilist(iel)
        jel = iel-1
        do while (jel.gt.0.and.ihlp.gt.ilist(jel))
          ilist(jel+1) = ilist(jel)
          jel = jel-1
        end do
        ilist(jel+1) = ihlp
      end do

      return
      end
*----------------------------------------------------------------------*
      subroutine isort_qs_asc(ilist,nel)
*----------------------------------------------------------------------*
*     quick sort with 3-element median (from numerical recepies)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     nstack = 50, msis = 10

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     ilist(nel)

      integer ::
     &     ii, rr, jj, kk, ll,
     &     jstack, ilp, ihlp
      integer, allocatable ::
     &     istack(:)
      

      allocate(istack(nstack))

      rr = nel
      ll = 1
      jstack = 1

      do
        if (rr-ll.le.msis) then
          ! straight insertion sort for small list:
          do ii = 2, nel
            ihlp = ilist(ii)
            jj = ii-1
            do while (jj.gt.0.and.ihlp.lt.ilist(jj))
              ilist(jj+1) = ilist(jj)
              jj = jj-1
            end do
            ilist(jj+1) = ihlp
          end do
          if (jstack.eq.1) exit
          rr=istack(jstack)
          jstack = jstack-1
          ll =istack(jstack)
          jstack = jstack-1
        else
          ! current sublist is ilist(rr:ll)
          ! take middle element kk
          kk = (ll+rr) / 2
          ! swap it to position ll+1
          ! and sort such that ilist(ll+1) .le. ilist(ll) .le. ilist(rr)
          ihlp = ilist(kk)
          ilist(kk) = ilist(ll+1)
          ilist(ll+1) = ihlp
          if (ilist(ll+1).gt.ilist(rr)) then
            ihlp = ilist(ll+1)
            ilist(ll+1) = ilist(rr)
            ilist(rr) = ihlp
          end if
          if (ilist(ll).gt.ilist(rr)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(rr)
            ilist(rr) = ihlp
          end if
          if (ilist(ll+1).gt.ilist(ll)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(ll+1)
            ilist(ll+1) = ihlp
          end if
          ii = ll+1
          jj = rr
          ilp = ilist(ll)
          ! inner loop of Q-sort:
          do
            ! scan array section from below for element .ge. 
            ! partioning element
            do
              ii = ii+1
              if (ilist(ii).ge.ilp) exit
            end do
            ! scan array section from above for element .le. 
            ! partioning element
            do 
              jj = jj-1
              if  (ilist(jj).le.ilp) exit
            end do
            ! if pointers crossed, we are done
            if (jj.lt.ii) exit
            ! else we exchange
            ihlp = ilist(ii)
            ilist(ii) = ilist(jj)
            ilist(jj) = ihlp
          end do
          ! insert partioning element
          ilist(ll) = ilist(jj)
          ilist(jj) = ilp
          
          jstack = jstack+2
          if (jstack > nstack) stop 'too small stack size in isort_qs'
          if (rr-ii+1 .ge. jj-ll) then
            istack(jstack) = rr
            istack(jstack-1) = ii
            rr = jj-1
          else
            istack(jstack) = jj-1
            istack(jstack-1) = ll
            ll = ii
          end if
            
        end if

      end do

      deallocate(istack)
      
      return
      end

*----------------------------------------------------------------------*
      subroutine isort_qs_dsc(ilist,nel)
*----------------------------------------------------------------------*
*     quick sort with 3-element median (from numerical recepies)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     nstack = 50, msis = 10

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     ilist(nel)

      integer ::
     &     ii, rr, jj, kk, ll,
     &     jstack, ilp, ihlp
      integer, allocatable ::
     &     istack(:)
      

      allocate(istack(nstack))

      rr = nel
      ll = 1
      jstack = 1

      do
        if (rr-ll.le.msis) then
          ! straight insertion sort for small list:
          do ii = 2, nel
            ihlp = ilist(ii)
            jj = ii-1
            do while (jj.gt.0.and.ihlp.gt.ilist(jj))
              ilist(jj+1) = ilist(jj)
              jj = jj-1
            end do
            ilist(jj+1) = ihlp
          end do
          if (jstack.eq.1) exit
          rr=istack(jstack)
          jstack = jstack-1
          ll =istack(jstack)
          jstack = jstack-1
        else
          ! current sublist is ilist(rr:ll)
          ! take middle element kk
          kk = (ll+rr) / 2
          ! swap it to position ll+1
          ! and sort such that ilist(ll+1) .ge. ilist(ll) .ge. ilist(rr)
          ihlp = ilist(kk)
          ilist(kk) = ilist(ll+1)
          ilist(ll+1) = ihlp
          if (ilist(ll+1).lt.ilist(rr)) then
            ihlp = ilist(ll+1)
            ilist(ll+1) = ilist(rr)
            ilist(rr) = ihlp
          end if
          if (ilist(ll).lt.ilist(rr)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(rr)
            ilist(rr) = ihlp
          end if
          if (ilist(ll+1).lt.ilist(ll)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(ll+1)
            ilist(ll+1) = ihlp
          end if
          ii = ll+1
          jj = rr
          ilp = ilist(ll)
          ! inner loop of Q-sort:
          do
            ! scan array section from below for element .ge. 
            ! partioning element
            do
              ii = ii+1
              if (ilist(ii).le.ilp) exit
            end do
            ! scan array section from above for element .le. 
            ! partioning element
            do 
              jj = jj-1
              if  (ilist(jj).ge.ilp) exit
            end do
            ! if pointers crossed, we are done
            if (jj.lt.ii) exit
            ! else we exchange
            ihlp = ilist(ii)
            ilist(ii) = ilist(jj)
            ilist(jj) = ihlp
          end do
          ! insert partioning element
          ilist(ll) = ilist(jj)
          ilist(jj) = ilp
          
          jstack = jstack+2
          if (jstack > nstack) stop 'too small stack size in isort_qs'
          if (rr-ii+1 .ge. jj-ll) then
            istack(jstack) = rr
            istack(jstack-1) = ii
            rr = jj-1
          else
            istack(jstack) = jj-1
            istack(jstack-1) = ll
            ll = ii
          end if
            
        end if

      end do

      deallocate(istack)
      
      return
      end
