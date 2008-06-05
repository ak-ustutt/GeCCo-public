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
      subroutine idxsort(ilist,ireo,nel,mode)
*----------------------------------------------------------------------*
*     as isort, but providing a reordering array
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
     &     ilist(nel), ireo(nel)
 
      integer ::
     &     iel

c      ! init reordering array
c      do iel = 1, nel
c        ireo(iel) = iel
c      end do

      if (mode.gt.0) then
        if (nel.lt.nswitch) then
          call idxsort_is_asc(ilist,ireo,nel)
        else
          call idxsort_qs_asc(ilist,ireo,nel)
        end if
      else if (mode.lt.0) then
        if (nel.lt.nswitch) then
          call idxsort_is_dsc(ilist,ireo,nel)
        else
          ! not checked yet
          call quit(0,'isort','check descending quicksort!')
c          call idxsort_qs_dsc(ilist,nel)
        end if
      else
        call quit(0,'isort','mode == 0')
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine idxsort8(ilist,ireo,nel,mode)
*----------------------------------------------------------------------*
*     as isort, but providing a reordering array, for integer(8) list
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
     &     ireo(nel)
      integer(8), intent(inout) ::
     &     ilist(nel)
 
      integer ::
     &     iel

c      ! init reordering array
c      do iel = 1, nel
c        ireo(iel) = iel
c      end do

      if (mode.gt.0) then
        if (nel.lt.nswitch) then
          call idxsort8_is_asc(ilist,ireo,nel)
        else
          call idxsort8_qs_asc(ilist,ireo,nel)
        end if
      else if (mode.lt.0) then
        if (nel.lt.nswitch) then
          call idxsort8_is_dsc(ilist,ireo,nel)
        else
          ! not checked yet
          call quit(0,'isort','check descending quicksort!')
c          call idxsort8_qs_dsc(ilist,nel)
        end if
      else
        call quit(0,'isort','mode == 0')
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine idxsortx(xlist,ireo,nel,mode)
*----------------------------------------------------------------------*
*     as idxsort, but reordering a list of real(8) values
*     mode > 0 -- ascending order
*     mode < 0 -- descending order
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'

      integer, parameter ::
     &     nswitch = 10

      integer, intent(in) ::
     &     mode, nel
      real(8), intent(inout) ::
     &     xlist(nel)
      integer, intent(inout) ::
     &     ireo(nel)
 
      integer ::
     &     iel

      if (mode.gt.0) then
        if (nel.lt.nswitch) then
          call idxsortx_is_asc(xlist,ireo,nel)
        else
          call idxsortx_qs_asc(xlist,ireo,nel)
        end if
      else if (mode.lt.0) then
        if (nel.lt.nswitch) then
          call idxsortx_is_dsc(xlist,ireo,nel)
        else
          ! not checked yet
          call quit(0,'isort','check descending quicksort!')
c          call idxsortx_qs_dsc(xlist,nel)
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

*----------------------------------------------------------------------*
      subroutine idxsort_is_asc(ilist,ireo,nel)
*----------------------------------------------------------------------*
*     insertion sort (ascending order)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     ilist(nel), ireo(nel)

      integer ::
     &     iel, jel, ihlp, jhlp

      do iel = 2, nel
        ihlp = ilist(iel)
        jhlp = ireo(iel)
        jel = iel-1
        do while (jel.gt.0.and.ihlp.lt.ilist(jel))
          ilist(jel+1) = ilist(jel)
          ireo(jel+1)  = ireo(jel)
          jel = jel-1
        end do
        ilist(jel+1) = ihlp
        ireo(jel+1) = jhlp
      end do

      return
      end
*----------------------------------------------------------------------*
      subroutine idxsort_is_dsc(ilist,ireo,nel)
*----------------------------------------------------------------------*
*     insertion sort (descending order)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     ilist(nel), ireo(nel)

      integer ::
     &     iel, jel, ihlp, jhlp

      do iel = 2, nel
        ihlp = ilist(iel)
        jhlp = ireo(iel)
        jel = iel-1
        do while (jel.gt.0.and.ihlp.gt.ilist(jel))
          ilist(jel+1) = ilist(jel)
          ireo(jel+1)  = ireo(jel)
          jel = jel-1
        end do
        ilist(jel+1) = ihlp
        ireo(jel+1)  = jhlp
      end do

      return
      end
*----------------------------------------------------------------------*
      subroutine idxsort_qs_asc(ilist,ireo,nel)
*----------------------------------------------------------------------*
*     quick sort with 3-element median (from numerical recepies)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     nstack = 50, msis = 10

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     ilist(nel), ireo(nel)

      integer ::
     &     ii, rr, jj, kk, ll,
     &     jstack, ilp, ihlp, jhlp
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
            jhlp = ireo(ii)
            jj = ii-1
            do while (jj.gt.0.and.ihlp.lt.ilist(jj))
              ilist(jj+1) = ilist(jj)
              ireo(jj+1)  = ireo(jj)
              jj = jj-1
            end do
            ilist(jj+1) = ihlp
            ireo(jj+1)  = jhlp
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
          ihlp = ireo(kk)
          ireo(kk) = ireo(ll+1)
          ireo(ll+1) = ihlp
          if (ilist(ll+1).gt.ilist(rr)) then
            ihlp = ilist(ll+1)
            ilist(ll+1) = ilist(rr)
            ilist(rr) = ihlp
            ihlp = ireo(ll+1)
            ireo(ll+1) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (ilist(ll).gt.ilist(rr)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(rr)
            ilist(rr) = ihlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (ilist(ll+1).gt.ilist(ll)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(ll+1)
            ilist(ll+1) = ihlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(ll+1)
            ireo(ll+1) = ihlp
          end if
          ii = ll+1
          jj = rr
          ilp = ilist(ll)
          jhlp = ireo(ll)
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
            ihlp = ireo(ii)
            ireo(ii) = ireo(jj)
            ireo(jj) = ihlp
          end do
          ! insert partioning element
          ilist(ll) = ilist(jj)
          ilist(jj) = ilp
          ireo(ll) = ireo(jj)
          ireo(jj) = jhlp
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
      subroutine idxsort_qs_dsc(ilist,ireo,nel)
*----------------------------------------------------------------------*
*     quick sort with 3-element median (from numerical recepies)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     nstack = 50, msis = 10

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     ilist(nel), ireo(nel)

      integer ::
     &     ii, rr, jj, kk, ll,
     &     jstack, ilp, ihlp, jhlp
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
            jhlp = ireo(ii)
            jj = ii-1
            do while (jj.gt.0.and.ihlp.gt.ilist(jj))
              ilist(jj+1) = ilist(jj)
              ireo(jj+1) = ireo(jj)
              jj = jj-1
            end do
            ilist(jj+1) = ihlp
            ireo(jj+1) = jhlp
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
          ihlp = ireo(kk)
          ireo(kk) = ireo(ll+1)
          ireo(ll+1) = ihlp
          if (ilist(ll+1).lt.ilist(rr)) then
            ihlp = ilist(ll+1)
            ilist(ll+1) = ilist(rr)
            ilist(rr) = ihlp
            ihlp = ireo(ll+1)
            ireo(ll+1) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (ilist(ll).lt.ilist(rr)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(rr)
            ilist(rr) = ihlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (ilist(ll+1).lt.ilist(ll)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(ll+1)
            ilist(ll+1) = ihlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(ll+1)
            ireo(ll+1) = ihlp
          end if
          ii = ll+1
          jj = rr
          ilp = ilist(ll)
          jhlp = ireo(ll)
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
            ihlp = ireo(ii)
            ireo(ii) = ireo(jj)
            ireo(jj) = ihlp
          end do
          ! insert partioning element
          ilist(ll) = ilist(jj)
          ilist(jj) = ilp
          ireo(ll) = ireo(jj)
          ireo(jj) = jhlp
          
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
      subroutine idxsort8_is_asc(ilist,ireo,nel)
*----------------------------------------------------------------------*
*     insertion sort (ascending order)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      integer(8), intent(inout) ::
     &     ilist(nel)
      integer, intent(inout) ::
     &     ireo(nel)

      integer(8) ::
     &     ihlp
      integer ::
     &     iel, jel, jhlp

      do iel = 2, nel
        ihlp = ilist(iel)
        jhlp = ireo(iel)
        jel = iel-1
        do while (jel.gt.0.and.ihlp.lt.ilist(jel))
          ilist(jel+1) = ilist(jel)
          ireo(jel+1)  = ireo(jel)
          jel = jel-1
        end do
        ilist(jel+1) = ihlp
        ireo(jel+1) = jhlp
      end do

      return
      end
*----------------------------------------------------------------------*
      subroutine idxsort8_is_dsc(ilist,ireo,nel)
*----------------------------------------------------------------------*
*     insertion sort (descending order)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      integer(8), intent(inout) ::
     &     ilist(nel)
      integer, intent(inout) ::
     &     ireo(nel)

      integer(8) ::
     &     ihlp
      integer ::
     &     iel, jel, jhlp

      do iel = 2, nel
        ihlp = ilist(iel)
        jhlp = ireo(iel)
        jel = iel-1
        do while (jel.gt.0.and.ihlp.gt.ilist(jel))
          ilist(jel+1) = ilist(jel)
          ireo(jel+1)  = ireo(jel)
          jel = jel-1
        end do
        ilist(jel+1) = ihlp
        ireo(jel+1)  = jhlp
      end do

      return
      end
*----------------------------------------------------------------------*
      subroutine idxsort8_qs_asc(ilist,ireo,nel)
*----------------------------------------------------------------------*
*     quick sort with 3-element median (from numerical recepies)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     nstack = 50, msis = 10

      integer, intent(in) ::
     &     nel
      integer(8), intent(inout) ::
     &     ilist(nel)
      integer, intent(inout) ::
     &     ireo(nel)

      integer(8) ::
     &     ihlp
      integer ::
     &     ii, rr, jj, kk, ll,
     &     jstack, ilp, jhlp
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
            jhlp = ireo(ii)
            jj = ii-1
            do while (jj.gt.0.and.ihlp.lt.ilist(jj))
              ilist(jj+1) = ilist(jj)
              ireo(jj+1)  = ireo(jj)
              jj = jj-1
            end do
            ilist(jj+1) = ihlp
            ireo(jj+1)  = jhlp
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
          ihlp = ireo(kk)
          ireo(kk) = ireo(ll+1)
          ireo(ll+1) = ihlp
          if (ilist(ll+1).gt.ilist(rr)) then
            ihlp = ilist(ll+1)
            ilist(ll+1) = ilist(rr)
            ilist(rr) = ihlp
            ihlp = ireo(ll+1)
            ireo(ll+1) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (ilist(ll).gt.ilist(rr)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(rr)
            ilist(rr) = ihlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (ilist(ll+1).gt.ilist(ll)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(ll+1)
            ilist(ll+1) = ihlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(ll+1)
            ireo(ll+1) = ihlp
          end if
          ii = ll+1
          jj = rr
          ilp = ilist(ll)
          jhlp = ireo(ll)
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
            ihlp = ireo(ii)
            ireo(ii) = ireo(jj)
            ireo(jj) = ihlp
          end do
          ! insert partioning element
          ilist(ll) = ilist(jj)
          ilist(jj) = ilp
          ireo(ll) = ireo(jj)
          ireo(jj) = jhlp
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
      subroutine idxsort8_qs_dsc(ilist,ireo,nel)
*----------------------------------------------------------------------*
*     quick sort with 3-element median (from numerical recepies)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     nstack = 50, msis = 10

      integer, intent(in) ::
     &     nel
      integer(8), intent(inout) ::
     &     ilist(nel)
      integer, intent(inout) ::
     &     ireo(nel)

      integer(8) ::
     &     ihlp
      integer ::
     &     ii, rr, jj, kk, ll,
     &     jstack, ilp, jhlp
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
            jhlp = ireo(ii)
            jj = ii-1
            do while (jj.gt.0.and.ihlp.gt.ilist(jj))
              ilist(jj+1) = ilist(jj)
              ireo(jj+1) = ireo(jj)
              jj = jj-1
            end do
            ilist(jj+1) = ihlp
            ireo(jj+1) = jhlp
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
          ihlp = ireo(kk)
          ireo(kk) = ireo(ll+1)
          ireo(ll+1) = ihlp
          if (ilist(ll+1).lt.ilist(rr)) then
            ihlp = ilist(ll+1)
            ilist(ll+1) = ilist(rr)
            ilist(rr) = ihlp
            ihlp = ireo(ll+1)
            ireo(ll+1) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (ilist(ll).lt.ilist(rr)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(rr)
            ilist(rr) = ihlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (ilist(ll+1).lt.ilist(ll)) then
            ihlp = ilist(ll)
            ilist(ll) = ilist(ll+1)
            ilist(ll+1) = ihlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(ll+1)
            ireo(ll+1) = ihlp
          end if
          ii = ll+1
          jj = rr
          ilp = ilist(ll)
          jhlp = ireo(ll)
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
            ihlp = ireo(ii)
            ireo(ii) = ireo(jj)
            ireo(jj) = ihlp
          end do
          ! insert partioning element
          ilist(ll) = ilist(jj)
          ilist(jj) = ilp
          ireo(ll) = ireo(jj)
          ireo(jj) = jhlp
          
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
      subroutine idxsortx_is_asc(xlist,ireo,nel)
*----------------------------------------------------------------------*
*     insertion sort (ascending order)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      real(8), intent(inout) ::
     &     xlist(nel)
      integer, intent(inout) ::
     &     ireo(nel)

      real(8) ::
     &     xhlp
      integer ::
     &     iel, jel, jhlp

      do iel = 2, nel
        xhlp = xlist(iel)
        jhlp = ireo(iel)
        jel = iel-1
        do while (jel.gt.0.and.xhlp.lt.xlist(jel))
          xlist(jel+1) = xlist(jel)
          ireo(jel+1)  = ireo(jel)
          jel = jel-1
        end do
        xlist(jel+1) = xhlp
        ireo(jel+1) = jhlp
      end do

      return
      end
*----------------------------------------------------------------------*
      subroutine idxsortx_is_dsc(xlist,ireo,nel)
*----------------------------------------------------------------------*
*     insertion sort (descending order)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      real(8), intent(inout) ::
     &     xlist(nel)
      integer, intent(inout) ::
     &     ireo(nel)

      real(8) ::
     &     xhlp
      integer ::
     &     iel, jel, jhlp

      do iel = 2, nel
        xhlp = xlist(iel)
        jhlp = ireo(iel)
        jel = iel-1
        do while (jel.gt.0.and.xhlp.gt.xlist(jel))
          xlist(jel+1) = xlist(jel)
          ireo(jel+1)  = ireo(jel)
          jel = jel-1
        end do
        xlist(jel+1) = xhlp
        ireo(jel+1)  = jhlp
      end do

      return
      end
*----------------------------------------------------------------------*
      subroutine idxsortx_qs_asc(xlist,ireo,nel)
*----------------------------------------------------------------------*
*     quick sort with 3-element median (from numerical recepies)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     nstack = 50, msis = 10

      integer, intent(in) ::
     &     nel
      real(8), intent(inout) ::
     &     xlist(nel)
      integer, intent(inout) ::
     &     ireo(nel)

      real(8) ::
     &     xhlp, xlp
      integer ::
     &     ii, rr, jj, kk, ll,
     &     jstack, ihlp, jhlp
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
            xhlp = xlist(ii)
            jhlp = ireo(ii)
            jj = ii-1
            do while (jj.gt.0.and.xhlp.lt.xlist(jj))
              xlist(jj+1) = xlist(jj)
              ireo(jj+1)  = ireo(jj)
              jj = jj-1
            end do
            xlist(jj+1) = xhlp
            ireo(jj+1)  = jhlp
          end do
          if (jstack.eq.1) exit
          rr=istack(jstack)
          jstack = jstack-1
          ll =istack(jstack)
          jstack = jstack-1
        else
          ! current sublist is xlist(rr:ll)
          ! take middle element kk
          kk = (ll+rr) / 2
          ! swap it to position ll+1
          ! and sort such that xlist(ll+1) .le. xlist(ll) .le. xlist(rr)
          xhlp = xlist(kk)
          xlist(kk) = xlist(ll+1)
          xlist(ll+1) = xhlp
          ihlp = ireo(kk)
          ireo(kk) = ireo(ll+1)
          ireo(ll+1) = ihlp
          if (xlist(ll+1).gt.xlist(rr)) then
            xhlp = xlist(ll+1)
            xlist(ll+1) = xlist(rr)
            xlist(rr) = xhlp
            ihlp = ireo(ll+1)
            ireo(ll+1) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (xlist(ll).gt.xlist(rr)) then
            xhlp = xlist(ll)
            xlist(ll) = xlist(rr)
            xlist(rr) = xhlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (xlist(ll+1).gt.xlist(ll)) then
            xhlp = xlist(ll)
            xlist(ll) = xlist(ll+1)
            xlist(ll+1) = xhlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(ll+1)
            ireo(ll+1) = ihlp
          end if
          ii = ll+1
          jj = rr
          xlp = xlist(ll)
          jhlp = ireo(ll)
          ! inner loop of Q-sort:
          do
            ! scan array section from below for element .ge. 
            ! partioning element
            do
              ii = ii+1
              if (xlist(ii).ge.xlp) exit
            end do
            ! scan array section from above for element .le. 
            ! partioning element
            do 
              jj = jj-1
              if  (xlist(jj).le.xlp) exit
            end do
            ! if pointers crossed, we are done
            if (jj.lt.ii) exit
            ! else we exchange
            xhlp = xlist(ii)
            xlist(ii) = xlist(jj)
            xlist(jj) = xhlp
            ihlp = ireo(ii)
            ireo(ii) = ireo(jj)
            ireo(jj) = ihlp
          end do
          ! insert partioning element
          xlist(ll) = xlist(jj)
          xlist(jj) = xlp
          ireo(ll) = ireo(jj)
          ireo(jj) = jhlp
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
      subroutine idxsortx_qs_dsc(xlist,ireo,nel)
*----------------------------------------------------------------------*
*     quick sort with 3-element median (from numerical recepies)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     nstack = 50, msis = 10

      integer, intent(in) ::
     &     nel
      real(8), intent(inout) ::
     &     xlist(nel)
      integer, intent(inout) ::
     &     ireo(nel)

      real(8) ::
     &     xhlp, xlp
      integer ::
     &     ii, rr, jj, kk, ll,
     &     jstack, ihlp, jhlp
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
            xhlp = xlist(ii)
            jhlp = ireo(ii)
            jj = ii-1
            do while (jj.gt.0.and.xhlp.gt.xlist(jj))
              xlist(jj+1) = xlist(jj)
              ireo(jj+1) = ireo(jj)
              jj = jj-1
            end do
            xlist(jj+1) = xhlp
            ireo(jj+1) = jhlp
          end do
          if (jstack.eq.1) exit
          rr=istack(jstack)
          jstack = jstack-1
          ll =istack(jstack)
          jstack = jstack-1
        else
          ! current sublist is xlist(rr:ll)
          ! take middle element kk
          kk = (ll+rr) / 2
          ! swap it to position ll+1
          ! and sort such that xlist(ll+1) .ge. xlist(ll) .ge. xlist(rr)
          xhlp = xlist(kk)
          xlist(kk) = xlist(ll+1)
          xlist(ll+1) = xhlp
          ihlp = ireo(kk)
          ireo(kk) = ireo(ll+1)
          ireo(ll+1) = ihlp
          if (xlist(ll+1).lt.xlist(rr)) then
            xhlp = xlist(ll+1)
            xlist(ll+1) = xlist(rr)
            xlist(rr) = xhlp
            ihlp = ireo(ll+1)
            ireo(ll+1) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (xlist(ll).lt.xlist(rr)) then
            xhlp = xlist(ll)
            xlist(ll) = xlist(rr)
            xlist(rr) = xhlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(rr)
            ireo(rr) = ihlp
          end if
          if (xlist(ll+1).lt.xlist(ll)) then
            xhlp = xlist(ll)
            xlist(ll) = xlist(ll+1)
            xlist(ll+1) = xhlp
            ihlp = ireo(ll)
            ireo(ll) = ireo(ll+1)
            ireo(ll+1) = ihlp
          end if
          ii = ll+1
          jj = rr
          xlp = xlist(ll)
          jhlp = ireo(ll)
          ! inner loop of Q-sort:
          do
            ! scan array section from below for element .ge. 
            ! partioning element
            do
              ii = ii+1
              if (xlist(ii).le.xlp) exit
            end do
            ! scan array section from above for element .le. 
            ! partioning element
            do 
              jj = jj-1
              if  (xlist(jj).ge.xlp) exit
            end do
            ! if pointers crossed, we are done
            if (jj.lt.ii) exit
            ! else we exchange
            xhlp = xlist(ii)
            xlist(ii) = xlist(jj)
            xlist(jj) = xhlp
            ihlp = ireo(ii)
            ireo(ii) = ireo(jj)
            ireo(jj) = ihlp
          end do
          ! insert partioning element
          xlist(ll) = xlist(jj)
          xlist(jj) = xlp
          ireo(ll) = ireo(jj)
          ireo(jj) = jhlp
          
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

