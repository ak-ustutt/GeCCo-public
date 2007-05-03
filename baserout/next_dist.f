*----------------------------------------------------------------------*
      logical function next_dist(idist,nel,imnmx,inc)
*----------------------------------------------------------------------*
*     increment a distribution of nel numbers of idist(nel) with
*     first index running fastest within the bounds given by 
*     imnmx(2,nel)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     inc, nel, imnmx(2,nel)
      integer, intent(inout) ::
     &     idist(nel)

      integer ::
     &     iel

      next_dist = .false.
      if (inc.gt.0) then
        do iel = 1, nel
          if (idist(iel).ge.imnmx(2,iel)) then
            idist(iel) = imnmx(1,iel)
          else
            idist(iel) = idist(iel)+inc
            next_dist = .true.
            exit
          end if
        end do
      else if (inc.lt.0) then
        do iel = 1, nel
          if (idist(iel).le.imnmx(1,iel)) then
            idist(iel) = imnmx(2,iel)
          else
            idist(iel) = idist(iel)+inc
            next_dist = .true.
            exit
          end if
        end do
      end if
        
      return
      end
