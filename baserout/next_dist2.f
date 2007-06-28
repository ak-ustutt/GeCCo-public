*----------------------------------------------------------------------*
      logical function next_dist2(idist,nel,imn,imx,inc)
*----------------------------------------------------------------------*
*     increment a distribution of nel numbers of idist(nel) with
*     first index running fastest within the bounds given by 
*     imn(nel),imx(nel)
*     clone of next_dist
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     inc, nel, imn(nel) ,imx(nel)
      integer, intent(inout) ::
     &     idist(nel)

      integer ::
     &     iel

      next_dist2 = .false.
      if (inc.gt.0) then
        do iel = 1, nel
          if (idist(iel).ge.imx(iel)) then
            idist(iel) = imn(iel)
          else
            idist(iel) = idist(iel)+inc
            next_dist2 = .true.
            exit
          end if
        end do
      else if (inc.lt.0) then
        do iel = 1, nel
          if (idist(iel).le.imn(iel)) then
            idist(iel) = imx(iel)
          else
            idist(iel) = idist(iel)+inc
            next_dist2 = .true.
            exit
          end if
        end do
      end if
        
      return
      end
