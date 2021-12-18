*----------------------------------------------------------------------*
      integer function idxssd(idspc,iyssg,nelmax,nspc)
*----------------------------------------------------------------------*
*     get index of subspace distribution:
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     nelmax, nspc,
     &     idspc(nelmax), iyssg(nelmax,nspc)

      integer ::
     &     iel, idx

      idx = 0
      ! lexical index of previous subspace distribution
      do iel = 1, nelmax
        idx = idx + iyssg(iel,idspc(iel))
      end do

      idxssd = idx+1
      
      return
      end
