*----------------------------------------------------------------------*
      integer function idxssg(nel,idspc,iyssg,iwssg,nelmax,nspc)
*----------------------------------------------------------------------*
*     get index of subspace graph:
*     idspc(1:nel) contains the subspace distribution of the
*     previous subspaces. Obtain the lexical index of that distribution
*     and add a further offset to be calculated from the weights
*     of nodes from which previous branch-offs to the current subspace
*     were possible
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nel, nelmax, nspc,
     &     idspc(nel+1), iyssg(nelmax,nspc), iwssg(0:nelmax,nspc)

      integer ::
     &     idx, iel, icspc, ispc, nelend

      if (ntest.ge.100) then
      ! error traps for debugging
        if (nel.lt.0) call quit(1,'idxssg','nel.lt.0')
        if (nel.ge.nelmax) call quit(1,'idxssg','nel.ge.nelmax')
        if (idspc(nel+1).le.0) call quit(1,'idxssg','idspc(nel+1).le.0')
        if (idspc(nel+1).gt.nspc) call quit(1,'idxssg',
     &       'idspc(nel+1).gt.nspc')
      end if

      idx = 0
      ! lexical index of previous subspace distribution
      do iel = 1, nel
        idx = idx + iyssg(iel,idspc(iel))
      end do

      ! current space
      icspc = idspc(nel+1)

      ! sum of weights of previous nodes of prev. space
      ! + those from which we could have
      ! branched to the current subspace give an additional offset
      do ispc = 2, icspc
        if (ispc.lt.icspc) then
          nelend=nelmax-1
        else
          nelend=nel-1
        end if
        do iel = 0, nelend
          ! a 0 indicates that no branching was possible
          if (iwssg(iel,icspc).ne.0) idx = idx+iwssg(iel,icspc-1)
        end do
      end do
      ! extra offset because of subspace 1 graph
      if (icspc.gt.1) idx = idx+1
      idxssg = idx+1
      
      return
      end
