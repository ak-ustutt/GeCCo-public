*----------------------------------------------------------------------*
      integer function rank_ivec(rank,ivec,nel)
*----------------------------------------------------------------------*
*     return on rank the rank of the entries on ivec and the 
*     index of the permutation as function result (ordering conforms
*     with the order produced by next_perm())
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel, ivec(nel)
      integer, intent(out) ::
     &     rank(nel)

      integer ::
     &     idx, jdx, value, rk

      rank_ivec = 0
      rank(1:nel) = 0
      do idx = 1, nel-1
        value = ivec(idx)
        rk = 0
        do jdx = idx+1, nel
          if (     (value.gt.ivec(jdx)) ) rk = rk+1
          if (.not.(value.gt.ivec(jdx)) ) rank(jdx) = rank(jdx)+1
        end do
        rank(idx) = rank(idx)+rk
        rank_ivec = rank_ivec*(nel-idx+1) + rk
      end do

      return
      end
