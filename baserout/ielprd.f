
      integer function ielprd(ivec,nel)
      implicit none

      integer nel, ivec(nel)
      integer iel

      ielprd = 1
      do iel = 1, nel
        ielprd = ielprd*ivec(iel)
      end do
      
      return
      end
