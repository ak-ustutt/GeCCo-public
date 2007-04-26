*----------------------------------------------------------------------*
      integer function imltlist(inum,ilist,nel,inc)
*----------------------------------------------------------------------*
*     return multiplicity of integer inum in list ilist
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     inum, nel, inc, ilist(nel)

      integer ::
     &     idx, ipos

      ipos = 1
      imltlist = 0
      do idx = 1, nel
        if (ilist(ipos).eq.inum) imltlist = imltlist+1
        ipos = ipos+inc
      end do

      return
      end

