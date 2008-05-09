*----------------------------------------------------------------------*
      integer function i8mltlist(inum,ilist,nel,inc)
*----------------------------------------------------------------------*
*     return multiplicity of integer inum in list ilist
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     inum, nel, inc
      integer(8), intent(in) ::
     &     ilist(nel)

      integer ::
     &     idx, ipos

      ipos = 1
      i8mltlist = 0
      do idx = 1, nel
        if (ilist(ipos).eq.inum) i8mltlist = i8mltlist+1
        ipos = ipos+inc
      end do

      return
      end

