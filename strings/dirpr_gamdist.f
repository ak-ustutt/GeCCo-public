*----------------------------------------------------------------------*
      subroutine dirpr_gamdist(idgamout,idgamin1,idgamin2,nel)
*----------------------------------------------------------------------*
*     element-wise direct product of IRREP distributions with nel 
*     elements
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     nel, idgamin1(nel), idgamin2(nel)

      integer, intent(out) ::
     &     idgamout(nel)

      integer ::
     &     idx

      do idx = 1, nel
        idgamout(idx) = multd2h(idgamin1(idx),idgamin2(idx))
      end do

      return
      end
