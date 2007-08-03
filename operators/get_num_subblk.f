*----------------------------------------------------------------------*
      subroutine get_num_subblk(csub,asub,occ,njoined)
*----------------------------------------------------------------------*
*     count number of non-zero occupations in occ
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'

      integer, intent(out) ::
     &     csub, asub
      integer, intent(in) ::
     &     njoined, occ(ngastp,2,njoined)

      integer ::
     &     ijoin, hpvx

      csub = 0
      asub = 0
      do ijoin = 1, njoined
        do hpvx = 1, ngastp
          if (occ(hpvx,1,ijoin).gt.0) csub = csub+1
          if (occ(hpvx,2,ijoin).gt.0) asub = asub+1
        end do
      end do
      
      return
      end
*----------------------------------------------------------------------*
