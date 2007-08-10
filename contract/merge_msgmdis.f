*----------------------------------------------------------------------*
      subroutine merge_msgmdis(msdis_r,gamdis_r,
     &                         nel_r,
     &                         msdis_1,gamdis_1,
     &                         msdis_2,gamdis_2,
     &                         map_inf)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     msdis_1(*), msdis_2(*),
     &     gamdis_1(*), gamdis_2(*),
     &     map_inf(*), nel_r
      integer, intent(out) ::
     &     msdis_r(*), gamdis_r(*)

      integer ::
     &     idx, jdx, nel, iel
      
c dbg
c      print *,'nel_r = ',nel_r
c      print *,'map_inf= ',map_inf(1:4)
c dbg
      idx = 0
      do jdx = 1, nel_r
        idx = idx+1
        nel = map_inf(idx)
        msdis_r(jdx)  = 0
        gamdis_r(jdx) = 1
        do iel = 1, nel
          idx = idx+1
c dbg
c          print *,'1 jdx,idx, map, ms, gam: ',jdx,idx,map_inf(idx)
c          print *,'                       ',
c     &         msdis_1(map_inf(idx)), gamdis_1(map_inf(idx))
c dbg
          msdis_r(jdx) = msdis_r(jdx) +
     &         msdis_1(map_inf(idx))
          gamdis_r(jdx) = multd2h(gamdis_r(jdx),
     &         gamdis_1(map_inf(idx)))
        end do
        idx = idx+1
        nel = map_inf(idx)
        do iel = 1, nel
          idx = idx+1
c dbg
c          print *,'2 jdx,idx, map, ms, gam: ',jdx,idx,map_inf(idx)
c          print *,'                       ',
c     &         msdis_2(map_inf(idx)), gamdis_2(map_inf(idx))
c dbg
          msdis_r(jdx) = msdis_r(jdx) + msdis_2(map_inf(idx))
          gamdis_r(jdx) = multd2h(gamdis_r(jdx),
     &         gamdis_2(map_inf(idx)))
        end do
      end do
      
      return
      end
