*----------------------------------------------------------------------*
      subroutine contr_clean_arcs(arc,narc)
*----------------------------------------------------------------------*
*     clean up arc and xarc information:
*      look for arcs with identical links and sum up
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(inout) ::
     &     narc
      type(cntr_arc) ::
     &     arc(narc)

      integer ::
     &     iarc, ivtx1, ivtx2, jarc, jvtx1, jvtx2, narc_count

      do iarc = 1, narc
        ivtx1 = arc(iarc)%link(1)
        ivtx2 = arc(iarc)%link(2)
        if (ivtx1.eq.0.and.ivtx2.eq.0) cycle
        if (arc(iarc)%occ_cnt(1,1).lt.0) cycle ! ignore proto-arcs
        do jarc = iarc+1, narc
          jvtx1 = arc(jarc)%link(1)
          jvtx2 = arc(jarc)%link(2)
          if (ivtx1.ne.jvtx1.or.ivtx2.ne.jvtx2) cycle
          if (arc(jarc)%occ_cnt(1,1).lt.0) cycle ! ignore proto-arcs
          ! equal arcs:
          ! sum up
          arc(iarc)%occ_cnt =
     &         arc(iarc)%occ_cnt + arc(jarc)%occ_cnt
          ! mark as deleted
          arc(jarc)%link = 0
        end do
      end do

      narc_count = 0
      ! remove deleted arcs
      do iarc = 1, narc
        ivtx1 = arc(iarc)%link(1)
        ivtx2 = arc(iarc)%link(2)
        if (ivtx1.eq.0.and.ivtx2.eq.0) cycle
        narc_count = narc_count+1
        if (iarc.gt.narc_count)
     &       arc(narc_count) = arc(iarc)
      end do

      narc = narc_count
      
      end
