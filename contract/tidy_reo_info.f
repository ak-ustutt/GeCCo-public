      subroutine tidy_reo_info(reo_info)
*----------------------------------------------------------------------*
*     fix routine: check whether the reo's initiated by
*     reorder_supvtx and reorder_supvtx_x are compatible and
*     merge reo's between identical nodes ...
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_reorder_info.h'

      type(reorder_info), intent(inout) ::
     &     reo_info

      integer ::
     &     ireo, jreo, isvtx, jsvtx, ito, jto, ifrom, jfrom,
     &     nreo, nreo_new

      if (reo_info%nreo.le.1) return

      nreo = reo_info%nreo
      nreo_new = nreo
      do ireo = 1, nreo
        isvtx = reo_info%reo(ireo)%idxsuper
        if (isvtx.le.0) cycle

        ito   = reo_info%reo(ireo)%to
        ifrom = reo_info%reo(ireo)%from

        do jreo = ireo+1, nreo
          jsvtx = reo_info%reo(jreo)%idxsuper
          if (isvtx.ne.jsvtx) cycle

          jto   = reo_info%reo(jreo)%to
          if (ito.ne.jto) cycle

          jfrom = reo_info%reo(jreo)%from
          if (ifrom.ne.jfrom) cycle
          
          reo_info%reo(ireo)%occ_shift =
     &         reo_info%reo(ireo)%occ_shift +
     &         reo_info%reo(jreo)%occ_shift
          nreo_new = nreo_new-1
          reo_info%reo(jreo)%idxsuper = -1 ! mark as deleted
c dbg
          print *,'merging: ',ireo,jreo
c dbg

        end do

      end do

      ! remove deleted items

      jreo = 0
      do ireo = 1, nreo
        if (reo_info%reo(ireo)%idxsuper.le.0) cycle

        jreo = jreo+1

        if (jreo.lt.ireo)
     &       reo_info%reo(jreo) = reo_info%reo(ireo)

      end do

      reo_info%nreo = nreo_new

c dbg
      print *,'tidy_reo: ',nreo,' -> ',nreo_new
c dbg

      return
      end
