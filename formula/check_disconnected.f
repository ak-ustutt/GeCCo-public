*----------------------------------------------------------------------*
      subroutine check_disconnected(contr)
*----------------------------------------------------------------------*
*     check for disconnected vertices
*     if present, add zero contractions to all other vertices
*
*     only for use during evaluation of contraction
*       i.e. fact_costX(), frm_schedX()
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(inout) ::
     &     contr

      logical ::
     &     connected(contr%nvtx)
      integer ::
     &     nvtx, narc, iarc, ivtx, jvtx, nvtx_nc, narc_new
      type(cntr_arc), pointer ::
     &     arc(:)

c dbg
c      print *,'in check_disconnected'
c dbg
      if (contr%nsupvtx.le.1) return

      nvtx = contr%nvtx
      narc = contr%narc

      arc => contr%arc

      connected(1:nvtx) = .false.
      do iarc = 1, narc 
        connected(arc(iarc)%link(1)) = .true.
        connected(arc(iarc)%link(2)) = .true.
      end do

      ! do not consider vertices that are part of a supervertex
      do ivtx = 1, nvtx
        do jvtx = 1,nvtx
          if (connected(jvtx)) cycle
          if (contr%svertex(ivtx).eq.contr%svertex(jvtx)
     &         .and.connected(ivtx))
     &         connected(jvtx) = .true.
        end do
      end do
c dbg
c      print *,'connected: ',connected(1:nvtx)
c dbg

      nvtx_nc = 0
      narc_new = 0
      do ivtx = 1, nvtx
        if (.not.connected(ivtx)) then
          nvtx_nc = nvtx_nc+1
          narc_new = narc_new + nvtx-nvtx_nc
          ! exclude dummy arcs within supervertex
          do jvtx = 1, ivtx-1
            if (contr%svertex(ivtx).eq.contr%svertex(jvtx)) 
     &           narc_new = narc_new - 1
          end do
        end if
      end do
c dbg
c      print *,'adding arcs:: ',narc_new
c dbg

      if (narc_new.eq.0) return

      call resize_contr(contr,nvtx,narc+narc_new,0,contr%nfac)

      ! set link again !!
      arc => contr%arc

      contr%narc = narc+narc_new

      iarc = narc
      do ivtx = 1, nvtx
        if (.not.connected(ivtx)) then

          do jvtx = 1, nvtx
            if (ivtx.eq.jvtx) cycle
            if (jvtx.lt.ivtx.and..not.connected(jvtx)) cycle
            ! exclude dummy arcs within supervertex
            if (contr%svertex(ivtx).eq.contr%svertex(jvtx)) cycle
            iarc = iarc+1
c dbg
c            if (iarc.gt.narc+narc_new) stop 'cd: oha!'
c dbg
            arc(iarc)%link(1) = min(ivtx,jvtx)
            arc(iarc)%link(2) = max(ivtx,jvtx)
            arc(iarc)%occ_cnt = 0
          end do
        end if
      end do

      return
      end
