      subroutine store_xarcs(xarcs,ixarc,
     &      ivtxder,maxarc,arc)
      ! slave routine for contr_deriv2
      
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, intent(inout) ::
     &     ixarc
      integer, intent(in) ::
     &     ivtxder,maxarc
      type(cntr_arc), intent(in) ::
     &     arc
      type(cntr_arc), intent(inout) ::
     &     xarcs(maxarc)

      logical ::
     &     dag_occ
      integer ::
     &     il
      integer ::
     &     occ_scr(ngastp,2)

      if (ixarc+2.gt.maxarc)
     &     call quit(1,'store_xarcs','maxarc seems too small!')
      ixarc = ixarc+1
      do il = 1, 2
        if (arc%link(il).eq.ivtxder) then
          ! still not decided to which vertex of the result
          ! we will contribute, so we set a zero:
          xarcs(ixarc)%link(2) = 0
          dag_occ = il.eq.1     ! transpose occupation
        else
          xarcs(ixarc)%link(1) =
     &         arc%link(il)
          ! shift by -1 as vertex ivtxder will be deleted
          if (xarcs(ixarc)%link(1).gt.ivtxder)
     &         xarcs(ixarc)%link(1) = xarcs(ixarc)%link(1)-1
        end if
      end do
      if (.not.dag_occ) then
        occ_scr = arc%occ_cnt
      else
        occ_scr = iocc_dagger(arc%occ_cnt)
      end if

      xarcs(ixarc)%occ_cnt(1:ngastp,1:2) = occ_scr
c                                ! store excitation part ...
c      xarcs(ixarc)%occ_cnt(1:ngastp,1:2) =
c     &     iocc_xdn(1,occ_scr)
c      ixarc = ixarc+1
c      xarcs(ixarc)%link = xarcs(ixarc-1)%link
c                                ! ... and de-excitation part
c      xarcs(ixarc)%occ_cnt(1:ngastp,1:2) =
c     &     iocc_xdn(2,occ_scr)
            

      return
      end
