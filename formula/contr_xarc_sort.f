*----------------------------------------------------------------------*
      subroutine contr_xarc_sort(contr)
*----------------------------------------------------------------------*
*     if possible, reorder vertices such that they contribute to xarcs
*     in ascending order
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(inout), target ::
     &     contr

      integer ::
     &     nvtx, narc, nxarc, nj, ivtx, iarc, hlp
      integer ::
     &     svertex_hlp(contr%nvtx)
      type(cntr_arc), pointer ::
     &     xarc(:), arc(:)
      type(cntr_vtx), pointer ::
     &     vtx_hlp(:)
      logical ::
     &     changed, swap
      integer, pointer ::
     &     svmap(:,:), xl4vtx(:), ireo(:), ireo2(:)
      
      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,'contr_xarc_sort')


      nxarc = contr%nxarc

      ! if there are no xarcs, we can quit
      if (nxarc.eq.0) return

      nvtx = contr%nvtx

      ! only one vertex: nothing to do
      if (nvtx.le.1) return

      nj = 0
      xarc => contr%xarc
      do iarc = 1, nxarc
        nj = max(nj,xarc(iarc)%link(2))
      end do

      ! also this is an uninteresting case
      if (nj.le.1) return

      narc = contr%narc
      arc => contr%arc
      

      allocate(xl4vtx(nvtx),svmap(nvtx,nj),ireo(nvtx),ireo2(nvtx))

      xl4vtx(1:nvtx) = 0
      svmap(1:nvtx,1:nj) = 0
      do iarc = 1, nxarc
        ivtx = xarc(iarc)%link(1)
        xl4vtx(ivtx) = xl4vtx(ivtx)+1
        svmap(ivtx,xl4vtx(ivtx)) = xarc(iarc)%link(2)
      end do

      do ivtx = 1, nvtx
        if (xl4vtx(ivtx).gt.2)
     &   call quit(1,'contr_xarc_sort','xl4vtx(ivtx).gt.2: Test this!')
        if (xl4vtx(ivtx).eq.2) then
          if (svmap(ivtx,2).le.svmap(ivtx,1))
     &   call quit(1,'contr_xarc_sort','not ascending svmap!')
        end if
      end do

      do ivtx = 1, nvtx
        ireo(ivtx) = ivtx
      end do

      changed = .true.
      do while (changed)
!     dbg
!        write(lulog,*) 'current ireo: ',ireo
!     dbg
        
        changed = .false.
        do ivtx = 2, nvtx
          if (svmap(ireo(ivtx-1),1).eq.svmap(ireo(ivtx),1)
     &         .and.svmap(ireo(ivtx),1)==0) cycle
          swap = svmap(ireo(ivtx-1),1).gt.svmap(ireo(ivtx),1)
          if(.not.swap.and.svmap(ireo(ivtx-1),1).eq.svmap(ireo(ivtx),1))
     &         then
            swap = svmap(ireo(ivtx-1),2).gt.svmap(ireo(ivtx),2)
          end if
!     dbg
!          if (swap) write(lulog,*) 'request swap of ',
!     &         ireo(ivtx-1),ireo(ivtx)
!     dbg
          if (swap.and.may_commute(ireo(ivtx-1),ireo(ivtx))) then
            changed = .true.
!     dbg
!            write(lulog,*) 'request accepted'
!     dbg
            
            hlp = ireo(ivtx)
            ireo(ivtx)=ireo(ivtx-1)
            ireo(ivtx-1) = hlp
          end if
        end do
      end do

      ! note: ireo(idx) = oldidx
      !   set ireo2(oldidx) = idx
      do ivtx = 1, nvtx
        ireo2(ireo(ivtx)) = ivtx
      end do
!     dbg
!      write(lulog,*) 'final reo:  ',ireo
!      write(lulog,*) 'final reo2: ',ireo2
!     dbg

      allocate(vtx_hlp(nvtx))
      vtx_hlp = contr%vertex

      do ivtx = 1, nvtx
        contr%vertex(ivtx) = vtx_hlp(ireo(ivtx))
      end do

      deallocate(vtx_hlp)

      svertex_hlp = contr%svertex
      do ivtx = 1, nvtx
        contr%svertex(ivtx) = svertex_hlp(ireo(ivtx))
      end do

      do iarc = 1, narc
        contr%arc(iarc)%link(1) = ireo2(contr%arc(iarc)%link(1))
        contr%arc(iarc)%link(2) = ireo2(contr%arc(iarc)%link(2))
      end do

      do iarc = 1, nxarc
        contr%xarc(iarc)%link(1) = ireo2(contr%xarc(iarc)%link(1))
      end do

      call update_svtx4contr(contr)

      deallocate(svmap,xl4vtx,ireo,ireo2)
      
      return

      contains

      logical function may_commute(ivtx1,ivtx2)

      implicit none
      integer ivtx1, ivtx2

      integer iarc
      
      may_commute = .true.
      do iarc = 1, narc
!     dbg
!        write(lulog,*) 'checking: ',iarc
!     dbg
        
        if (arc(iarc)%link(1)==ivtx1 .and.
     &       arc(iarc)%link(2)==ivtx2) then
          may_commute = .false.
          exit
        end if
      end do

      return

      end function
      
      end
