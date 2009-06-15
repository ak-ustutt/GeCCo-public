*----------------------------------------------------------------------*
      subroutine topo_remove_arcs(topo,nvtx,list,nlist)
*----------------------------------------------------------------------*
*     carry out contraction: remove all arcs from the list
*       (i.e. zero the respective entries in topo)
*       at the same time, get the sign associated with this
*       event
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     nvtx, nlist,
     &     list(2,nlist)
      integer(8), intent(inout) ::
     &     topo(nvtx,nvtx)

      integer ::
     &     idx, ivtx, jvtx

      do idx = 1, nlist
        ivtx = list(1,idx)
        jvtx = list(2,idx)
        topo(ivtx,jvtx) = 0
        topo(jvtx,ivtx) = 0
      end do

      end 
