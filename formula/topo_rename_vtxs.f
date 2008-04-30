      subroutine topo_rename_vtxs(svertex,vtx,
     &     idx_svtx_new,idx_vtx_new,
     &     vtx_list,nvtx,nlist)

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     idx_svtx_new, idx_vtx_new, nvtx, nlist, vtx_list(nlist)
      integer, intent(out) ::
     &     svertex(nvtx)
      integer(8), intent(out) ::
     &     vtx(nvtx)

      integer ::
     &     idx, ivtx

      do idx = 1, nlist
        ivtx = vtx_list(idx)
        svertex(ivtx) = idx_svtx_new
        vtx(ivtx) =
     &       sign(abs(idx_vtx_new)*pack_base*pack_base+idx,idx_vtx_new)
      end do

      return
      end
