      subroutine topo_rename_vtxs(svertex,vtx,
     &     idx_svtx_new,idx_vtx_new,idx_pack_new,mode,
     &     vtx_list,nvtx,nlist)
*
*     mode =  0: use idx_vtx_new
*     mode != 0: use idx_pack_new

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     idx_svtx_new, idx_vtx_new, nvtx, nlist, vtx_list(nlist),
     &     mode
      integer, intent(out) ::
     &     svertex(nvtx)
      integer(8), intent(in) ::
     &     idx_pack_new(nlist)
      integer(8), intent(out) ::
     &     vtx(nvtx)

      integer ::
     &     idx, ivtx

      integer(8), external ::
     &     pack_vtx

      if (mode.eq.0) then
        do idx = 1, nlist
          ivtx = vtx_list(idx)
          svertex(ivtx) = idx_svtx_new
          vtx(ivtx) = pack_vtx(idx_vtx_new,idx,.false.)
        end do
      else
        do idx = 1, nlist
          ivtx = vtx_list(idx)
          svertex(ivtx) = idx_svtx_new
          vtx(ivtx) = idx_pack_new(idx)
        end do
      end if

      return
      end
