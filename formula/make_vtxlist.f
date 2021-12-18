      subroutine make_vtxlist(vtx_list,nvtx_cnt,svertex,nvtx)

      implicit none

      integer, intent(in) ::
     &     nvtx, svertex(nvtx)
      integer, intent(inout) ::
     &     nvtx_cnt, vtx_list(*)

      logical ::
     &     new
      integer ::
     &     nvtx_cnt_new, idx, kdx, ivtx, jvtx, kvtx, isvtx

      nvtx_cnt_new = nvtx_cnt
      do idx = 1, nvtx_cnt
        ivtx = vtx_list(idx)
        isvtx = svertex(ivtx)
        ! scan svertex and look for joined vertices
        jvtx = -1
        do kvtx = 1, nvtx
          if (kvtx.ne.ivtx.and.svertex(kvtx).eq.isvtx) then
            ! already present in vtx_list?
            new = .true.
            do kdx = 1, nvtx_cnt_new
              new = new.and.vtx_list(kdx).ne.kvtx
            end do
            if (new) then
              nvtx_cnt_new = nvtx_cnt_new+1
              vtx_list(nvtx_cnt_new) = kvtx
            end if
          end if
        end do
      end do

c dbg
c      if (nvtx_cnt.ne.nvtx_cnt_new)
c     &     print *,' it DID happen ...'
c dbg
      nvtx_cnt = nvtx_cnt_new

      return
      end
