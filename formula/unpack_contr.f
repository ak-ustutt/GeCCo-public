      subroutine unpack_contr(contr,svertex,vtx,topo,xlines,nvtx,nj_res)

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, intent(in) ::
     &     nj_res, nvtx
      type(contraction), intent(inout) ::
     &     contr
      integer, intent(in) ::
     &     svertex(nvtx)
      integer(8), intent(in) ::
     &     vtx(nvtx), topo(nvtx,nvtx),
     &     xlines(nvtx,nj_res)

      integer ::
     &     narc, nxarc, iarc, ivtx, jvtx, ij,
     &     idx_op, iblk_op, icnt, iadj
      integer(8) ::
     &     avtx, base

      integer ::
     &     occ(ngastp,2)

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(cntr_arc), pointer ::
     &     arc(:), xarc(:)

      integer, external ::
     &     int8_expand
      
      base = pack_base

      ! count arcs
      narc = 0
      do jvtx = 1, nvtx
        do ivtx = 1, jvtx
          if (topo(ivtx,jvtx).gt.0) narc = narc+1
        end do
      end do
      ! count xarcs
      nxarc = 0
      do ij = 1, nj_res
        do ivtx = 1, nvtx
          if (xlines(ivtx,ij).gt.0) nxarc = nxarc+1
        end do
      end do

      call resize_contr(contr,nvtx,narc,nxarc,0)

      vertex => contr%vertex
      arc    => contr%arc
      xarc   => contr%xarc
      
      contr%nvtx  = nvtx
      contr%narc  = narc
      contr%nxarc = nxarc

      ! unpack vertices
      do ivtx = 1, nvtx
        avtx = abs(vtx(ivtx))
        iadj    = avtx/(base**6)
        avtx = mod(avtx,(base**6))
        idx_op  = sign(avtx/(base*base),vtx(ivtx))
        iblk_op = mod(avtx,base*base)
        vertex(ivtx)%idx_op = idx_op
        vertex(ivtx)%iblk_op = iblk_op
        vertex(ivtx)%dagger  = iadj.eq.1
      end do

      ! unpack arcs
      iarc = 0
      do jvtx = 1, nvtx
        do ivtx = 1, jvtx
          if (topo(ivtx,jvtx).gt.0) then
            iarc = iarc+1
            occ = 0
            icnt = int8_expand(topo(ivtx,jvtx),base,occ)
            arc(iarc)%link(1) = ivtx
            arc(iarc)%link(2) = jvtx
            arc(iarc)%occ_cnt = occ
          end if
        end do
      end do
      
      ! unpack xarcs
      iarc = 0
      do ij = 1, nj_res
        do ivtx = 1, nvtx
          if (xlines(ivtx,ij).gt.0) then
            iarc = iarc+1
            occ = 0
            icnt = int8_expand(xlines(ivtx,ij),base,occ)
            xarc(iarc)%link(1) = ivtx
            xarc(iarc)%link(2) = ij
            xarc(iarc)%occ_cnt = occ
          end if
        end do
      end do

      ! store svertex
      do ivtx = 1, nvtx
        contr%svertex(ivtx) = svertex(ivtx)
      end do

      ! update other svertex info
      call update_svtx4contr(contr)
      
      return
      end
      
