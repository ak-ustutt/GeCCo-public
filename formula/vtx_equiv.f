*----------------------------------------------------------------------*
      logical function vtx_equiv(ivtx,jvtx,vtx,svtx,topo,xlines,nvtx,nj)
*----------------------------------------------------------------------*
*     are vertices ivtx, jvtx equivalent?
*     we call them equivalent, if after removing either of them
*     the same diagram results
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ivtx, jvtx, nvtx, nj, svtx(nvtx)
      integer(8), intent(inout) ::
     &     vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,nj)

      integer(8) ::
     &     vtx1(nvtx-1), topo1(nvtx-1,nvtx-1), xlines1(nvtx-1,nj+1),
     &     vtx2(nvtx-1), topo2(nvtx-1,nvtx-1), xlines2(nvtx-1,nj+1)

      integer ::
     &     idx, jdx, iidx, jjdx, ireo(nvtx-1),
     &     svtx1(nvtx-1), svtx2(nvtx-1)

      integer, external ::
     &     i8list_cmp

      ! probably do this test on before calling me:
      vtx_equiv = .false.
      if (vtx(ivtx).ne.vtx(jvtx)) return

      ! remove ivtx
      iidx = 0
      do idx = 1, nvtx
        if (idx.eq.ivtx) cycle
        iidx = iidx+1
        vtx1(iidx) = vtx(idx)
        if (svtx(idx).gt.svtx(ivtx)) then
          svtx1(iidx) = svtx(idx) - 1
        else
          svtx1(iidx) = svtx(idx)
        end if
        xlines1(iidx,1:nj) = xlines(iidx,1:nj)
        xlines1(iidx,nj+1) = topo(iidx,ivtx)
        jjdx = 0
        do jdx = 1, nvtx
          if (idx.eq.ivtx) cycle
          jjdx = jjdx+1
          topo1(jjdx,iidx) = topo(jdx,idx)
        end do
      end do

      ! remove jvtx
      iidx = 0
      do idx = 1, nvtx
        if (idx.eq.jvtx) cycle
        iidx = iidx+1
        vtx2(iidx) = vtx(idx)
        if (svtx(idx).gt.svtx(jvtx)) then
          svtx2(iidx) = svtx(idx) - 1
        else
          svtx2(iidx) = svtx(idx)
        end if
        xlines2(iidx,1:nj) = xlines(iidx,1:nj)
        xlines2(iidx,nj+1) = topo(iidx,jvtx)
        jjdx = 0
        do jdx = 1, nvtx
          if (idx.eq.jvtx) cycle
          jjdx = jjdx+1
          topo2(jjdx,iidx) = topo(jdx,idx)
        end do
      end do

      call topo_make_unique2(ireo,vtx1,svtx1,topo1,xlines1,nvtx-1,nj+1)
      call topo_make_unique2(ireo,vtx2,svtx2,topo2,xlines2,nvtx-1,nj+1)

      vtx_equiv = i8list_cmp(vtx1,vtx2,nvtx-1).eq.0
      if (ntest.ge.100) write(luout,*) 'vtx_equiv > (1): ',vtx_equiv
      vtx_equiv = vtx_equiv.and.
     &            i8list_cmp(xlines1,xlines2,(nvtx-1)*(nj+1)).eq.0
      if (ntest.ge.100) write(luout,*) 'vtx_equiv > (2): ',vtx_equiv
      vtx_equiv = vtx_equiv.and.
     &            i8list_cmp(topo1,topo2,(nvtx-1)*(nvtx-1)).eq.0
      if (ntest.ge.100) write(luout,*) 'vtx_equiv > (3): ',vtx_equiv

      return
      end
