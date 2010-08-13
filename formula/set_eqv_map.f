*----------------------------------------------------------------------*
      subroutine set_eqv_map(neqv,idx_eqv,
     &                       vtx,svertex,topo,xlines,nvtx,nj)
*----------------------------------------------------------------------*
*     generate equivalence map needed for topology analysis
*
*     neqv(nvtx) :          number of equivalent vertices (i.e. same
*                           op+blk and commuting) for first occurrence
*                           -1 for second etc occurrence of vertex
*
*     idx_eqv(nvtx,nvtx):   idx_eqv(1:neqv(ivtx),ivtx) contains
*                           the equivalent vertices (ivtx is first
*                           occurrence of vertex)
*
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nvtx, nj, svertex(nvtx)
      integer(8), intent(in) ::
     &     topo(nvtx,nvtx),xlines(nvtx,nj),vtx(nvtx)
      integer, intent(out) ::
     &     neqv(nvtx), idx_eqv(nvtx,nvtx)

      integer ::
     &     ivtx, jvtx
      integer, pointer ::
     &     occ_cnt(:,:)

      logical, external ::
     &     vtx_equiv

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'this is set_eqv_map')
      end if
      
      neqv(1:nvtx) = 1
      do ivtx = 1, nvtx
        idx_eqv(1,ivtx) = ivtx
        idx_eqv(2:nvtx,ivtx) = 0
      end do
      
      do ivtx = 1, nvtx
        do jvtx = 1, ivtx-1
          if (neqv(jvtx).lt.0) cycle
          ! vertices equivalent?
          if (vtx(ivtx).eq.vtx(jvtx).and. ! <-avoid too many calls to vtx_equiv
     &        vtx_equiv(ivtx,jvtx,vtx,svertex,topo,xlines,nvtx,nj)) then
            neqv(jvtx) = neqv(jvtx)+1
            neqv(ivtx) = -1
            idx_eqv(neqv(jvtx),jvtx) = ivtx
          end if          
        end do
      end do

      if (ntest.ge.100) then
        write(luout,*) 'neqv:'
        call iwrtma(neqv,1,nvtx,1,nvtx)
        write(luout,*) 'idx_eqv:'
        call iwrtma(idx_eqv,nvtx,nvtx,nvtx,nvtx)
      end if

      return
      end
