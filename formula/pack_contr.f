*----------------------------------------------------------------------*
      subroutine pack_contr(svertex,vtx,topo,xlines,contr,nj_res)
*----------------------------------------------------------------------*
*          
*     transform contraction info on contr in into:
*         svertex               - which super-vertex the vertex belongs to
*         vtx                   - condensed operator/block/adj info
*         topo                  - the contraction pattern
*         xlines (if nj_res>0)  - 
*          
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, intent(in) ::
     &     nj_res
      type(contraction), intent(in) ::
     &     contr
      integer, intent(out) ::
     &     svertex(contr%nvtx)
      integer(8), intent(out) ::
     &     vtx(contr%nvtx), topo(contr%nvtx,contr%nvtx),
     &     xlines(contr%nvtx,nj_res)

      integer ::
     &     nvtx, narc, nxarc, ivtx, jvtx, iarc, iblk_op, idx_op,
     &     iadj
      integer(8) ::
     &     base

      integer ::
     &     occ(ngastp,2)

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(cntr_arc), pointer ::
     &     arc(:), xarc(:)

      integer(8), external ::
     &     int8_pack

      base =  pack_base
      
      nvtx  = contr%nvtx
      narc  = contr%narc
      nxarc = contr%nxarc

      vertex => contr%vertex
      arc    => contr%arc
      xarc   => contr%xarc

      ! -----------------------
      ! make packed vertex list
      ! -----------------------
      do ivtx = 1, nvtx
        idx_op = vertex(ivtx)%idx_op
        iblk_op = vertex(ivtx)%iblk_op
        iadj = 0
        if (vertex(ivtx)%dagger) iadj = 1
        vtx(ivtx) = sign(iadj*(base**4)+abs(idx_op)*(base**2)+iblk_op,
     &                   idx_op)
        if (iblk_op.gt.base*base) then
          write(luout,*) 'base    = ',base*base
          write(luout,*) 'iblk_op = ',iblk_op
          call quit(1,'contr2topo',
     &         'incredibly large iblk_op encountered')
        end if
      end do

      ! -----------------------------
      ! make packed contraction table
      ! -----------------------------
      topo = 0
      do iarc = 1, narc
        ivtx = arc(iarc)%link(1)
        jvtx = arc(iarc)%link(2)
        occ  = arc(iarc)%occ_cnt
        topo(ivtx,jvtx) = int8_pack(occ,ngastp*2,base)
        occ  = iocc_dagger(occ)
        topo(jvtx,ivtx) = int8_pack(occ,ngastp*2,base)
      end do

      svertex(1:nvtx) = contr%svertex(1:nvtx)

      if (nj_res.le.0) return

      ! -------------------------------
      ! make packed external line table
      ! -------------------------------
      xlines = 0
      do iarc = 1, nxarc
        ivtx = xarc(iarc)%link(1)
        jvtx = xarc(iarc)%link(2)
        occ  = xarc(iarc)%occ_cnt
        xlines(ivtx,jvtx) = int8_pack(occ,ngastp*2,base)
      end do

      return
      end
