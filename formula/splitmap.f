*----------------------------------------------------------------------*
      subroutine splitmap(splmap,nvtx_rem,
     &                    vtxmap,vtxinf,topomap,nvtx)
*----------------------------------------------------------------------*
*     on entry:
*       vtxmap(nvtx) contains -1 if the vertex remains after splitting
*             off a part of the contraction and non-0 (number of super
*             vertex of the intermediate that was split off) if this
*             vertex is removed
*     vtxinf(nvtx)  unique op+block identifier (orig. order)
*     topomap(nvtx,nvtx)  topology of contraction
*
*     splmap(nvtx) contains the place in the fragment contraction
*             after splitting; a negative number is the position of
*             an empty vertex; several vertices of the split-off part
*             may be merged into a single vertex (dep. on type of
*             intermediate.
*
*     the intent of this routine is to group the remaining vertices
*     correctly around the empty vertices such that in later steps
*     the correct connections can be made
*
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     nvtx
      integer, intent(out) ::
     &     splmap(nvtx), nvtx_rem
      integer, intent(in) ::
     &     vtxmap(nvtx), vtxinf(nvtx), topomap(nvtx,nvtx)

      integer ::
     &     ivtx_rem, ivtx_spl_last, isupervtx_spl_last,
     &     ivtx, jvtx, kvtx, lvtx, nsuper_spl,  isuper, ncontr
      integer ::
     &     ireo(nvtx), topomap_cp(nvtx,nvtx), vtxmap_cp(nvtx)

      integer, external ::
     &     imltlist

c dbg
c      print *,'in splitmap'
c      print *,'vtxmap: ',vtxmap(1:nvtx)
c dbg

c      ! we start by setting up splmap for the given ordering
      ivtx_rem = 0
      ivtx_spl_last = -1
      isupervtx_spl_last = -1
      do ivtx = 1, nvtx
        if (vtxmap(ivtx).ge.0) then
          if (vtxmap(ivtx).gt.isupervtx_spl_last) then
            ivtx_rem = ivtx_rem + 1
            isupervtx_spl_last = vtxmap(ivtx)
            ivtx_spl_last = ivtx_rem
          end if
        else
          ivtx_rem = ivtx_rem + 1
        end if        
      end do

      nvtx_rem = ivtx_rem

      nsuper_spl = isupervtx_spl_last

      ! next we set up a reordering array such that for each
      ! super-vertex (which is represented by a single empty vertex
      ! in the remainder contraction) excitation-connections go only
      ! to vertices above and deexcitaion connections go only to 
      ! vertices below

      do ivtx = 1, nvtx
        ireo(ivtx) = ivtx
      end do

      topomap_cp = topomap
      vtxmap_cp  = vtxmap

c dbg
c      call iwrtma(vtxmap_cp,nvtx,1,nvtx,1)
c      call iwrtma(topomap_cp,nvtx,nvtx,nvtx,nvtx)
c dbg

      do isuper = 1, nsuper_spl
        ! how may vertices (with external connections) contribute?
        ncontr = imltlist(isuper,vtxmap,nvtx,1)
        if (ncontr.le.1) cycle ! no problem

c dbg
c        print *,'caring for isuper = ',isuper,' ncontr = ',ncontr
c dbg
        
        ! loop over contributing vertices
        do ivtx = 1, nvtx
          if (vtxmap_cp(ivtx).ne.isuper) cycle
c dbg
c          print *,'contributing vertex: ',ivtx
c dbg
          ! vertices above: shift up as much as possible
          do jvtx = 1, ivtx-1
            if (vtxmap_cp(jvtx).ge.0) cycle ! only connections to other vtx's
            if (topomap_cp(jvtx,ivtx).ne.0) then
c dbg
c              print *,'1> connected to ',jvtx
c dbg
              lvtx = jvtx
              do kvtx = jvtx-1, 1, -1
c dbg
c                print *,'kvtx,jvtx,topo: ',
c     &               kvtx,jvtx,topomap_cp(kvtx,jvtx)
c dbg
                if (topomap_cp(kvtx,jvtx).ne.0) then
c                  lvtx = kvtx+1
                  exit
                end if
                lvtx = lvtx-1
              end do
              call shift_vtx(jvtx,lvtx)
            end if
          end do
          ! vertices below: shift down as much as possible
          do jvtx = nvtx, ivtx+1, -1
            if (vtxmap_cp(jvtx).ge.0) cycle ! only connections to other vtx's
            if (topomap_cp(jvtx,ivtx).ne.0) then
c dbg
c              print *,'2> connected to ',jvtx
c dbg
              lvtx = jvtx
              do kvtx = jvtx+1, nvtx
                if (topomap_cp(kvtx,jvtx).ne.0) then
c                  lvtx = kvtx-1
                  exit
                end if
                lvtx = lvtx+1
              end do
              call shift_vtx(jvtx,lvtx)
            end if
          end do
        end do

      end do

c dbg
c      print *,'suggested reordering'
c      call iwrtma(ireo,nvtx,1,nvtx,1)
c dbg

      ! we finish by setting up splmap for the given ordering
      ivtx_rem = 0
      ivtx_spl_last = -1
      isupervtx_spl_last = -1
      do ivtx = 1, nvtx
        jvtx = ireo(ivtx)
        if (vtxmap(jvtx).ge.0) then
          if (vtxmap(jvtx).gt.isupervtx_spl_last) then
            ivtx_rem = ivtx_rem + 1
            splmap(jvtx) = -ivtx_rem
            isupervtx_spl_last = vtxmap(jvtx)
            ivtx_spl_last = ivtx_rem
          else
            splmap(jvtx) = -ivtx_spl_last
          end if
        else
          ivtx_rem = ivtx_rem + 1
          splmap(jvtx) = ivtx_rem
        end if        
      end do
c dbg
c      print *,'splmap: ',splmap
c dbg

      return

      contains
      
      subroutine shift_vtx(ivtx_old,ivtx_new)

      implicit none

      integer ::
     &     ivtx_old, ivtx_new
      integer ::
     &     inc, idx, ihelp, ivhelp(nvtx)

      if (ivtx_old.eq.ivtx_new) return
c dbg
c      if (kvtx.ne.lvtx) print *,'shifting: ',ivtx_old,ivtx_new
c dbg

      inc = +1
      if (ivtx_old.gt.ivtx_new) inc = -1

      ! shift in reo array
      ihelp = ireo(ivtx_old)
      do idx = ivtx_old, ivtx_new-inc, inc
        ireo(idx) = ireo(idx+inc)        
      end do
      ireo(ivtx_new) = ihelp

      ! shift in vtxmap_cp
      ihelp = vtxmap_cp(ivtx_old)
      do idx = ivtx_old, ivtx_new-inc, inc
        vtxmap_cp(idx) = vtxmap_cp(idx+inc)        
      end do
      vtxmap_cp(ivtx_new) = ihelp

      ! shift in topomap_cp
      ! rows
      ivhelp(1:nvtx) = topomap_cp(1:nvtx,ivtx_old)
      do idx = ivtx_old, ivtx_new-inc, inc
        topomap_cp(1:nvtx,idx) = topomap_cp(1:nvtx,idx+inc)
      end do
      topomap_cp(1:nvtx,ivtx_new) = ivhelp(1:nvtx)

      ! cols
      ivhelp(1:nvtx) = topomap_cp(ivtx_old,1:nvtx)
      do idx = ivtx_old, ivtx_new-inc, inc
        topomap_cp(idx,1:nvtx) = topomap_cp(idx+inc,1:nvtx)
      end do
      topomap_cp(ivtx_new,1:nvtx) = ivhelp(1:nvtx)

c dbg
c      print *,'new ireo:'
c      call iwrtma(ireo,nvtx,1,nvtx,1)
c      print *,'new vtxmap:'
c      call iwrtma(vtxmap_cp,nvtx,1,nvtx,1)
c      print *,'new topomap:'
c      call iwrtma(topomap_cp,nvtx,nvtx,nvtx,nvtx)
c dbg

      return
      end subroutine

      end
