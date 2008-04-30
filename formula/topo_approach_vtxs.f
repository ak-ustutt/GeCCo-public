      subroutine topo_approach_vtxs(ireo,
     &     svertex,vtx,topo,xlines,
     &     nvtx,nj,vtx_list,nlist)
      
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     nvtx, nj, nlist,
     &     vtx_list(nlist)
      integer, intent(inout) ::
     &     svertex(nvtx)
      integer(8), intent(inout) ::
     &     vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,nj)
      integer, intent(out) ::
     &     ireo(nvtx)
      
      integer ::
     &     iord(nvtx)
      integer ::
     &     idx, jdx, ivtxr, ivtxa, ivtxb, ivtxt, jvtx

      ! check list
      do idx = 1, nlist
        do jdx = 1, idx-1
          if (vtx_list(jdx).ge.vtx_list(idx)) then
            write(luout,*) '>',vtx_list(1:nlist)
            call quit(1,'topo_approach_vtxs',
     &           'a unique, ascending list was expected')

          end if
        end do
      end do
      
      do idx = 1, nvtx
        ireo(idx) = idx
        iord(idx) = idx
      end do

      ! loop over list
      do idx = 1, nlist
        ! reference point
        ivtxr = ireo(vtx_list(idx))
c dbg
c        print *,'ref = ',ivtxr,' <- ', vtx_list(idx)
c dbg
        ! vertices above reference
        do jdx = 1, idx-1
          ivtxa = ireo(vtx_list(jdx))
c dbg
c          print *,'a = ',ivtxa,' <- ', vtx_list(jdx)
c dbg
          if (ivtxr-ivtxa.gt.1) then
            ! try to move down
            ivtxt = ivtxa
            do jvtx = ivtxa+1, ivtxr-1
              ivtxt = jvtx
c              if (topo(jvtx,ivtxa).ne.0) then
              if (.not.may_commute(jvtx,ivtxa)) then
                ivtxt = jvtx-1
                exit
              end if
            end do
            call shift_vtx(ivtxa,ivtxt)
          end if
        end do
        ! vertices below reference
        do jdx = idx+1, nlist
          ivtxb = ireo(vtx_list(jdx))
c dbg
c          print *,'b = ',ivtxb,' <- ', vtx_list(jdx)
c dbg
          if (ivtxb-ivtxr.gt.1) then
            ! try to move up
            ivtxt = ivtxb
            do jvtx = ivtxb-1, ivtxr+1, -1
c dbg
c              print *,'jvtx,topo: ',jvtx,topo(jvtx,ivtxb)
c dbg
              ivtxt = jvtx
c              if (topo(jvtx,ivtxb).ne.0) then
              if (.not.may_commute(jvtx,ivtxb)) then
                ivtxt = jvtx+1
                exit
              end if
            end do
            call shift_vtx(ivtxb,ivtxt)
          end if
        end do

      end do

      contains

      logical function may_commute(ivtx1,ivtx2)

      implicit none

      integer, intent(in) ::
     &     ivtx1, ivtx2

      integer(8) ::
     &     scr(nj)

      integer, external ::
     &     i8mltlist

      may_commute = topo(ivtx1,ivtx2).eq.0 .and.     
     &              svertex(ivtx1).ne.svertex(ivtx2)
      call i8list_ovl(scr,xlines(ivtx1,1),xlines(ivtx2,1),nj,nvtx)
      may_commute = may_commute .and.
     &     i8mltlist(0,scr,nj,1).le.1

      return
      end function

      subroutine shift_vtx(ivtx_old,ivtx_new)

      implicit none

      integer ::
     &     ivtx_old, ivtx_new
      integer ::
     &     inc, idx, ij, ihelp
      integer(8) ::
     &     ivhelp(nvtx)

      if (ivtx_old.eq.ivtx_new) return
c dbg
c      if (ivtx_old.ne.ivtx_new) print *,'shifting: ',ivtx_old,ivtx_new
c dbg

      ! shift in ord array
      call shift_ivec(iord,ivtx_old,ivtx_new,nvtx)

      do idx = 1, nvtx
        ireo(iord(idx)) = idx
      end do

      ! shift in svertex
      call shift_ivec(svertex,ivtx_old,ivtx_new,nvtx)
      ! shift in vtx
      call shift_i8vec(vtx,ivtx_old,ivtx_new,nvtx)
      ! shift in xlines
      do ij = 1, nj
        call shift_i8vec(xlines(1:nvtx,ij),ivtx_old,ivtx_new,nvtx)
      end do

      inc = +1
      if (ivtx_old.gt.ivtx_new) inc = -1
      ! shift in topo
      ! rows
      ivhelp(1:nvtx) = topo(1:nvtx,ivtx_old)
      do idx = ivtx_old, ivtx_new-inc, inc
        topo(1:nvtx,idx) = topo(1:nvtx,idx+inc)
      end do
      topo(1:nvtx,ivtx_new) = ivhelp(1:nvtx)

      ! cols
      ivhelp(1:nvtx) = topo(ivtx_old,1:nvtx)
      do idx = ivtx_old, ivtx_new-inc, inc
        topo(idx,1:nvtx) = topo(idx+inc,1:nvtx)
      end do
      topo(ivtx_new,1:nvtx) = ivhelp(1:nvtx)

c dbg
c      print *,'new iord: ',iord(1:nvtx)
c      print *,'new ireo: ',ireo(1:nvtx)
c      print *,'new topomap:'
c      do idx = 1, nvtx
c        print '(10i9.8)',topo(idx,1:nvtx)
c      end do
c dbg

      return
      end subroutine

      subroutine shift_ivec(vec,idx1,idx2,len)

      implicit none

      integer, intent(in) ::
     &     len, idx1, idx2
      integer, intent(inout) ::
     &     vec(len)

      integer ::
     &     ihelp, idx, inc

      inc = +1
      if (idx1.gt.idx2) inc = -1

      ihelp = vec(idx1)
      do idx = idx1, idx2-inc, inc
        vec(idx) = vec(idx+inc)
      end do
      vec(idx2) = ihelp

      return
      end subroutine

      subroutine shift_i8vec(vec,idx1,idx2,len)

      implicit none

      integer, intent(in) ::
     &     len, idx1, idx2
      integer(8), intent(inout) ::
     &     vec(len)

      integer(8) ::
     &     ihelp
      integer ::
     &     idx,inc

      inc = +1
      if (idx1.gt.idx2) inc = -1

      ihelp = vec(idx1)
      do idx = idx1, idx2-inc, inc
        vec(idx) = vec(idx+inc)
      end do
      vec(idx2) = ihelp

      return
      end subroutine

      end
