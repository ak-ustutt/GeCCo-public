*----------------------------------------------------------------------*
      subroutine topo_approach_vtxs(ireo,sh_sign,
     &     svertex,vtx,topo,xlines,
     &     nvtx,nj,vtx_list,nlist)
*----------------------------------------------------------------------*
*     given a contraction in matrix form and a list of vertices:
*     try to move the vertices as close together as commutativity
*     with the intervening operators allows
*     this is the step to take prior to calling topo_merge_vtxs
*     we also calculate a possible sign change due to the reordering
*     of the operators
*----------------------------------------------------------------------*
      
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
     &     ireo(nvtx), sh_sign
      
      integer ::
     &     iord(nvtx)
      integer ::
     &     idx, jdx, ivtxr, ivtxrm1, ivtxrp1, ivtx, jvtx

      sh_sign = 1

      if (nlist.eq.1) then
        do idx = 1, nvtx
          ireo(idx) = idx
        end do
        return
      end if

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
c dbg
c      print *,'ireo (initial): ',ireo(1:nvtx)
c dbg
      ! loop over list
      ! a) first vertex
      ivtxr = ireo(vtx_list(1))
      ! position of second vertex
      ivtxrp1 = ireo(vtx_list(2))
      ! if vertices lie in between ...
c dbg
c        print *,'(1) ivtxr, ivtxrp1: ',ivtxr,ivtxrp1
c        call prt_contr_p(luout,svertex,vtx,topo,
c     &       xlines,nvtx,nj)
c dbg
      if (ivtxrp1-ivtxr.gt.1) then
        ! ... try to shift those up
        do ivtx = ivtxr+1, ivtxrp1-1
          jvtx = ivtx-1
          do while(jvtx.ge.1)
            if (.not.may_commute(jvtx,ivtx)) exit
            jvtx = jvtx-1
          end do
          call shift_vtx(ivtx,jvtx+1)
        end do
      end if
      ! update ivtxr position (ireo has changed)
      ivtxr = ireo(vtx_list(1))
c dbg
c        print *,'(2) ivtxr, ivtxrp1: ',ivtxr,ivtxrp1
c        call prt_contr_p(luout,svertex,vtx,topo,
c     &       xlines,nvtx,nj)
c dbg
      ! if still someone stands between us: push that vertex 
      ! below position of next vertex on list (if possible)
      if (ivtxrp1-ivtxr.gt.1) then
        do ivtx = ivtxrp1-1, ivtxr+1, -1
          jvtx = ivtx+1
          do while(jvtx.le.ivtxrp1)
            if (.not.may_commute(jvtx,ivtx)) exit
            jvtx = jvtx+1
          end do
          call shift_vtx(ivtx,jvtx-1)
        end do
      end if
c dbg
c        print *,'(3)'
c        call prt_contr_p(luout,svertex,vtx,topo,
c     &       xlines,nvtx,nj)
c dbg

      ! process further vertices on list
      do idx = 2, nlist-1
        ivtxr = ireo(vtx_list(idx))
        ivtxrp1 = ireo(vtx_list(idx+1))
        ivtxrm1 = ireo(vtx_list(idx-1))
c dbg
c        print *,'(4) ivtxr, ivtxrp1, ivtxrm1 ',ivtxr,ivtxrp1,ivtxrm1
c        call prt_contr_p(luout,svertex,vtx,topo,
c     &       xlines,nvtx,nj)
c dbg
        ! pushing up only, if upper vertex was separated anyway
        if (ivtxrp1-ivtxr.gt.1 .and. ivtxr-ivtxrm1.gt.1) then
          do ivtx = ivtxr+1, ivtxrp1-1
            jvtx = ivtx-1
            do while(jvtx.ge.1)
              if (.not.may_commute(jvtx,ivtx)) exit
              jvtx = jvtx-1
            end do
            call shift_vtx(ivtx,jvtx+1)
          end do
        end if
        ivtxr = ireo(vtx_list(idx)) ! update
c dbg
c        print *,'(5) ivtxr, ivtxrp1, ivtxrm1 ',ivtxr,ivtxrp1,ivtxrm1
c        call prt_contr_p(luout,svertex,vtx,topo,
c     &       xlines,nvtx,nj)
c dbg
        ! else we try pushing down ...
        if (ivtxrp1-ivtxr.gt.1) then
          do ivtx = ivtxrp1-1, ivtxr+1, -1
            jvtx = ivtx+1
            do while(jvtx.le.ivtxrp1)
              if (.not.may_commute(jvtx,ivtx)) exit
              jvtx = jvtx+1
            end do
            call shift_vtx(ivtx,jvtx-1)
          end do
        end if

      end do

c dbg
c      print *,'ireo (final): ',ireo(1:nvtx)
c      if (sh_sign.ne.1) print *,'sh_sign = ',sh_sign
c dbg
      return

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
     &     inc, idx, ij, ihelp, n_shift, n_pass
      integer(8) ::
     &     ivhelp(nvtx), i8occ

      integer, external ::
     &     nca_i8occ

      if (ivtx_old.eq.ivtx_new) return
c dbg
c      if (ivtx_old.ne.ivtx_new) print *,'shifting: ',ivtx_old,ivtx_new
c dbg
      i8occ = sum(xlines(ivtx_old,1:nj))+sum(topo(ivtx_old,1:nvtx))
      n_shift = nca_i8occ(i8occ)
c dbg
c      print *,'i8occ,n_shift: ',i8occ,n_shift
c dbg
      ! count vertex occupations that we pass while shifting
      n_pass  = 0
      do idx = ivtx_old+1, ivtx_new
        i8occ = sum(xlines(idx,1:nj))+sum(topo(idx,1:nvtx))
        n_pass = n_pass + nca_i8occ(i8occ)
      end do
c dbg
c      print *,'n_pass: ',i8occ,n_pass
c dbg
      if (mod(n_shift*n_pass,2).ne.0) sh_sign = -1*sh_sign

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
