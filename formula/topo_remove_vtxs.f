*----------------------------------------------------------------------*
      subroutine topo_remove_vtxs(     
     &           svertex_new,vtx_new,topo_new,xlines_new,
     &           svertex,    vtx,    topo,    xlines,
     &           vtx_list,nlist,
     &           nvtx,nvtx_new,nj,nj_new)
*----------------------------------------------------------------------*
*     remove a list of vertices and generate new vtx, topo, xlines ...
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 000

      integer, intent(in) ::
     &     nlist,nvtx,nvtx_new,nj,nj_new,
     &     vtx_list(nlist), svertex(nvtx)
      integer(8), intent(in) ::
     &     vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,nj)
      integer(8), intent(out) ::
     &     vtx_new(nvtx_new), topo_new(nvtx_new,nvtx_new),
     &     xlines_new(nvtx_new,nj_new)
      integer, intent(out) ::
     &     svertex_new(nvtx_new)

      integer ::
     &     dim_error, iidx, jjdx, idx, jdx, ilist, nvtx1

      integer(8) ::
     &     xlines_scr_u(nvtx_new,nlist),
     &     xlines_scr_l(nvtx_new,nlist)

      integer, external ::
     &     imltlist
      logical, external ::
     &     zero_i8vec

      if (nvtx_new.eq.0) return

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'topo_remove_vtxs at work')
        write(lulog,*) 'on entry:'
        call prt_contr_p(lulog,svertex,vtx,topo,xlines,nvtx,nj)        
        write(lulog,*) 'vtx_list: ',vtx_list(1:nlist)
      end if

      xlines_scr_u = 0
      xlines_scr_l = 0

      xlines_new = 0

      ! store on new arrays, with the vertices on the list
      ! removed
      dim_error = 0
      iidx = 0
      the_loop: do idx = 1, nvtx
        if (imltlist(idx,vtx_list,nlist,1).gt.0) cycle the_loop
        iidx = iidx+1
        if (iidx.gt.nvtx_new) then
          dim_error = 1
          exit the_loop
        end if
        svertex_new(iidx) = svertex(idx)
        vtx_new(iidx) = vtx(idx)
        xlines_new(iidx,1:nj) = xlines(idx,1:nj)
        ! collect dangling connections
        do ilist = 1, nlist
          if (vtx_list(ilist).lt.idx) then
            ! connections to lower vertices
            xlines_scr_l(iidx,ilist) = topo(idx,vtx_list(ilist))
          else
            ! connections to upper vertices
            xlines_scr_u(iidx,ilist) = topo(idx,vtx_list(ilist))
          end if
        end do
        jjdx = 0
        do jdx = 1, nvtx
          if (imltlist(jdx,vtx_list,nlist,1).gt.0) cycle
          jjdx = jjdx+1
          if (jjdx.gt.nvtx_new) then
            dim_error = 2
            exit the_loop
          end if
          topo_new(jjdx,iidx) = topo(jdx,idx)
        end do
      end do the_loop

      if (dim_error.gt.0) then
        write(lulog,*) 'position = ',dim_error
        call quit(1,'topo_remove_vtxs','dimension error (see above)')
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'xlines_scr_u:'
        do idx = 1, nvtx_new
          write(lulog,'(1x,8i8.8)') xlines_scr_u(idx,1:nlist)
        end do
        write(lulog,*) 'xlines_scr_l:'
        do idx = 1, nvtx_new
          write(lulog,'(1x,8i8.8)') xlines_scr_l(idx,1:nlist)
        end do
      end if
  
      ! what to do with new xlines?
      ! two simple cases
      if (zero_i8vec(xlines,nvtx*nj,1).and.nj.eq.1) then
        if (nlist.eq.1.and.nj_new.eq.2) then
          xlines_new(1:nvtx_new,1) = xlines_scr_u(1:nvtx_new,1)
          xlines_new(1:nvtx_new,2) = xlines_scr_l(1:nvtx_new,1)
        else if (zero_i8vec(xlines_scr_u,nvtx_new*nlist,1)) then
          if (nj_new.eq.nlist) then
            do idx = 1, nj_new
              ! reverse sequence as result should be conjugate
              ! of derivative operator
              jdx = nj_new+1-idx
              xlines_new(1:nvtx_new,idx) =
     &             xlines_scr_l(1:nvtx_new,jdx)
            end do
          else
            dim_error = 3
          end if
        else if (zero_i8vec(xlines_scr_l,nvtx_new*nlist,1)) then
          if (nj_new.eq.nlist) then
            do idx = 1, nj_new
              ! cf. above
              jdx = nj_new+1-idx
              xlines_new(1:nvtx_new,idx) =
     &             xlines_scr_u(1:nvtx_new,jdx)
            end do
          else
            dim_error = 4
          end if
        else
          call quit(1,'topo_remove_vtxs',
     &         'new case occurred (for nj=1)')
        end if
      ! allow 2nd derivations
      ! if operand has open lines, they will be deleted
      ! => reduces dimension of operator!
      ! if you want to keep the open lines, insert a unit operator
      ! prior to the first differentiation!
      else if (nj.eq.2.and.nlist.eq.1.and.nj_new.eq.3) then

        ! choose lowest number of vertices belonging to first supervertex
        ! somewhat arbitrary, could be more sophisticated
        do nvtx1 = nvtx,1,-1
          if (xlines(nvtx1,1).ne.0) exit
        end do
        if (.not.all(xlines(1:nvtx1,2).eq.0).or.
     &      .not.all(xlines(nvtx1+1:nvtx,1).eq.0))
     &    call quit(1,'topo_remove_vtxs',
     &       'new case occurred (complicated supervertex structure)')
c        if (nvtx1.gt.vtx_list(1)) then
c          do idx = 3,2,-1
c            xlines_new(1:nvtx_new,idx) = xlines_new(1:nvtx_new,idx-1)
c          end do
c          xlines_new(1:nvtx_new,1) = 0
c          xlines_new(1:nvtx_new,1) = xlines_scr_u(1:nvtx_new,1)
c          xlines_new(1:nvtx_new,2) = xlines_new(1:nvtx_new,2)
c     &                             + xlines_scr_l(1:nvtx_new,1)
c        else if (nvtx1+1.lt.vtx_list(1).or.nvtx1.eq.0) then
c          if (all(xlines(nvtx1+1:vtx_list(1)-1,1:nj).eq.0).and.
c     &        nvtx1.gt.0) call quit(1,'topo_remove_vtxs',
c     &       'new case occurred (multiple possibilities(2))')
c          xlines_new(1:nvtx_new,2) = xlines_new(1:nvtx_new,2)
c     &                             + xlines_scr_u(1:nvtx_new,1)
c          xlines_new(1:nvtx_new,3) = xlines_scr_l(1:nvtx_new,1)
c        else
c          call quit(1,'topo_remove_vtxs',
c     &       'new case occurred (multiple possibilities(1))')
          !preliminary: assume that derivation wrt lower vtx is done last
          xlines_new(1:nvtx_new,2) = xlines_new(1:nvtx_new,2)
     &                             + xlines_scr_u(1:nvtx_new,1)
          xlines_new(1:nvtx_new,3) = xlines_scr_l(1:nvtx_new,1)
c        end if

      else if (nj.eq.2.and.nlist.eq.1.and.nj_new.eq.2
     &         .and.all(xlines_scr_u(1:nvtx_new,1).eq.0)) then
        ! 2nd deriv. where operand is directly below the gap
        ! -> no middle vertex required
        xlines_new(1:nvtx_new,2) = xlines_new(1:nvtx_new,2)
     &                           + xlines_scr_l(1:nvtx_new,1)

      else if (nj.eq.2.and.nlist.eq.1.and.nj_new.eq.2) then
        ! here it is unclear which result vtx the upper lines belong to
        ! put upper lines on upper result vertex if possible (no reo.)
        if (vtx_list(1).eq.1.or.all(xlines_new(1:vtx_list(1)-1,2).eq.0))
     &     then
          xlines_new(1:nvtx_new,1) = xlines_new(1:nvtx_new,1)
     &                             + xlines_scr_u(1:nvtx_new,1)
        else ! has to be put to lower vertex
          xlines_new(1:nvtx_new,2) = xlines_new(1:nvtx_new,2)
     &                             + xlines_scr_u(1:nvtx_new,1)

        end if
        ! put lower lines on lower vertex in any case
        xlines_new(1:nvtx_new,2) = xlines_new(1:nvtx_new,2)
     &                           + xlines_scr_l(1:nvtx_new,1)

      else
c dbg
c        ! e.g. for effective Hamiltonian through double differentiation
c        xlines_new = xlines_new + xlines_scr_u
c dbgend
          call quit(1,'topo_remove_vtxs',
     &       'new case occurred')
      end if

      if (dim_error.gt.0) then
        write(lulog,*) 'position = ',dim_error
        call quit(1,'topo_remove_vtxs','dimension error (see above)')
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'finally:'
        call prt_contr_p(lulog,svertex_new,vtx_new,topo_new,
     &       xlines_new,nvtx_new,nj_new)        
      end if

      return
      end
