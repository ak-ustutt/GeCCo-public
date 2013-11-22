*----------------------------------------------------------------------*
      subroutine topo_rank_change(allowed,xlines,topo,
     &     occ,vtx_list,nj,nvtx,nj_res)
*----------------------------------------------------------------------*
*     check whether contraction is still valid after changing the rank
*     of the specified operator (modifies xlines, if allowed)
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      
      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nj, nvtx, nj_res
      logical, intent(out) ::
     &     allowed
      integer(8), intent(inout) ::
     &     xlines(nvtx,nj_res)
      integer(8), intent(in) ::
     &     topo(nvtx,nvtx), occ(nj)
      integer, intent(in) ::
     &     vtx_list(nj)

      integer ::
     &     ivtx, jvtx, ij, ij_res
      integer(8) ::
     &     ovl, occ_rem

      integer(8), external ::
     &     occ_overlap_p

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'here speaks topo_rank_change')
      end if

      if (nj_res.ne.1)
     &     call quit(1,'topo_rank_change','nj_res>1 not yet allowed')
      ! in this case we have to possibly care for multiple possibilities
      ! of xline-deletions (if the deleted index is distributed among several
      ! joined vertices); for the sake of caution, we exit here!

      ! loop over relevant vertices
      nj_loop: do ij = 1, nj
        ivtx = vtx_list(ij)

        ! check overlap with sum over topo fields: 
        ! must be equal to sum over topo, as we may not delete
        ! any contracted lines
        allowed = .true.
        occ_rem = 0
        do jvtx = 1, nvtx
          occ_rem = occ_rem + topo(ivtx,jvtx)
        end do
        ovl = occ_overlap_p(occ_rem,occ(ij))
        if (ntest.ge.100) then
          write(lulog,*) 'sum(topo),new_op,ovl:'
          write(lulog,'(3(x,i8.8))') occ_rem,occ(ij),ovl
        end if
        allowed = ovl.eq.occ_rem
        if (.not.allowed) exit

        occ_rem = occ(ij)-occ_rem

        if (ntest.ge.100) then
          write(lulog,*) 'check OK!! now the xlines:'
        end if

        ! get overlap of xlines (here we need to actually adapt for nj_res>1)
        do ij_res = 1, nj_res
          ovl = occ_overlap_p(xlines(ivtx,ij_res),occ_rem)
          if (ntest.ge.100) then
            write(lulog,*) 'xline,remainder,ovl'
            write(lulog,'(3(x,i8.8))') xlines(ivtx,ij_res),occ_rem,ovl
          end if
          xlines(ivtx,ij_res) = ovl
        end do

      end do nj_loop

      if (ntest.eq.100) then
        write(lulog,*) 'final xlines:'
        write(lulog,'(4(x,i8.8))') xlines(1:nvtx,1)
      end if

      return
      end
