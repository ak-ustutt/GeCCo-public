*----------------------------------------------------------------------*
      subroutine joinmap4contr(vtxmap,contr,nvtx_all,
     &                         iop_intm,ipos_intm,
     &                         svmap_intm,nvtx_intm,njoined_intm)
*----------------------------------------------------------------------*
*     make a map: 
*     for each vertex in expanded contraction --
*       positive number == vertex in original contraction contr
*       negative number == vertex of expanded intermediate
*     the intermediate is either identified by iop_intm (if >=0) or
*     the positions are given by ipos_intm(1:njoined)
*     svmap_intm: cf. svmap4contr()
*     nvtx_intm: number of vertices of expanded intermediate
*     vtxmap must have the size  contr%nvtx-njoined_intm+nvtx_intm
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     nvtx_intm, njoined_intm, nvtx_all,
     &     iop_intm, ipos_intm(njoined_intm), svmap_intm(nvtx_intm)
      integer, intent(out) ::
     &     vtxmap(nvtx_all)

      logical ::
     &     use_ipos, is_intm, active, error
      integer ::
     &     ivtx_total, ijoin, ivtx_intm, ivtx_st, ivtx, jvtx, ii

      ! use this information only of activated; must be positive numbers
      use_ipos = ipos_intm(1).ne.-1
      error = .false.

      ivtx_total = 0
      ijoin = 1
      ivtx_intm = 1
      do ivtx = 1, contr%nvtx
        ! check whether this is the intermediate to be replaced
        is_intm = contr%vertex(ivtx)%idx_op.eq.iop_intm
        ! only replace the active one (given by ipos_intm)
        if (use_ipos.and.is_intm) then
          active=.false.
          do ii = 1, njoined_intm
            active = active.or.ipos_intm(ii).eq.ivtx
          end do
          is_intm = is_intm.and.active
        end if

        if (.not.is_intm) then
          ! keep this vertex
          error = ivtx_total.eq.nvtx_all
          ivtx_total = min(nvtx_all,ivtx_total+1)
          vtxmap(ivtx_total) = ivtx
        else
          ! insert vertices of intermediate definition, which
          ! contribute to current primitive vertex of intermediate
          ivtx_st = ivtx_intm
          do jvtx = ivtx_st, nvtx_intm
            if (svmap_intm(jvtx).le.ijoin) then
              error = ivtx_total.eq.nvtx_all
              ivtx_total = min(nvtx_all,ivtx_total+1)
              vtxmap(ivtx_total) = -jvtx
              ivtx_intm = jvtx+1 ! next start would be from next vtx
            else
              ivtx_intm = jvtx ! start here next time
              ijoin = ijoin+1
              exit
            end if
          end do
        end if
      end do

      if (error) then
        write(lulog,*) 'svmap_intm =    ',svmap_intm
        if (use_ipos) write(lulog,*) 'ipos_vtx = ',ipos_intm
        write(lulog,*) 'vtxmap = ',vtxmap
        call quit(1,'joinmap4contr','more vertices than expected?')
      end if
      
      return
      end
