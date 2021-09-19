*----------------------------------------------------------------------*
      subroutine pseudo_svmap2(svmap,contr,njoined_in)
*----------------------------------------------------------------------*
*     get pseudo svmap in which each vertex is assigned to some
*     result operator vertex except vertices without open lines (=>0).
*     If the assignment is not unique, we try to find an assignment
*     in which each result vertex is represented.
*     Needed for correct joining of complicated intermediates with
*     surrounding vertices
*
*     matthias, dec 2009
*     modified version: andreas, sept 2021
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(in), target ::
     &     contr
      integer, intent(in) ::
     &     njoined_in
      integer, intent(out) ::
     &     svmap(contr%nvtx)

      integer ::
     &     nvtx, ivtx, ij, idx1, idx2, icnt, nij, mxij, base
      integer ::
     &     occ_x(ngastp,2)
      integer, pointer ::
     &     svtx(:)
      logical ::
     &     is_ex(njoined_in), is_dx(njoined_in)
      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:)

      integer, external ::
     &     int8_expand
      
      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,'pseudo_svmap2')

      nvtx = contr%nvtx

      base = pack_base

      allocate(svtx(nvtx),vtx(nvtx),topo(nvtx,nvtx),
     &         xlines(nvtx,njoined_in))
      call pack_contr(svtx,vtx,topo,xlines,contr,njoined_in)
      if (ntest.ge.100) then
        write(lulog,*) 'no unique map! Approximate map using xlines:'
        call prt_contr_p(lulog,svtx,vtx,topo,xlines,nvtx,njoined_in)
      end if
      svmap = 0
      do ivtx = 1, nvtx
        nij = 0  ! collect number of contributions that this operator makes
        mxij = 0 ! the max ij that it contributes to
        do ij = 1, njoined_in
          if (xlines(ivtx,ij).gt.0) then
            nij = nij+1
            mxij = ij
            occ_x = 0
            icnt = int8_expand(xlines(ivtx,ij),base,occ_x)
            is_ex(ij) = occ_x(ihole,2).gt.0 .or.
     &              occ_x(ipart,1).gt.0 .or.
     &              occ_x(ivale,1).gt.0 .or.
     &              occ_x(iextr,1).gt.0
            is_dx(ij) = occ_x(ihole,1).gt.0 .or.
     &              occ_x(ipart,2).gt.0 .or.
     &              occ_x(ivale,2).gt.0 .or.
     &              occ_x(iextr,2).gt.0
          end if            
        end do
        if (nij.le.1) then
          svmap(ivtx) = mxij
        else if (nij.eq.2) then
          ! get the two positions
          idx1 = 0; idx2 = 0
          do ij = 1, njoined_in
            if (xlines(ivtx,ij).gt.0) then
              if (idx1==0) then
                idx1 = ij
              else
                idx2 = ij
              end if
            end if
          end do
          if (ntest.ge.100) then
            write(lulog,'(1x,"is_ex: ",10l2)') is_ex
            write(lulog,'(1x,"is_dx: ",10l2)') is_dx
            write(lulog,'(1x,"idx1,idx2: ",2i4)') idx1, idx2
          end if
          if (is_dx(idx1)) then
! if first of the two is linked down, we have to assign the op to this svertex of the result operator
            svmap(ivtx) = idx1
          else if (is_ex(idx2)) then
! if second is linked up, (and if the first one has no down-link) we assign the second one
            svmap(ivtx) = idx2
          else
            svmap(ivtx) = idx1  ! in all other cases: put it here (but no guarantee this will give a proper diagram)
          end if
        else
          call quit(1,'pseudo_svmap2','3 vertices: extend me')
        end if
      end do
      deallocate(svtx,vtx,topo,xlines)

      if (ntest.ge.100)
     &     write(lulog,'(x,a,10i5)') 'svmap: ',svmap(1:nvtx)

      return
      end
