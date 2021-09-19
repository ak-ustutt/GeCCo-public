*----------------------------------------------------------------------*
      subroutine pseudo_svmap(svmap,contr,njoined_in)
*----------------------------------------------------------------------*
*     get pseudo svmap in which each vertex is assigned to some
*     result operator vertex except vertices without open lines (=>0).
*     If the assignment is not unique, we try to find an assignment
*     in which each result vertex is represented.
*     Needed for correct joining of complicated intermediates with
*     surrounding vertices
*
*     matthias, dec 2009
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
     &     nvtx, ivtx, ij
      integer, pointer ::
     &     svtx(:)
      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:)
      
      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,'pseudo_svmap')

      nvtx = contr%nvtx

      allocate(svtx(nvtx),vtx(nvtx),topo(nvtx,nvtx),
     &         xlines(nvtx,njoined_in))
      call pack_contr(svtx,vtx,topo,xlines,contr,njoined_in)
      if (ntest.ge.100) then
        write(lulog,*) 'no unique map! Approximate map using xlines:'
        call prt_contr_p(lulog,svtx,vtx,topo,xlines,nvtx,njoined_in)
      end if
      svmap = 0
      do ivtx = 1, nvtx
        do ij = 1, njoined_in
          if (xlines(ivtx,ij).gt.0) then
            svmap(ivtx) = ij
            if (ij.ge.ivtx) exit
          else if (ij.eq.ivtx.and.svmap(ivtx).eq.0) then
            svmap(ivtx) = ij ! account for zero vertices
          end if
        end do
      end do
      deallocate(svtx,vtx,topo,xlines)

      if (ntest.ge.100)
     &     write(lulog,'(x,a,10i5)') 'svmap: ',svmap(1:nvtx)

      return
      end
