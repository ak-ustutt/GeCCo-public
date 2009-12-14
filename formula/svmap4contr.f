*----------------------------------------------------------------------*
      subroutine svmap4contr(svmap,contr,occ_vtx_in,njoined_in,unique)
*----------------------------------------------------------------------*
*     assign each vertex in contr a number which tells us to which
*     super-vertex node this vertex contributes an external line
*     (0 if completely contracted).
*     occ_vtx is obtained from occvtx4contr with mode==0
*     (i.e. the first njoined entries describe the result occupation)
*     if njoined.eq.-1, interpret the result occupation as density
*     storage, i.e. divide into EX and DX part
*
*     TO BE WORKED OVER: does probably not work in the general case
*
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
     &     occ_vtx_in(ngastp,2,contr%nvtx+njoined_in), njoined_in
      integer, intent(out) ::
     &     svmap(contr%nvtx)
      logical, intent(out) ::
     &     unique

      integer ::
     &     nvtx, narc, ivtx, idx1, idx2, iarc, isuper, njoined, ij
      type(cntr_arc), pointer ::
     &     arc(:)
      integer, pointer ::
     &     occ_vtx(:,:,:), occ_res(:,:,:), occ_cnt(:,:),
     &     svtx(:)
      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:)
      
      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'svmap4contr')

      nvtx = contr%nvtx

      if (ntest.ge.100) then
        write(luout,*) 'occupations on input:'
        call wrt_occ_n(luout,occ_vtx_in,abs(njoined_in))
        call wrt_occ_n(luout,occ_vtx_in(1,1,abs(njoined_in)+1),nvtx)
      end if

      njoined = njoined_in

      if (njoined .gt. 0) then
      ! get a copy of the vertex occupations on occ_vtx and occ_res
        allocate(occ_vtx(ngastp,2,nvtx),occ_res(ngastp,2,njoined))
        occ_vtx = occ_vtx_in(1:ngastp,1:2,njoined+1:njoined+nvtx)
        occ_res = occ_vtx_in(1:ngastp,1:2,1:njoined)
      else
        njoined = 2
        allocate(occ_vtx(ngastp,2,nvtx),occ_res(ngastp,2,2))
        occ_vtx = occ_vtx_in(1:ngastp,1:2,2:nvtx+1)
        occ_res(1:ngastp,1:2,1) = iocc_xdn(2,occ_vtx_in(1:ngastp,1:2,1))
        occ_res(1:ngastp,1:2,2) = iocc_xdn(1,occ_vtx_in(1:ngastp,1:2,1))
        if (ntest.ge.100) then
          write(luout,*) 'converted result to:'
          call wrt_occ_n(luout,occ_res,2)
        end if
      end if

      narc = contr%narc
      arc => contr%arc
      ! remove all contracted indices
      do iarc = 1, narc
        idx1 = arc(iarc)%link(1)
        idx2 = arc(iarc)%link(2)
        occ_cnt => arc(iarc)%occ_cnt
        occ_vtx(1:ngastp,1:2,idx1) = occ_vtx(1:ngastp,1:2,idx1)
     &       - occ_cnt
        occ_vtx(1:ngastp,1:2,idx2) = occ_vtx(1:ngastp,1:2,idx2)
     &       - iocc_dagger(occ_cnt)
      end do

      if (ntest.ge.100) then
        write(luout,*) 'remainder occupations:'
        call wrt_occ_n(luout,occ_vtx,nvtx)
      end if

      ! loop over remainder an try to assign to components of super-vertex
      unique = .true.
      do ivtx = 1, nvtx
        if (iocc_zero(occ_vtx(1:ngastp,1:2,ivtx))) then
          svmap(ivtx) = 0
        else
          do isuper = 1, njoined+1
            ! a small trap:
            if (isuper.eq.njoined+1) then
              unique = .false. ! use information from xlines instead (below)
              exit
cmh              write(luout,*) 'current contraction:'
cmh              call prt_contr3(luout,contr,occ_vtx_in(1,1,1+njoined_in))
cmh              write(luout,*) 'result, vertices, vertices reduced:'
cmh              call wrt_occ_n(luout,occ_vtx_in,njoined)
cmh              call wrt_occ_n(luout,occ_vtx_in(1,1,njoined+1),nvtx)
cmh              call wrt_occ_n(luout,occ_vtx,nvtx)
cmh              call quit(1,'svmap4contr','something is strange')
            end if
            if (iocc_zero(occ_res(1:ngastp,1:2,isuper))) cycle
            if (iocc_bound('<=',occ_vtx(1:ngastp,1:2,ivtx),.false.,
     &           occ_res(1:ngastp,1:2,isuper),.false.))then
              svmap(ivtx) = isuper
c ??
              occ_res(1:ngastp,1:2,isuper) =
     &             occ_res(1:ngastp,1:2,isuper)-
     &             occ_vtx(1:ngastp,1:2,ivtx)
c ??
              exit
            end if
          end do

        end if

      end do

      ! quick fix: middle zero vertex would not be accounted for
      if (njoined_in.eq.3.and.svmap(3).eq.3.and.svmap(2).eq.0)
     &        unique = .false.

      ! if unique assignment fails, do approximate assignment using xlines
      if (.not.unique) then
        allocate(svtx(nvtx),vtx(nvtx),topo(nvtx,nvtx),
     &           xlines(nvtx,njoined_in))
        call pack_contr(svtx,vtx,topo,xlines,contr,njoined_in)
        if (ntest.ge.100) then
          write(luout,*) 'no unique map! Approximate map using xlines:'
          call prt_contr_p(luout,svtx,vtx,topo,xlines,nvtx,njoined_in)
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
      end if

      if (ntest.ge.100)
     &     write(luout,'(x,a,10i5)') 'svmap: ',svmap(1:nvtx)

      deallocate(occ_vtx,occ_res)

      return
      end
