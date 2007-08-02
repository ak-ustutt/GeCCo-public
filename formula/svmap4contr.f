*----------------------------------------------------------------------*
      subroutine svmap4contr(svmap,contr,occ_vtx_in,njoined_in)
*----------------------------------------------------------------------*
*     assign each vertex in contr a number which tells us to which
*     super-vertex node this vertex contributes an external line
*     (0 if completely contracted).
*     occ_vtx is obtained from occvtx4contr with mode==0
*     (i.e. the first njoined entries describe the result occupation)
*     if njoined.eq.-1, interpret the result occupation as density
*     storage, i.e. divide into EX and DX part
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

      integer ::
     &     nvtx, narc, ivtx, idx1, idx2, iarc, isuper, njoined
      type(cntr_arc), pointer ::
     &     arc(:)
      integer, pointer ::
     &     occ_vtx(:,:,:), occ_res(:,:,:), occ_cnt(:,:)
      
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
      isuper = 1
      do ivtx = 1, nvtx
        if (iocc_zero(occ_vtx(1:ngastp,1:2,ivtx))) then
          svmap(ivtx) = 0
        else
          if (iocc_zero(occ_res(1:ngastp,1:2,isuper))) then
            isuper = isuper + 1
          end if
          svmap(ivtx) = isuper
          if (iocc_bound('<=',occ_vtx(1:ngastp,1:2,ivtx),.false.,
     &                        occ_res(1:ngastp,1:2,isuper),.false.))then
            occ_res(1:ngastp,1:2,isuper) =
     &           occ_res(1:ngastp,1:2,isuper)-occ_vtx(1:ngastp,1:2,ivtx)
          else
            write(luout,*) 'current contraction:'
            call prt_contr3(contr,occ_vtx)
            write(luout,*) 'result, vertices, vertices reduced:'
            call wrt_occ_n(luout,occ_vtx_in,njoined)
            call wrt_occ_n(luout,occ_vtx_in(1,1,njoined+1),nvtx)
            call wrt_occ_n(luout,occ_vtx,nvtx)
            call quit(1,'svmap4contr','something is strange')
          end if
        end if

      end do

      if (ntest.ge.100)
     &     write(luout,'(x,a,10i5)') 'svmap: ',svmap(1:nvtx)

      deallocate(occ_vtx,occ_res)

      return
      end
