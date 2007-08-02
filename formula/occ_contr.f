*----------------------------------------------------------------------*
      subroutine occ_contr(occ,contr,occ_vtx_in,njoined)
*----------------------------------------------------------------------*
*     get resulting occupation of contr
*     occ_vtx is obtained from occvtx4contr with mode==1
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
     &     occ_vtx_in(ngastp,2,contr%nvtx), njoined
      integer, intent(out) ::
     &     occ(ngastp,2,njoined)

      integer ::
     &     nvtx, narc, idx1, idx2, iarc, type, type_last
      type(cntr_arc), pointer ::
     &     arc(:)
      integer, pointer ::
     &     occ_vtx(:,:,:), occ_cnt(:,:)
      
      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'occ_contr')

      nvtx = contr%nvtx

      if (ntest.ge.100) then
        write(luout,*) 'occupations on input:'
        call wrt_occ_n(luout,occ_vtx_in,nvtx)
      end if

      ! get a copy of occ_vtx
      allocate(occ_vtx(ngastp,2,nvtx))
      occ_vtx = occ_vtx_in

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

      ! loop over remainder-occupations
      type_last = -1
      idx2 = 0
      do idx1 = 1, nvtx        
        ! get type of remainder: 0, ex, dx, general
        if (iocc_zero(occ_vtx(1:ngastp,1:2,idx1))) then
          type = 0
        else if (iocc_zero(iocc_xdn(2,occ_vtx(1:ngastp,1:2,idx1)))) then
          type = 1
        else if (iocc_zero(iocc_xdn(1,occ_vtx(1:ngastp,1:2,idx1)))) then
          type = 2
        else
          type = 3
        end if
        !  0: remove -- do nothing
        !  ex/dx: if previous vertex was of same type, add up
        if (type.eq.1.or.type.eq.2) then
          if (type.eq.type_last) then
            occ_vtx(1:ngastp,1:2,idx2) = occ_vtx(1:ngastp,1:2,idx2) +
     &                                   occ_vtx(1:ngastp,1:2,idx1)
          else 
            idx2 = idx2+1
            if (idx2.lt.idx1)
     &           occ_vtx(1:ngastp,1:2,idx2) = occ_vtx(1:ngastp,1:2,idx1)
          end if
          type_last = type
        else if (type.eq.3) then
          !  general: remains as is
          idx2 = idx2+1
          if (idx2.lt.idx1)
     &         occ_vtx(1:ngastp,1:2,idx2) = occ_vtx(1:ngastp,1:2,idx1)
          type_last = type
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'after joining vertices (remaining:',idx2,')'
        call wrt_occ_n(luout,occ_vtx,idx2)
      end if

      ! exception: density storage of (pure DX|pure EX)
      if (idx2.eq.2.and.njoined.eq.1) then
        if (iocc_nonzero(iocc_xdn(1,occ_vtx(1:ngastp,1:2,1))).or.
     &      iocc_nonzero(iocc_xdn(2,occ_vtx(1:ngastp,1:2,2)))) then
          call wrt_occ_n(luout,occ_vtx,2)
          call quit(1,'occ_contr',
     &         'cannot express supervertex as density')
        end if
        occ(1:ngastp,1:2,1) = occ_vtx(1:ngastp,1:2,1) +
     &                    occ_vtx(1:ngastp,1:2,2)
      else if (idx2.ne.njoined) then
        write(luout,*) 'expected super-vertices: ',njoined
        write(luout,*) 'but found the following final occupation:'
        call wrt_occ_n(luout,occ_vtx,idx2)
        call quit(1,'occ_contr',
     &         'resulting supervertex is incompatible with expected
     &       number of joined vertices')
      else
        occ(1:ngastp,1:2,1:njoined) = occ_vtx(1:ngastp,1:2,1:njoined)
      end if

      deallocate(occ_vtx)

      return
      end
