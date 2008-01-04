*----------------------------------------------------------------------*
      subroutine occ_contr(occ,ierr,contr,occ_vtx_in,njoined)
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
     &     ierr,
     &     occ(ngastp,2,njoined)

      integer ::
     &     nvtx, narc, nins, idx1, idx2, iarc, type, type_last
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

      ierr = 0

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
      else if (idx2.gt.njoined) then
        ierr = -1-idx2
c        write(luout,*) 'expected super-vertices: ',njoined
c        write(luout,*) 'but found the following final occupation:'
c        call wrt_occ_n(luout,occ_vtx,idx2)
c        call quit(1,'occ_contr',
c     &         'resulting supervertex is incompatible with expected
c     &       number of joined vertices')
      else if (idx2.lt.njoined) then
        ! try to fix: insert 0-occ after DX, before EX
c dbg
c        print *,'fixing 0-occ'
c dbg
        nvtx = idx2
        nins = njoined-nvtx
        idx1 = 1
        do idx2 = 1, nvtx
          ! no de-excitation follows?
          if (iocc_zero(iocc_xdn(2,occ_vtx(1:ngastp,1:2,idx2)))) then
            occ(1:ngastp,1:2,idx1:idx1-1+nins) = 0
            occ(1:ngastp,1:2,idx1+nins) = occ_vtx(1:ngastp,1:2,idx2)
            idx1 = idx1+nins+1
          else
            occ(1:ngastp,1:2,idx1) = occ_vtx(1:ngastp,1:2,idx2)
            idx1 = idx1+1
          end if
        end do
        if (idx1.le.njoined) then
          occ(1:ngastp,1:2,idx1:njoined) = 0
        end if
      else
        occ(1:ngastp,1:2,1:njoined) = occ_vtx(1:ngastp,1:2,1:njoined)
      end if

      deallocate(occ_vtx)

      return
      end
