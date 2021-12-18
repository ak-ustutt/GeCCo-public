*----------------------------------------------------------------------*
      subroutine occ_contr2(occ,ierr,contr,njoined)
*----------------------------------------------------------------------*
*     get resulting occupation of contr
*     version with contr containing xarc info -> easy game
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
     &     njoined
      integer, intent(out) ::
     &     ierr,
     &     occ(ngastp,2,njoined)

      integer ::
     &     nxarc, ixarc, ivtx
      type(cntr_arc), pointer ::
     &     xarc(:)
      
      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,'occ_contr2')

      nxarc = contr%nxarc
      xarc  => contr%xarc

      ierr = 0
      occ = 0
      do ixarc = 1, nxarc
        ivtx = xarc(ixarc)%link(2)
        if (ivtx.le.0.or.ivtx.gt.njoined) then
          ierr = min(ierr,-1-ivtx)
          cycle
        end if
        occ(1:ngastp,1:2,ivtx) = occ(1:ngastp,1:2,ivtx)
     &                         + xarc(ixarc)%occ_cnt
      end do

      if (ntest.ge.100) then
        if (ierr.eq.0) then
          write(lulog,*) 'resulting occ:'
          call wrt_occ_n(lulog,occ,njoined)
        else
          write(lulog,*) 'error code: ',ierr
        end if
      end if

      return
      end
