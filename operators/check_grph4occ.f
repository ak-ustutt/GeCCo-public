*----------------------------------------------------------------------*
      logical function check_grph4occ(iocc,irst,
     &             str_info,ihpvgas,ngas,njoined)
*----------------------------------------------------------------------*
*     check whether all parts of the operator can be addressed by
*     the available graphs
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      integer, parameter ::
     &     ntest = 00

      type(strinf), intent(in) ::
     &     str_info
      integer, intent(in) ::
     &     ngas, ihpvgas(ngas), njoined,
     &     iocc(ngastp,2,njoined), irst(2,ngas,2,2,njoined)

      integer ::
     &     idx_gr(ngastp,2,njoined)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_func,'check_grph4occ')
        call wrt_occ_n(luout,iocc,njoined)
      end if

      call get_grph4occ(idx_gr,iocc,irst,
     &       str_info,ihpvgas,ngas,njoined,.false.)

      if (ntest.ge.100) then
        call wrt_occ_n(luout,idx_gr,njoined)        
      end if

      check_grph4occ = idx_gr(1,1,1).ge.0

      return
      end
