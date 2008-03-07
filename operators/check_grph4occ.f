*----------------------------------------------------------------------*
      logical function check_grph4occ(iocc,irst,njoined,
     &             str_info,orb_info)
*----------------------------------------------------------------------*
*     check whether all parts of the operator can be addressed by
*     the available graphs
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     njoined,
     &     iocc(ngastp,2,njoined),
     &     irst(2,orb_info%ngas,2,2,orb_info%nspin,njoined)

      integer ::
     &     ijoin,
     &     idx_gr(ngastp,2,njoined)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_func,'check_grph4occ')
        call wrt_occ_n(luout,iocc,njoined)
        do ijoin = 1, njoined
          call wrt_rstr(luout,irst(1,1,1,1,1,ijoin),orb_info%ngas)
        end do
      end if

      call get_grph4occ(idx_gr,iocc,irst,njoined,
     &       str_info,orb_info,.false.)

      if (ntest.ge.100) then
        call wrt_occ_n(luout,idx_gr,njoined)        
      end if

      check_grph4occ = idx_gr(1,1,1).ge.0

      return
      end
