*----------------------------------------------------------------------*
      subroutine set_ps_op(oper,name,iocc,irst,njoined,orb_info)
*----------------------------------------------------------------------*
*     set up a pseudo-operator (only single block) 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_orbinf.h'

      type(operator), intent(out) ::
     &     oper
      character, intent(in) ::
     &     name*(*)
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     njoined,
     &     iocc(ngastp,2,njoined), irst(2,orb_info%ngas,2,2,njoined)

      integer ::
     &     ica, ihpv, idx, ijoin, ngas

      integer, external ::
     &     ielsum

      ngas = orb_info%ngas
      
      oper%dagger = .false.

      oper%name = name

      oper%n_occ_cls = 1
      oper%njoined = njoined

      ! initial allocation
      call init_operator(oper,orb_info)

      oper%ihpvca_occ(1:ngastp,1:2,1:njoined) =
     &           iocc(1:ngastp,1:2,1:njoined)
      oper%formal = .false.
      oper%formal_blk(1) = .false.

      oper%ica_occ(1:2,1) = 0
      do ijoin = 1, njoined
        oper%ica_occ(1,1) = oper%ica_occ(1,1)+
     &       ielsum(iocc(1,1,ijoin),ngastp)
        oper%ica_occ(2,1) = oper%ica_occ(2,1)+
     &       ielsum(iocc(1,2,ijoin),ngastp)
      end do

      oper%igasca_restr(1:2,1:ngas,1:2,1:2,1:njoined) =
     &             irst(1:2,1:ngas,1:2,1:2,1:njoined)

      return
      end
