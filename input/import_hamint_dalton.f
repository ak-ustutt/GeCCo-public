*----------------------------------------------------------------------*
      subroutine import_hamint_dalton(hop,ffham,str_info,orb_info)
*----------------------------------------------------------------------*
*     import one- and two-electron densities from DALTON
*     environment
*     we need:
*      SIRIFC for fock matrix
*      MOTWOINT for two-electron integrals
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_operator.h'

      type(operator), intent(in) ::
     &     hop
      type(filinf), intent(inout) ::
     &     ffham
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      ! read reference energy and fock matrix from SIRIFC
      ! and sort fock matrix into operator file
      call import_fock_dalton(ffham,hop,str_info,orb_info)

      ! get 2-electron integrals and sort them as well
      call import_h2_dalton(ffham,hop,str_info,orb_info)

      return
      end
