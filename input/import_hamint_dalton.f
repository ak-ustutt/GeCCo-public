*----------------------------------------------------------------------*
      subroutine import_hamint_dalton(hlist,str_info,orb_info)
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
      include 'def_me_list.h'
      
      type(me_list), intent(inout) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      ! read reference energy and fock matrix from SIRIFC
      ! and sort fock matrix into operator file
      call import_fock_dalton(hlist,str_info,orb_info,.false.,'-')

      ! get 2-electron integrals and sort them as well
      call import_h2_dalton(hlist,str_info,orb_info)

      return
      end
