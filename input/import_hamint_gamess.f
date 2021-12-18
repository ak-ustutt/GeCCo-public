*----------------------------------------------------------------------*
      subroutine import_hamint_gamess(hlist,str_info,orb_info)
*----------------------------------------------------------------------*
*     import one- and two-electron integrals from GAMESS
*     environment
*     we need:
*      MOINTS for fock matrix
*      MOINTS for two-electron integrals
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'par_gamess.h'
      
      type(me_list), intent(inout) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      ! read reference energy and fock matrix from MOINTS
      ! and sort fock matrix into operator file
      call import_fock_dalton(hlist,str_info,orb_info,.true.,moints)

      ! get 2-electron integrals and sort them as well
      call import_h2_gamess(hlist,str_info,orb_info)

      return
      end
