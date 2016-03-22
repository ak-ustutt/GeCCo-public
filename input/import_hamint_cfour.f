*----------------------------------------------------------------------*
      subroutine import_hamint_cfour(hlist,str_info,orb_info)
*----------------------------------------------------------------------*
*     import one- and two-electron integrals from CFOUR
*     environment
*     we need:
*      fockin.dat for fock matrix
*      HF2 for two-electron integrals (or whatever name passed in orb_info)
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

      ! read reference energy and fock matrix from fockin.dat
      ! and sort fock matrix into operator file
      call import_fock_cfour(hlist,str_info,orb_info)

      ! get 2-electron integrals and sort them as well
      call import_h2_cfour(hlist,str_info,orb_info)

      return
      end
