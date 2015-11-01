*---------------------------------------------------------------------*
!>     wrapper for set_user_op ( very old)
!>
!>     set up user-provided operator
!>     irestr is chosen appropriately
!>     @param[inout] op the operator struct
!>     @param[in] name name of the operator
!>     @param[in] dagger should the transposed of the indicated current operator be created?
!>     @param[in] occ_def definition of the ca stings 
!>     @param[in] ndef length of occ_def (indirect?) 
!>     @param[in] orb_info information about the orbital spaces
*----------------------------------------------------------------------*
      subroutine set_uop(op,name,dagger,
     &     occ_def,ndef,orb_info)
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(operator), intent(inout) ::
     &     op
      character, intent(in) ::
     &     name*(*)
      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     occ_def, ndef

      type(orbinf) ::
     &     orb_info
      integer ::
     &     ncadiff,
     &     irestr(2,orb_info%ngas,2,2)

      ncadiff = 0
      call set_restr_for_uop()

      call set_user_op(op,name,optyp_operator,
     &     dagger,
     &     occ_def,ndef,irestr,orb_info)

      return
      
      contains

*----------------------------------------------------------------------*
      subroutine set_restr_for_uop()
*----------------------------------------------------------------------*

      implicit none

      integer ::
     &     ica, igastp, igas, max_rank

      integer, external ::
     &     ifndmax

      max_rank = ifndmax(occ_def,1,ngastp*2*ndef,1)
      irestr(1:2,1:orb_info%ngas,1:2,1:2) = 0
      do ica = 1, 2
        do igas = 1, orb_info%ngas
          irestr(1,igas,ica,1) = 0
          irestr(2,igas,ica,1) = max_rank
        end do
      end do

      return
      end subroutine set_restr_for_uop

      end
