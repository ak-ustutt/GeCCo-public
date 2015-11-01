*----------------------------------------------------------------------*
!>     wrapper for set_user_op
!>
!>     set up user-provided operator 
!>     irestr is chosen appropriately
!>     @param[inout] op the operator struct
!>     @param[in] name name of the operator
!>     @param[in] occ_def
!>     @param[in] ndef 
!>     @param[in] njoined number if joined vertices
!>     @param[in] occ_def definition of the ca stings 
!>     @param[in] freeze ???
!>     @param[in] min_formal lower end of the formal blocks
!>     @param[in] orb_info information about the orbital spaces
*----------------------------------------------------------------------*
      subroutine set_uop2(op,name,
     &     occ_def,ndef,njoined,freeze,min_formal,orb_info)
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
      integer, intent(in) ::
     &     occ_def(ngastp,2,ndef*njoined), ndef, njoined, min_formal
      integer, intent(in) ::
     &     freeze(2,njoined)

      type(orbinf) ::
     &     orb_info
      integer ::
     &     ncadiff, 
     &     irestr(2,orb_info%ngas,2,2)

      integer ::
     &     ndef_loc
      integer, pointer ::
     &     occ_def_loc(:,:,:)

      ndef_loc = ndef
      allocate(occ_def_loc(ngastp,2,ndef*njoined))
      occ_def_loc = occ_def

      call remove_inactive_occ(occ_def_loc,ndef_loc,freeze,
     &                                             njoined,orb_info)

      ncadiff = 0
      call set_restr_for_uop()

      call set_user_op2(op,name,optyp_operator,
     &     occ_def_loc,ndef_loc,njoined,irestr,freeze,
     &                              min_formal,orb_info)

      deallocate(occ_def_loc)
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

      max_rank = ifndmax(occ_def,1,ngastp*2*ndef*njoined,1)
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
