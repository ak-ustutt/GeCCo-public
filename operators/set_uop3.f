*----------------------------------------------------------------------*
!>     wrapper for set_user_op2
!>
!>     set up user-provided operator 
!>     irestr is chosen appropriately
!>     @param[inout] op the operator struct
!>     @param[in] name name of the operator
!>     @param[in] occ_descr descriptor string for op (see process_occ_descr)
!>     @param[in] njoined number if joined vertices
!>     @param[in] freeze ???
!>     @param[in] min_formal lower end of the formal blocks
!>     @param[in] orb_info information about the orbital spaces
*----------------------------------------------------------------------*
      subroutine set_uop3(op,name,
     &     occ_descr,njoined,freeze,min_formal,orb_info)
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(operator), intent(inout) ::
     &     op
      character(len=*), intent(in) ::
     &     name, occ_descr
      integer, intent(in) ::
     &     njoined, min_formal
c      logical, intent(in) ::
c     &     freeze(2)
      integer, intent(in) ::
     &     freeze(2,njoined)

      type(orbinf) ::
     &     orb_info
      integer ::
     &     irestr(2,orb_info%ngas,2,2)
      integer ::
     &     maxlist,ndef
   
      integer, parameter ::
     &     maxdef = 512
      integer ::
     &     occ_def(ngastp,2,maxdef)


      ! process descriptor
      maxlist = maxdef/njoined
      call process_occ_descr(occ_def,ndef,
     &                       occ_descr,njoined,maxlist)

      call remove_inactive_occ(occ_def,ndef,freeze,njoined,orb_info)

      call set_restr_for_uop()

      call set_user_op2(op,name,optyp_operator,
     &     occ_def,ndef,njoined,irestr,freeze,min_formal,orb_info)

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
