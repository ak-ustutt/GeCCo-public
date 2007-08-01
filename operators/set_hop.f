*----------------------------------------------------------------------*
      subroutine set_hop(op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,orb_info)
*----------------------------------------------------------------------*
*     wrapper for set_genop
*     set up hamiltonian-like operator (minrank to maxrank)
*     hpvx_mnmx,irestr are chosen appropriately
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 100

      type(operator), intent(inout) ::
     &     op
      character, intent(in) ::
     &     name*(*)
      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     absym, casym, gamma, s2, ms,
     &     min_rank, max_rank

      type(orbinf) ::
     &     orb_info
      integer ::
     &     ncadiff,
     &     hpvx_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2)

      ncadiff = 0
      call set_hpvx_and_restr_for_hop()

      call set_genop(op,name,optyp_operator,
     &     dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,hpvx_mnmx,irestr,orb_info)

      return
      
      contains

*----------------------------------------------------------------------*
      subroutine set_hpvx_and_restr_for_hop()
*----------------------------------------------------------------------*

      implicit none

      integer ::
     &     ica, igastp, igas

      do ica = 1, 2
        do igastp = 1, ngastp
          if (orb_info%nactt_hpv(igastp).gt.0) then
            hpvx_mnmx(1,igastp,ica) = 0
            hpvx_mnmx(2,igastp,ica) = max_rank
          else
            hpvx_mnmx(1,igastp,ica) = 0
            hpvx_mnmx(2,igastp,ica) = 0
          end if
        end do
      end do
      irestr(1:2,1:orb_info%ngas,1:2,1:2) = 0
      do ica = 1, 2
        do igas = 1, orb_info%ngas
          irestr(1,igas,ica,1) = 0
          irestr(2,igas,ica,1) = max_rank
        end do
      end do

      return
      end subroutine set_hpvx_and_restr_for_hop

      end
