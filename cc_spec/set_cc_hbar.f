*----------------------------------------------------------------------*
      subroutine set_cc_hbar(op,name,
     &     max_rank_t,orb_info)
*----------------------------------------------------------------------*
*     wrapper for set_genop2
*     set up blocks of Hbar which should be precomputed
*     using some empirical rules
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
      integer, intent(in) ::
     &     max_rank_t

      type(orbinf) ::
     &     orb_info

      integer ::
     &     ncadiff, min_xrank, max_xrank, min_rank, max_rank,
     &     hpvx_mnmx(2,ngastp),
     &     hpvxca_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2)

      ncadiff = 0
      min_rank = 1
      max_rank = max(2,max_rank_t-1)
      min_xrank = -max_rank
      max_xrank = max_rank-1
      hpvx_mnmx = 0
      hpvx_mnmx(2,IHOLE) = 2*max_rank_t
      hpvx_mnmx(2,IPART) = max_rank_t
      call set_hpvx_and_restr_for_hbar()

      call set_genop2(op,name,optyp_operator,
     &     min_rank,max_rank,ncadiff,
     &     min_xrank,max_xrank,
     &     hpvx_mnmx,hpvxca_mnmx,
     &     irestr,1,orb_info)

      return
      
      contains

*----------------------------------------------------------------------*
      subroutine set_hpvx_and_restr_for_hbar()
*----------------------------------------------------------------------*

      implicit none

      integer ::
     &     ica, igastp, igas

      do ica = 1, 2
        do igastp = 1, ngastp
          if (orb_info%nactt_hpv(igastp).gt.0)then
            hpvxca_mnmx(1,igastp,ica) = 0
            hpvxca_mnmx(2,igastp,ica) = max_rank
          else
            hpvxca_mnmx(1,igastp,ica) = 0
            hpvxca_mnmx(2,igastp,ica) = 0
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
      end subroutine set_hpvx_and_restr_for_hbar

      end
