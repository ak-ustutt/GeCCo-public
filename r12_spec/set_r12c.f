*----------------------------------------------------------------------*
      subroutine set_r12c(op,name,dagger,
     &     min_rank,max_rank,ncadiff,iformal,orb_info)
*----------------------------------------------------------------------*
*     wrapper for set_genop
*     set up coefficient vector for R12
*     hpvx_mnmx,irestr are chosen appropriately
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 000

      type(operator), intent(inout) ::
     &     op
      character, intent(in) ::
     &     name*(*)
      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     min_rank, max_rank, ncadiff,iformal

      type(orbinf) ::
     &     orb_info
      integer ::
     &     hpvx_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2)

      call set_hpvx_and_restr_for_c()

      call set_genop(op,name,optyp_operator,
     &     dagger,
     &     min_rank,max_rank,ncadiff,hpvx_mnmx,irestr,iformal,
     &     orb_info)

      return
      
      contains

*----------------------------------------------------------------------*
*     Subroutine to set the restrictions for the coefficient operator of
*     the R12 operators.
*----------------------------------------------------------------------*
      subroutine set_hpvx_and_restr_for_c()
      
      implicit none 

      integer ::
     &     ica, igas

      hpvx_mnmx(1:2,1:ngastp,1:2)=0
      hpvx_mnmx(1,ihole,1)=min_rank
      hpvx_mnmx(2,ihole,1)=min(2,max_rank)
      hpvx_mnmx(1,ihole,2)=min_rank
      hpvx_mnmx(2,ihole,2)=max_rank
      hpvx_mnmx(2,ipart,1)=max(0,max_rank-2)

      irestr(1:2,1:orb_info%ngas,1:2,1:2)=0
      do ica = 1, 2
        do igas = 1, orb_info%ngas
          if (orb_info%ihpvgas(igas,1).eq.ihole) then
            irestr(1,igas,ica,1) = 0
            if (ica.eq.1) then
              irestr(2,igas,ica,1) = min(2,max_rank)
            else
              irestr(2,igas,ica,1) = max_rank
            end if
          end if
          if (orb_info%ihpvgas(igas,1).eq.ipart) then
            irestr(1,igas,ica,1) = 0
            irestr(2,igas,ica,1) = max(0,max_rank-2)
          end if
        end do
      end do
      
      return
      end subroutine set_hpvx_and_restr_for_c

      end
