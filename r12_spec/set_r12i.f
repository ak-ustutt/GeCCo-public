*----------------------------------------------------------------------*
      subroutine set_r12i(op,name,n_ap,
     &     min_rank,max_rank,ncadiff,iformal,orb_info)
*----------------------------------------------------------------------*
*     wrapper for set_genop
*     set up integrals for R12
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
      integer, intent(in) ::
     &     n_ap
      integer, intent(in) ::
     &     min_rank, max_rank, ncadiff,iformal

      type(orbinf) ::
     &     orb_info
      integer ::
     &     max_xrank, min_xrank, hpvx_mnmx(2,ngastp),
     &     hpvxca_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2)
      logical ::
     &     freeze(2)

      call set_hpvx_and_restr_for_int()

      min_xrank = -max_rank
      max_xrank = max_rank
      freeze(1) = .false.
      freeze(2) = .true.
      call set_genop2(op,name,optyp_operator,
     &     min_rank,max_rank,ncadiff,
     &     min_xrank,max_xrank,
     &     hpvx_mnmx,hpvxca_mnmx,
     &     irestr,iformal,freeze,
     &     orb_info)

      return
      
      contains

c-----------------------------------------------------------------------
      subroutine set_hpvx_and_restr_for_int()
c-----------------------------------------------------------------------
      implicit none

      integer::
     &     igastp,igas

      hpvxca_mnmx(1:2,1:ngastp,1:2)=0

      if (n_ap.eq.0) then
        hpvxca_mnmx(1:2,IHOLE,2)=(/0,2/)
c        hpvxca_mnmx(1:2,IHOLE,2)=2
      else
        hpvxca_mnmx(1:2,IHOLE,2)=(/0,2/)
        hpvxca_mnmx(1:2,IPART,2)=(/0,n_ap/)
      end if

      if (orb_info%nactt_hpv(IVALE).gt.0) then
        hpvxca_mnmx(1:2,IVALE,2)=(/0,2/)
      end if
      
      do igastp=1,ngastp
        if(igastp.eq.IEXTR)then
          hpvxca_mnmx(1,igastp,1)=0
          if (orb_info%caborb.gt.0) hpvxca_mnmx(2,igastp,1)=iformal-1
c dbg
          ! Quick fix.
          if(trim(name).eq.'G')then
            hpvxca_mnmx(1,igastp,2)=0
            if (orb_info%caborb.gt.0) hpvxca_mnmx(2,igastp,2)=iformal-1
          endif
c dbg
        else
          if(igastp.ne.IVALE.or.orb_info%nactt_hpv(IVALE).gt.0)then
            hpvxca_mnmx(1,igastp,1)=0
            hpvxca_mnmx(2,igastp,1)=max_rank
          endif
        endif  
      enddo

      hpvx_mnmx = 0
      hpvx_mnmx(2,IHOLE) = 4
      if (n_ap.eq.0) then
        hpvx_mnmx(2,IPART) = 2
      else
        hpvx_mnmx(2,IPART) = 2+n_ap
      end if
      hpvx_mnmx(2,IEXTR) = hpvxca_mnmx(2,IEXTR,1)
      if (orb_info%nactt_hpv(IVALE).gt.0) 
     &   hpvx_mnmx(2,IVALE) = 4

      irestr(1,1:orb_info%ngas,1:2,1:2)=0
      irestr(2,1:orb_info%ngas,1:2,1:2)=max_rank
      
c      do igas=1,orb_info%ngas
c        
c        irestr(1:2,ihole,2,1)=max_rank
c      
c        if(igastp.ne.ivale)then
c          irestr(1,igas,1,1)=0
c          irestr(2,igas,1,1)=max_rank
c        endif  
c      enddo

      return
      end subroutine set_hpvx_and_restr_for_int

      end
