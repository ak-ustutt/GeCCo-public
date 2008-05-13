*----------------------------------------------------------------------*
      subroutine set_r12gem(op,name,n_ap,
     &     min_rank,max_rank,ansatz,orb_info)
*----------------------------------------------------------------------*
*     wrapper for set_genop
*     set up r12 geminal factor for the given ansatz
*     min_rank might be 1 or 2
*     max_rank>2 used for generation of formal C*R12 operator
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
     &     min_rank, max_rank, ansatz, n_ap

      type(orbinf), intent(in) ::
     &     orb_info
      integer ::
     &     iformal, ncadiff,
     &     min_x_rank, max_x_rank,
     &     min_p_rank, max_p_rank,
     &     min_h_rank, max_h_rank,
     &     hpvx_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2)

c dbg
      print *,'ansatz = ',ansatz
c dbg
      if(ansatz.eq.1)then
        min_x_rank=2
        max_x_rank=2
        min_p_rank=0
        max_p_rank=max_rank-2
      elseif(ansatz.eq.2)then
        min_x_rank=0
        max_x_rank=2
        min_p_rank=0
        max_p_rank=max_rank
      elseif(ansatz.eq.3)then
        min_x_rank=1
        max_x_rank=2
        min_p_rank=0
        max_p_rank=max_rank-1
      endif  
      min_h_rank=min_rank
      max_h_rank=max_rank

      iformal = 0
      ncadiff = 0

      call set_hpvx_and_restr_for_r()

      call set_genop(op,name,optyp_operator,
     &     .false.,
     &     min_rank,max_rank,ncadiff,hpvx_mnmx,irestr,iformal,
     &     orb_info)

      return
      
      contains

c-----------------------------------------------------------------------
      subroutine set_hpvx_and_restr_for_r()
c-----------------------------------------------------------------------
      implicit none

      integer::
     &     ica,igastp,igas

      ! Constraints on the operator are made depending on which ansatz
      ! is being used.
      hpvx_mnmx(1:2,1:ngastp,1:2)=0
      
      hpvx_mnmx(1:2,IPART,1) = (/min_p_rank,max_p_rank/)
      hpvx_mnmx(1:2,IEXTR,1) = (/min_x_rank,max_x_rank/)
      
      hpvx_mnmx(1:2,IHOLE,2) = (/0,2/)
      hpvx_mnmx(1:2,IPART,2) = (/0,n_ap/)

c      do ica=1,2
c        do igastp=1,ngastp
c          if(orb_info%nactt_hpv(igastp).gt.0.or.igastp.eq.iextr)then
c            if(ica.eq.2.and.igastp.eq.IHOLE)then
c              hpvx_mnmx(1,igastp,ica)=min_h_rank
c              hpvx_mnmx(2,igastp,ica)=max_h_rank
c            elseif(ica.eq.1)then
c              if(igastp.eq.IPART)then
c                hpvx_mnmx(1,igastp,ica)=min_p_rank
c                hpvx_mnmx(2,igastp,ica)=max_p_rank
c              elseif(igastp.eq.iextr)then
c                hpvx_mnmx(1,igastp,ica)=min_x_rank
c                hpvx_mnmx(2,igastp,ica)=max_x_rank 
c              endif
c            endif  
c          else
c            hpvx_mnmx(1,igastp,ica)=0
c            hpvx_mnmx(2,igastp,ica)=0
c          endif
c        enddo
c      enddo  

      irestr(1:2,1:orb_info%ngas,1:2,1:2)=0
      do ica=1,2
        do igas=1,orb_info%ngas
          irestr(1,igas,ica,1)=0
          irestr(2,igas,ica,1)=max_rank
        enddo
      enddo  
  
      return
      end subroutine set_hpvx_and_restr_for_r

      end
