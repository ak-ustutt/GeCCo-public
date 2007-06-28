*----------------------------------------------------------------------*
      subroutine set_r12_operators(op_list,nops,orb_info)
*----------------------------------------------------------------------*
*     hard-wired set-up of the basic operators needed for the
*     r12 excitation operators
*     Modified version of set_cc_operators. (GWR March 2007)
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'ifc_input.h'
      include 'par_opnames_gen.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'explicit.h'

      type(operator_list), intent(inout), target ::
     &     op_list
      integer, intent(inout) ::
     &     nops
      type(orbinf) ::
     &     orb_info
      type(operator_list), pointer ::
     &     list_pnt
      logical ::
     &     dagger
      character ::
     &     name*(len_opname)
      integer ::
     &     absym, casym, s2, ms, min_rank, max_rank, ncadiff,
     &     gamma, iarr(1),  min_h_rank, max_h_rank,
     &     min_p_rank, max_p_rank, min_x_rank, max_x_rank, iformal
      integer ::
     &     ihpv_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2)

c      write(luout,'(/"R12 operator definition subroutine.")')

      list_pnt=>op_list
      do while(associated(list_pnt%next))
        list_pnt=>list_pnt%next
      enddo
      ! The above loop gets us to the end of the list. 
      if (associated(list_pnt%op))then
        allocate(list_pnt%next)
        list_pnt%next%prev=>list_pnt
        list_pnt=>list_pnt%next
        nullify(list_pnt%next)
      endif
      allocate(list_pnt%op)  

      call get_argument_value('method.R12','ansatz',ival=ansatze)
      if(ansatze.eq.2.or.ansatze.lt.1)then
        write(luout,'("Error: Undefined R12 ansatz requested.")')
        stop
      endif

      call get_argument_value('method.R12','triples',ival=trir12)

      if(ansatze.eq.1)then
        min_x_rank=2
        max_x_rank=2
        min_p_rank=0
        max_p_rank=0
      elseif(ansatze.eq.2)then
        min_x_rank=0
        max_x_rank=2
        min_p_rank=0
        max_p_rank=2
      elseif(ansatze.eq.3)then
        min_x_rank=1
        max_x_rank=2
        min_p_rank=0
        max_p_rank=1
      endif  
      min_h_rank=2
      max_h_rank=2

      nops=nops+1
      ! New entry: linear R12 operator.
      name=op_r12
      dagger=.false.
      absym=0
      casym=0
      gamma=1
      s2=0
      ms=0
      min_rank=2
      max_rank=2
      ncadiff=0
      iformal=max_x_rank
  
      call set_hpvx_and_restr_for_r()

      call set_genop(list_pnt%op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info%iad_gas,orb_info%ihpvgas,orb_info%ngas)  
      
      ! New entry: variable coefficient operator associated with R12.
      nops=nops+1
      allocate(list_pnt%next)
      list_pnt%next%prev=>list_pnt
      list_pnt=>list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%op)
      
      name = op_c12
      dagger=.false.
      absym=0
      casym=0
      gamma=1
      s2=0
      ms=0
      ncadiff=0
      min_rank=2
      if(trir12.eq.1)then
        max_rank=3
      else  
        max_rank=2
      endif
      iformal=max_rank+1
      call set_hpvx_and_restr_for_c()
      
      call set_genop(list_pnt%op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info%iad_gas,orb_info%ihpvgas,orb_info%ngas)

      ! New entry: adjoint of linear R12 operator.
      nops=nops+1
      allocate(list_pnt%next)
      list_pnt%next%prev=>list_pnt
      list_pnt=>list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%op)

      name = op_rba
      ! Define an excitation operator to ensure the same storage
      ! sequence as for R, however, consider its adjoint.
      dagger=.true.
      absym=0
      casym=0
      gamma=1
      s2=0
      ms=0
      ncadiff=0
      min_rank=2
      max_rank=2
      iformal=max_x_rank

      call set_hpvx_and_restr_for_r()

      call set_genop(list_pnt%op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info%iad_gas,orb_info%ihpvgas,orb_info%ngas)

      ! New entry: variable coefficient operator associated with R12bar.
      nops=nops+1
      allocate(list_pnt%next)
      list_pnt%next%prev=>list_pnt
      list_pnt=>list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%op)
      
      name = op_cba
      dagger=.true.
      absym=0
      casym=0
      gamma=1
      s2=0
      ms=0
      ncadiff=0
      min_rank=2
      if(trir12.eq.1)then
        max_rank=3
      else  
        max_rank=2
      endif
      iformal=max_rank+1

      call set_hpvx_and_restr_for_c()
      
      call set_genop(list_pnt%op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info%iad_gas,orb_info%ihpvgas,orb_info%ngas)

      return

      contains
c-----------------------------------------------------------------------
      subroutine set_hpvx_and_restr_for_r()
c-----------------------------------------------------------------------
      implicit none

      integer::
     &     ica,igastp,igas

c      if(.not.(orb_info%nactt_hpv(iextr).gt.0))
c     &     stop 'Stop: Attempt to perform R12 calc. with no external
c     &           orbitals.'


      ! Constraints on the operator are made depending on which ansatz
      ! is being used.
      orb_info%nactt_hpv(ngastp)=2
      ihpv_mnmx(1:2,1:ngastp,1:2)=0
      do ica=1,2
        do igastp=1,ngastp
          if(orb_info%nactt_hpv(igastp).gt.0)then
            if(ica.eq.2.and.igastp.eq.ihole)then
              ihpv_mnmx(1,igastp,ica)=min_h_rank
              ihpv_mnmx(2,igastp,ica)=max_h_rank
            elseif(ica.eq.1)then
              if(igastp.eq.ipart)then
                ihpv_mnmx(1,igastp,ica)=min_p_rank
                ihpv_mnmx(2,igastp,ica)=max_p_rank
              elseif(igastp.eq.iextr)then
                ihpv_mnmx(1,igastp,ica)=min_x_rank
                ihpv_mnmx(2,igastp,ica)=max_x_rank 
              endif
            endif  
          else
            ihpv_mnmx(1,igastp,ica)=0
            ihpv_mnmx(2,igastp,ica)=0
          endif
        enddo
      enddo  

      irestr(1:2,1:orb_info%ngas,1:2,1:2)=0
      do ica=1,2
        do igas=1,orb_info%ngas
          irestr(1,igas,ica,1)=min_rank
          irestr(2,igas,ica,1)=max_rank
        enddo
      enddo  
  
      return
      end subroutine set_hpvx_and_restr_for_r

*----------------------------------------------------------------------*
*     Subroutine to set the restrictions for the coefficient operator of
*     the R12 operators.
*----------------------------------------------------------------------*
      subroutine set_hpvx_and_restr_for_c()
      
      implicit none 

      orb_info%nactt_hpv(ngastp)=2

      ihpv_mnmx(1:2,1:ngastp,1:2)=0
      ihpv_mnmx(1,ihole,1:2)=min_rank
      ihpv_mnmx(2,ihole,1:2)=2
      if(trir12.eq.1)then
        ihpv_mnmx(2,ihole,2)=3
        ihpv_mnmx(2,ipart,1)=1
      endif        

      irestr(1:2,1:orb_info%ngas,1:2,1:2)=0
      irestr(1,ihole,1:2,1)=min_rank
      irestr(2,ihole,1:2,1)=max_rank
      if(trir12.eq.1)then
        irestr(2,ipart,1:2,1)=1
      endif  
      
      return
      end subroutine set_hpvx_and_restr_for_c

      end
