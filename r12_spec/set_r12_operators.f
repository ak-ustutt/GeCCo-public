*----------------------------------------------------------------------*
      subroutine set_r12_operators(op_info,orb_info)
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
      include 'mdef_operator_info.h'
      include 'explicit.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      type(orbinf) ::
     &     orb_info
      type(operator), pointer ::
     &     op_pnt, r12_pnt, c12_pnt, cba_pnt, rint_pnt, ecc_pnt,
     &     v_pnt, b_pnt
      type(operator_array),pointer ::
     &     ops_array(:)
      logical ::
     &     dagger
      character ::
     &     name*(len_opname)
      integer ::
     &     absym, casym, s2, ms, min_rank, max_rank, ncadiff,
     &     gamma, iarr(1),  min_h_rank, max_h_rank,
     &     min_p_rank, max_p_rank, min_x_rank, max_x_rank, iformal,
     &     tkmax, idx, idx_top, idx_c12, idx_ecc
      integer ::
     &     ihpv_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2),
     &     occ_def(ngastp,2,1)

      integer, external ::
     &     idx_oplist2

c      write(luout,'(/"R12 operator definition subroutine.")')

c      call get_argument_value('method.R12','ansatz',ival=ansatze)
c      if(ansatze.gt.3.or.ansatze.lt.1)then
c        write(luout,'("Error: Undefined R12 ansatz requested.")')
c        stop
c      endif
c
c      call get_argument_value('method.R12','triples',ival=trir12)

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

      call add_operator(op_r12,op_info)
      idx = idx_oplist2(op_r12,op_info)
      op_pnt => op_info%op_arr(idx)%op
      r12_pnt => op_pnt

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
      iformal=0
  
      call set_hpvx_and_restr_for_r()

      call set_genop(op_pnt,name,optyp_operator,
     &     dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info)  
      
      ! New entry: variable coefficient operator associated with R12.
      call add_operator(op_c12,op_info)
      idx_c12 = idx_oplist2(op_c12,op_info)
      op_pnt => op_info%op_arr(idx_c12)%op
      c12_pnt => op_pnt
      
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
      
      call set_genop(op_pnt,name,optyp_operator,
     &     dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info)

      ! New entry: adjoint of linear R12 operator.
      call add_operator(op_rba,op_info)
      idx = idx_oplist2(op_rba,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,r12_pnt,orb_info)
      ! we define an excitation operator to ensure same
      ! storage sequence as for T
      op_pnt%dagger = .true.  ! but we consider the conjugate

      ! New entry: variable coefficient operator associated with R12bar.
      call add_operator(op_cba,op_info)
      idx = idx_oplist2(op_cba,op_info)
      op_pnt => op_info%op_arr(idx)%op
      cba_pnt=>op_pnt
      call clone_operator(op_pnt,c12_pnt,orb_info)
      ! we define an excitation operator to ensure same
      ! storage sequence as for T
      op_pnt%dagger = .true.  ! but we consider the conjugate

c      ! Form the total residual for the CC-R12 function: OMG=C+T
c      call add_operator(op_omg_candt,op_info)
c      idx = idx_oplist2(op_omg_candt,op_info)
c      op_pnt => op_info%op_arr(idx)%op

c      idx_top=idx_oplist2(op_top,op_info)
c      call clone_operator(op_pnt,op_info%op_arr(idx_top)%op,orb_info)
c      call join_operator(op_pnt,op_info%op_arr(idx_c12)%op,orb_info)

      ! New entry: the R12 integrals (<ab|r12|cd>)
      call add_operator(op_rint,op_info)
      idx = idx_oplist2(op_rint,op_info)
      op_pnt => op_info%op_arr(idx)%op
      rint_pnt => op_pnt

      name=op_rint
      dagger=.false.
      absym=0
      casym=0
      gamma=1
      s2=0
      ms=0
      ncadiff=0
      min_rank=2
      max_rank=2
      call set_hpvx_and_restr_for_int()
      iformal=2

      call set_genop(op_pnt,name,optyp_operator,
     &     dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info)

      ! New entry: the R12 integrals' adjoint (<ab|r12|cd>)
      call add_operator(op_rinba,op_info)
      idx = idx_oplist2(op_rinba,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,rint_pnt,orb_info)
      op_pnt%dagger = .true.

      ! New entry: commutator integrals <ab|[T1+T2,r12]|cd>
      call add_operator(op_ttr,op_info)
      idx=idx_oplist2(op_ttr,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,rint_pnt,orb_info)

      ! <ab|[T1+T2,r12]|cd>+
      call add_operator(op_ttr_bar,op_info)
      idx=idx_oplist2(op_ttr_bar,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,rint_pnt,orb_info)
      op_pnt%dagger = .true.

      ! Definition of intermediates (operators where there are 
      ! uncontracted arcs between constituent operator lines).
      ! Get indices of required operators.
      idx_ecc=idx_oplist2(op_ccen,op_info)
      ecc_pnt => op_info%op_arr(idx_ecc)%op

      ! Anti-symmetric V^{ij}_{kl}.
      call add_operator(op_v_inter,op_info)
      idx=idx_oplist2(op_v_inter,op_info)
      op_pnt => op_info%op_arr(idx)%op
      v_pnt => op_pnt
      allocate(ops_array(2))
      ops_array(1)%op=>ecc_pnt
      ops_array(2)%op=>c12_pnt
      call set_gen_intermediate(op_pnt,op_v_inter,
     &     ops_array,2,orb_info)
      deallocate(ops_array)

      ! Adjoint of V.
      call add_operator(op_vbar_inter,op_info)
      idx=idx_oplist2(op_vbar_inter,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,v_pnt,orb_info)
      ! This operator should be daggered, but this causes inconsistencies
      ! later on. Need to rectify this at some point. GWR 14/08/2007
c      op_pnt%dagger=.true.

      ! New entry: the CC-R12 residual. Same as V.
      call add_operator(op_omgr12,op_info)
      idx = idx_oplist2(op_omgr12,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,v_pnt,orb_info)

      ! Delta integrals.
      call add_operator(op_del_inter,op_info)
      idx=idx_oplist2(op_del_inter,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,v_pnt,orb_info)

      ! Anti-symmetric B^{ij}_{kl}.
      call add_operator(op_b_inter,op_info)
      idx=idx_oplist2(op_b_inter,op_info)
      op_pnt => op_info%op_arr(idx)%op
      b_pnt => op_pnt
      allocate(ops_array(2))
      ops_array(1)%op=>ecc_pnt
      ops_array(2)%op=>cba_pnt
c      ops_array(3)%op=>c12_pnt
      call set_gen_intermediate(op_pnt,op_b_inter,
     &     ops_array,2,orb_info)
      deallocate(ops_array)

      ! Symmetrised B (for MP2 approx 1A).
      call add_operator(op_b_symm,op_info)
      idx=idx_oplist2(op_b_symm,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,b_pnt,orb_info)

      ! Inverse of the diagonal of B for MP2-R12 approx 1A.
      call add_operator(op_b_inv,op_info)
      idx=idx_oplist2(op_b_inv,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,b_pnt,orb_info)

      ! New entry: the CC-R12 diagonal. Same as B (at least for approx A).
      call add_operator(op_diar12,op_info)
      idx = idx_oplist2(op_diar12,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,b_pnt,orb_info)

c      ! Adjoint of B.
c      call add_operator(op_bbar_inter,op_info)
c      idx=idx_oplist2(op_bbar_inter,op_info)
c      op_pnt => op_info%op_arr(idx)%op
c      call clone_operator(op_pnt,b_pnt,orb_info)
c      op_pnt%dagger=.true.

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

c-----------------------------------------------------------------------
      subroutine set_hpvx_and_restr_for_s()
c-----------------------------------------------------------------------
      implicit none

      integer::
     &     ica,igastp,igas

      ! Constraints on the operator are made depending on which ansatz
      ! is being used.
      do ica = 1, 2
        do igastp = 1, ngastp
          if (orb_info%nactt_hpv(igastp).gt.0.and.
     &        ((ica.eq.1.and.(igastp.eq.ipart.or.igastp.eq.ivale)).or.
     &         (ica.eq.2.and.(igastp.eq.ihole.or.igastp.eq.ivale)) ))
     &           then
            ihpv_mnmx(1,igastp,ica) = 0
            ihpv_mnmx(2,igastp,ica) = max_rank
          else
            ihpv_mnmx(1,igastp,ica) = 0
            ihpv_mnmx(2,igastp,ica) = 0
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
      end subroutine set_hpvx_and_restr_for_s

c-----------------------------------------------------------------------
      subroutine set_hpvx_and_restr_for_int()
c-----------------------------------------------------------------------
      implicit none

      integer::
     &     igastp,igas

      ihpv_mnmx(1:2,1:orb_info%ngas,1:2)=0
      
      ihpv_mnmx(1:2,ihole,2)=2
      
      do igastp=1,ngastp
        if(igastp.ne.ivale)then
          ihpv_mnmx(1,igastp,1)=0
          ihpv_mnmx(2,igastp,1)=max_rank
        endif  
      enddo

      irestr(1:2,1:orb_info%ngas,1:2,1:2)=0
      
      irestr(1:2,ihole,2,1)=max_rank
      
      do igas=1,orb_info%ngas
        if(igastp.ne.ivale)then
          irestr(1,igas,1,1)=0
          irestr(2,igas,1,1)=max_rank
        endif  
      enddo

      return
      end subroutine set_hpvx_and_restr_for_int

      end
