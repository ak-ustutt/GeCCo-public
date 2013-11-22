*----------------------------------------------------------------------*
      subroutine set_r12intm(op,name,dagger,
     &     min_rank,max_rank,ncadiff,iformal,op_info,orb_info)
*----------------------------------------------------------------------*
*     set up intermediate for R12 with open lines within
*     settings analoguous to R12coeff
*     hpvx_mnmx,irestr are chosen appropriately
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'opdim.h'
c      include 'def_operator.h'
c      include 'def_operator_array.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'par_opnames_gen.h'
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 00

      type(operator), intent(inout) ::
     &     op
      character, intent(in) ::
     &     name*(*)
      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     min_rank, max_rank, ncadiff, iformal

      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(inout) ::
     &     orb_info

      integer ::
     &     idx, idx_dum, i, j
      integer ::
     &     occ0(ngastp,2)
      logical ::
c     &     r12fix
     &     cclog
      type(operator), pointer ::
     &     dum1_pnt, dum2_pnt, dum3_pnt
      type(operator_array), pointer ::
     &     opscr(:)

      character(5) ::
     &     op_dummy1*5 = 'DUM-1',
     &     op_dummy2*5 = 'DUM-2',
     &     op_dummy3*5 = 'DUM-3'
     

      integer, external ::
     &     idx_oplist2

c      ! Do we want amplitudeless intermediates?
c      call get_argument_value('method','CC',lval=cclog)
      cclog = is_keyword_set('method.CC').gt.0

      if(ntest.ge.100)then
        write(lulog,*) '-------------'
        write(lulog,*) ' set_r12intm '
        write(lulog,*) '-------------'
        write(lulog,*) 'Intermediate: ',trim(name)
c        if(r12fix) write(lulog,*)'Fully contracted intermediate'
      endif

c      if(.not.r12fix)then
        allocate(opscr(2))
        allocate(opscr(1)%op)
        allocate(opscr(2)%op)

        ! no external lines, so:
        occ0 = 0
        call set_uop(opscr(1)%op,'scr1',.false.,
     &       occ0,1,orb_info)

        ! the internal lines are defined by a coefficient-like 
        ! operator, so we use ...
        call set_r12c(opscr(2)%op,'scr2',dagger,
     &       min_rank,max_rank,ncadiff,iformal,orb_info)
        
        ! ... and define the intermediate:
        call set_gen_intermediate(op,name,
     &       opscr,2,orb_info)

        call dealloc_operator(opscr(1)%op)
        call dealloc_operator(opscr(2)%op)
        deallocate(opscr(1)%op)
        deallocate(opscr(2)%op)
        deallocate(opscr)
c      else
c        ! no external lines, so:
c        occ0 = 0
c        if(trim(name).ne.op_x_inter)then
c          call set_uop(op,name,.false.,
c     &         occ0,1,orb_info)
c        else
c          allocate(opscr(2))
c          allocate(opscr(1)%op)
c          allocate(opscr(2)%op)
c          call set_uop(opscr(1)%op,'scr1',.false.,
c     &         occ0,1,orb_info)
c
c          occ0(1,1:2) = 1
c          call set_uop(opscr(2)%op,'scr2',.false.,
c     &         occ0,1,orb_info)
c
c          call set_gen_intermediate(op,name,
c     &         opscr,2,orb_info)
c
c          call dealloc_operator(opscr(1)%op)
c          call dealloc_operator(opscr(2)%op)
c          deallocate(opscr(1)%op)
c          deallocate(opscr(2)%op)
c          deallocate(opscr)
c        endif
        do idx = 1, op%n_occ_cls
          op%formal_blk(idx) = .false.
        enddo

c      endif

      op%formal = .true.
      do idx = 1, op%n_occ_cls
        op%formal = op%formal.and.op%formal_blk(idx)
      end do

c      ! Extra terms required for CC-R12 intermediates.
c      ! Initially only for V.
c      if((trim(name).eq.op_v_inter.or.trim(name).eq.op_vbar_inter)
c     &     .and.cclog)then
c        if(trim(name).eq.op_v_inter)then
c          i=1
c          j=2
c        else
c          i=2
c          j=1
c        endif
c        allocate(opscr(2))
c        do idx = 1, 2
c          ! Shape of the operator.
c          occ0 = 0
c          call add_operator(op_dummy1,op_info)
c          idx_dum = idx_oplist2(op_dummy1,op_info)
c          dum1_pnt => op_info%op_arr(idx_dum)%op
c          occ0(2,i) = idx
c          call set_uop(dum1_pnt,op_dummy1,.false.,
c     &         occ0,1,orb_info)
c
cc          if(.not.r12fix)then
c            ! Spacer operator.
c            occ0 = 0
c            call add_operator(op_dummy2,op_info)
c            idx_dum = idx_oplist2(op_dummy2,op_info)
c            dum2_pnt => op_info%op_arr(idx_dum)%op
c            occ0(1,i) = 2
c            occ0(1,j) = 2-idx
c            call set_uop(dum2_pnt,op_dummy2,.false.,
c     &           occ0,1,orb_info)
c        
c            opscr(1)%op => dum1_pnt
c            opscr(2)%op => dum2_pnt
c
c            ! Temporary intermediate.
c            call add_operator(op_dummy3,op_info)
c            idx_dum = idx_oplist2(op_dummy3,op_info)
c            dum3_pnt => op_info%op_arr(idx_dum)%op
c          
c            call set_gen_intermediate(dum3_pnt,op_dummy3,
c     &           opscr,2,orb_info)
c
c            call join_operator(op,dum3_pnt,orb_info)
c
c            call del_operator(op_dummy3,op_info)
c            call del_operator(op_dummy2,op_info)
cc          else
cc            call join_operator(op,dum1_pnt,orb_info)
cc          endif
c          call del_operator(op_dummy1,op_info)
c        enddo
c        deallocate(opscr)
c      endif

      if (ntest.ge.100) then
        write(lulog,*) 'generated: '
        call print_op_occ(lulog,op)
        write(lulog,*) 'formal: ',op%formal_blk(1:op%n_occ_cls)
      end if

      return
      
      end
