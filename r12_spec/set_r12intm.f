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
        write(luout,*) '-------------'
        write(luout,*) ' set_r12intm '
        write(luout,*) '-------------'
        write(luout,*) 'Intermediate: ',trim(name)
c        if(r12fix) write(luout,*)'Fully contracted intermediate'
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

      ! Extra terms required for CC-R12 intermediates.
      ! Initially only for V.
      if((trim(name).eq.op_v_inter.or.trim(name).eq.op_vbar_inter)
     &     .and.cclog)then
        if(trim(name).eq.op_v_inter)then
          i=1
          j=2
        else
          i=2
          j=1
        endif
        allocate(opscr(2))
        do idx = 1, 2
          ! Shape of the operator.
          occ0 = 0
          call add_operator(op_dummy1,op_info)
          idx_dum = idx_oplist2(op_dummy1,op_info)
          dum1_pnt => op_info%op_arr(idx_dum)%op
          occ0(2,i) = idx
          call set_uop(dum1_pnt,op_dummy1,.false.,
     &         occ0,1,orb_info)

c          if(.not.r12fix)then
            ! Spacer operator.
            occ0 = 0
            call add_operator(op_dummy2,op_info)
            idx_dum = idx_oplist2(op_dummy2,op_info)
            dum2_pnt => op_info%op_arr(idx_dum)%op
            occ0(1,i) = 2
            occ0(1,j) = 2-idx
            call set_uop(dum2_pnt,op_dummy2,.false.,
     &           occ0,1,orb_info)
        
            opscr(1)%op => dum1_pnt
            opscr(2)%op => dum2_pnt

            ! Temporary intermediate.
            call add_operator(op_dummy3,op_info)
            idx_dum = idx_oplist2(op_dummy3,op_info)
            dum3_pnt => op_info%op_arr(idx_dum)%op
          
            call set_gen_intermediate(dum3_pnt,op_dummy3,
     &           opscr,2,orb_info)

            call join_operator(op,dum3_pnt,orb_info)

            call del_operator(op_dummy3,op_info)
            call del_operator(op_dummy2,op_info)
c          else
c            call join_operator(op,dum1_pnt,orb_info)
c          endif
          call del_operator(op_dummy1,op_info)
        enddo
        deallocate(opscr)
      endif

      if (ntest.ge.100) then
        write(luout,*) 'generated: '
        call print_op_occ(luout,op)
        write(luout,*) 'formal: ',op%formal_blk(1:op%n_occ_cls)
      end if

      return
      
      end
