*----------------------------------------------------------------------*
      subroutine find_contr_w_intm(success,fpa_found,contr_rpl,
     &                             fl_tgt,fpl_intm,nterms,
     &                             op_info)
*----------------------------------------------------------------------*
*
*     given: a formula list starting at fl_tgt
*            a pointer list to terms of an intermediate (nterms terms)
*
*     check whether current item on formula list contains a term I_i of the
*     intermediate I
*     if yes, factorize  T = T0*I_i and look for the other terms 
*     originating from T0*I_j (i/=j)
*     
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'
      include 'def_formula_item_array.h'

      integer, parameter ::
     &     ntest = 100
      
      logical, intent(out) ::
     &     success
      integer, intent(in) ::
     &     nterms
      type(formula_item), target, intent(in) ::
     &     fl_tgt
      type(formula_item_array), intent(out) ::
     &     fpa_found(nterms)
      type(contraction), intent(inout) ::
     &     contr_rpl
      type(formula_item_list), target, intent(in) ::
     &     fpl_intm
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     len_i_max, idxop_tgt, iblk_tgt, iterm, nfound,
     &     idxop_current, iblk_current, njoined, idx
      type(formula_item_list), pointer ::
     &     fpl_intm_pnt
      logical ::
     &     assigned(nterms), success1, success2, dagger
      type(contraction) ::
     &     contr_t0, contr_int
      type(contraction), pointer ::
     &     contr_tgt(:), contr_i
      type(formula_item), pointer ::
     &     fl_tgt_pnt

      logical, external ::
     &     contr_in_contr, cmp_contr

      if (ntest.ge.100) then
        write(luout,*) '==========================='
        write(luout,*) ' entered find_contr_w_intm'
        write(luout,*) '==========================='
      end if

      call init_contr(contr_t0)

      if (fl_tgt%command.ne.command_add_contribution)
     &     call quit(1,'find_contr_w_intm','[ADD] expected on entry')
      ! get current operator, blk
      idxop_tgt = fl_tgt%contr%idx_res
      iblk_tgt = fl_tgt%contr%iblk_res
c      iocc_tgt = op_info%op_arr(idxop_tgt)%
c     &           op%ihpvca_occ(1:ngastp,1:2,iblk_tgt)

      success1 = .false.
      fpl_intm_pnt => fpl_intm
      ! loop intermediates and check whether I_i is contained in T
      len_i_max = 0
      iterm = 0
      do
        iterm = iterm+1
        if (contr_in_contr(fpl_intm_pnt%item%contr,
     &                     fl_tgt%contr,           op_info)) then
          ! take I_i with longest contraction
          if (fpl_intm_pnt%item%contr%nvtx.gt.len_i_max) then
            len_i_max = fpl_intm_pnt%item%contr%nvtx
            ! let first element of fpa_found point to T
            fpa_found(1)%item => fl_tgt
            ! remember I_i
            contr_i => fpl_intm_pnt%item%contr
            success1 = .true.
            ! mark this term as "assigned"
            assigned(1:nterms) = .false.
            assigned(iterm) = .true.
          end if
        end if
        if (.not.associated(fpl_intm_pnt%next)) exit
        fpl_intm_pnt => fpl_intm_pnt%next
      end do

c dbg
c      print *,'success 1',success1
      print *,'nterms',nterms
c      print *,'assoc',associated(fl_tgt%next)
c dbg

      if (success1) then
        ! get factor, vertices and arcs associated with T_0
c        call split_contr(contr_t0,contr_i,fl_tgt%contr,op_info)
        call split_contr2(.true.,contr_t0,contr_i,fl_tgt%contr,op_info)
        if (ntest.ge.100) then
          write(luout,*) 'considering contraction:'
          call prt_contr2(luout,fl_tgt%contr,op_info)
          write(luout,*) 'split into T0 '
          call prt_contr2(luout,contr_t0,op_info)
          write(luout,*) 'and I'
          call prt_contr2(luout,contr_i,op_info)
        end if
      end if

      ! more than one term (which should be the usual case)?
      ! look for the remaining terms
      if (success1.and.nterms.gt.1.and.associated(fl_tgt%next)) then

        allocate(contr_tgt(nterms))
        do iterm = 1, nterms
          call init_contr(contr_tgt(iterm))
        end do

        success2 = .false.
        fpl_intm_pnt => fpl_intm
        iterm = 0
        do
          iterm = iterm+1
          ! make target contractions that we need to find
          if (.not.assigned(iterm))then
            call join_contr2(contr_tgt(iterm),
     &           contr_t0,fpl_intm_pnt%item%contr,
     &           idxop_tgt,iblk_tgt,op_info)
c dbg
            print *,'iterm',iterm
            print *,'targeted'
            call prt_contr2(luout,contr_tgt(iterm),op_info)
c dbg
          endif    
          if (.not.associated(fpl_intm_pnt%next)) exit
          fpl_intm_pnt => fpl_intm_pnt%next
        end do

        nfound = 1 ! the first term is already on fpa_found list
        ! continue with next record of fl_tgt
        ! increment to %next will be at top of tgt_loop
        fl_tgt_pnt => fl_tgt
        idxop_current = idxop_tgt

        ! loop items of formula
        tgt_loop: do
          ! go to next record
          if (.not.associated(fl_tgt_pnt%next)) exit tgt_loop
          fl_tgt_pnt => fl_tgt_pnt%next

          ! result block OK?
          select case (fl_tgt_pnt%command)
            case(command_end_of_formula)
              exit tgt_loop
            case(command_set_target_init,command_set_target_update)
              idxop_current = fl_tgt_pnt%target
              cycle tgt_loop
            case(command_add_contribution)
              ! we are looking for same operator/block
              if (idxop_current.ne.idxop_tgt) cycle tgt_loop
              iblk_current = fl_tgt_pnt%contr%iblk_res
              if (iblk_current.ne.iblk_tgt) cycle tgt_loop
            case default
              write(luout,*) 'command = ',fl_tgt_pnt%command
              call quit(1,'find_contr_w_intm',
     &             'not prepared for that command (see above)')
          end select

          ! compare with generated target contractions
          term_loop: do iterm = 1, nterms
            if (assigned(iterm)) cycle
c dbg
            print *,'comparing: iterm = ',iterm
            print *,'assigned: ',assigned(1:nterms)
            call prt_contr2(6,fl_tgt_pnt%contr,op_info)
            call prt_contr2(6,contr_tgt(iterm),op_info)
c dbg
            if (cmp_contr(fl_tgt_pnt%contr,
     &                    contr_tgt(iterm),.false.)) then
c dbg
c              print *,'hurra'
c dbg
              assigned(iterm) = .true.
              nfound = nfound+1
c dbg
              print *,'nfound',nfound
c dbg
              fpa_found(nfound)%item => fl_tgt_pnt
              ! all terms found? let's go
              success2 =  nfound.eq.nterms
              if (success2) exit tgt_loop
              exit term_loop
            end if
          end do term_loop

        end do tgt_loop

        deallocate(contr_tgt)

      else if (nterms.eq.1) then
        success2 = .true.
      end if

      if (ntest.ge.100) then
        write(luout,*) 'at the end: ',success1,success2
      end if

      success = success1.and.success2

      ! provide contraction with intermediate
      if (success) then
        ! contr_t0 is still in memory
        ! set up contr_i, only one super-vertex
        call init_contr(contr_int)
        ! result is operator itself
        contr_int%idx_res = fpl_intm%item%contr%idx_res
        contr_int%iblk_res = fpl_intm%item%contr%iblk_res
        contr_int%dagger = fpl_intm%item%contr%dagger
        njoined = op_info%op_arr(contr_int%idx_res)%op%njoined
        dagger  = contr_int%dagger
        call resize_contr(contr_int,njoined,0,0,0)
        contr_int%nvtx = njoined
        contr_int%nsupvtx = 1
        contr_int%svertex(1:njoined) = 1
        call update_svtx4contr(contr_int)
        contr_int%narc = 0
        contr_int%nfac = 0
        contr_int%fac = 1d0
        contr_int%vertex(1:njoined)%idx_op = contr_int%idx_res
        contr_int%vertex(1:njoined)%dagger = dagger
        if (.not.dagger) then
          do idx = 1, njoined
            contr_int%vertex(idx)%iblk_op =
     &           (contr_int%iblk_res-1)*njoined+idx
          end do
        else
          do idx = 1, njoined
            contr_int%vertex(idx)%iblk_op =
     &           (contr_int%iblk_res)*njoined+1-idx
          end do
        end if

        ! make new contraction
c dbg
c        print *,'the interesting call to join_contr2'
c dbg
        call join_contr2(contr_rpl,
     &                  contr_t0,contr_int,
     &                  idxop_tgt,iblk_tgt,op_info)

        call dealloc_contr(contr_int)

        if (ntest.ge.100) then
          write(luout,*) 'generated term:'
          call prt_contr2(luout,contr_rpl,op_info)
        end if

      end if

      call dealloc_contr(contr_t0)

      return
      end
