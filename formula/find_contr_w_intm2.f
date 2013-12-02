*----------------------------------------------------------------------*
      subroutine find_contr_w_intm2(success,fpl_found,contr_rpl,
     &                             fl_tgt,fpl_intm,iposs,
     &                             nmod_max,nmod,imod,xmod,
     &                             op_info)
*----------------------------------------------------------------------*
*
*     given: a formula list starting at fl_tgt
*            a pointer list to terms of an intermediate
*
*     check whether current item on formula list contains a term I_i of the
*     intermediate I
*     if yes, factorize  T = T0*I_i and look for the other terms 
*     originating from T0*I_j (i/=j)
*     
*     slightly improved version: allows factorization of intermediates
*       which imply symmetrization of external lines in their defining
*       formula (like Ttilde(iajb) = T(iajb) + T(ia)T(jb))
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
     &     ntest = 00
      
      logical, intent(out) ::
     &     success
      integer, intent(in) ::
     &     iposs, nmod_max
      integer, intent(out) ::
     &     nmod, imod(nmod_max)
      real(8), intent(out) ::
     &     xmod(nmod_max)
      type(formula_item), target, intent(in) ::
     &     fl_tgt
      type(formula_item_list), intent(out), target ::
     &     fpl_found
      type(contraction), intent(inout) ::
     &     contr_rpl
      type(formula_item_list), target, intent(in) ::
     &     fpl_intm
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idxop_tgt, iblk_tgt, iterm, nfound,
     &     idxop_current, iblk_current, njoined, idx, nterms_gen
      integer, pointer ::
     &     len_list(:), iterm_list(:)
      type(formula_item_list), pointer ::
     &     fpl_intm_pnt, fpl_found_pnt
      logical ::
     &     success1, success2, dagger
      logical, pointer ::
     &     assigned(:)
      type(contraction) ::
     &     contr_t0, contr_int
      type(contraction), pointer ::
     &     contr_i
      type(formula_item), pointer ::
     &     fl_tgt_pnt, fl_t0_i, fl_t0_i_pnt

      logical, external ::
     &     contr_in_contr, cmp_contr

      if (ntest.ge.100) then
        write(lulog,*) '============================='
        write(lulog,*) ' entered find_contr_w_intm 2'
        write(lulog,*) '============================='
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
      success2 = .false.
      fpl_intm_pnt => fpl_intm
      ! loop intermediates and check whether I_i is contained in T
      iterm = 0
      allocate(len_list(iposs),iterm_list(iposs))
      len_list = 0
      iterm_list = 0
      nmod = 0
      do
        iterm = iterm+1
        if (contr_in_contr(fpl_intm_pnt%item%contr,
     &                     fl_tgt%contr,           op_info)) then
c          ! take I_i with longest contraction
c          if (fpl_intm_pnt%item%contr%nvtx.gt.len_i_max) then
c            len_i_max = fpl_intm_pnt%item%contr%nvtx
c            ! remember I_i
c            contr_i => fpl_intm_pnt%item%contr
c            success1 = .true.
c          end if
          ! create list of iposs terms with decreasing nvtx
          do idx = 1, iposs
            if (fpl_intm_pnt%item%contr%nvtx.gt.len_list(idx)) then
              if (idx.gt.1) then
                len_list(idx-1) = len_list(idx)
                iterm_list(idx-1) = iterm_list(idx)
              end if
              if (idx.eq.iposs) then
                len_list(idx) = fpl_intm_pnt%item%contr%nvtx
                iterm_list(idx) = iterm
              end if
            else if (idx.gt.1) then
              len_list(idx-1) = fpl_intm_pnt%item%contr%nvtx
              iterm_list(idx-1) = iterm
            end if
          end do
        end if
        if (.not.associated(fpl_intm_pnt%next)) exit
        fpl_intm_pnt => fpl_intm_pnt%next
      end do

      ! take I_i with iposs-th largest number of vertices
      fpl_intm_pnt => fpl_intm
      iterm = 0
      do
        iterm = iterm + 1
        if (iterm.eq.iterm_list(1)) then
          ! remember I_i
          contr_i => fpl_intm_pnt%item%contr
          success1 = .true.
        end if
        if (.not.associated(fpl_intm_pnt%next)) exit
        fpl_intm_pnt => fpl_intm_pnt%next
      end do
      deallocate(len_list,iterm_list)

      if (success1) then
        ! get factor, vertices and arcs associated with T_0
c        if (fl_tgt%contr%nvtx.le.4) then
c         call split_contr2(.true.,contr_t0,contr_i,fl_tgt%contr,op_info)
c        else
          call split_contr3(contr_t0,contr_i,fl_tgt%contr,op_info,
     &                      success1)
c        end if
c        call split_contr2(.true.,contr_t0,contr_i,fl_tgt%contr,op_info)

        if (ntest.ge.100.and.success1) then
          write(lulog,*) 'considering contraction:'
          call prt_contr2(lulog,fl_tgt%contr,op_info)
          write(lulog,*) 'split into T0 '
          call prt_contr2(lulog,contr_t0,op_info)
          write(lulog,*) 'and I'
          call prt_contr2(lulog,contr_i,op_info)
        end if
      end if

      success2 = .true.
      ! more than one term (which should be the usual case)?
      ! look for the remaining terms
      if (success1.and.associated(fl_tgt%next)) then

        allocate(fl_t0_i)
        call init_formula(fl_t0_i)
        fl_t0_i_pnt => fl_t0_i

        success2 = .false.
        fpl_intm_pnt => fpl_intm
        iterm = 0
        do
          ! make target contractions that we need to find
          call join_contr2a(1,fl_t0_i_pnt,nterms_gen,contr_rpl,
     &           contr_t0,fpl_intm_pnt%item%contr,
     &           idxop_tgt,iblk_tgt,op_info)
          iterm = iterm+nterms_gen
          do while(associated(fl_t0_i_pnt%next))
            fl_t0_i_pnt => fl_t0_i_pnt%next
          end do
          if (.not.associated(fpl_intm_pnt%next)) exit
          fpl_intm_pnt => fpl_intm_pnt%next
        end do

        nterms_gen = iterm
        if (ntest.ge.100) then
          write(lulog,*) 'looking for ',nterms_gen,' terms'
        end if

        allocate(assigned(nterms_gen))
        assigned(1:nterms_gen) = .false.

        nfound = 0
        ! continue with present record of fl_tgt
        ! increment to %next will be at top of tgt_loop
        fl_tgt_pnt => fl_tgt
        idxop_current = idxop_tgt

        ! loop items of formula
        tgt_loop: do
          ! result block OK?
          select case (fl_tgt_pnt%command)
          case(command_end_of_formula)
            exit tgt_loop
          case(command_set_target_init,command_set_target_update)
            idxop_current = fl_tgt_pnt%target
            fl_tgt_pnt => fl_tgt_pnt%next
            cycle tgt_loop
          case(command_add_contribution)
            ! we are looking for same operator/block
            if (idxop_current.ne.idxop_tgt) then
              fl_tgt_pnt => fl_tgt_pnt%next
              cycle tgt_loop
            end if
            iblk_current = fl_tgt_pnt%contr%iblk_res
            if (iblk_current.ne.iblk_tgt) then
              fl_tgt_pnt => fl_tgt_pnt%next
              cycle tgt_loop
            end if
          case default
            write(lulog,*) 'command = ',fl_tgt_pnt%command
            call quit(1,'find_contr_w_intm',
     &             'not prepared for that command (see above)')
          end select

          ! compare with generated target contractions
          iterm = 0
          fl_t0_i_pnt => fl_t0_i
          term_loop: do
            iterm = iterm+1
       
            if (.not.associated(fl_t0_i_pnt%contr))
     &           call quit(1,'find_contr_w_intm2',
     &                       'this should not happen')

            if (.not.assigned(iterm)) then
c dbg
c              print *,'comparing: iterm = ',iterm
c              print *,'assigned: ',assigned(1:nterms_gen)
c              call prt_contr2(6,fl_tgt_pnt%contr,op_info)
c              call prt_contr2(6,fl_t0_i_pnt%contr,op_info)
c dbg
              if (cmp_contr(fl_tgt_pnt%contr,
     &                      fl_t0_i_pnt%contr,.false.)) then

c      If you wish to factor out two identical intermediates in a term,
c      e.g. 1/2*(a+b)(a+b), comment out the previous two lines
c      and comment in the following three lines
c      and the paragraph titled "would factor have to be changed?"
c      Then factor out (a+b) twice. This should lead to
c      1/2*aa+ab+1/2*bb --> 1/2*a(a+b)+1/2*b(a+b) --> 1/2*(a+b)(a+b)
c              if (cmp_contr(fl_tgt_pnt%contr,
c     &                      fl_t0_i_pnt%contr,.true.)
c     &            .and.nmod.lt.nmod_max) then

c dbg
c                print *,'OK!'
c dbg
                assigned(iterm) = .true.
                nfound = nfound+1
                if (nfound.eq.1) then
                  fpl_found_pnt => fpl_found
                else
                  call new_formula_plist_entry(fpl_found_pnt)
                  fpl_found_pnt => fpl_found_pnt%next
                end if
c      Please don't delete, could be of use (see comment above)
c                ! would factor have to be changed?
c                if (abs(fl_tgt_pnt%contr%fac-fl_t0_i_pnt%contr%fac)
c     &              .ge.1d-12) then
c                  if (nfound.eq.1) call quit(1,'find_contr_w_intm2',
c     &                 'factor change of leading term needed')
c                  nmod = nmod + 1
c                  if (nmod.eq.nmod_max) call warn('find_contr_w_intm2',
c     &            'consider increasing nmod_max in factor_out_subexpr2')
c                  imod(nmod) = nfound
c                  xmod(nmod) = fl_tgt_pnt%contr%fac
c     &                         - fl_t0_i_pnt%contr%fac
c                end if
                fpl_found_pnt%item => fl_tgt_pnt
                ! all terms found? let's go
                success2 =  nfound.eq.nterms_gen
                if (success2) exit tgt_loop
                exit term_loop
              end if
            end if

            if (.not.associated(fl_t0_i_pnt%next).or.
     &               fl_t0_i_pnt%next%command.eq.command_end_of_formula)
     &           exit term_loop
            fl_t0_i_pnt => fl_t0_i_pnt%next

          end do term_loop

          ! go to next record
          if (.not.associated(fl_tgt_pnt%next)) exit tgt_loop
          fl_tgt_pnt => fl_tgt_pnt%next


        end do tgt_loop

        deallocate(assigned)
        call dealloc_formula_list(fl_t0_i)
        deallocate(fl_t0_i)

      end if

      if (ntest.ge.100) then
        write(lulog,*) 'at the end: ',success1,success2
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

        allocate(fl_t0_i)
        call init_formula(fl_t0_i)
        fl_t0_i_pnt => fl_t0_i

        ! make new contraction
        call join_contr2a(0,fl_t0_i_pnt,nterms_gen,contr_rpl,
     &                  contr_t0,contr_int,
     &                  idxop_tgt,iblk_tgt,op_info)

        call dealloc_formula_list(fl_t0_i)
        deallocate(fl_t0_i)

        call dealloc_contr(contr_int)

        if (ntest.ge.100) then
          write(lulog,*) 'generated term:'
          call prt_contr2(lulog,contr_rpl,op_info)
        end if

      end if

      call dealloc_contr(contr_t0)

      return
      end
