*----------------------------------------------------------------------*
      subroutine optimize_kernel(pass,lti_cnt,success,fl_opt,
     &     op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     automatic optimization of the sequence of binary contractions
*     on fl_opt
*
*     andreas, oct. 2012
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'ioparam.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'
      include 'ifc_fl_aux.h'

      integer, parameter ::
     &     ntest = 000
 
      integer, parameter ::
     &    min_replace = 10,
     &    min_scal = 6
     
      character(len=15), parameter ::
     &     i_am = 'optimize_kernel'

      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(inout) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(formula_item), intent(inout), target ::
     &     fl_opt
      integer, intent(in) ::
     &     pass
      integer, intent(inout) ::
     &     lti_cnt
      logical, intent(out) ::
     &     success

      character(len=len_opname) ::
     &     label_old, label_new, label_op
      logical ::
     &     ok, tra_op, ldummy, skip
      integer ::
     &     idxfl, ngas, nspin, nlist, ilist, nlist2, rank, nj,
     &     cost_intm, size_intm, cost_ref, idx_rpl, iblk_op,
     &     idummy, icmp
      real(8) ::
     &     fact, fact_itf, 
     &     cpu, sys, wall, cpu0, sys0, wall0
      type(formula_item), pointer ::
     &     fl_pnt_o, fl_pnt_i, 
     &     fl_pnt_mark1, fl_pnt_mark2, fl_pnt_mark3, fl_pnt_mark4
c dbg
     &     ,fl_pnt_watch
c dbg
      type(formula_item_list), target ::
     &     fpl_marks, fpl_marks2
      type(formula_item_list), pointer ::
     &     fpl_marks_pnt, fpl_marks2_pnt
      type(operator), pointer ::
     &     intm_new, op

      integer, pointer ::
     &     occ(:,:,:), rst(:,:,:,:,:,:), svdummy(:)
      
      logical, external ::
     &     list_cmp
      integer, external ::
     &     cmp_bcontr_spc, ielsum, len_blk_est
c dbg
      logical watch
c
c
      watch = .false.
c dbg


      if (ntest.gt.0) call write_title(lulog,wst_dbg_subr,i_am)
      if (ntest.gt.0) write(lulog,*) ' pass ',pass
c dbg
c      print *,'check at entry'
c      call check_formula_list(fl_opt,op_info)
c dbg

      ngas = orb_info%ngas
      nspin = orb_info%nspin

      success = .false.
      fl_pnt_o => fl_opt
      idxfl = 0

      ! reference for cost estimate
      ! n^6 with at least two virtual or ext, or n^8 in active
      cost_ref = 3000
      cost_ref = max(cost_ref,
     &           len_blk_est((/2,1,0,0,2,1,0,0/),
     &                       1,orb_info))
      cost_ref = max(cost_ref,
     &           len_blk_est((/1,1,1,0,1,1,1,0/),
     &                       1,orb_info))
      cost_ref = max(cost_ref,
     &           len_blk_est((/0,1,2,0,0,1,2,0/),
     &                       1,orb_info))
      cost_ref = max(cost_ref,
     &           len_blk_est((/0,1,3,0,0,1,3,0/),
     &                       1,orb_info))

      outer_loop: do

C        if (ntest.ge.1000) 
C     &      call print_form_item(lulog,idxfl,fl_pnt_o,op_info)

        if (fl_pnt_o%command.eq.command_end_of_formula.or.
     &      fl_pnt_o%command.eq.command_set_target_init) 
     &    exit outer_loop

        ! for pass=0: look for [NEW INTERMEDIATE]
        if (pass.eq.0) then
          if (fl_pnt_o%command.ne.command_new_intermediate) then
            fl_pnt_o => fl_pnt_o%next
            cycle outer_loop
          end if

          ! none of the constituents must be a _STIN
          if (fl_pnt_o%parent1(1:5).eq.'_STIN'.or.
     &        fl_pnt_o%parent2(1:5).eq.'_STIN') then
            fl_pnt_o => fl_pnt_o%next
            cycle outer_loop
          end if

          ! if the next entry is not a [CONTRACT] we skip here
          ! should be adapted for [REO]
          ! we also skip for traces
          if (fl_pnt_o%next%command.ne.command_bc.and.
     &        fl_pnt_o%next%command.ne.command_bc_reo) then
            fl_pnt_o => fl_pnt_o%next
            cycle outer_loop
          end if
          ! we also skip for traces
          if (fl_pnt_o%next%bcontr%n_operands.ne.2) then
            fl_pnt_o => fl_pnt_o%next
            cycle outer_loop
          end if

          if (ntest.ge.100)  then
            write(lulog,*) 'now examining: '
            call print_form_item(lulog,idxfl,fl_pnt_o,op_info)
            call print_form_item(lulog,idxfl,fl_pnt_o%next,op_info)
          end if

          ! mark this entry and proceed to next entry (should be
          ! the contraction that has the new intm. as result)
          fl_pnt_mark1 => fl_pnt_o
          label_old = fl_pnt_mark1%interm%name
          nj   = fl_pnt_mark1%interm%njoined
          rank = ielsum(fl_pnt_mark1%interm%ica_occ,2*nj)
          fl_pnt_o => fl_pnt_o%next
          fl_pnt_mark2 => fl_pnt_o

          ! check this assumption
          if (trim(fl_pnt_mark2%bcontr%label_res).ne.label_old) 
     &      goto 1001
c dbg
c          print *,'OK'
c dbg

        else

          ! other passes: look for contractions
          ! we exclude reorderings for the moment
          if (fl_pnt_o%command.ne.command_bc.and.
     &        fl_pnt_o%command.ne.command_add_bc)!.and.
!     &        fl_pnt_o%command.ne.command_bc_reo.and.
!     &        fl_pnt_o%command.ne.command_add_bc_reo) then
     &     then
            fl_pnt_o => fl_pnt_o%next
            cycle outer_loop          
          end if

          ! should be a binary contr.
          ! and result must not be an _STIN
          if (fl_pnt_o%bcontr%n_operands.ne.2.or.
     &        fl_pnt_o%bcontr%label_res(1:5).eq.'_STIN') then
            fl_pnt_o => fl_pnt_o%next
            cycle outer_loop
          end if

          if (ntest.ge.100)  then
c            write(lulog,*) 'now examining: '
c            call print_form_item(lulog,idxfl,fl_pnt_o,op_info)
          end if

          ! mark this entry
          fl_pnt_mark1 => fl_pnt_o

        end if

        fl_pnt_i => fl_pnt_o   

        call init_formula_plist(fpl_marks)
        fpl_marks_pnt => fpl_marks
        nlist = 0 
        if (pass.gt.0) then
          call init_formula_plist(fpl_marks2)
          fpl_marks2_pnt => fpl_marks2
          nlist2 = 0
        end if

c dbg
c        print *,'starting inner loop '
c        call print_form_item(lulog,idxfl,fl_pnt_i%next,op_info)
c        print *,'form until end:'
c        call print_form_list(lulog,fl_pnt_i%next,op_info)
c dbg

        inner_loop: do

          if (.not.associated(fl_pnt_i%next)) goto 1003
          fl_pnt_i => fl_pnt_i%next
          if (fl_pnt_i%command.eq.command_end_of_formula.or.
     &        fl_pnt_i%command.eq.command_set_target_init) 
     &      exit inner_loop

          if (pass.eq.0) then
            ! pass=0: look for [NEW INTERMEDIATE]
            if (fl_pnt_i%command.ne.command_new_intermediate) 
     &         cycle inner_loop
            ! interm. must have the same parents:
            if (trim(fl_pnt_i%parent1).ne.trim(fl_pnt_mark1%parent1).or.
     &          trim(fl_pnt_i%parent2).ne.trim(fl_pnt_mark1%parent2))
     &         cycle inner_loop
            ! transposition of parents must be the same
            if (fl_pnt_i%tra1.neqv.fl_pnt_mark1%tra1.or.
     &          fl_pnt_i%tra2.neqv.fl_pnt_mark1%tra2)
     &         cycle inner_loop
            ! the occupation of the intermediates must match
            if (fl_pnt_i%interm%njoined.ne.fl_pnt_mark1%interm%njoined)
     &         cycle inner_loop
            if (.not.list_cmp(fl_pnt_i%interm%ihpvca_occ,
     &                    fl_pnt_mark1%interm%ihpvca_occ,
     &                        fl_pnt_i%interm%njoined*ngastp*2))
     &         cycle inner_loop
            if (.not.list_cmp(fl_pnt_i%interm%igasca_restr,
     &                    fl_pnt_mark1%interm%igasca_restr,
     &                    fl_pnt_mark1%interm%njoined*8*ngas*nspin)) 
     &         cycle inner_loop
            ! mark entry
            fl_pnt_mark3 => fl_pnt_i
            label_old = fl_pnt_mark3%interm%name
            fl_pnt_i => fl_pnt_i%next
            ! the next entry should be the corresponding contraction
            if (fl_pnt_i%command.ne.command_bc.and.
     &          fl_pnt_i%command.ne.command_bc_reo) goto 1001
            if (trim(fl_pnt_mark3%interm%name).ne.
     &                                 trim(fl_pnt_i%bcontr%label_res)) 
     &          goto 1001
            ! last check: the contractions must match
            if (fl_pnt_i%command.ne.fl_pnt_mark2%command)
     &         cycle inner_loop
            ! must check that reorderings are equal ... not yet done so:
            if (fl_pnt_i%command.eq.command_bc_reo)
     &         cycle inner_loop
            if (fl_pnt_i%bcontr%n_cnt.ne.fl_pnt_mark2%bcontr%n_cnt)
     &         cycle inner_loop
            if (abs(fl_pnt_i%bcontr%fact-fl_pnt_mark2%bcontr%fact)
     &                                                       .gt.1d-12)
     &         cycle inner_loop
            if (fl_pnt_i%bcontr%iblk_op1.ne.
     &                                    fl_pnt_mark2%bcontr%iblk_op1)
     &         cycle inner_loop
            if (fl_pnt_i%bcontr%iblk_op2.ne.
     &                                    fl_pnt_mark2%bcontr%iblk_op2)
     &         cycle inner_loop
            if (fl_pnt_i%bcontr%tra_op1.neqv.
     &                                     fl_pnt_mark2%bcontr%tra_op1)
     &         cycle inner_loop
            if (fl_pnt_i%bcontr%tra_op2.neqv.
     &                                     fl_pnt_mark2%bcontr%tra_op2)
     &         cycle inner_loop
            if (fl_pnt_i%bcontr%tra_res.neqv.
     &                                     fl_pnt_mark2%bcontr%tra_res)
     &         cycle inner_loop
          
            if (.not.list_cmp(fl_pnt_i%bcontr%occ_cnt,
     &                    fl_pnt_mark2%bcontr%occ_cnt,
     &                    fl_pnt_mark2%bcontr%n_cnt*ngastp*2 ))
     &         cycle inner_loop
            if (.not.list_cmp(fl_pnt_i%bcontr%rst_cnt,
     &                    fl_pnt_mark2%bcontr%rst_cnt,
     &                    fl_pnt_mark2%bcontr%n_cnt*8*ngas*nspin)) 
     &         cycle inner_loop
            ! match: put the two entries to list
            if (ntest.ge.100) then
              write(lulog,*) 'found matching intermediate'
              call print_form_item(lulog,idxfl,fl_pnt_mark3,op_info)
              call print_form_item(lulog,idxfl,fl_pnt_i,op_info)
            end if
            nlist = nlist+1
            fpl_marks_pnt%item => fl_pnt_mark3
            call new_formula_plist_entry(fpl_marks_pnt)
            fpl_marks_pnt => fpl_marks_pnt%next
            fpl_marks_pnt%item => fl_pnt_i
            call new_formula_plist_entry(fpl_marks_pnt)
            fpl_marks_pnt => fpl_marks_pnt%next

          else
            ! other passes:
            ! compare binary contraction, must have common result op.
c            call quit(1,i_am,'baustelle')

            if (fl_pnt_i%command.ne.command_bc.and.
     &        fl_pnt_i%command.ne.command_add_bc) !.and.
     &       ! fl_pnt_i%command.ne.command_bc_reo.and.
     &       ! fl_pnt_i%command.ne.command_add_bc_reo) 
     &      then
              cycle inner_loop     
            end if

            ! TO DO: if we allow REO here, check that these are the same!
            icmp = cmp_bcontr_spc(fl_pnt_mark1%bcontr,fl_pnt_i%bcontr)

            ! the above will only be 1 or 2 if 1 or 2 was not
            ! a short-time intermediate
            ! we must make sure, however, that the other operator
            ! (if it is a short-time intermediate) has the same
            ! occupation
            if (icmp.eq.1) then
              if (.not.list_cmp(fl_pnt_mark1%bcontr%occ_op2,
     &                          fl_pnt_i%bcontr%occ_op2,
     &                          2*ngastp*fl_pnt_i%bcontr%nj_op2).or.
     &            .not.list_cmp(fl_pnt_mark1%bcontr%rst_op2,
     &                          fl_pnt_i%bcontr%rst_op2,
     &                          8*fl_pnt_i%bcontr%ngas*
     &                            fl_pnt_i%bcontr%nspin*
     &                            fl_pnt_i%bcontr%nj_op2)) icmp = 0
            else if (icmp.eq.2) then
              if (.not.list_cmp(fl_pnt_mark1%bcontr%occ_op1,
     &                          fl_pnt_i%bcontr%occ_op1,
     &                          2*ngastp*fl_pnt_i%bcontr%nj_op1).or.
     &            .not.list_cmp(fl_pnt_mark1%bcontr%rst_op1,
     &                          fl_pnt_i%bcontr%rst_op1,
     &                          8*fl_pnt_i%bcontr%ngas*
     &                            fl_pnt_i%bcontr%nspin*
     &                            fl_pnt_i%bcontr%nj_op1)) icmp = 0
            end if
c dbg
c            if (fl_pnt_mark1%bcontr%nj_op1.eq.2.and.
c     &          trim(fl_pnt_mark1%bcontr%label_op1).eq.'_STIN0001') then
c              print *,'comparing'
c              call print_form_item(lulog,idxfl,fl_pnt_mark1,op_info)
c              call print_form_item(lulog,idxfl,fl_pnt_i,op_info)
c              print *,'icmp = ',icmp
c            end if
c dbg
            ! same second op: put to list_m1
            if (icmp.eq.2) then
              nlist = nlist+1
              fpl_marks_pnt%item => fl_pnt_i
              call new_formula_plist_entry(fpl_marks_pnt)
              fpl_marks_pnt => fpl_marks_pnt%next   
            end if
            ! same first op: put to list_m2
            if (icmp.eq.1) then
              nlist2 = nlist2+1
              fpl_marks2_pnt%item => fl_pnt_i
              call new_formula_plist_entry(fpl_marks2_pnt)
              fpl_marks2_pnt => fpl_marks2_pnt%next   
            end if
          end if

        end do inner_loop

        if (ntest.ge.100) write(lulog,*)'found matching entries: ',nlist
        if (ntest.ge.100.and.pass.gt.0)
     &                   write(lulog,*)'found also:             ',nlist2

        ! pass == 0
        if (pass.eq.0.and.nlist.gt.0) then
          size_intm = len_blk_est(fl_pnt_mark1%interm%ihpvca_occ,
     &                            fl_pnt_mark1%interm%njoined,
     &                            orb_info)
          ! TO BE FIXED FOR N_OPERAND==1
          cost_intm = len_blk_est(fl_pnt_mark2%bcontr%occ_ex1,
     &                            fl_pnt_mark2%bcontr%nj_op1,
     &                            orb_info)
     &              * len_blk_est(fl_pnt_mark2%bcontr%occ_ex2,
     &                            fl_pnt_mark2%bcontr%nj_op2,
     &                            orb_info)
     &              * len_blk_est(fl_pnt_mark2%bcontr%occ_cnt,
     &                            fl_pnt_mark2%bcontr%n_cnt,
     &                            orb_info)

          if (nlist.ge.100.or.cost_intm.gt.cost_ref) then
          ! if the list is non-empty
          ! rename first entry to _LTIN
          lti_cnt = lti_cnt+1
          if (lti_cnt.ge.10000) 
     &          call quit(1,i_am,'LTIN counter exceeded')
          write(label_new,'("_LTIN",i4.4)') lti_cnt

          if (ntest.ge.100) write(lulog,*)'new interm: ',trim(label_new)
c dbg
c          print *,trim(label_new),' -> ',nlist,'  size ',size_intm,
c     &          cost_intm
c dbg
          label_old = fl_pnt_mark1%interm%name
          ! allocate new
          call fl_item_rename_intm(fl_pnt_mark1,label_new,
     &                             size_intm.le.lblk_da,orb_info)
C          allocate(intm_new)
C          intm_new%name = label_new
C          call clone_operator(intm_new,fl_pnt_mark1%interm,.false.,
C     &                                                         orb_info)
C          call dealloc_operator(fl_pnt_mark1%interm)
C          deallocate(fl_pnt_mark1%interm)
C          fl_pnt_mark1%interm => intm_new
C          if (size_intm.le.lblk_da) fl_pnt_mark1%incore = 1 ! keep these incore
          ! look for next appearance of _STIN and rename
          fl_pnt_mark2%bcontr%label_res = label_new

          fpl_marks_pnt => fpl_marks
          do ilist = 0, nlist
c dbg
c            print *,'ilist = ',ilist
c dbg
            if (ilist.gt.0) then
              ! the first entry is the definition
              label_old = fpl_marks_pnt%item%interm%name  ! the _STIN name may differ 
              fl_pnt_i => fpl_marks_pnt%item%prev ! start search from here
c dbg
c              print *,'original start: '
c              call print_form_list(lulog,fl_pnt_i,op_info)
c dbg
              ! delete this entry from the formula list
c dbg
c              print *,'deleting 1'
c dbg
              call delete_fl_node(fpl_marks_pnt%item)
              deallocate(fpl_marks_pnt%item)
              fpl_marks_pnt => fpl_marks_pnt%next
              ! the second entry is the contraction
              ! remember the pointer to the formula item after that contr.
              !fl_pnt_i => fpl_marks_pnt%item%next
              ! and delete
c dbg
c              print *,'deleting 2'
c dbg
              call delete_fl_node(fpl_marks_pnt%item)
              deallocate(fpl_marks_pnt%item)
              fpl_marks_pnt => fpl_marks_pnt%next
            else 
              ! first time: just ...
              !fl_pnt_i => fl_pnt_mark2
              ! better use:
              fl_pnt_i => fl_pnt_mark2
            end if
            ! now look for the first (and only) occurrence of the _STIN
            ! in a contraction
c dbg
c              print *,'replacing ',trim(label_old),' by ',
c     &                             trim(label_new)
c              print *,'start from '
c              !call print_form_item(lulog,idxfl,fl_pnt_i,op_info)
c              call print_form_list(lulog,fl_pnt_i,op_info)
c dbg
            replace_loop: do
              fl_pnt_i => fl_pnt_i%next
c dbg
c              print *,'now at'
c              call print_form_item(lulog,idxfl,fl_pnt_i,op_info)
c dbg
              if (fl_pnt_i%command.eq.command_end_of_formula.or.
     &            fl_pnt_i%command.eq.command_set_target_init) 
     &              goto 1002
              ! possibly, the intermediate occurred in def of further
              ! intermediate, so check for that:
              if (fl_pnt_i%command.eq.command_new_intermediate) then
                if (trim(fl_pnt_i%parent1).eq.trim(label_old)) then
                  fl_pnt_i%parent1 = label_new
                end if
                if (trim(fl_pnt_i%parent2).eq.trim(label_old)) then
                  fl_pnt_i%parent2 = label_new
                end if
              end if
              if (fl_pnt_i%command.eq.command_add_bc.or.
     &            fl_pnt_i%command.eq.command_add_bc_reo.or.
     &            fl_pnt_i%command.eq.command_bc.or.
     &            fl_pnt_i%command.eq.command_bc_reo) then
                if (trim(fl_pnt_i%bcontr%label_op1)
     &                          .eq.trim(label_old)) then
                  fl_pnt_i%bcontr%label_op1 = label_new
                  fl_pnt_mark3 => fl_pnt_i
                  exit replace_loop
                end if
                if (trim(fl_pnt_i%bcontr%label_op2)
     &                          .eq.trim(label_old)) then
                  fl_pnt_i%bcontr%label_op2 = label_new
                  fl_pnt_mark3 => fl_pnt_i
                  exit replace_loop
                end if
              end if
            end do replace_loop
          end do ! nlist
          ! -----------------------------------------------
          ! now, the fl_pnt_mark3 should point to the last
          ! formula item in which the intermediate was used
          ! -----------------------------------------------

          ! insert a [DELETE INTERMEDIATE]
c dbg
c          print *,'now inserting a DEL for '//trim(label_new)
c dbg
          call insert_fl_item(fl_pnt_mark3,command_del_intermediate,-1)
          fl_pnt_i%next%label = label_new

          end if
        ! pass > 0
        else if (pass.gt.0.and.(nlist.gt.0.or.nlist2.gt.0)) then
          ! at this line, fl_pnt_mark1 points to the item
          ! where the first contraction is
          ! we will now define the new intermediate, with
          ! fl_pnt_mark2 pointing at the most recent entry
          ! fl_pnt_mark2 should by construction preceed fl_pnt_mark1
          !    (in the list) 

          ! if the list is non-empty
          ! process the longer list of list_m1, list_m2
c dbg
c          if (nlist.ge.nlist2) print *,'chose list for op1'
c          if (nlist.lt.nlist2) print *,'chose list for op2'
c dbg
          if (nlist.ge.nlist2) then
            idx_rpl = 1
            fpl_marks_pnt => fpl_marks
          else
            nlist = nlist2
            idx_rpl = 2
            fpl_marks_pnt => fpl_marks2
          end if

c dbg
c          print *,'idx_rpl = ',idx_rpl
c dbg
          ! define new intermediate
          lti_cnt = lti_cnt+1
          if (lti_cnt.gt.10000) 
     &          call quit(1,i_am,'LTIN counter exceeded')

          write(label_new,'("_LTIN",i4.4)') lti_cnt
          if (ntest.ge.100) write(lulog,*)'new interm: ',trim(label_new)
          
c dbg
          !watch = lti_cnt.eq.10
          if (watch) print *,'watching!'
c          print *,'intermediate def'
c dbg
          if (idx_rpl.eq.1) label_op = fl_pnt_mark1%bcontr%label_op1
          if (idx_rpl.eq.2) label_op = fl_pnt_mark1%bcontr%label_op2
          if (idx_rpl.eq.1) iblk_op = fl_pnt_mark1%bcontr%iblk_op1
          if (idx_rpl.eq.2) iblk_op = fl_pnt_mark1%bcontr%iblk_op2
          if (idx_rpl.eq.1) tra_op = fl_pnt_mark1%bcontr%tra_op1
          if (idx_rpl.eq.2) tra_op = fl_pnt_mark1%bcontr%tra_op2
          if (idx_rpl.eq.1) nj = fl_pnt_mark1%bcontr%nj_op1
          if (idx_rpl.eq.2) nj = fl_pnt_mark1%bcontr%nj_op2
          if (idx_rpl.eq.1) occ => fl_pnt_mark1%bcontr%occ_op1
          if (idx_rpl.eq.2) occ => fl_pnt_mark1%bcontr%occ_op2
          if (idx_rpl.eq.1) rst => fl_pnt_mark1%bcontr%rst_op1
          if (idx_rpl.eq.2) rst => fl_pnt_mark1%bcontr%rst_op2
          skip = .false.
          if (label_op(1:5).eq.'_STIN') then
            ! just look for definition and replace
c dbg
            if (watch) print *,'  case STIN first --- '//trim(label_op)
            if (watch) print *,'  look backward from: '
            if (watch) 
     &         call print_form_item(lulog,idxfl,fl_pnt_mark1,op_info)
            if (watch) 
     &      call print_form_item(lulog,idxfl,fl_pnt_mark1%prev,op_info)
c dbg
            fl_pnt_mark2 => find_fl_item(fl_pnt_mark1,
     &                                command=command_new_intermediate,
     &                                label=label_op,
     &                                backward=.true.)
            if (.not.associated(fl_pnt_mark2)) 
     &            call quit(1,i_am,'problem')
c dbg
            fl_pnt_watch => fl_pnt_mark2
c dbg
            ! rename here
            label_old = fl_pnt_mark2%interm%name
            ! if one of the parents is a _STIN, we should use a different
            ! operator as parent
            if (fl_pnt_mark2%parent1(1:5).eq.'_STIN'.or.
     &          fl_pnt_mark2%parent2(1:5).eq.'_STIN') then
c dbg
c              print *,'need better parents: '
c dbg
              ! loop over the other terms in the list and take a
              ! better choice
              fpl_marks2_pnt => fpl_marks_pnt
              do
                if (idx_rpl.eq.1) then
                  label_op = fpl_marks2_pnt%item%bcontr%label_op1
                  tra_op = fpl_marks2_pnt%item%bcontr%tra_op1
                else if (idx_rpl.eq.2) then
                  label_op = fpl_marks2_pnt%item%bcontr%label_op2
                  tra_op = fpl_marks2_pnt%item%bcontr%tra_op2
                end if
                ! a real operator? great! take it!
                if (label_op(1:5).ne.'_STIN') then
                  fl_pnt_mark2%parent1 = label_op
                  fl_pnt_mark2%parent2 = '---'
                  fl_pnt_mark2%tra1 = tra_op
                  fl_pnt_mark2%tra2 = .false.
                  exit
                else
c dbg
c                  print *,'resolving another STIN ...'
c dbg
                  ! another STIN ??? find its definition
                  fl_pnt_mark3 => find_fl_item(fpl_marks2_pnt%item,
     &                                command=command_new_intermediate,
     &                                label=label_op,
     &                                backward=.true.)
                  if (.not.associated(fl_pnt_mark3))
     &                    call quit(1,i_am,'problem (2)')
                  ! hopefully this one is defined w/o _STINs
                  if (fl_pnt_mark3%parent1(1:5).ne.'_STIN'.and.
     &                fl_pnt_mark3%parent2(1:5).ne.'_STIN') then
                    fl_pnt_mark2%parent1 = fl_pnt_mark3%parent1 
                    fl_pnt_mark2%parent2 = fl_pnt_mark3%parent2 
                    fl_pnt_mark2%tra1 = fl_pnt_mark3%tra1
                    fl_pnt_mark2%tra2 = fl_pnt_mark3%tra2
                    exit
                  end if
                end if
                ! the list has always one more entry ...
                if (.not.associated(fpl_marks2_pnt%next%next)) then
c                  call quit(1,i_am,'damn! trapped in _STINs ...')
                   skip = .true.
                   lti_cnt = lti_cnt-1 ! decrease counter again
                   exit
                end if
                fpl_marks2_pnt => fpl_marks2_pnt%next
              end do
            end if
            if (.not.skip) then
c dbg
c              print *,'   renaming'
c dbg
              call fl_item_rename_intm(fl_pnt_mark2,label_new,.false.,
     &                               orb_info)
            else
              nlist = -1
c dbg
c              print *,'giving up ...'
c dbg
            end if


          else 
c dbg
            if (watch) print *,'  case OP first --- '//trim(label_op)
c dbg
            ! insert before first occurrence:
            fl_pnt_mark2 => fl_pnt_mark1%prev
            call insert_fl_item(fl_pnt_mark2,
     &                          command_new_intermediate,-1)
            ! store info
            fl_pnt_mark2 => fl_pnt_mark2%next
            call store_def_intm(fl_pnt_mark2,
     &            label_new,occ,rst,nj,1,
     &            label_op,'---',.false.,tra_op,.false.,
     &            orb_info)
c dbg
            fl_pnt_watch => fl_pnt_mark2
c dbg

          end if
          ! define as sum of second/first operators
c dbg
c          if (nlist.gt.-1) print *,'now processing list of terms'
c dbg
          list_loop: do ilist = 0, nlist
c dbg
            if (watch) print *,'  ilist = ',ilist,'/',nlist
c dbg
            if (ilist.eq.0) then
              fl_pnt_mark3 => fl_pnt_mark1
            else
              fl_pnt_mark3 => fpl_marks_pnt%item
              fpl_marks_pnt => fpl_marks_pnt%next
            end if
            if (idx_rpl.eq.1) label_op = fl_pnt_mark3%bcontr%label_op1
            if (idx_rpl.eq.2) label_op = fl_pnt_mark3%bcontr%label_op2
            if (idx_rpl.eq.1) iblk_op = fl_pnt_mark3%bcontr%iblk_op1
            if (idx_rpl.eq.2) iblk_op = fl_pnt_mark3%bcontr%iblk_op2
            if (idx_rpl.eq.1) tra_op = fl_pnt_mark3%bcontr%tra_op1
            if (idx_rpl.eq.2) tra_op = fl_pnt_mark3%bcontr%tra_op2
            if (idx_rpl.eq.1) nj = fl_pnt_mark3%bcontr%nj_op1
            if (idx_rpl.eq.2) nj = fl_pnt_mark3%bcontr%nj_op2
            if (idx_rpl.eq.1) occ => fl_pnt_mark3%bcontr%occ_op1
            if (idx_rpl.eq.2) occ => fl_pnt_mark3%bcontr%occ_op2
            if (idx_rpl.eq.1) rst => fl_pnt_mark3%bcontr%rst_op1
            if (idx_rpl.eq.2) rst => fl_pnt_mark3%bcontr%rst_op2

            fact = fl_pnt_mark3%bcontr%fact
            fact_itf = fl_pnt_mark3%bcontr%fact_itf
c dbg
c            print *,'  fact = ',fact
c dbg
            ! operator: add to intermediate
            if (label_op(1:1).ne.'_') then
c dbg
              if (watch) print *,'  case: OP --- '//trim(label_op)
c dbg
              if (ilist.eq.0) then
                call insert_fl_item(fl_pnt_mark2,
     &                          command_cp_intm,-1)
              else
                call insert_fl_item(fl_pnt_mark2,
     &                          command_add_intm,-1)
              end if
c dbg
              if (watch) print *,'item inserted'
c dbg
              ! dummy map for itf
              allocate(svdummy(3))
              svdummy(1:3) = 0
              fl_pnt_mark2 => fl_pnt_mark2%next
              call store_bc(fl_pnt_mark2,
     &                  fact,fact_itf,
     &                  label_new,label_op,'---',
     &                  1,iblk_op,idummy,
     &                  .false.,tra_op,ldummy,
     &                  nj,nj,0,
     &                  occ,occ,idummy,
     &                  rst,rst,idummy,
     &                  idummy,idummy,idummy,
     &                  idummy,idummy,idummy,0,
     &                  idummy,idummy,
     &                  idummy,idummy,
     &                  svdummy,
     &                  orb_info)

              deallocate(svdummy)
c dbg
              if (watch) print *,' ... and written'
c dbg
            else if (label_op(1:2).eq.'_S') then
c dbg
              if (watch) print *,'  case: _STIN'
c
c dbg
              if (ilist.gt.0) then
                ! _STI: find definition and remove
                fl_pnt_mark4 => find_fl_item(fl_pnt_mark3,
     &                                command=command_new_intermediate,
     &                                label=label_op,
     &                                backward=.true. )
                if (.not.associated(fl_pnt_mark4))
     &             call quit(1,i_am,'mistmist')
                ! delete the definition
                call delete_fl_node(fl_pnt_mark4)
                deallocate(fl_pnt_mark4)
              end if
              ! find the contraction and relink
              fl_pnt_mark4 => find_fl_item(fl_pnt_mark3,
     &                     command_list=(/command_bc,command_bc_reo,
     &                                    command_reorder/),
     &                                nlist = 3,
     &                                label_res=label_op,
     &                                backward=.true. )
              if (.not.associated(fl_pnt_mark4)) then
                 write(lulog,*) 'Did not find defining [CONTR] for '//
     &                       trim(label_op)
                 write(lulog,*) 'I started from ... (search backwards)'
                 call print_form_item(lulog,idxfl,fl_pnt_mark3,op_info)
                 write(lulog,*) 'Preceeding item ...'
                 call print_form_item(lulog,idxfl,
     &                        fl_pnt_mark3%prev,op_info)
                 call quit(1,i_am,'mistmistmist')
              end if

              ! insert contraction at mark 4 behind mark 2
              ! take care to  move all required terms  
              call optk_move_bcontr(fl_pnt_mark4,fl_pnt_mark2)
              ! the result must be renamed and the fact must be adapted
              fl_pnt_mark4%bcontr%fact = fl_pnt_mark4%bcontr%fact*fact
              fl_pnt_mark4%bcontr%fact_itf 
     %                           = fl_pnt_mark4%bcontr%fact_itf*fact_itf
              fl_pnt_mark4%bcontr%label_res = label_new
              fl_pnt_mark4%bcontr%iblk_res = 1
              ! if not the first term: contract and add:
              if (ilist.gt.0.and.fl_pnt_mark4%command.eq.command_bc)
     &            fl_pnt_mark4%command = command_add_bc
              if (ilist.gt.0.and.fl_pnt_mark4%command.eq.command_bc_reo)
     &            fl_pnt_mark4%command = command_add_bc_reo
              ! if this is a REO rename here, too
              if (fl_pnt_mark4%command.eq.command_bc_reo.or.
     &            fl_pnt_mark4%command.eq.command_add_bc_reo) then
                fl_pnt_mark4%reo%label_in  = label_new
                fl_pnt_mark4%reo%label_out = label_new
              end if

              ! shift mark2 pointer such that further items are appended 
              ! after this contraction
c dbg
c             print *,'after inserting: formula starting from old mark 2'
c             call print_form_list(lulog,fl_pnt_mark2,op_info)
c             print *,'--- finito ---'
c dbg
              fl_pnt_mark2 => fl_pnt_mark4
            else if (label_op(1:2).eq.'_L') then
              ! _LTI: find definition and add, make sure that defined before (relink)
               call warn(i_am,'avoiding LTINs presently ...')
               call quit(1,i_am,'There is a bug here ... try '//
     &                          '"routes auto_opt=F"')
               cycle list_loop
            else
               call quit(1,i_am,'illegal operator name: '//
     &                   trim(label_op))
            end if
            ! delete the node
c dbg
            if (watch) then
              print *,'Formula starting from beg. of new inter:'
              call print_form_list(lulog,fl_pnt_watch,op_info)
              print *,'--- finito ---'
            end if
c dbg
            ! the item for ilist==0 is needed below for defining
            ! the new contraction with the present intermediate ...!!
            if (ilist.gt.0) then
c dbg
            if (watch) print *,'  deleting '//
     &                 trim(fl_pnt_mark3%bcontr%label_op1)//'*'//
     &                 trim(fl_pnt_mark3%bcontr%label_op2)
c dbg
              call delete_fl_node(fl_pnt_mark3)
              deallocate(fl_pnt_mark3)
            end if
          end do list_loop

          if (.not.skip) then
c dbg
c          print *,'defining new bin. c. '
c dbg
            ! rename in binary contr.
            if (idx_rpl.eq.1)
     &          fl_pnt_mark1%bcontr%label_op1 = label_new
            if (idx_rpl.eq.2)
     &          fl_pnt_mark1%bcontr%label_op2 = label_new
            if (idx_rpl.eq.1)
     &          fl_pnt_mark1%bcontr%iblk_op1 = 1
            if (idx_rpl.eq.2)
     &          fl_pnt_mark1%bcontr%iblk_op2 = 1
            ! the intermediate is never transposed;
            ! future improvement: it might sometimes be better to 
            ! define the transpose of the intermediate, instead.
            if (idx_rpl.eq.1)
     &          fl_pnt_mark1%bcontr%tra_op1 = .false.
            if (idx_rpl.eq.2)
     &          fl_pnt_mark1%bcontr%tra_op2 = .false.
            ! the factors were absorbed into the intermediate
            fl_pnt_mark1%bcontr%fact = 1d0
            fl_pnt_mark1%bcontr%fact_itf = 1d0
            ! delete intm. thereafter
            call insert_fl_item(fl_pnt_mark1,
     &                          command_del_intermediate,-1)
            fl_pnt_mark1%next%label=label_new
            success = .true.
         end if
c dbg
c         print *,'processing done for this term ...'
c dbg

        end if

        call dealloc_formula_plist(fpl_marks)
        if (pass.gt.0) call dealloc_formula_plist(fpl_marks2)

        fl_pnt_o => fl_pnt_o%next

      end do outer_loop

      if (ntest.gt.0) write(lulog,*) 'end of '//trim(i_am)

      return
      
 1001 write(lulog,*) 'I am still assuming that '//
     &          '[NEW INTERMEDIATE] is followed by its definition'
      write(lulog,*) 'Expected result: '//trim(label_old)
      write(lulog,*) 'But I found: '
      call print_form_item(lulog,idxfl,fl_pnt_o,op_info)
      write(lulog,*) 'the index is possibly wrong, ignore that'
      call quit(1,i_am,'assumption not true? buggy formula?') 
 1002 write(lulog,*) '[NEW TARGET] or [END] occured while '//
     &     'trying to replace '//trim(label_old)
      call quit(1,i_am,'buggy formula?')
 1003 write(lulog,*) 'formula ended unexpectedly (without [END])'
      write(lulog,*) 'last item was:'
      call print_form_item(lulog,idxfl,fl_pnt_i,op_info)
      call quit(1,i_am,'buggy formula?')
      end

