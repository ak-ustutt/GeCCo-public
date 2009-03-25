*----------------------------------------------------------------------*
      subroutine sum_hermite(fl_tgt,sum,strict,op_info)
*----------------------------------------------------------------------*
*     find terms that are hermitian conjugates and:
*     * sum==.false.:
*     remove terms that are hermitian conjugates of previous ones
*     applies a factor 1/2 to self-adjoined terms and put a
*     [SYMMETRIZE] statement at the end of the formula list
*     * sum==.true.:
*     sum up the hermitian pair (self-adjoint: do nothing)
*     for scalar results only!
*     * strict==.true.: we quit if no adjoint was found
*
*     NOTE: use sum_terms before calling this routine, as we only
*           look for a single adjoint term
*           If used with sum==.true. and strict==.false.,
*           it will reduce the number of terms, also in the presence
*           of non-hermitian terms
*           If used with sum==.false. and strict==.false.,
*           it will effectively symmetrize the expression.
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_array.h'
      include 'def_formula_item_list.h'
      
      type(formula_item), target, intent(inout) ::
     &     fl_tgt      
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     sum, strict

      logical ::
     &     diag, found
      integer ::
     &     idxop_tgt, iblk_tgt, iblk_tgt_trp, iterm, jterm,
     &     njoined_tgt

      type(contraction) ::
     &     contr_scr

      type(formula_item), pointer ::
     &     fl_tgt_pnt, fl_tgt_pnt_next, fl_tgt_current
      type(operator), pointer ::
     &     op_tgt
      integer, pointer ::
     &     occ_tgt(:,:,:)

      logical, external ::
     &     cmp_contr
      integer, external ::
     &     vtx_type, iblk_occ

      if (ntest.ge.100) then
        write(luout,*) '=========================='
        write(luout,*) ' info from sum_hermite'
        write(luout,*) '=========================='
      end if

      if (.not.sum) ! but it likely works
     &     call quit(1,'sum_hermite',
     &                 '.not.sum -- have not checked this option')

      if (fl_tgt%command.ne.command_set_target_init) then
        call quit(1,'sum_hermite',
     &         'must start with [INIT]')
      end if

      iterm = 0
      fl_tgt_current => fl_tgt
      ! loop over target items
      tgt_loop: do

        ! new operator target ?
        if (fl_tgt_current%command.eq.command_set_target_init) then
          idxop_tgt = fl_tgt_current%target
          if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'sum_hermite',
     &         'unexpected end of list (target)')
          if (ntest.ge.100) then
            write(luout,'(70("="))')
            write(luout,*) 'New operator target: ',idxop_tgt
            write(luout,'(70("="))')
          end if
          fl_tgt_current => fl_tgt_current%next
          
          ! make sure that target is either a scalar or
          ! has the "hermitian" attribute
          op_tgt => op_info%op_arr(idxop_tgt)%op

          if (.not.op_tgt%hermitian.eq.1.and.
     &        .not.vtx_type(op_tgt).eq.0) then
            ! anti-hermitian is not yet checked
            call quit(1,'sum_hermite',
     &           'operator is not declared Hermitian: "'//
     &           trim(op_tgt%name)//'"')
          end if

          occ_tgt => op_tgt%ihpvca_occ
          njoined_tgt = op_tgt%njoined

        end if

        iblk_tgt = fl_tgt_current%contr%iblk_res
        iblk_tgt_trp =
     &       iblk_occ(occ_tgt(1:,1:,(iblk_tgt-1)*njoined_tgt+1),
     &                          .true.,op_tgt)
        ! is this a diagonal block?
        diag = iblk_tgt .eq. iblk_tgt_trp

        iterm = iterm+1
        jterm = iterm+1
        if (ntest.ge.100) then
          write(luout,*) 'current term: # ',iterm
          call prt_contr2(luout,fl_tgt_current%contr,op_info)
        end if

        found = .false.

        ! if diagonal, strict for self-adjoint first
        if (diag) then
          call init_contr(contr_scr)
          call copy_contr(fl_tgt_current%contr,contr_scr)
          call transpose_contr(contr_scr,op_info)
          
          if (cmp_contr(contr_scr,
     &         fl_tgt_current%contr,.true.)) then
            if (.not.sum)
     &           fl_tgt_current%contr%fac =
     &           fl_tgt_current%contr%fac*0.5d0
            found = .true.
          end if
        end if

        if (.not.found) then
          if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'sum_hermite',
     &         'unexpected end of list (target)')
          fl_tgt_pnt => fl_tgt_current%next
          search_loop: do
            ! we search within a single domain only
            if (fl_tgt_pnt%command.ne.command_add_contribution)
     &         exit search_loop

            ! as the node might get deleted, better save the next pointer
            if (.not.associated(fl_tgt_pnt%next))
     &         call quit(1,'sum_hermite',
     &         'unexpected end of list (target, inner loop)')
            fl_tgt_pnt_next => fl_tgt_pnt%next
          
            if (fl_tgt_pnt%contr%iblk_res.eq.iblk_tgt_trp) then
              if (fl_tgt_pnt%contr%idx_res.ne.idxop_tgt)
     &           call quit(1,'sum_hermite',
     &           'suspicious change of operator target')

              call init_contr(contr_scr)
              call copy_contr(fl_tgt_pnt%contr,contr_scr)
              call transpose_contr(contr_scr,op_info)

              if (cmp_contr(contr_scr,
     &                    fl_tgt_current%contr,.true.)) then
                found = .true.
                if (ntest.ge.100) then
                  write(luout,*) 'found equal term: # ',jterm
                  call prt_contr2(luout,fl_tgt_pnt%contr,
     &                 op_info)
                  write(luout,*) 'now deleting'
                end if
                if (sum)
     &               fl_tgt_current%contr%fac =
     &               fl_tgt_current%contr%fac + contr_scr%fac
                call delete_fl_node(fl_tgt_pnt)
                deallocate(fl_tgt_pnt)
                exit search_loop ! we look only for one (see header)
              else if (ntest.ge.1000) then              
                write(luout,*) 'non-equiv term: # ',jterm
                call prt_contr2(luout,fl_tgt_pnt%contr,
     &               op_info)
                write(luout,*) 'transposed:'
                call prt_contr2(luout,contr_scr,
     &               op_info)
              end if
            end if

            fl_tgt_pnt => fl_tgt_pnt_next
            jterm = jterm+1
          end do search_loop

        end if

        if (strict.and..not.found) then
          write(luout,*) 'no adjungate found for term:'
          call prt_contr2(luout,fl_tgt_current%contr,op_info)
          call quit(1,'sum_hermite','check failed')
        end if
        ! advance to next item
        ! end of formula list?
        if (.not.associated(fl_tgt_current%next)) exit tgt_loop
        ! go to next item
        fl_tgt_current => fl_tgt_current%next
        if (fl_tgt_current%command.eq.command_end_of_formula)
     &                                             exit tgt_loop

      end do tgt_loop
        
      return
      end
