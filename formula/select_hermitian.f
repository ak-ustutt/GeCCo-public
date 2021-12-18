*----------------------------------------------------------------------*
      subroutine select_hermitian(fl_tgt,op_info)
*----------------------------------------------------------------------*
*     delete all but self-adjoint terms
*
*     adapted from sum_hermite by matthias, fall 2009
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
      include 'ifc_formula.h'
      
      type(formula_item), target, intent(inout) ::
     &     fl_tgt      
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     diag, found
      integer ::
     &     idxop_tgt, iblk_tgt, iblk_tgt_trp, iterm,
     &     njoined_tgt

      type(contraction) ::
     &     contr_scr

      type(formula_item), pointer ::
     &     fl_tgt_pnt, fl_tgt_pnt_next, fl_tgt_current, fl_tgt_next
      type(operator), pointer ::
     &     op_tgt
      integer, pointer ::
     &     occ_tgt(:,:,:)

      logical, external ::
     &     cmp_contr
      integer, external ::
     &     vtx_type, iblk_occ

      if (ntest.ge.100) then
        write(lulog,*) '=========================='
        write(lulog,*) ' info from select_hermitian'
        write(lulog,*) '=========================='
      end if

      if (fl_tgt%command.ne.command_set_target_init) then
        call quit(1,'select_hermitian',
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
     &         call quit(1,'select_hermitian',
     &         'unexpected end of list (target)')
          if (ntest.ge.100) then
            write(lulog,'(70("="))')
            write(lulog,*) 'New operator target: ',idxop_tgt
            write(lulog,'(70("="))')
          end if
          fl_tgt_current => fl_tgt_current%next
          
          ! make sure that target is either a scalar or
          ! has the "hermitian" attribute
          op_tgt => op_info%op_arr(idxop_tgt)%op

          if (.not.op_tgt%hermitian.eq.1.and.
     &        .not.vtx_type(op_tgt).eq.0) then
            ! anti-hermitian is not yet checked
            call quit(1,'select_hermitian',
     &           'operator is not declared Hermitian: "'//
     &           trim(op_tgt%name)//'"')
          end if

          occ_tgt => op_tgt%ihpvca_occ
          njoined_tgt = op_tgt%njoined

        end if

        fl_tgt_next => fl_tgt_current%next

        iblk_tgt = fl_tgt_current%contr%iblk_res
        iblk_tgt_trp =
     &       iblk_occ(occ_tgt(1:,1:,(iblk_tgt-1)*njoined_tgt+1),
     &                          .true.,op_tgt,
     &                op_tgt%blk_version(iblk_tgt))
        ! is this a diagonal block?
        diag = iblk_tgt .eq. iblk_tgt_trp

        iterm = iterm+1
        if (ntest.ge.100) then
          write(lulog,*) 'current term: # ',iterm
          call prt_contr2(lulog,fl_tgt_current%contr,op_info)
        end if

        found = .false.

        ! if diagonal, strict for self-adjoint first
        if (diag) then
          call init_contr(contr_scr)
          call copy_contr(fl_tgt_current%contr,contr_scr)
          call transpose_contr(contr_scr,op_info)
          found = cmp_contr(contr_scr,fl_tgt_current%contr,.true.)
        end if

        ! delete the node if not self-adjoint
        if (.not.found) call delete_fl_node(fl_tgt_current)

        ! advance to next item
        ! end of formula list?
        if (.not.associated(fl_tgt_next)) exit tgt_loop
        ! go to next item
        fl_tgt_current => fl_tgt_next
        if (fl_tgt_current%command.eq.command_end_of_formula)
     &                                             exit tgt_loop

      end do tgt_loop
        
      return
      end
