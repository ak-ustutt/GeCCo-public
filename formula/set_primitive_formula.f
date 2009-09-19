*----------------------------------------------------------------------*
      subroutine set_primitive_formula(form,idx_op,fac,idx_opsh,init,
     &     op_info)
*----------------------------------------------------------------------*
*     store blocks of operator in formula format, corresponds to
*       F = Op(block1) + Op(block2) + ....
*     we may select blocks by giving a shape operator which must 
*     contain a subset of blocks of op
*     init == .true. :  set [INIT]
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'par_opnames_gen.h'

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     init
      type(formula_item), intent(inout), target ::
     &     form
      integer, intent(in) ::
     &     idx_op, idx_opsh
      real(8), intent(in) ::
     &     fac
      type(operator_info), intent(in), target ::
     &     op_info

      logical ::
     &     blk_select, unit_exception
      type(operator), pointer ::
     &     op, opsh
      type(formula_item), pointer ::
     &     form_pnt
      integer ::
     &     iblk_op, iblk_opsh, iblk_off, iblk_off_sh, ijoin, njoined

      integer, external ::
     &     iblk_occ
      logical, external ::
     &     occ_is_diag_blk
      
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_primitive_formula')        
      end if

      blk_select = idx_op.ne.idx_opsh

      op   => op_info%op_arr(idx_op)%op
      unit_exception = trim(op%name).eq.trim(op_unity)
      opsh => op_info%op_arr(idx_opsh)%op

      if (ntest.ge.100) then
        write(luout,*) 'op:   ',trim(op%name)
        call print_op_occ(luout,op)
        write(luout,*) 'opsh: ',trim(opsh%name)
        call print_op_occ(luout,opsh)
      end if

      if (opsh%name.eq.op_unity)
     &     call quit(1,'set_primitive_formula',
     &     'Unit operator cannot be the shape defining operator!')
      form_pnt => form

      if (form_pnt%command.ne.command_end_of_formula)
     &     call quit(1,'set_primitive_formula',
     &     'formula seems positioned incorrectly on input')
      if (associated(form_pnt%contr))
     &     call quit(1,'set_primitive_formula',
     &     'contr is already associated on entry')

      if (init) then
        call new_formula_item(form,command_set_target_init,idx_opsh)
        form_pnt => form_pnt%next
      end if

      njoined = opsh%njoined
      do iblk_opsh = 1, opsh%n_occ_cls
        iblk_off_sh = (iblk_opsh-1)*njoined
        if (unit_exception) then
          if (.not.occ_is_diag_blk(
     &         opsh%ihpvca_occ(1,1,iblk_off_sh+1),njoined))
     &         cycle
          iblk_op = 1
        else if (blk_select) then          
c dbg
c          print *,'iblk_opsh = ',iblk_opsh
c          call wrt_occ_n(6,opsh%ihpvca_occ(1,1,iblk_off_sh+1),njoined)
c          print *,'opsh%dagger = ',opsh%dagger
c          print *,'op%dagger = ',opsh%dagger
c dbg
          iblk_op = iblk_occ(opsh%ihpvca_occ(1,1,iblk_off_sh+1),
     &                       opsh%dagger,op,
     &                       opsh%blk_version(iblk_opsh))
c dbg
c          print *,'result: iblk_op = ',iblk_op
c          print *,'op blks: ',op%blk_version(:)
c          print *,'version: ',opsh%blk_version(iblk_opsh)
c dbg
          if (iblk_op .lt. 1) cycle
        else
          iblk_op = iblk_opsh
        end if
        iblk_off = (iblk_op-1)*njoined

        ! new entry
        call new_formula_item
     &       (form_pnt,command_add_contribution,idx_opsh)

        ! set "contraction"
        call resize_contr(form_pnt%contr,njoined,0,njoined,0)
        form_pnt%contr%fac = fac
        form_pnt%contr%idx_res = idx_opsh
c??        form_pnt%contr%iblk_res = iblk_off_sh+1
        form_pnt%contr%iblk_res = iblk_opsh
        form_pnt%contr%nvtx = njoined
        form_pnt%contr%svertex(1:njoined) = 1
        form_pnt%contr%nxarc = njoined
        form_pnt%contr%vertex(1:njoined)%idx_op = idx_op
        do ijoin = 1, njoined
          if (.not.unit_exception) then
            form_pnt%contr%vertex(ijoin)%iblk_op = iblk_off+ijoin
          else
            form_pnt%contr%vertex(ijoin)%iblk_op = iblk_off+1
          end if
          form_pnt%contr%xarc(ijoin)%link(1) = ijoin
          form_pnt%contr%xarc(ijoin)%link(2) = ijoin
          if (.not.unit_exception) then
            if (.not.op%dagger) then
              form_pnt%contr%xarc(ijoin)%occ_cnt =
     &           op%ihpvca_occ(1:ngastp,1:2,iblk_off+ijoin)
            else
              form_pnt%contr%xarc(ijoin)%occ_cnt(1:ngastp,1) =
     &             op%ihpvca_occ(1:ngastp,2,iblk_off_sh+njoined+1-ijoin)
              form_pnt%contr%xarc(ijoin)%occ_cnt(1:ngastp,2) =
     &             op%ihpvca_occ(1:ngastp,1,iblk_off_sh+njoined+1-ijoin)
            end if
          else 
            form_pnt%contr%xarc(ijoin)%occ_cnt =
     &         opsh%ihpvca_occ(1:ngastp,1:2,iblk_off_sh+ijoin)
          end if
        end do

        ! update svertex info
        call update_svtx4contr(form_pnt%contr)
        
        form_pnt => form_pnt%next
      end do

      if (ntest.ge.100) then
        write(luout,*) 'result:'
        call print_form_list(luout,form,op_info)
      end if

      return
      end
