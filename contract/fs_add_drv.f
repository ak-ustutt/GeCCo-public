*----------------------------------------------------------------------*
      subroutine fs_add_drv(xret_blk,type_xret,idx_tgt,me_tgt,
     &         fl_item,update,
     &         op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     driver for [ADD] operation in frm_schedX() (X>=2)
*     ... and [COPY] as well ...
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'par_opnames_gen.h'

      real(8), intent(inout), target ::
     &     xret_blk(*)
      logical, intent(in) ::
     &     update
      integer, intent(in) ::
     &     type_xret, idx_tgt
      type(me_list), intent(inout), target ::
     &     me_tgt
      type(formula_item), intent(in) ::
     &     fl_item

      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(inout) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      character(maxlen_bc_label) ::
     &     label_res, label_op
      integer ::
     &     idx_res, idx_op, iblk_res, iblk_op, idxme, type_xret_loc
      real(8) ::
     &     fact
      logical ::
     &     tra_res, tra_op, copy_only
      real(8), target ::
     &     xret_dummy(1)

      integer, pointer ::
     &     op2list(:)
      type(me_list_array), pointer ::
     &     mel_arr(:)
      type(me_list), pointer ::
     &     me_op, me_res
      type(binary_contr), pointer ::
     &     add_info
      real(8), pointer ::
     &     xret_pnt(:)

      integer, external ::
     &     idx_oplist2


      mel_arr => op_info%mel_arr
      op2list => op_info%op2list

      add_info => fl_item%bcontr

      copy_only = fl_item%command.eq.command_cp_intm

      if (add_info%n_operands.ne.1.or.add_info%n_cnt.ne.0) then
        write(lulog,*) '[ADD]: something is wrong:'
        call prt_bcontr(lulog,add_info)
        call quit(1,'fs_add_drv','something is wrong')
      end if

      label_res = add_info%label_res
      label_op  = add_info%label_op1

      idx_res = idx_oplist2(label_res,op_info)
      if (update.and.idx_res.ne.idx_tgt) then
        write(lulog,*) 'formula target: ',
     &       trim(op_info%op_arr(idx_tgt)%op%name)
        write(lulog,*) '[ADD] result: ',trim(label_res)
        call quit(1,'fs_add_drv','inconsistency: target<>result')
      end if
      idx_op = idx_oplist2(label_op,op_info)

      iblk_res = add_info%iblk_res
      iblk_op  = add_info%iblk_op1
      tra_res = add_info%tra_res
      tra_op  = add_info%tra_op1

      if (update) then
        me_res => me_tgt
        type_xret_loc = type_xret
        xret_pnt => xret_blk(iblk_res:iblk_res)
      else
        idxme = op2list(idx_res)
        me_res => mel_arr(idxme)%mel
        if (me_res%fhand%unit.lt.0)
     &       call file_open(me_res%fhand)
        type_xret_loc = 0
        xret_pnt => xret_dummy
      end if

      fact = add_info%fact

      ! special: unit operator
      if (label_op.eq.op_unity) then
        call quit(1,'fs_add_drv',
     &       'add_unity still misses update of xret')
        call add_unity(fact,1,me_res,iblk_res,orb_info,str_info)
      else
        idxme = op2list(idx_op)
        me_op => mel_arr(idxme)%mel 
        if (me_op%fhand%unit.le.0)
     &       call file_open(me_op%fhand)
c?        njoined = add_info%nj_op1
c?        ! iblkop fix:
c?        iblkop  = (iblkop-1)/njoined + 1
        if (.not.tra_op.and.     tra_res .or.
     &           tra_op.and..not.tra_res) then
          call add_opblk_transp(xret_pnt,type_xret_loc,fact,
     &             me_op,me_res,tra_op,tra_res,
     &             iblk_op,iblk_res,
     &             op_info,str_info,orb_info,copy_only)
        else
          call add_opblk(xret_pnt,type_xret_loc,fact,
     &             me_op,me_res,
     &             iblk_op,iblk_res,orb_info,copy_only)
        end if

      end if

      return
      end
