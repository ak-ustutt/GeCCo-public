*----------------------------------------------------------------------*
      subroutine fs_reo_drv(xret_blk,type_xret,idx_tgt,me_tgt,
     &         fl_item,update,add,
     &         op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     driver for [REORDER] and [REORDER][ADD] of list
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
      include 'def_reorder_info.h'
      include 'def_formula_item.h'
      include 'par_opnames_gen.h'

      integer, parameter ::
     &     ntest = 00

      real(8), intent(inout), target ::
     &     xret_blk(*)
      integer, intent(in) ::
     &     type_xret, idx_tgt
      logical, intent(in) ::
     &     update, add
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

      integer ::
     &     ngas, nspin, 
     &     n_cnt, idx_opin, idx_opout,
     &     iblk_opin, iblk_opout, iblk_tmp,
     &     nj_opin, nj_opout, nj_tmp, idxme,
     &     type_xret_loc,
     &     idoff_opin, idoff_opout
      logical ::
     &     tra_opin, tra_opout
      type(reorder_info) ::
     &     reo_info
      real(8), target ::
     &     xret_dummy(1)
      real(8) ::
     &     fact

      type(operator), pointer ::
     &     op_tmp
      type(me_list), pointer ::
     &     me_tmp

      integer, pointer ::
     &     op2list(:)
      type(me_list_array), pointer ::
     &     mel_arr(:)
      type(me_list), pointer ::
     &     me_opin, me_opout
      type(reorder), pointer ::
     &     reo_inf0
      type(filinf), pointer ::
     &     ffop
      integer, pointer ::
     &     iocc_opin(:,:,:), iocc_opout(:,:,:),
     &     irst_opin(:,:,:,:,:,:),irst_opout(:,:,:,:,:,:)
      real(8), pointer ::
     &     xret_pnt(:)

      integer, external ::
     &     idx_oplist2
      

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'fs_reo_drv at work')
      end if

      ngas = orb_info%ngas
      nspin = orb_info%nspin

      mel_arr => op_info%mel_arr
      op2list => op_info%op2list

      reo_inf0 => fl_item%reo
      fact = 1d0
      if (add) fact = fl_item%bcontr%fact

      idx_opout = idx_oplist2(reo_inf0%label_out,op_info)

      if (update.and.idx_opout.ne.idx_tgt) then
        write(luout,*) 'formula target: ',
     &       trim(op_info%op_arr(idx_tgt)%op%name)
        write(luout,*) '[REO] result: ',trim(reo_inf0%label_out)
        call quit(1,'fs_reo_drv','inconsistency: target<>result')
      end if

      idx_opin = idx_oplist2(reo_inf0%label_in,op_info)
      idx_opout = idx_oplist2(reo_inf0%label_out,op_info)

      iblk_opin = reo_inf0%iblk_in
      iblk_opout = reo_inf0%iblk_out

      nj_opin = reo_inf0%nj_in
      nj_opout = reo_inf0%nj_out

      tra_opin =  reo_inf0%tra_in
      tra_opout = reo_inf0%tra_out

      if (.not.update) then
        idxme = op2list(idx_opout)
        me_opout => mel_arr(idxme)%mel
        xret_pnt => xret_dummy
        type_xret_loc = 0
      else
        me_opout => me_tgt
        xret_pnt => xret_blk(iblk_opout:iblk_opout)
        type_xret_loc = type_xret
      end if

      idxme = op2list(idx_opin)
      me_opin => mel_arr(idxme)%mel

c      igamt_opin = me_opin%gamt
c      mst_opin   = me_opin%mst

c      idxme = op2list(idx_opout)
c      me_opout => mel_arr(idxme)%mel

c      igamt_opout = me_opout%gamt
c      mst_opout   = me_opout%mst

      iocc_opout => reo_inf0%occ_opout
      iocc_opin => reo_inf0%occ_opin

      irst_opout => reo_inf0%rst_opout
      irst_opin => reo_inf0%rst_opin

c      iocc_opout => reo_inf0%occ_opout
c      irst_opout => reo_inf0%rst_opout
 
      call init_reo_info(reo_info)
      call interface_reo_info(reo_info,reo_inf0,
     &     str_info,orb_info)

      ! translate records to offset in file:
      ffop => me_opin%fhand
      idoff_opin = ffop%length_of_record*(ffop%current_record-1)
      if (ffop%unit.le.0) call file_open(ffop)
      ffop => me_opout%fhand
      idoff_opout = ffop%length_of_record*(ffop%current_record-1)
      if (ffop%unit.le.0) call file_open(ffop)

      call reo_op_wmaps_c(fact,
     &     update,xret_pnt,type_xret_loc,
     &     me_opin,me_opout,
     &     tra_opin, tra_opout,
     &     iblk_opin,iblk_opout,
     &     idoff_opin,idoff_opout,
     &     reo_info,
     &     str_info,strmap_info,orb_info)

      call dealloc_reo_info(reo_info)

      return
      end

