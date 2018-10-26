*----------------------------------------------------------------------*
      subroutine fs_contr_drv(xret_blk,type_xret,idx_tgt,me_tgt,
     &         fl_item,update,add,reo,
     &         op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     driver for [CONTR] operation in frm_schedX() (X>=2)
*     if (update): [ADD] to target
*     if (add):    [ADD] to target or other iterm
*     if (reo):    [REORDER] result (on-the-fly in contr_op1op2)
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
     &     ntest = 1000

      real(8), intent(inout), target ::
     &     xret_blk(*)
      integer, intent(in) ::
     &     type_xret, idx_tgt
      logical, intent(in) ::
     &     update, reo, add
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
     &     n_cnt, idx_res, idx_op1, idx_op2,
     &     iblk_res, iblk_op1, iblk_op2, iblk_tmp,
     &     nj_res, nj_op1, nj_op2, nj_tmp, idxme,
     &     igamt_res, mst_res, igamt_op1, igamt_op2, mst_op1, mst_op2,
     &     type_xret_loc,
     &     idoff_op1, idoff_op2, idoff_res
      logical ::
     &     self, tra_res, tra_op1, tra_op2
      real(8) ::
     &     fact
      type(reorder_info) ::
     &     reo_info
      real(8), target ::
     &     xret_dummy(1)

      type(operator), pointer ::
     &     op_tmp
      type(me_list), pointer ::
     &     me_tmp

      integer, pointer ::
     &     op2list(:)
      type(me_list_array), pointer ::
     &     mel_arr(:)
      type(binary_contr), pointer ::
     &     bc_info
      type(me_list), pointer ::
     &     me_res, me_op1, me_op2
      type(reorder), pointer ::
     &     reo_inf0
      type(filinf), pointer ::
     &     ffop
      integer, pointer ::
     &     merge_op1(:), merge_op2(:), merge_op1op2(:), merge_op2op1(:),
     &     iocc_op1(:,:,:), iocc_op2(:,:,:), iocc_res(:,:,:),
     &     irst_op1(:,:,:,:,:,:),irst_op2(:,:,:,:,:,:),
     &     irst_res(:,:,:,:,:,:),irst_tmp(:,:,:,:,:,:),
     &     irst_cnt(:,:,:,:,:,:),
     &     iocc_ex1(:,:,:), iocc_ex2(:,:,:),
     &     iocc_cnt(:,:,:), iocc_tmp(:,:,:),
     &     irst_ex1(:,:,:,:,:,:),irst_ex2(:,:,:,:,:,:)
      real(8), pointer ::
     &     xret_pnt(:)
      integer, target ::
     &     idummy(1)

      integer, external ::
     &     idx_oplist2
      

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'fs_contr_drv at work')
        write(lulog,*) 'update, add, reo: ',update, add, reo
      end if

      ngas = orb_info%ngas
      nspin = orb_info%nspin

      mel_arr => op_info%mel_arr
      op2list => op_info%op2list

      bc_info  => fl_item%bcontr
      reo_inf0 => fl_item%reo

      n_cnt = bc_info%n_cnt
      
      if (bc_info%n_operands.gt.2.or.bc_info%n_cnt.eq.0) then
        write(lulog,*) '[CONTR]: something is wrong:'
        call prt_bcontr(lulog,bc_info)
        call quit(1,'fs_contr_drv','something is wrong')
      end if
      self = bc_info%n_operands.eq.1


      idx_res = idx_oplist2(bc_info%label_res,op_info)
      if (update.and.idx_res.ne.idx_tgt) then
        write(lulog,*) 'formula target: ',
     &       trim(op_info%op_arr(idx_tgt)%op%name)
        write(lulog,*) '[CONTR] result: ',trim(bc_info%label_res)
        call quit(1,'fs_contr_drv','inconsistency: target<>result')
      end if

      idx_op1 = idx_oplist2(bc_info%label_op1,op_info)
      if (.not.self) idx_op2 = idx_oplist2(bc_info%label_op2,op_info)

      iblk_res = bc_info%iblk_res
      iblk_op1 = bc_info%iblk_op1
      if (.not.self) iblk_op2 = bc_info%iblk_op2

      nj_res = bc_info%nj_res
      nj_op1 = bc_info%nj_op1
      if (.not.self) nj_op2 = bc_info%nj_op2
      if (self) nj_op2 = 0

      tra_res = bc_info%tra_res
      tra_op1 = bc_info%tra_op1
      if (.not.self) tra_op2 = bc_info%tra_op2

      if (.not.update) then
        idxme = op2list(idx_res)
        me_res => mel_arr(idxme)%mel
        xret_pnt => xret_dummy
        type_xret_loc = 0
      else
        me_res => me_tgt
        xret_pnt => xret_blk(iblk_res:iblk_res)
        type_xret_loc = type_xret
      end if

      igamt_res = me_res%gamt
      mst_res = me_res%mst

      idxme = op2list(idx_op1)
      me_op1 => mel_arr(idxme)%mel

      igamt_op1 = me_op1%gamt
      mst_op1   = me_op1%mst

      if (.not.self) then
        idxme = op2list(idx_op2)
        me_op2 => mel_arr(idxme)%mel

        igamt_op2 = me_op2%gamt
        mst_op2   = me_op2%mst
      else
        igamt_op2 = 1
        mst_op2   = 0
      end if
c dbg
c      write(lulog,*) "DBG info: fs_contr_drv"
c      write(lulog,*) "ME op1 and record:",
c     &     trim(me_op1%label),me_op1%fhand%current_record
c      if (.not.self)
c     &     write(lulog,*) "ME op2 and record:",
c     &     trim(me_op2%label),me_op2%fhand%current_record
c dbgend

      iocc_res => bc_info%occ_res
      iocc_op1 => bc_info%occ_op1
      iocc_ex1 => bc_info%occ_ex1

      irst_res => bc_info%rst_res
      irst_op1 => bc_info%rst_op1
      irst_ex1 => bc_info%rst_ex1

      if (.not.self) then
        iocc_op2 => bc_info%occ_op2
        irst_op2 => bc_info%rst_op2
        iocc_ex2 => bc_info%occ_ex2
        irst_ex2 => bc_info%rst_ex2
      else 
        allocate(iocc_op2(ngastp,2,1),
     &           irst_op2(2,ngas,2,2,nspin,1),
     &           iocc_ex2(ngastp,2,1),
     &           irst_ex2(2,ngas,2,2,nspin,1))
        iocc_op2 = 0
        irst_op2 = 0
        iocc_ex2 = 0
        irst_ex2 = 0
      end if
      iocc_cnt => bc_info%occ_cnt
      irst_cnt => bc_info%rst_cnt

      merge_op1 => bc_info%merge_op1
      if (.not.self) then
        merge_op2 => bc_info%merge_op2
      else
        merge_op2 => idummy
      end if
      merge_op1op2 => bc_info%merge_op1op2
      merge_op2op1 => bc_info%merge_op2op1

      fact = bc_info%fact

      call init_reo_info(reo_info)
      if (reo) then
        allocate(op_tmp,me_tmp)
        call interface_reo_info(reo_info,reo_inf0,
     &       str_info,orb_info)
        iocc_tmp => reo_inf0%occ_opin
        irst_tmp => reo_inf0%rst_opin
        iblk_tmp = 1
        nj_tmp = reo_inf0%nj_in
        iocc_res => reo_inf0%occ_opout
        irst_res => reo_inf0%rst_opout
        if (reo_inf0%nj_out.ne.nj_res) then
          write(lulog,*) reo_inf0%nj_out,nj_res
          call quit(1,'fs_contr_drv','reo_inf0%nj_out.ne.nj_res!')
        end if
        call set_ps_op(op_tmp,'_OP_TMP',
     &       iocc_tmp,irst_tmp,nj_tmp,
     &       orb_info)
        me_tmp%op => op_tmp
        ! Should Ms be fixed?
        me_tmp%fix_vertex_ms = me_res%fix_vertex_ms
        call set_ps_list(me_tmp,'_OP_TMP',
     &           0,0,mst_res,igamt_res,0,
     &           str_info,strmap_info,orb_info)
      else
        iblk_tmp = iblk_res
        iocc_tmp => bc_info%occ_res
        irst_tmp => bc_info%rst_res
        me_tmp => me_res
      end if

      ! translate records to offset in file:
      ! (makes life easier in case we once decide to use
      ! one scratch file only: no changes to contr_op1op2 necessary)
      ffop => me_op1%fhand
      idoff_op1 = ffop%length_of_record*(ffop%current_record-1)
      if (ffop%unit.le.0) call file_open(ffop)
      if (.not.self) ffop => me_op2%fhand
      if (.not.self)
     &     idoff_op2 = ffop%length_of_record*(ffop%current_record-1)
      if (.not.self.and.ffop%unit.le.0) call file_open(ffop)
      ffop => me_res%fhand
      idoff_res = ffop%length_of_record*(ffop%current_record-1)
      if (ffop%unit.le.0) call file_open(ffop)

      if (ntest.ge.100)
     &         write(lulog,*) 'calling contraction kernel'
          ! do the contraction
      call contr_op1op2(fact,1d0,
     &       add,self,xret_pnt,type_xret_loc,
     &       me_op1,me_op2,me_res, me_tmp,
     &       tra_op1, tra_op2, tra_res,
     &       iblk_op1,iblk_op2,iblk_res,iblk_tmp,
     &       idoff_op1,idoff_op2,idoff_res,
     &       iocc_ex1,iocc_ex2,iocc_cnt,
     &       irst_ex1,irst_ex2,irst_cnt,2,
     &       iocc_op1, iocc_op2, iocc_res, iocc_tmp,
     &       irst_op1,irst_op2,irst_res, irst_tmp,
     &       merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &       nj_op1, nj_op2, nj_res, n_cnt,
     &       mst_op1,mst_op2,mst_res,
     &       igamt_op1,igamt_op2,igamt_res,
     &       reo_info,
     &       str_info,strmap_info,orb_info)
      if (ntest.ge.100)
     &         write(lulog,*) 'returned from contraction kernel'

      if (self) deallocate(iocc_op2,irst_op2,iocc_ex2)

      if (reo) then
        call dealloc_me_list(me_tmp)
        call dealloc_operator(op_tmp)
        call dealloc_reo_info(reo_info)
        deallocate(op_tmp,me_tmp)
      end if

      return
      end

