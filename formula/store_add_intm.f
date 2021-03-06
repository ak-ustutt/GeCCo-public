      subroutine store_add_intm(fl_item,
     &     contr,op_info,orb_info)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_orbinf.h'

      type(formula_item), intent(inout) ::
     &     fl_item
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      character(len=len_opname) ::
     &     label_res, label_op
      integer ::
     &     iblk_res, iblk_op, nj_res, nj_op, idummy,
     &     idxblk, nidx
      logical ::
     &     tra_res, tra_op, ldummy
      integer, pointer ::
     &     occ_res(:,:,:), occ_op(:,:,:),
     &     rst_res(:,:,:,:,:,:), rst_op(:,:,:,:,:,:),
     &     itf_index_info(:)
      type(operator), pointer ::
     &     op_res, op_add

      real(8) ::
     &     fact, fact_itf

      integer, external ::
     &     get_nidx4contr

c dbg
c      print *,'call to store_add_intm'
c dbg      

      op_res => op_info%op_arr(contr%idx_res)%op
      op_add => op_info%op_arr(contr%vertex(1)%idx_op)%op
      
      idummy = -1
      ldummy = .false.

      label_res = op_res%name
      label_op  = op_add%name
      nj_res = op_res%njoined
      nj_op  = op_add%njoined
      iblk_res  = contr%iblk_res
c      iblk_op   = contr%vertex(1)%iblk_op
      ! still the strange storage in contr:
      iblk_op   = (contr%vertex(1)%iblk_op-1)/nj_op + 1
      tra_res  = contr%dagger
      tra_op   = contr%vertex(1)%dagger

      idxblk = (iblk_res-1)*nj_res+1
      occ_res => op_res%ihpvca_occ(1:,1:,idxblk:idxblk-1+nj_res)
      rst_res => op_res%igasca_restr(1:,1:,1:,1:,1:,
     &                                   idxblk:idxblk-1+nj_res)
      idxblk = (iblk_op-1)*nj_op+1
      occ_op => op_add%ihpvca_occ(1:,1:,idxblk:idxblk-1+nj_op)
      rst_op => op_add%igasca_restr(1:,1:,1:,1:,1:,
     &                                   idxblk:idxblk-1+nj_op)
      ! dummy for itf
      nidx = get_nidx4contr(contr)
      allocate(itf_index_info(3+2*nidx)) ! need to store the index info twice

      fact = contr%fac
      if (contr%total_sign==0)
     &     call quit(1,'store_add_intm','Something is wicked !')
      fact_itf = dble(contr%total_sign)*contr%fac
      
      call itf_set_index_info(itf_index_info,
     &     contr,contr,.true.,1,1,1,contr%nvtx) ! last two arguments are dummy

      call store_bc(fl_item,
     &     fact,fact_itf,
     &     label_res,label_op,'---',
     &     iblk_res,iblk_op,idummy,
     &     tra_res,tra_op,ldummy,
     &     nj_res,nj_op,0,
     &     occ_res,occ_op,idummy,
     &     rst_res,rst_op,idummy,
     &     idummy,idummy,idummy,
     &     idummy,idummy,idummy,0,
     &     idummy,idummy,
     &     idummy,idummy,
     &     itf_index_info,
     &     orb_info)
     
      deallocate(itf_index_info)

      return
      end
