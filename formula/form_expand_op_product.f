*----------------------------------------------------------------------*
      subroutine form_expand_op_product(form_res,
     &     title,label_res,label,nlabels,idx_sv,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     Driver routine to set up the operator product
*
*     OP(label_res) = OP(label(1)) * OP(label(2)) ... OP(label(n))
*
*     initial version, not allowing for any special desires ...
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 1000

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, intent(in) ::
     &     nlabels, idx_sv(nlabels)
      character(*), intent(in) ::
     &     label_res, label(nlabels), title

      type(formula), intent(inout), target ::
     &     form_res

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist_res
      type(formula_item), pointer ::
     &     fl_pnt
      integer ::
     &     idxop(nlabels), idx_sv_cp(nlabels)
      integer ::
     &     nterms, ilabel, idx, idxres, nvtx, nops

      integer, external::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.ge.100) then
        write(luout,*) '===================================='
        write(luout,*) ' output from form_expand_op_product'
        write(luout,*) '===================================='
        write(luout,*) 'nlabels = ',nlabels
        write(luout,*) 'idx_sv  = ',idx_sv(1:nlabels)
      end if

      call atim_csw(cpu0,sys0,wall0)

      idxres = idx_oplist2(label_res,op_info)
      if (idxres.le.0)
     &     call quit(1,'form_expand_op_product',
     &     'label not on list: '//trim(label_res))
      do ilabel = 1, nlabels
        idxop(ilabel) = idx_oplist2(label(ilabel),op_info)
        if (idxop(ilabel).le.0)
     &       call quit(1,'form_expand_op_product',
     &       'label not on list: '//trim(label(ilabel)))
      end do

      nvtx = nlabels
      idx_sv_cp = idx_sv
      nops = nlabels
      call unique_list(idx_sv_cp,nops)

      ! init formula
      call init_formula(flist_res)
      fl_pnt => flist_res
      call new_formula_item(flist_res,
     &     command_set_target_init,idxres)
      fl_pnt => fl_pnt%next
      call expand_op_product2(fl_pnt,idxres,
     &     1d0,nvtx,nops,
     &     idxop,idx_sv,
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     0,0,
     &     op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'Generated formula')
        call print_form_list(luout,flist_res,op_info)
      endif

      ! Assign canonical name and comment.
      form_res%comment = title
      ! Write to disc.
      write(name,'(a,".fml")') trim(form_res%label)
      call file_init(form_res%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_res%fhand,flist_res,title)

      ! Delete linked list
      call dealloc_formula_list(flist_res)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'Operator product',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end


