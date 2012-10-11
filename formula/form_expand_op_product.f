*----------------------------------------------------------------------*
      subroutine form_expand_op_product(init,form_res,fac,
     &     title,label_res,label,nlabels,
     &     idx_sv,iblkmin,iblkmax,
     &     connect,nconnect,
     &     avoid,navoid,
     &     inproj,ninproj,
     &     descr,ndescr,
     &     fix_in,op_info,orb_info)
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
     &     ntest = 0

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, intent(in) ::
     &     nlabels, nconnect, navoid, ninproj, ndescr,
     &     idx_sv(nlabels), iblkmin(nlabels), iblkmax(nlabels),
     &     connect(nconnect*2), avoid(navoid*2), inproj(ninproj*4)
      character(len=*) ::
     &     descr(ndescr)
      character(len=*), intent(in) ::
     &     label_res, label(nlabels), title
      logical, intent(in) ::
     &     init, fix_in
      real(8), intent(in) ::
     &     fac

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
     &     nterms, ilabel, idx, idxres, nvtx, nops, len
      logical ::
     &     transpose

      integer, external::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.ge.100) then
        write(luout,*) '===================================='
        write(luout,*) ' output from form_expand_op_product'
        write(luout,*) '===================================='
        if (.not.init) write(luout,*) '(adding to existing formula)'
        write(luout,*) 'result op: ',trim(label_res)
        write(luout,*) 'nlabels = ',nlabels
        do ilabel = 1, nlabels
          write(luout,'(1x,i4,1x,a)') ilabel, trim(label(ilabel))
        end do
        write(luout,*) 'idx_sv  = ',idx_sv(1:nlabels)
        write(luout,*) 'iblkmin  = ',iblkmin(1:nlabels)
        write(luout,*) 'iblkmax  = ',iblkmax(1:nlabels)
        write(luout,*) 'fac  = ',fac
      end if

      ! check for conflicts
      if (ninproj.gt.0.and.ndescr.gt.0) then
         write(luout,*) 'ninproj = ',ninproj
         write(luout,*) 'ndescr  = ',ndescr
         call quit(1,'form_expand_op_product',
     &       'cannot use old "inproj" and new "descr" simultaneously!')
      end if

      call atim_csw(cpu0,sys0,wall0)

      idxres = idx_oplist2(label_res,op_info)
      if (idxres.le.0)
     &     call quit(1,'form_expand_op_product',
     &     'label not on list: '//trim(label_res))
      do ilabel = 1, nlabels
        len = len_trim(label(ilabel))
        transpose = (label(ilabel)(len-1:len).eq.'^+') 
        if (transpose) len = len-2
        idxop(ilabel) = idx_oplist2(label(ilabel)(1:len),op_info)
        if (idxop(ilabel).le.0)
     &       call quit(1,'form_expand_op_product',
     &       'label not on list: '//label(ilabel)(1:len))
        if (transpose) idxop(ilabel) = -idxop(ilabel)
      end do

      nvtx = nlabels
      idx_sv_cp = idx_sv
      nops = nlabels
      call unique_list(idx_sv_cp,nops)

      ! init formula
      call init_formula(flist_res)
      fl_pnt => flist_res
      if (init) then
        call new_formula_item(flist_res,
     &       command_set_target_init,idxres)
        fl_pnt => fl_pnt%next
      else
        call read_form_list(form_res%fhand,flist_res,.true.)
        do while(fl_pnt%command.ne.command_end_of_formula)
          fl_pnt => fl_pnt%next
        end do
      end if

      if (ndescr.eq.0) then
        ! old routine
        call expand_op_product2(fl_pnt,idxres,
     &     fac,nvtx,nops,
     &     idxop,idx_sv,
     &     iblkmin,iblkmax,
     &     connect,nconnect,
     &     avoid,navoid,
     &     inproj,ninproj,
     &     fix_in,op_info)
      else
        ! new routine
        call expand_op_product3(fl_pnt,idxres,
     &     fac,nvtx,nops,
     &     idxop,idx_sv,
     &     iblkmin,iblkmax,
     &     connect,nconnect,
     &     avoid,navoid,
     &     descr,ndescr,
     &     op_info)
      end if

      if(ntest.ge.10)then
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


