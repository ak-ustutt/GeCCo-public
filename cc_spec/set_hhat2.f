*----------------------------------------------------------------------*
      subroutine set_hhat2(formula_hhat,
     &     title,name_hhat,name_h,name_t,
     &     iblkmax_t,op_info)
*----------------------------------------------------------------------*
*     generate the formula for Hhat = e^{-T1}He^{T1}
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
c      include 'par_formnames_gen.h'
      include 'ifc_input.h'

      type(formula), intent(inout), target ::
     &     formula_hhat

      character*(*), intent(in) ::
     &     name_hhat, name_t, name_h, title

      integer, intent(in) ::
     &     iblkmax_t

      type(operator_info), intent(in) ::
     &     op_info

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_hhat
      type(formula_item), pointer ::
     &     form_pnt

      integer ::
     &     nterms,
     &     idxhhat,idxham,idxtop,trunc_type

      logical ::
     &     truncate

      integer, external ::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '======================'
        write(luout,*) ' output from set_hhat'
        write(luout,*) '======================'
      end if

      call atim_csw(cpu0,sys0,wall0)

      idxhhat = idx_oplist2(name_hhat,op_info)
      idxham  = idx_oplist2(name_h,op_info)
      idxtop  = idx_oplist2(name_t,op_info)
      if (idxham.lt.0.or.idxhhat.lt.0.or.idxtop.lt.0)
     &     call quit(1,'set_hhat2',
     &     'required operators are not yet defined')

      ! initialize formula
      call init_formula(form_hhat)
      form_pnt => form_hhat
      ! put [INIT] at the beginning
      call new_formula_item(form_pnt,command_set_target_init,idxhhat)
      form_pnt => form_pnt%next

      ! expand e^{-T1} H e^T1 
      call expand_op_bch(form_pnt,4,idxhhat,
     &     1d0,-1,idxham,1d0,idxtop,1,iblkmax_t,op_info)

      if (iblkmax_t.eq.2) then
        call get_argument_value('method.R12','trunc',ival=trunc_type)
        truncate = trunc_type.ge.0
        if (is_keyword_set('method.truncate').gt.0) then
          truncate = is_keyword_set('method.truncate').gt.0
          if(truncate)
     &       call get_argument_value('method.truncate','trunc_type',
     &                                ival=trunc_type)
        end if
        if (truncate)
     &       call r12_truncation(form_pnt,trunc_type,
     &       0,idxham,0,idxtop,op_info)
      end if

      ! reorder
      call reorder_formula(form_hhat,op_info)
c dbg
c      print *,'generated (1):'
c      call print_form_list(luout,form_hhat,op_info)
c      stop 'testing'
c dbg      
      ! write to disc
c      formula_hhat%label = label_cchhat
      formula_hhat%comment = title
      write(name,'(a,".fml")') trim(formula_hhat%label)
      call file_init(formula_hhat%fhand,name,ftyp_sq_unf,0)
      call write_form_list(formula_hhat%fhand,form_hhat,title)

      call dealloc_formula_list(form_hhat)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'CC Hhat',cpu-cpu0,sys-sys0,wall-wall0)

      end
