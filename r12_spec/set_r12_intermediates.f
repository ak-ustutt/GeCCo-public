*----------------------------------------------------------------------*
c      subroutine set_r12_intermediates(formula_vint,op_info,orb_info)
      subroutine set_r12_intermediates(form_info,op_info,orb_info)
*----------------------------------------------------------------------*
*     Routine which calls the routines required to generate the formulae 
*     for R12 intermediates.
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'par_opnames_gen.h'
c      include 'def_formula.h'
      include 'mdef_formula_info.h'
      include 'def_formula_item.h'
      include 'par_formnames_gen.h'

c      type(formula), intent(inout), target ::
c     &     formula_vint

      type(formula_info), intent(inout) ::
     &     form_info

      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula), pointer ::
     &     form_pnt

      integer ::
     &     idx

      integer,external ::
     &     idx_oplist2, idx_formlist

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '==================================='
        write(luout,*) ' output from set_r12_intermediates'
        write(luout,*) '==================================='
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! Add the V-intermediate.
      call add_formula(form_info,label_r12_vint)
      idx = idx_formlist(label_r12_vint,form_info)
      form_pnt => form_info%form_arr(idx)%form

      call set_v_intermediate(form_pnt,op_info,orb_info)

      ! Add the B-intermediate.
      call add_formula(form_info,label_r12_bint)
      idx = idx_formlist(label_r12_bint,form_info)
      form_pnt => form_info%form_arr(idx)%form

      call set_b_intermediate(form_pnt,op_info,orb_info)


      call atim_csw(cpu,sys,wall)
      call prtim(luout,'Total R12 interm.',cpu-cpu0,sys-sys0,wall-wall0)

      end
