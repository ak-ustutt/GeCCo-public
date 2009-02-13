*----------------------------------------------------------------------*
      subroutine set_experimental_formula(form_exp,
     &     title,label,nlabels,ansatz,approx,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     template for setting up formulae
*     ansatz, approx can be freely used
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'ifc_input.h'

      type(formula), intent(inout), target ::
     &     form_exp
      integer, intent(in) ::
     &     ansatz, nlabels
      character(*), intent(in) ::
     &     label(nlabels), title, approx

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info

#if defined (EXPERIMENTAL)
      call set_my_experimental_formula(form_exp,
     &     title,label,nlabels,ansatz,approx,
     &     op_info,orb_info)
#endif

      return
      end
*----------------------------------------------------------------------*
      subroutine set_experimental_keywords()
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'

#if defined (EXPERIMENTAL)
      call set_my_experimental_keywords()
#endif

      return
      end
*----------------------------------------------------------------------*
      subroutine set_experimental_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'
      include 'opdim.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      character(*), intent(in) ::
     &     env_type

#if defined(EXPERIMENTAL)
      call set_my_experimental_targets(tgt_info,orb_info,env_type)
#endif

      return
      end

