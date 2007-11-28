      subroutine fac_r12_inter(fform,form_info,op_info,
     &     orb_info)
*----------------------------------------------------------------------*
*     Wrapper which calls the routines to factorise out the R12 
*     intermediate terms.
*     GWR November 2007
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'par_opnames_gen.h'
      include 'mdef_formula_info.h'
      include 'def_formula_item.h'
      include 'par_formnames_gen.h'

      type(formula), intent(inout) ::
     &     fform
      type(formula_info), intent(inout) ::
     &     form_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      ! Local variables.
      type(formula_item) ::
     &     form_lag

      character ::
     &     name*(form_maxlen_label*2)
      character, allocatable ::
     &     ops*64(:)
      integer ::
     &     nops

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0


      if (ntest.eq.100) then
        write(luout,*) '==========================='
        write(luout,*) ' output from fac_r12_inter '
        write(luout,*) '==========================='
        write(luout,*) 'called for ',trim(fform%label)
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! Extract the input formula.
      call init_formula(form_lag)
      call read_form_list(fform%fhand,form_lag)

      ! Factor out the V-intermediate.
      call factor_v(form_lag,op_info,orb_info)

      ! Factor out the V+-intermediate.
      call factor_vbar(form_lag,op_info,orb_info)

      ! Factor out B-intermediate.
      call factor_b(form_lag,op_info,orb_info)

      ! Delete the nodes in the formula which have not been factorised
      ! with the preceding intermediates.

c dbg ! Deleting of B.
      nops = 3
      allocate(ops(nops))
      ops(1) = op_r12
      ops(2) = op_rba
      ops(3) = op_b_inter
      call delete_non_fact(ops,nops,form_lag,op_info)
      deallocate(ops)

      ! Save the factored formula to the old formula file.
      call file_delete(fform%fhand)
      write(name,'(a,".fml")') trim(fform%label)
      call file_init(fform%fhand,name,ftyp_sq_unf,0)
      call write_form_list(fform%fhand,form_lag,fform%comment)

      call atim_csw(cpu,sys,wall)
      write(luout,'(/a)')trim(fform%label)
      call prtim(luout,'Factor R12 interm.',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

