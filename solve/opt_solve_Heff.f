*----------------------------------------------------------------------*
      subroutine opt_solve_Heff(n_states,
     &     op_info,form_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     Solve the EVP problem for the effective Hamiltonian
*
*     It is intended to be used by the solve_nleq to update
*     the eigenvectors and eigenvalues of Heff in each iteration
*
*     This subroutine is "dirty": it explicitly uses labels for 
*     the ME lists, assuming that these where correctly defined. It
*     should become obsolete with a more general solver, able to 
*     evaluate targets.
*
*     Yuri, Ago 2015
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'mdef_formula_info.h'
      include 'ifc_adv.h'
      include 'mdef_target_info.h'
      include 'ifc_input.h'

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer ::
     &     n_states

      type(formula_info) ::
     &     form_info
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info

      integer ::
     &     req_state, idxmel, idoff, idx,
     &     rec_Heff, rec_C_MS, rec_C0_bar,
     &     i_state, i_state_tr, j_state

      type(me_list), pointer ::
     &     mel_Heff, mel_C_MS, mel_E_MS, mel_C0_bar, mel_C_ji

      type(filinf), pointer ::
     &     ffop, ffop2

      character(len_target_name) ::
     &     c_st_i, c_st_j

      real(8) ::
     &     coef_Heff(n_states), C_i, ratio_C_ji,
     &     Heff_tr(n_states*n_states),Heff_s(n_states*n_states)
      real(8), parameter ::
     &     small_den = 0.01d0

      logical ::
     &     closeit, closeit2, Heff_symm

      integer, external ::
     &     idx_mel_list, get_mel_record, idx_formlist
      character(len_target_name), external ::
     &     state_label

      real(8) ::
     &    xdum

      type(formula), pointer ::
     &     f_eval
      type(formula_item) ::
     &     fl_eval

      call get_argument_value('method.MRCC','req_state',
     &     ival=req_state)
      call get_argument_value('method.MRCC','Heff_symm',
     &     lval=Heff_symm)
      
      idxmel = idx_mel_list("ME_C0_bar",op_info)
      mel_C0_bar => op_info%mel_arr(idxmel)%mel

      idxmel = idx_mel_list('ME_pack_Heff_MS',op_info)
      mel_Heff => op_info%mel_arr(idxmel)%mel

      idxmel = idx_mel_list('ME_C_MS',op_info)
      mel_C_MS => op_info%mel_arr(idxmel)%mel

      idxmel = idx_mel_list('ME_E_MS',op_info)
      mel_E_MS => op_info%mel_arr(idxmel)%mel
      
      rec_C0_bar = get_mel_record( mel_C0_bar)
      rec_Heff = get_mel_record( mel_Heff)
      rec_C_MS = get_mel_record( mel_C_MS)

      idx = idx_formlist('FOPT_pack_Heff_MS',form_info)
      if (idx.le.0)
     &     call quit(1,'opt_get_C0bar',
     &     'did not find formula FOPT_pack_P_S')
      f_eval => form_info%form_arr(idx)%form
      call read_form_list(f_eval%fhand,fl_eval,.true.)
      do i_state = 1, n_states*n_states
       call evaluate2(fl_eval,.true.,.true.,
     &      op_info,str_info,strmap_info,orb_info,xdum,.false.)
       call mel_adv_state(mel_Heff,n_states*n_states)
       call mel_adv_state(mel_C0_bar,n_states)
       if ( MOD(i_state, n_states) .EQ. 0) then
        call op_adv_state(["C0_1"],1,n_states,op_info,.true.)
        call op_adv_state(["T"],1,n_states,op_info,.false.)
       end if
      end do

!     Diagonalize the symmetrized Heff
      if (Heff_symm) then
       do i_state = 1, n_states*n_states
        i_state_tr = n_states*(MOD((i_state-1),n_states))
     &       + (i_state-1)/n_states + 1
        ffop => mel_Heff%fhand
        idoff = ffop%length_of_record*(i_state-1)
        call get_vec( ffop, Heff_tr(i_state_tr:),
     &       idoff+1, idoff+1)
        call get_vec( ffop, Heff_s(i_state:),
     &       idoff+1, idoff+1)
       end do
       do i_state = 1, n_states*n_states
        Heff_s(i_state) = (Heff_s(i_state) + Heff_tr(i_state))/2
        ffop => mel_Heff%fhand
        idoff = ffop%length_of_record*(i_state-1)
        call put_vec( ffop, [Heff_s(i_state)],
     &       idoff+1, idoff+1)
       end do
      end if

      call diag_packed_op(mel_Heff,mel_C_MS,mel_E_MS,
     &     n_states,verbose=.false.)
      call switch_mel_record( mel_E_MS, req_state)
      
!     Put coefficients in the right ME list:
      ffop => mel_C_MS%fhand
      if (.not.associated(ffop))
     &     call quit(1,'solve_nleq - solve Heff',
     &     'No file assigned to list: '//
     &     trim(mel_C_MS%label))
      if (ffop%unit.le.0) then
       call file_open(ffop)
       closeit = .true.
      else
       closeit = .false.
      end if
      
      idoff = ffop%length_of_record*(req_state-1)
      call get_vec( ffop, coef_Heff, idoff+1, idoff+n_states)
c     dbg
c      write(luout,*) "Heff eigenvec (from opt_solve_Heff): ", coef_Heff
c     dbg
      
      do i_state = 1, n_states
       c_st_i = state_label(i_state,.true.)

       C_i = coef_Heff(i_state)
       if (ABS(C_i).LT.small_den) then
        write(luout, '("WARNING: small coefficient: ", E10.5,'//
     &       '" Using ", F8.5, " as denominator.")') C_i, small_den
        C_i = sign(small_den, C_i)
       end if

       do j_state = 1, n_states
        c_st_j = state_label(j_state,.true.)
        ratio_C_ji = coef_Heff(j_state)/C_i
         
        idxmel = idx_mel_list(
     &       'ME_C_MS'//trim(c_st_j)//trim(c_st_i),
     &       op_info)
        mel_C_ji => op_info%mel_arr(idxmel)%mel
        ffop2 => mel_C_ji%fhand
        if (.not.associated(ffop2))
     &       call quit(1,'solve_nleq - solve Heff',
     &       'No file assigned to list: '//
     &       trim(mel_C_ji%label))
        if (ffop2%unit.le.0) then
         call file_open(ffop2)
         closeit2 = .true.
        else
         closeit2 = .false.
        end if

        idoff = ffop2%length_of_record*(req_state-1)
        call put_vec( ffop2, [ratio_C_ji],
     &       idoff+1, idoff+1)
        
        if (closeit2)
     &       call file_close_keep(ffop2)
        
        call switch_mel_record( mel_C_ji, req_state)
       end do
      end do
      
      if (closeit)
     &     call file_close_keep(ffop)
      
      call switch_mel_record( mel_C0_bar, rec_C0_bar)
      call switch_mel_record( mel_Heff, rec_Heff)
      call switch_mel_record( mel_C_MS, rec_C_MS)
            
      return
      end
