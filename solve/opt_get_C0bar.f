*----------------------------------------------------------------------*
      subroutine opt_get_C0bar(n_states,
     &     op_info,form_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     Calculate the S-1 matrix in the P-space and get C0_bar
*
*     It is intended to be used by the solve_nleq to update
*     the C0_bar ME-list
*
*     This subroutine is "dirty": it explicitly uses labels for 
*     the ME lists, assuming that these where correctly defined. It
*     should become obsolete with a more general solver, able to 
*     evaluate targets.
*
*     Yuri, Sep 2015
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
     &     idxmel, i_state, j_state

      type(me_list), pointer ::
     &     mel_C0, mel_C0_bar, mel_P_S, mel_P_Sinv, mel_pnt

      character(len_target_name) ::
     &     c_st

      logical ::
     &     assume_orth

      integer, external ::
     &     idx_mel_list, get_mel_record
      character(len_target_name), external ::
     &     state_label

      call get_argument_value('method.MRCC','assume_orth',
     &     lval=assume_orth)
      
      idxmel = idx_mel_list("ME_C0",op_info)
      mel_C0 => op_info%mel_arr(idxmel)%mel

      idxmel = idx_mel_list("ME_C0_bar",op_info)
      mel_C0_bar => op_info%mel_arr(idxmel)%mel

      idxmel = idx_mel_list("ME_pack_P_S",op_info)
      mel_P_S => op_info%mel_arr(idxmel)%mel

      idxmel = idx_mel_list("ME_pack_P_Sinv",op_info)
      mel_P_Sinv => op_info%mel_arr(idxmel)%mel

      if (.not.assume_orth) then
!     Simulate target 'EVAL_pack_P_S'
       do i_state = 1, n_states*n_states
        call evaluate('FOPT_pack_P_S',.true.,
     &       op_info,form_info,str_info,strmap_info,orb_info)
        call mel_adv_state(mel_P_S,n_states*n_states)
        call mel_adv_state(mel_C0,n_states)
        if ( MOD(i_state, n_states) .EQ. 0) then
         call op_adv_state(["C0_1"],1,n_states,op_info,.true.)
        end if
       end do
       
!     Simulate target 'invert_P_S'
       call inv_pack_op(mel_P_S,mel_P_Sinv,n_states)
       j_state = 1
       c_st = state_label(j_state,.true.)
       do i_state = 1, n_states*n_states
        idxmel = idx_mel_list("ME_P_Sinv"//trim(c_st),op_info)
        mel_pnt => op_info%mel_arr(idxmel)%mel
        call list_copy(mel_P_Sinv,mel_pnt,.false.)
        call mel_adv_state(mel_P_Sinv,n_states*n_states)
        call mel_adv_state(mel_pnt,n_states)
        if ( MOD(i_state, n_states) .EQ. 0) then
         j_state = j_state + 1
         c_st = state_label(j_state,.true.)
        end if
       end do
      end if

! Simulate target 'EVAL_C0_non_orth'
      do i_state = 1, n_states
       if (assume_orth) then
        call list_copy(mel_C0,mel_C0_bar,.false.)
       else
        call evaluate('FOPT_C0_bar',.true.,
     &       op_info,form_info,str_info,strmap_info,orb_info)
       end if
       call mel_adv_state(mel_C0_bar,n_states)
       if (assume_orth) then
        call mel_adv_state(mel_C0,n_states)
       else
        do j_state = 1, n_states
         c_st = state_label(j_state,.true.)
         idxmel = idx_mel_list("ME_P_Sinv"//trim(c_st),op_info)
         mel_pnt => op_info%mel_arr(idxmel)%mel
         call mel_adv_state(mel_pnt,n_states)
        end do
       end if
      end do

      do i_state = 1, n_states
       c_st = state_label(i_state,.true.)
       idxmel = idx_mel_list("ME_C0_bar"//trim(c_st),op_info)
       mel_pnt => op_info%mel_arr(idxmel)%mel
       call list_copy(mel_C0_bar,mel_pnt,.false.)
       call mel_adv_state(mel_C0_bar,n_states)
      end do
           
      return
      end
