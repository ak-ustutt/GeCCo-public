*----------------------------------------------------------------------*
      subroutine set_experimental_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets for response theory calculations with GeCCo
*
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer, parameter ::
     &     len_short = 20, len_long = 200

      integer ::
     &     min_rank, max_rank,
     &     isim, ncat, nint, icnt,
     &     isym, ms, msc, sym_arr(8),
     &     ord, op_par, len_op_exp, side, x_max_ord, maxord,
     &     freq_idx, digit, ilabels, ord2, op_par2, x_max_ord2

      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(2)

      character(len_short) ::
     &     opname, melname, defmelname, formname, formname2,
     &     opname2, melname2, defmelname2, hubname, optname,
     &     solvename

      character(len_long) ::
     &     opexp

      character ::
     &     op_name*4, op_parent*1, lagf_name*12, op_exp*50,
     &     pert_ord*1, lag_name*11, res_name*9, leftright*1,
     &     resx_name*10,res_lag_name*15,resl_lag_name*16,
     &     resr_lag_name*16, def_me_x_name*11, me_x_name*7,
     &     def_me_resx_name*14, me_resx_name*10, opt_name*8,
     &     def_me_xalt_name*11, lresp_name*8,
     &     def_me_l_name*15, me_l_name*11, solve_x_name*10,
     &     solve_xalt_name*10, eval_lag_name*11, opt_l_name*12,
     &     eval_props_name*13, pert*7

      logical ::
     &     comb_ok

      real(8),allocatable ::
     &     freq(:)

      logical, allocatable ::
     &     evaluate(:)

      integer, allocatable ::
     &     ifreq(:)

      logical, external ::
     &     next_comb

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting experimental targets ...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

*----------------------------------------------------------------------*
*     get maximum order and allocate arrays:
*----------------------------------------------------------------------*
      call get_argument_value('calculate.experimental','order',
     &                        ival=maxord)
      allocate(freq(maxord),evaluate(0:maxord))
      freq = 0
      evaluate = .false.
      
      ! to be improved:
      evaluate(maxord) = .true.
*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! define op_ham as H(0)
      call ord_parameters(-1,parameters,0,3,-1)
      call set_rule(op_ham,ttype_op,SET_ORDER,op_ham,
     &              1,1,parameters,1,tgt_info)

      ! define V(1)
      call add_target('V(1)',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,1,1,1,.false.)
      call set_rule('V(1)',ttype_op,DEF_HAMILTONIAN,'V(1)',
     &              1,1,parameters,1,tgt_info)
      call ord_parameters(-1,parameters,1,3,0)  ! defined as freq. "joker"
      call set_rule('V(1)',ttype_op,SET_ORDER,'V(1)',
     &              1,1,parameters,1,tgt_info)

c      ! define frequency components V(1)1, V(1)2, ...,V(1)maxord
c      opname(1:len_short) = ' '
c      opname(1:4) = 'V(1)'
c      do freq_idx = 1,maxord
c        write(opname(5:5),'(i1)') freq_idx
c        call add_target(trim(opname),ttype_op,.false.,tgt_info)
c        call hop_parameters(-1,parameters,1,1,1,.false.)
c        call set_rule(trim(opname),ttype_op,DEF_HAMILTONIAN,
c     &                trim(opname),
c     &                1,1,parameters,1,tgt_info)
c        call ord_parameters(-1,parameters,1,3,freq_idx)
c        call set_rule(trim(opname),ttype_op,SET_ORDER,trim(opname),
c     &                1,1,parameters,1,tgt_info)
c      end do

      ! Hnew will be used as sum of H(0) and V(1)
      call add_target('Hnew',ttype_op,.false.,tgt_info)
      call set_dependency('Hnew',op_ham,tgt_info)
      call cloneop_parameters(-1,parameters,op_ham,.true.)
      call set_rule('Hnew',ttype_op,CLONE_OP,'Hnew',1,1,
     &              parameters,1,tgt_info)

      ! define scalar response lagrangian
      call add_target('LRESP',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('LRESP',ttype_op,DEF_HAMILTONIAN,'LRESP',
     &              1,1,parameters,1,tgt_info)

      ! define excitation operator T
      call add_target('T',ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,.false.,1,2,0,1)
      call set_rule('T',ttype_op,DEF_EXCITATION,'T',
     &              1,1,parameters,1,tgt_info)
 
      ! define deexcitation operator L
      call add_target('L',ttype_op,.false.,tgt_info)
      call set_dependency('L','T',tgt_info)
      call cloneop_parameters(-1,parameters,'T',.true.)
      call set_rule('L',ttype_op,CLONE_OP,'L',1,1,
     &              parameters,1,tgt_info)

      ! define operators T(0), T(1), T(2), ...
      ! and L(0), L(1), L(2), ...
      ! following (2n+1) and (2n+2) rules regarding the maximum order
      do op_par = 1,2
        if (op_par.eq.1.) then
          op_parent = 'T'
          op_name = 'T(x)'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        else
          op_parent = 'L'
          op_name = 'L(x)'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        end if
        do ord=0,x_max_ord
          write(op_name(3:3),'(i1)') ord
          call add_target(op_name,ttype_op,.false.,tgt_info)
          call set_dependency(op_name,op_parent,tgt_info)
          call cloneop_parameters(-1,parameters,op_parent,.false.)
          call set_rule(op_name,ttype_op,CLONE_OP,op_name,1,1,
     &                  parameters,1,tgt_info)
          allocate(ifreq(ord))
          ifreq = -1
          call ord_parameters(-1,parameters,ord,op_par,ifreq)
          deallocate(ifreq)
          call set_rule(op_name,ttype_op,SET_ORDER,op_name,
     &                  1,1,parameters,1,tgt_info)
          ! define X(1)1, X(1)2, ..., X(2)12, X(2)13, ...
          if (ord.gt.0) then
            opname(1:len_short) = ' '
            opname(1:4) = op_name
            allocate(ifreq(ord))
            do digit = 1,ord
              ifreq(digit) = digit
            end do
            comb_ok = .true.
            do while (comb_ok)
              do digit = 1,ord
                write(opname(4+digit:4+digit),'(i1)') ifreq(digit)
              end do
              call add_target(trim(opname),ttype_op,.false.,tgt_info)
              call set_dependency(trim(opname),op_parent,tgt_info)
              call cloneop_parameters(-1,parameters,op_parent,.false.)
              call set_rule(trim(opname),ttype_op,CLONE_OP,
     &                      trim(opname),1,1,
     &                      parameters,1,tgt_info)
              call ord_parameters(-1,parameters,ord,op_par,ifreq)
              call set_rule(trim(opname),ttype_op,SET_ORDER,
     &                      trim(opname),
     &                      1,1,parameters,1,tgt_info)
              comb_ok = next_comb(ifreq,ord,maxord)
            end do
            deallocate(ifreq)
          end if
        end do
      end do

      ! define residuals O(n)_L(0) and O(n)_T(0)
      ! following (2n+1) and (2n+2) rules regarding the maximum order
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          res_name = 'O(x)_L'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        else
          op_parent = 'L'
          res_name = 'O(x)_T'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        end if
        do ord=0,x_max_ord
          write(res_name(3:3),'(i1)') ord
          call add_target(res_name,ttype_op,.false.,tgt_info)
          call set_dependency(res_name,op_parent,tgt_info)
          call cloneop_parameters(-1,parameters,op_parent,.false.)
          call set_rule(res_name,ttype_op,CLONE_OP,res_name,1,1,
     &                  parameters,1,tgt_info)
        end do
      end do

      ! define left and right residuals
      ! O(0)SX; O(1)SX1, O(1)SX2, ...; O(2)SX12, O(2)SX13, ...; ...
      ! following the (2n+1) and (2n+2) rules regarding the maximum order
      opname(1:len_short) = ' '
      opname(1:6) = 'O(n)SX'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'L'
          opname(6:6) = 'T'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        else
          op_parent = 'T'
          opname(6:6) = 'L'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        end if
        do side = 1,2
          if (side.eq.1) then
            leftright = 'L'
          else
            leftright = 'R'
          end if
          opname(5:5) = leftright
          do ord=op_par-1,x_max_ord
            write(opname(3:3),'(i1)') ord
            opname(7:len_short) = ' '

            if (ord.gt.0) then
              call add_target(trim(opname),ttype_op,.false.,tgt_info)
              call set_dependency(trim(opname),op_parent,tgt_info)
              call cloneop_parameters(-1,parameters,op_parent,.false.)
              call set_rule(trim(opname),ttype_op,CLONE_OP,
     &                      trim(opname),1,1,
     &                      parameters,1,tgt_info)
            end if

            allocate(ifreq(ord))
            do digit = 1,ord
              ifreq(digit) = digit
            end do
            comb_ok = .true.
            do while (comb_ok)
              do digit = 1,ord
                write(opname(6+digit:6+digit),'(i1)') ifreq(digit)
              end do
              call add_target(trim(opname),ttype_op,.false.,tgt_info)
              call set_dependency(trim(opname),op_parent,tgt_info)
              call cloneop_parameters(-1,parameters,op_parent,.false.)
              call set_rule(trim(opname),ttype_op,CLONE_OP,
     &                      trim(opname),1,1,
     &                      parameters,1,tgt_info)
              comb_ok = next_comb(ifreq,ord,maxord)
            end do
            deallocate(ifreq)
          end do
        end do
      end do

      ! define scalar response lagrangian LRESP(n) of order n
      lresp_name = 'LRESP(n)'
      do ord=0,maxord
        write(lresp_name(7:7),'(i1)') ord
        call add_target(lresp_name,ttype_op,.false.,tgt_info)
        call set_dependency(lresp_name,'LRESP',tgt_info)
        call cloneop_parameters(-1,parameters,'LRESP',.false.)
        call set_rule(lresp_name,ttype_op,CLONE_OP,lresp_name,1,1,
     &                parameters,1,tgt_info)
      end do

      ! Diagonal
      call add_target(op_dia,ttype_op,.false.,tgt_info)
      call set_dependency(op_dia,'T(0)',tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        'T(0)',.false.)
      call set_rule(op_dia,ttype_op,CLONE_OP,
     &              op_dia,1,1,
     &              parameters,1,tgt_info)

      ! T1 transformed Hamiltonian
      call add_target(op_hhat,ttype_op,.false.,tgt_info)
      call set_dependency(op_hhat,op_ham,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_ham,.false.)
      call set_rule(op_hhat,ttype_op,CLONE_OP,
     &              op_hhat,1,1,
     &              parameters,1,tgt_info)

      ! Hbar intermediate
      call add_target(op_hbar,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &                    .false.,1,2,0,1)
      call set_rule(op_hbar,ttype_op,DEF_CC_HBAR_OP,
     &              op_hbar,1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      ! define response lagrangian
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'RESP_LAG'
      labels(2) = 'LRESP'
      labels(3) = 'L'
      labels(4) = 'Hnew'
      labels(5) = 'T'
      call add_target('RESP_LAG',ttype_frm,.false.,tgt_info)
      call set_dependency('RESP_LAG','LRESP',tgt_info)
      call set_dependency('RESP_LAG','L',tgt_info)
      call set_dependency('RESP_LAG','Hnew',tgt_info)
      call set_dependency('RESP_LAG','T',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'response lagrange functional',0,'---')
      call set_rule('RESP_LAG',ttype_frm,DEF_EXP_FORMULA,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! perturbation expansion of H: Hnew=H+V(1)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'H_FORM'
      call add_target('H_FORM',ttype_frm,.false.,tgt_info)
      call set_dependency('H_FORM','Hnew',tgt_info)
      call set_dependency('H_FORM',op_ham,tgt_info)
      call set_dependency('H_FORM','V(1)',tgt_info)
      call def_form_parameters(-1,
     &     parameters,2,'Hnew='//op_ham//'+V(1)','pert exp of H')
      call set_rule('H_FORM',ttype_frm,DEF_FORMULA,
     &              labels,1,1,
     &              parameters,2,tgt_info)

      ! perturbation expansion of T and L: T=T(0)+T(1)+T(2)+.../L=...
      ! following the (2n+1) and (2n+2) rules regarding the maximum order
      do op_par = 1,2
        len_op_exp = 6
        if (op_par.eq.1) then
          op_parent = 'T'
          op_exp(1:len_op_exp) = 'T=T(0)'
          op_name = 'T(x)'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        else
          op_parent = 'L'
          op_exp(1:len_op_exp) = 'L=L(0)'
          op_name = 'L(x)'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        end if
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = op_parent//'_FORM'
        call add_target(op_parent//'_FORM',ttype_frm,.false.,tgt_info)
        call set_dependency(op_parent//'_FORM',op_parent,tgt_info)
        call set_dependency(op_parent//'_FORM',
     &                      op_parent//'(0)',tgt_info)
        do ord=1,x_max_ord
          write(op_name(3:3),'(i1)') ord
          call set_dependency(op_parent//'_FORM',op_name,tgt_info)
          len_op_exp = len_op_exp+5
          write(op_exp(len_op_exp-4:len_op_exp),'(a1,a4)') '+',op_name
        end do
        call def_form_parameters(-1,
     &       parameters,2,op_exp(1:len_op_exp),
     &       'pert exp of '//op_parent)
        call set_rule(op_parent//'_FORM',ttype_frm,DEF_FORMULA,
     &                labels,1,1,
     &                parameters,2,tgt_info)
      end do

c      ! frequency expansion of V(1): V(1)=V(1)1+V(1)2+...
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'V(1)_FORM'
c      call add_target('V(1)_FORM',ttype_frm,.false.,tgt_info)
c      opexp(1:len_long) = ' '
c      opexp(1:5) = 'V(1)='
c      len_op_exp = 5
c      opname(1:len_short) = ' '
c      opname(1:4) = 'V(1)'
c      do freq_idx = 1,maxord
c        write(opname(5:5),'(i1)') freq_idx
c        opexp(len_op_exp+1:len_op_exp+6) = trim(opname)//'+'
c        len_op_exp = len_op_exp + 6
c        call set_dependency('V(1)_FORM',trim(opname),tgt_info)
c      end do
c      opexp(len_op_exp:len_op_exp) = ' '
c      len_op_exp = len_op_exp - 1
c      call def_form_parameters(-1,
c     &     parameters,2,opexp(1:len_op_exp),
c     &     'freq exp of V(1)')
c      call set_rule('V(1)_FORM',ttype_frm,DEF_FORMULA,
c     &              labels,1,1,
c     &              parameters,2,tgt_info)

      ! frequency expansion of X(n), n>0: X(1) = X(1)1+X(1)2+..., ...
      ! following the (2n+1) and (2n+2) rules regarding the maximum order
      formname(1:len_short) = ' '
      formname(1:9) = 'X(n)_FORM'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_name = 'T(x)'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        else
          op_name = 'L(x)'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        end if
        do ord=1,x_max_ord
          write(op_name(3:3),'(i1)') ord
          formname(1:4) = op_name(1:4)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(formname)
          call add_target(trim(formname),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(formname),op_name,tgt_info)
          opexp(1:len_long) = ' '
          opexp(1:5) = op_name//'='
          len_op_exp = 5
          opname(1:len_short) = ' '
          opname(1:4) = op_name
          allocate(ifreq(ord))
          do digit = 1,ord
            ifreq(digit) = digit
          end do
          comb_ok = .true.
          do while (comb_ok)
            do digit = 1,ord
              write(opname(4+digit:4+digit),'(i1)') ifreq(digit)
            end do
            opexp(len_op_exp+1:len_op_exp+5+ord) = trim(opname)//'+'
            len_op_exp = len_op_exp + 5 + ord
            call set_dependency(trim(formname),trim(opname),tgt_info)
            comb_ok = next_comb(ifreq,ord,maxord)
          end do
          deallocate(ifreq)
          opexp(len_op_exp:len_op_exp) = ' '
          len_op_exp = len_op_exp - 1
          call def_form_parameters(-1,
     &         parameters,2,opexp(1:len_op_exp),
     &         'freq exp of '//op_name)
          call set_rule(trim(formname),ttype_frm,DEF_FORMULA,
     &                  labels,1,1,
     &                  parameters,2,tgt_info)
        end do
      end do

      ! expand response lagrangian with Hnew=H+V(1),T=T(0)+T(1)+...,L=...
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'RESP_LAGF'
      labels(2) = 'RESP_LAG'
      labels(3) = 'H_FORM'
      labels(4) = 'T_FORM'
      labels(5) = 'L_FORM'
      call add_target('RESP_LAGF',ttype_frm,.false.,tgt_info)
      call set_dependency('RESP_LAGF','RESP_LAG',tgt_info)
      call set_dependency('RESP_LAGF','H_FORM',tgt_info)
      call set_dependency('RESP_LAGF','T_FORM',tgt_info)
      call set_dependency('RESP_LAGF','L_FORM',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'full response lagrangian',3,'---')
      call set_rule('RESP_LAGF',ttype_frm,EXPAND,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! expand RESP_LAG(n) with X(1)=X(1)1+X(1)2+..., X(2)=X(2)12+X(2)13+...
      ! and V(1)=V(1)1+V(1)2+...
      formname(1:len_short) = ' '
      formname(1:8) = 'F_LAG(n)'
      lag_name = 'RESP_LAG(n)'
      do ord = 1,maxord
        write(formname(7:7),'(i1)') ord
        write(lag_name(10:10),'(i1)') ord
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = trim(formname)
        labels(2) = lag_name
c        labels(3) = 'V(1)_FORM'
        ilabels = 2
        call add_target(trim(formname),ttype_frm,.false.,tgt_info)
        call set_dependency(trim(formname),lag_name,tgt_info)
c        call set_dependency(trim(formname),'V(1)_FORM',tgt_info)
        formname2(1:len_short) = ' '
        formname2(1:9) = 'X(n)_FORM'
        do op_par = 1,2
          if (op_par.eq.1) then
            formname2(1:1) = 'T'
            x_max_ord = int((real(maxord)-1)/2+0.6)
          else
            formname2(1:1) = 'L'
            x_max_ord = int((real(maxord)-2)/2+0.6)
          end if
          do ord2=1,x_max_ord
            write(formname2(3:3),'(i1)') ord2
            ilabels = ilabels+1
            labels(ilabels) = trim(formname2)
            call set_dependency(trim(formname),trim(formname2),tgt_info)
          end do
        end do
        call form_parameters(-1,
     &       parameters,2,'full response lagrangian with freqs',
     &       ilabels-2,'---')
        call set_rule(trim(formname),ttype_frm,EXPAND,
     &                labels,ilabels,1,
     &                parameters,2,tgt_info)
      end do

      ! expand RESS_LAG(n)_X with X(1)=X(1)1+X(1)2+..., X(2)=X(2)12+X(2)13+...
      ! and V(1)=V(1)1+V(1)2+...
      ! following (2n+1) and (2n+2) rules
      formname(1:len_short) = ' '
      formname(1:14) = 'RESFS_LAG(n)_X'
      resl_lag_name = 'RESS_LAG(n)_X'
      do op_par2 = 1,2
        if (op_par2.eq.1) then
          formname(14:14) = 'T'
          resl_lag_name(13:13) = 'T'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        else
          formname(14:14) = 'L'
          resl_lag_name(13:13) = 'L'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        end if
        do ord = op_par2-1,x_max_ord
          write(formname(11:11),'(i1)') ord
          write(resl_lag_name(10:10),'(i1)') ord
          do side = 1,2
            if (side.eq.1) then
              formname(5:5) = 'L'
              resl_lag_name(4:4) = 'L'
            else
              formname(5:5) = 'R'
              resl_lag_name(4:4) = 'R'
            end if
            labels(1:20)(1:len_target_name) = ' '
            labels(1) = trim(formname)
            labels(2) = resl_lag_name
c            labels(3) = 'V(1)_FORM'
            ilabels = 2
            call add_target(trim(formname),ttype_frm,.false.,tgt_info)
            call set_dependency(trim(formname),resl_lag_name,tgt_info)
c            call set_dependency(trim(formname),'V(1)_FORM',tgt_info)
            formname2(1:len_short) = ' '
            formname2(1:9) = 'X(n)_FORM'
            do op_par = 1,2
              if (op_par.eq.1) then
                formname2(1:1) = 'T'
                x_max_ord2 = int((real(maxord)-1)/2+0.6)
              else
                formname2(1:1) = 'L'
                x_max_ord2 = int((real(maxord)-2)/2+0.6)
              end if
              do ord2=1,x_max_ord2
                write(formname2(3:3),'(i1)') ord2
                ilabels = ilabels+1
                labels(ilabels) = trim(formname2)
                call set_dependency(trim(formname),
     &                              trim(formname2),tgt_info)
              end do
            end do
            call form_parameters(-1,
     &           parameters,2,'left and right residuals with freqs',
     &           ilabels-2,'---')
            call set_rule(trim(formname),ttype_frm,EXPAND,
     &                    labels,ilabels,1,
     &                    parameters,2,tgt_info)
          end do
        end do
      end do

      ! extract lagrangian of orders 0 to maxord
      lagf_name = 'RESP_LAGF(n)'
      do ord = 0,maxord
        write(lagf_name(11:11),'(i1)') ord
        write(pert_ord,'(i1)') ord
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = lagf_name
        labels(2) = 'RESP_LAGF'
        labels(3) = 'LRESP'
        call add_target(lagf_name,ttype_frm,.false.,tgt_info)
        call set_dependency(lagf_name,'RESP_LAGF',tgt_info)
        call form_parameters(-1,
     &       parameters,2,'full lagrangian of pert order '//pert_ord,
     &       ord,'---')
        call set_rule(lagf_name,ttype_frm,EXTRACT_ORDER,
     &                labels,3,1,
     &                parameters,2,tgt_info)
      end do

      ! extract terms obeying (2n+1) and (2n+2) rules
      lagf_name = 'RESP_LAGF(n)'
      lag_name = 'RESP_LAG(n)'
      lresp_name = 'LRESP(n)'
      do ord = 0,maxord
        write(lagf_name(11:11),'(i1)') ord
        write(lag_name(10:10),'(i1)') ord
        write(pert_ord,'(i1)') ord
        write(lresp_name(7:7),'(i1)') ord
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = lag_name
        labels(2) = lagf_name
        labels(3) = lresp_name
        call add_target(lag_name,ttype_frm,.false.,tgt_info)
        call set_dependency(lag_name,lagf_name,tgt_info)
        call set_dependency(lag_name,lresp_name,tgt_info)
        call form_parameters(-1,
     &       parameters,2,'lagrangian of pert order '//pert_ord,
     &       -1,'---')
        call set_rule(lag_name,ttype_frm,EXTRACT_ORDER,
     &                labels,3,1,
     &                parameters,2,tgt_info)
      end do

      ! extract terms of F_LAG(n) with correct freq. pattern
      formname(1:len_short) = ' '
      formname(1:6) = 'LAG(n)'
      formname2(1:len_short) = ' '
      formname2(1:8) = 'F_LAG(n)'
      lresp_name = 'LRESP(n)'
      do ord = 1,maxord
        write(formname(5:5),'(i1)') ord
        write(formname2(7:7),'(i1)') ord
        write(lresp_name(7:7),'(i1)') ord
        write(pert_ord,'(i1)') ord
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = trim(formname)
        labels(2) = trim(formname2)
        labels(3) = lresp_name
        call add_target(trim(formname),ttype_frm,.false.,tgt_info)
        call set_dependency(trim(formname),trim(formname2),tgt_info)
        call set_dependency(trim(formname),lresp_name,tgt_info)
        allocate(ifreq(ord))
        do freq_idx = 1,ord
          ifreq(freq_idx) = freq_idx
        end do
        call form_parameters2(-1,
     &       parameters,2,'lagrangian of pert order '//pert_ord//
     &       ' with freqs',
     &       ord,ifreq)
        deallocate(ifreq)
        call set_rule(trim(formname),ttype_frm,EXTRACT_FREQ,
     &                labels,3,1,
     &                parameters,2,tgt_info)
c
c        call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c        call set_rule(trim(formname),ttype_frm,PRINT_FORMULA,
c     &                labels,3,1,parameters,2,tgt_info)
c
      end do

      ! extract terms of RESFS_LAG(n)_X with correct freq. pattern
      ! following (2n+1) and (2n+2) rules
      formname(1:len_short) = ' '
      formname(1:6) = 'R(n)SX'
      formname2(1:len_short) = ' '
      formname2(1:14) = 'RESFS_LAG(n)_X'
      opname(1:len_short) = ' '
      opname(1:6) = 'O(n)SX'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        else
          op_parent = 'L'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        end if
        formname(6:6) = op_parent
        formname2(14:14) = op_parent
        opname(6:6) = op_parent
        do ord = op_par-1,x_max_ord
          write(formname(3:3),'(i1)') ord
          write(formname2(11:11),'(i1)') ord
          write(opname(3:3),'(i1)') ord
          write(pert_ord,'(i1)') ord
          do side = 1,2
            if (side.eq.1) then
              leftright = 'L'
            else
              leftright = 'R'
            end if
            formname(5:5) = leftright
            formname2(5:5) = leftright
            opname(5:5) = leftright
            allocate(ifreq(ord))
            do digit = 1,ord
              ifreq(digit) = digit
            end do
            comb_ok = .true.
            do while (comb_ok)
              do digit = 1,ord
                write(formname(6+digit:6+digit),'(i1)') ifreq(digit)
                write(opname(6+digit:6+digit),'(i1)') ifreq(digit)
              end do
              labels(1:20)(1:len_target_name) = ' '
              labels(1) = trim(formname)
              labels(2) = trim(formname2)
              labels(3) = trim(opname)
              call add_target(trim(formname),ttype_frm,.false.,tgt_info)
              call set_dependency(trim(formname),trim(formname2),
     &                      tgt_info)
              call set_dependency(trim(formname),trim(opname),tgt_info)
              call form_parameters2(-1,
     &             parameters,2,'residual of pert order '//pert_ord//
     &             ' with freqs',
     &             ord,ifreq)
              call set_rule(trim(formname),ttype_frm,EXTRACT_FREQ,
     &                      labels,3,1,
     &                      parameters,2,tgt_info)
c
c        call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c        call set_rule(trim(formname),ttype_frm,PRINT_FORMULA,
c     &                labels,3,1,parameters,2,tgt_info)
c
              comb_ok = next_comb(ifreq,ord,maxord)
            end do
            deallocate(ifreq)
          end do
        end do
      end do

      ! derivatives dLAG(n)/dL(0) and dLAG(n)/dT(0)
      ! following (2n+1) and (2n+2) rules regarding the maximum order
      lagf_name = 'RESP_LAGF(n)'
      res_lag_name = 'RES_LAG(n)_X'
      res_name = 'O(n)_X'
      op_name = 'X(0)'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        else
          op_parent = 'L'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        end if
        res_lag_name(12:12) = op_parent
        res_name(6:6) = op_parent
        op_name(1:1) = op_parent
        do ord=0,x_max_ord
          write(lagf_name(11:11),'(i1)') ord
          write(res_lag_name(9:9),'(i1)') ord
          write(res_name(3:3),'(i1)') ord
          labels(1:20)(1:len_target_name)= ' '
          labels(1) = res_lag_name
          labels(2) = lagf_name
          labels(3) = res_name
          labels(4) = op_name
          labels(5) = ' '
          call add_target(res_lag_name,ttype_frm,.false.,tgt_info)
          call set_dependency(res_lag_name,lagf_name,tgt_info)
          call set_dependency(res_lag_name,res_name,tgt_info)
          call set_dependency(res_lag_name,op_name,tgt_info)
          call form_parameters(-1,
     &         parameters,2,'residual',1,'---')
          call set_rule(res_lag_name,ttype_frm,DERIVATIVE,
     &                  labels,5,1,
     &                  parameters,2,tgt_info)
        end do
      end do

      ! split residuals in left (containing Y(n)) and right part
      ! but not for Y(n)=T(0), because non-linear CC-equations
      ! following (2n+1) and (2n+2) rules regarding the maximum order
      resl_lag_name = 'RESL_LAG(n)_X'
      resr_lag_name = 'RESR_LAG(n)_X'
      res_lag_name = 'RES_LAG(n)_X'
      opname(1:len_short) = ' '
      opname(1:6) = 'O(0)LT'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          op_name = 'L(n)'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        else
          op_parent = 'L'
          op_name = 'T(n)'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        end if
        resl_lag_name(13:13) = op_parent
        resr_lag_name(13:13) = op_parent
        res_lag_name(12:12) = op_parent
        opname(6:6) = op_parent
        do ord = op_par-1,x_max_ord
          write(op_name(3:3),'(i1)') ord
          write(resl_lag_name(10:10),'(i1)') ord
          write(resr_lag_name(10:10),'(i1)') ord
          write(res_lag_name(9:9),'(i1)') ord
          write(opname(3:3),'(i1)') ord
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = resl_lag_name
          labels(2) = resr_lag_name
          labels(3) = res_lag_name
          labels(4) = trim(opname)
          opname(5:5) = 'R'
          labels(5) = opname
          labels(6) = op_name
          call add_target(resl_lag_name,ttype_frm,.false.,tgt_info)
          call add_target(resr_lag_name,ttype_frm,.false.,tgt_info)
          call set_joined_targets(resl_lag_name,resr_lag_name,tgt_info)
          call set_dependency(resl_lag_name,res_lag_name,tgt_info)
          call set_dependency(resl_lag_name,trim(opname),tgt_info)
          opname(5:5) = 'L'
          call set_dependency(resl_lag_name,trim(opname),tgt_info)
          call set_dependency(resl_lag_name,op_name,tgt_info)
          call set_rule(resl_lag_name,ttype_frm,LEQ_SPLIT,
     &                  labels,6,2,
     &                  'split residual',1,tgt_info)
        end do
      end do

      ! FORM_CCHHAT
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_cchhat
      labels(2) = op_hhat
      labels(3) = op_ham
      labels(4) = 'T(0)'
      call add_target(form_cchhat,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cchhat,op_hhat,tgt_info)
      call set_dependency(form_cchhat,op_ham,tgt_info)
      call set_dependency(form_cchhat,'T(0)',tgt_info)
      call set_rule(form_cchhat,ttype_frm,DEF_HHAT,
     &              labels,4,1,
     &              title_cchhat,1,tgt_info)

      ! FORM_CCHBAR
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_cchbar
      labels(2) = op_hbar
      labels(3) = op_ham
      labels(4) = 'T(0)'
      call add_target(form_cchbar,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cchbar,op_hbar,tgt_info)
      call set_dependency(form_cchbar,op_ham,tgt_info)
      call set_dependency(form_cchbar,'T(0)',tgt_info)
      call set_rule(form_cchbar,ttype_frm,DEF_CC_HBAR,
     &              labels,4,1,
     &              title_cchbar,1,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! OPT_Y(n)w: optimized formula for determination of Y(n)w
      ! depends on DEF_ME_H(0), DEF_ME_V(1)v (if n>0),
      !            DEF_ME_Y(n)w, DEF_ME_O(n)SXw (S=L,R)
      ! not for OPT_T(0) because of non-linear CC equations
      ! following (2n+1) and (2n+2) rules regarding the maximum order
      formname(1:len_short) = ' '
      formname(1:6) = 'R(n)LX'
      formname2(1:len_short) = ' '
      formname2(1:6) = 'R(n)RX'
      defmelname(1:len_short) = ' '
      defmelname(1:13) = 'DEF_ME_O(n)SX'
      defmelname2(1:len_short) = ' '
      optname(1:len_short) = ' '
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          optname(1:8) = 'OPT_L(n)'
          defmelname2(1:11) = 'DEF_ME_L(n)'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        else
          op_parent = 'L'
          optname(1:8) = 'OPT_T(n)'
          defmelname2(1:11) = 'DEF_ME_T(n)'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        end if
        formname(6:6) = op_parent
        formname2(6:6) = op_parent
        defmelname(13:13) = op_parent
        do ord = op_par-1,x_max_ord
          write(optname(7:7),'(i1)') ord
          optname(9:len_short) = ' '
          write(formname(3:3),'(i1)') ord
          formname(7:len_short) = ' '
          write(formname2(3:3),'(i1)') ord
          formname2(7:len_short) = ' '
          write(defmelname2(10:10),'(i1)') ord
          defmelname2(12:len_short) = ' '
          write(defmelname(10:10),'(i1)') ord
          defmelname(14:len_short) = ' '
          allocate(ifreq(ord))
          do digit = 1,ord
            ifreq(digit) = digit
          end do
          comb_ok = .true.
          do while (comb_ok)
            do digit = 1,ord
              write(optname(8+digit:8+digit),'(i1)') ifreq(digit)
              write(formname(6+digit:6+digit),'(i1)') ifreq(digit)
              write(formname2(6+digit:6+digit),'(i1)') ifreq(digit)
              write(defmelname(13+digit:13+digit),'(i1)') ifreq(digit)
              write(defmelname2(11+digit:11+digit),'(i1)') ifreq(digit)
            end do

            labels(1:20)(1:len_target_name) = ' '
            labels(1) = trim(optname)
            labels(2) = trim(formname)
            labels(3) = trim(formname2)
            ncat = 2  ! 2 formulae pasted into final formula
            nint = 0  ! no intermediate to factor out so far ...
            call add_target(trim(optname),ttype_frm,.false.,tgt_info)
            call get_argument_value('calculate.routes','simtraf',
     &                              ival=isim)
            if (isim.eq.1) then
              nint = 1
              call set_dependency(trim(optname),form_cchhat,tgt_info)
              call set_dependency(trim(optname),mel_hhatdef,tgt_info)
              labels(4) = form_cchhat
            else if (isim.eq.2) then
              nint = 1
              call set_dependency(trim(optname),form_cchbar,tgt_info)
              call set_dependency(trim(optname),meldef_hbar,tgt_info)
              labels(4) = form_cchbar
            end if
            call set_dependency(trim(optname),trim(formname),tgt_info)
            call set_dependency(trim(optname),trim(formname2),tgt_info)
            call set_dependency(trim(optname),mel_ham,tgt_info)
            if (ord.gt.0)
     &         call set_dependency(trim(optname),'DEF_ME_V(1)',
     &                             tgt_info)
            call set_dependency(trim(optname),trim(defmelname2),
     &                          tgt_info)
            defmelname(12:12) = 'L'
            call set_dependency(trim(optname),trim(defmelname),tgt_info)
            defmelname(12:12) = 'R'
            call set_dependency(trim(optname),trim(defmelname),tgt_info)
            call opt_parameters(-1,parameters,ncat,nint)
            call set_rule(trim(optname),ttype_frm,OPTIMIZE,
     &                    labels,ncat+nint+1,1,
     &                    parameters,1,tgt_info)

            comb_ok = next_comb(ifreq,ord,maxord)
          end do
          deallocate(ifreq)

        end do
      end do

      ! OPT_T(0): for non-linear CC-equations
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'OPT_T(0)'
      labels(2) = 'RES_LAG(0)_L'
      labels(3) = 'RESP_LAG(0)'
      ncat = 2  ! 2 formulae pasted into final formula
      nint = 0  ! no intermediate to factor out so far ...
      call add_target('OPT_T(0)',ttype_frm,.false.,tgt_info)
      call set_dependency('OPT_T(0)','RES_LAG(0)_L',tgt_info)
      call set_dependency('OPT_T(0)','RESP_LAG(0)',tgt_info)
      call set_dependency('OPT_T(0)','DEF_ME_LRESP(0)',tgt_info)
      call set_dependency('OPT_T(0)','DEF_ME_O(0)_L',tgt_info)
      call set_dependency('OPT_T(0)','DEF_ME_T(0)',tgt_info)
      call set_dependency('OPT_T(0)',mel_ham,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('OPT_T(0)',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! OPT_LRESP(n): for evaluation of LRESP(n)
      opt_l_name = 'OPT_LRESP(n)'
      def_me_l_name = 'DEF_ME_LRESP(n)'
      formname(1:len_short) = ' '
      formname(1:11) = 'RESP_LAG(0)'
      do ord = 0,maxord
        write(opt_l_name(11:11),'(i1)') ord
        write(def_me_l_name(14:14),'(i1)') ord
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = opt_l_name
        labels(2) = trim(formname)
        ncat = 1
        nint = 0
        call add_target(opt_l_name,ttype_frm,.false.,tgt_info)
        call set_dependency(opt_l_name,trim(formname),tgt_info)
        call set_dependency(opt_l_name,def_me_l_name,tgt_info)
        call opt_parameters(-1,parameters,ncat,nint)
        call set_rule(opt_l_name,ttype_frm,OPTIMIZE,
     &                labels,ncat+nint+1,1,
     &                parameters,1,tgt_info)
        formname(1:11) = 'LAG(n)     '
        write(formname(5:5),'(i1)') ord+1
      end do

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! for H(0) (op_ham), a ME-list is already defined and imported

      ! ME_V(1)1, ME_V(1)2,...
      call get_argument_value('calculate.experimental','pert',
     &                        str=pert)
      call get_argument_value('calculate.experimental','pert_sym',
     &                        ival=isym)
      opname(1:len_short) = ' '
      opname(1:4) = 'V(1)'
      melname(1:len_short) = ' '
      melname(1:7) = 'ME_V(1)'
      defmelname(1:len_short) = ' '
      defmelname(1:11) = 'DEF_ME_V(1)'
c      do freq_idx = 1,maxord
c        write(opname(5:5),'(i1)') freq_idx
c        write(melname(8:8),'(i1)') freq_idx
c        write(defmelname(12:12),'(i1)') freq_idx
        call add_target(trim(defmelname),ttype_opme,.false.,tgt_info)
        call set_dependency(trim(defmelname),trim(opname),tgt_info)
        ! (a) define
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = trim(melname)
        labels(2) = trim(opname)
        call me_list_parameters(-1,parameters,
     &       msc,0,isym,0,0,.false.)
        call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        ! (b) import
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = trim(melname)
        call import_parameters(-1,parameters,pert,'DALTON_SPECIAL')
        call set_rule(trim(defmelname),ttype_opme,IMPORT,
     &                labels,1,1,
     &                parameters,1,tgt_info)
c      end do

      ! ME_Y(n)1,..., ME_O(n)SX1,... (X=T,L; n=0,..,x_max_ord, S=L,R)
      ! not for ME_O(0)SL because of non-linear CC-equations
      defmelname(1:len_short) = ' '
      defmelname(1:11) = 'DEF_ME_Y(n)'
      melname(1:len_short) = ' '
      melname(1:7) = 'ME_Y(n)'
      opname(1:len_short) = ' '
      opname(1:4) = 'Y(n)'
      defmelname2(1:len_short) = ' '
      defmelname2(1:13) = 'DEF_ME_O(n)SX'
      melname2(1:len_short) = ' '
      melname2(1:9) = 'ME_O(n)SX'
      opname2(1:len_short) = ' '
      opname2(1:6) = 'O(n)SX'
      hubname(1:len_short) = ' '
      hubname(1:14) = 'HUB_DEFME_Y(n)'
      do op_par = 1,2
        if (op_par.eq.1) then
          opname2(6:6) = 'T'
          opname(1:1) = 'L'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        else
          opname2(6:6) = 'L'
          opname(1:1) = 'T'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        end if
        do ord = 0,x_max_ord
          write(opname(3:3),'(i1)') ord
          opname(5:len_short) = ' '
          write(opname2(3:3),'(i1)') ord
          opname2(7:len_short) = ' '
          allocate(ifreq(ord))
          do digit = 1,ord
            ifreq(digit) = digit
          end do
          comb_ok = .true.
          do while (comb_ok)
            do digit = 1,ord
              write(opname(4+digit:4+digit),'(i1)') ifreq(digit)
            end do
            melname(4:len_short) = opname(1:len_short-3)
            defmelname(8:len_short) = opname(1:len_short-7)
            call add_target(trim(defmelname),ttype_opme,.false.,
     &                      tgt_info)
            call set_dependency(trim(defmelname),trim(opname),tgt_info)
            ! (a) T(n) depends on T(n-1), L(n) depends on both L(n-1) and T(n)
            hubname(11:14) = opname(1:4)
            if (op_name(1:1).eq.'L') then
              hubname(11:11) = 'T'
              call set_dependency(trim(defmelname),trim(hubname),
     &                            tgt_info)
              hubname(11:14) = opname(1:4)
            end if
            if (ord.gt.0) then
              write(hubname(13:13),'(i1)') ord-1
              call set_dependency(trim(defmelname),trim(hubname),
     &                            tgt_info)
            end if
            ! (b) ME_X(n)
            labels(1:20)(1:len_target_name) = ' '
            labels(1) = trim(melname)
            labels(2) = trim(opname)
            call me_list_parameters(-1,parameters,
     &           msc,0,1,0,0,.false.)
            call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                    labels,2,1,
     &                    parameters,1,tgt_info)
            comb_ok = next_comb(ifreq,ord,maxord)
          end do
          deallocate(ifreq)
          ! (c) ME_OS(n)_X (S=L,R; X=T,L; n=0,...,maxord)
          do side = 1,2
            if (side.eq.1) then
              leftright = 'L'
            else
              leftright = 'R'
            end if
            opname2(5:5) = leftright
            allocate(ifreq(ord))
            do digit = 1,ord
              ifreq(digit) = digit
            end do
            comb_ok = .true.
            do while (comb_ok)
              do digit = 1,ord
                write(opname2(6+digit:6+digit),'(i1)') ifreq(digit)
              end do
              melname2(4:len_short) = opname2(1:len_short-3)
              defmelname2(8:len_short) = opname2(1:len_short-7)
              if (ord.ge.op_par-1) then
                call add_target(trim(defmelname2),
     &                          ttype_opme,.false.,tgt_info)
                call set_dependency(trim(defmelname2),trim(opname2),
     &                              tgt_info)
                labels(1:20)(1:len_target_name) = ' '
                labels(1) = trim(melname2)
                labels(2) = trim(opname2)
                call me_list_parameters(-1,parameters,
     &               msc,0,1,0,0,.false.)
                call set_rule(trim(defmelname2),ttype_opme,DEF_ME_LIST,
     &                        labels,2,1,
     &                        parameters,1,tgt_info)
              end if
              comb_ok = next_comb(ifreq,ord,maxord)
            end do
            deallocate(ifreq)
          end do
        end do
      end do

      ! ME_O(0)_L
      call add_target('DEF_ME_O(0)_L',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_O(0)_L','O(0)_L',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_O(0)_L'
      labels(2) = 'O(0)_L'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF_ME_O(0)_L',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! ME_LRESP(n)
      def_me_l_name = 'DEF_ME_LRESP(n)'
      me_l_name = 'ME_LRESP(n)'
      lresp_name = 'LRESP(n)'
      do ord = 0,maxord
        write(def_me_l_name(14:14),'(i1)') ord
        write(me_l_name(10:10),'(i1)') ord
        write(lresp_name(7:7),'(i1)') ord
        call add_target(def_me_l_name,ttype_opme,.false.,tgt_info)
        call set_dependency(def_me_l_name,lresp_name,tgt_info)
        if (ord.gt.0) call set_dependency(def_me_l_name,
     &                     'DEF_ME_V(1)',tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = me_l_name
        labels(2) = lresp_name
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(def_me_l_name,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
      end do

      ! Preconditioner
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
      call add_target(mel_dia1,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_dia1,mel_ham,tgt_info)
      call set_dependency(mel_dia1,op_dia,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &     'T(0)',.false.)
c      call set_rule('DIAG',ttype_op,CLONE_OP,
c     &     'D2',1,1,
c     &     parameters,1,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_dia1
      labels(2) = op_dia
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_dia1,ttype_opme,DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
      labels(1) = mel_dia1
      labels(2) = mel_ham
      call set_rule(mel_dia1,ttype_opme,PRECONDITIONER,
     &              labels,2,1,
     &              parameters,0,tgt_info)

      ! HHat definition
      call add_target(mel_hhatdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_hhatdef,op_hhat,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_hhat
      labels(2) = op_hhat
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_hhatdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! Hbar definition
      call add_target(meldef_hbar,ttype_opme,.false.,tgt_info)
      call set_dependency(meldef_hbar,op_hbar,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_hbar
      labels(2) = op_hbar
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(meldef_hbar,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*
  
      ! solve CC-equations 
      call add_target('SOLVE_T(0)',ttype_gen,.false.,tgt_info)
      call set_dependency('SOLVE_T(0)','OPT_T(0)',tgt_info)
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
      call set_dependency('SOLVE_T(0)',mel_dia1,tgt_info)
      call solve_parameters(-1,parameters,2, 1,1,'DIA')
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_T(0)'
      labels(2) = 'ME_O(0)_L'
      labels(3) = mel_dia1
      labels(4) = 'ME_LRESP(0)'
      labels(5) = 'OPT_T(0)'
      call set_rule('SOLVE_T(0)',ttype_opme,SOLVENLEQ,
     &     labels,5,2,
     &     parameters,2,tgt_info)

      ! solve linear equations for X(n) (X=T,L, n=0,...,x_max_ord, not for T(0))
      if (maxord.gt.0) then
        freq(1) = 0d0
        call get_argument_value('calculate.experimental','freq',
     &                          xarr=freq(1:maxord-1))
        freq(maxord) = -sum(freq)
      end if
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
      call solve_parameters(-1,parameters,2, 1,1,'DIA')
      solvename(1:len_short) = ' '
      defmelname(1:len_short) = ' '
      opname(1:len_short) = ' '
      melname(1:len_short) = ' '
      optname(1:len_short) = ' '
      do op_par = 1,2
        if (op_par.eq.1) then
          defmelname(1:11) = 'DEF_ME_L(n)'
          melname(1:7) = 'ME_L(n)'
          opname(1:6) = 'O(n)LT'
          optname(1:8) = 'OPT_L(n)'
          solvename(1:10) = 'SOLVE_L(n)'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        else
          defmelname(1:11) = 'DEF_ME_T(n)'
          melname(1:7) = 'ME_T(n)'
          opname(1:6) = 'O(n)LL'
          optname(1:8) = 'OPT_T(n)'
          solvename(1:10) = 'SOLVE_T(n)'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        end if
        do ord = op_par-1,x_max_ord
          write(solvename(9:9),'(i1)') ord
          write(defmelname(10:10),'(i1)') ord
          write(melname(6:6),'(i1)') ord
          write(opname(3:3),'(i1)') ord
          write(optname(7:7),'(i1)') ord
          allocate(ifreq(ord))
          do digit = 1,ord
            ifreq(digit) = digit
          end do
          comb_ok = .true.
          do while (comb_ok)
            do digit = 1,ord
              write(solvename(10+digit:10+digit),'(i1)') ifreq(digit)
              write(defmelname(11+digit:11+digit),'(i1)') ifreq(digit)
              write(optname(8+digit:8+digit),'(i1)') ifreq(digit)
              write(opname(6+digit:6+digit),'(i1)') ifreq(digit)
              write(melname(7+digit:7+digit),'(i1)') ifreq(digit)
            end do
            call add_target(trim(solvename),ttype_gen,.false.,tgt_info)
            call set_dependency(trim(solvename),trim(defmelname),
     &                          tgt_info)
            call set_dependency(trim(solvename),trim(optname),tgt_info)
            call set_dependency(trim(solvename),mel_dia1,tgt_info)
            ! SOLVE_T(n) depends on SOLVE_T(n-1),
            ! SOLVE_L(n) depends on SOLVE_L(n-1) and SOLVE_T(n)
            hubname(1:14) = 'HUB_SOLVE_Y(n)'
            hubname(5:14) = solvename(1:10)
            if (op_par.eq.1) then
              hubname(11:11) = 'T'
              call set_dependency(trim(solvename),hubname(1:14),
     &                            tgt_info)
              hubname(5:14) = solvename(1:10)
            end if
            if (ord.gt.0) then
              write(hubname(13:13),'(i1)') ord-1
              call set_dependency(trim(solvename),hubname(1:14),
     &                            tgt_info)
            end if
            ! SOLVE_Y(n)
            labels(1:20)(1:len_target_name) = ' '
            labels(1) = trim(melname)
            labels(2) = mel_dia1
            labels(3) = trim(opname)
            opname(5:5) = 'R'
            labels(4) = trim(opname)
            opname(5:5) = 'L'
            labels(5) = trim(optname)
            call set_rule(trim(solvename),ttype_opme,SOLVELEQ,
     &           labels,5,1,
     &           parameters,2,tgt_info)
            comb_ok = next_comb(ifreq,ord,maxord)
          end do
          deallocate(ifreq)
        end do
      end do

      ! evaluate RESP_LAG(n)
      opt_l_name = 'OPT_LRESP(n)'
      eval_lag_name = 'EVAL_LAG(n)'
      hubname(1:14) = 'HUB_SOLVE_X(n)'
      def_me_l_name = 'DEF_ME_LRESP(n)'
      do ord = 0,maxord
        write(opt_l_name(11:11),'(i1)') ord
        write(eval_lag_name(10:10),'(i1)') ord
        write(def_me_l_name(14:14),'(i1)') ord
        call add_target(eval_lag_name,ttype_gen,evaluate(ord),tgt_info)
        ! first solve T(k), L(k) equations according to (2n+1) and (2n+2) rules
        hubname(11:11) = 'T'
        x_max_ord = int((real(ord)-1)/2+0.6)
        write(hubname(13:13),'(i1)') x_max_ord
        call set_dependency(eval_lag_name,hubname(1:14),tgt_info)
        hubname(11:11) = 'L'
        x_max_ord = int((real(ord)-2)/2+0.6)
        write(hubname(13:13),'(i1)') x_max_ord
        call set_dependency(eval_lag_name,hubname(1:14),tgt_info)
        ! evaluate
        call set_dependency(eval_lag_name,opt_l_name,tgt_info)
        call set_dependency(eval_lag_name,def_me_l_name,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = opt_l_name
        call set_rule(eval_lag_name,ttype_opme,EVAL,
     &       labels,1,0,
     &       parameters,0,tgt_info)
      end do

      ! hub for DEF_ME_Y(n)-dependencies
      defmelname(1:len_short) = ' '
      defmelname(1:11) = 'DEF_ME_Y(n)'
      hubname(1:len_short) = ' '
      hubname(1:14) = 'HUB_DEFME_Y(n)'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        else
          op_parent = 'L'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        end if
        defmelname(8:8) = op_parent
        hubname(11:11) = op_parent
        do ord = 0,x_max_ord
          write(defmelname(10:10),'(i1)') ord
          defmelname(12:len_short) = ' '
          write(hubname(13:13),'(i1)') ord
          call add_target(trim(hubname),ttype_gen,.false.,tgt_info)
          allocate(ifreq(ord))
          do digit = 1,ord
            ifreq(digit) = digit
          end do
          comb_ok = .true.
          do while (comb_ok)
            do digit = 1,ord
              write(defmelname(11+digit:11+digit),'(i1)') ifreq(digit)
            end do
            call set_dependency(trim(hubname),trim(defmelname),tgt_info)
            comb_ok = next_comb(ifreq,ord,maxord)
          end do
          deallocate(ifreq)
        end do
      end do   

c      ! hub for DEF_ME_V(1)-dependencies
c      defmelname(1:len_short) = ' '
c      defmelname(1:11) = 'DEF_ME_V(1)'
c      call add_target('HUB_DEFME_V(1)',ttype_gen,.false.,tgt_info)
c      do freq_idx = 1,maxord
c        write(defmelname(12:12),'(i1)') freq_idx
c        call set_dependency('HUB_DEFME_V(1)',trim(defmelname),tgt_info)
c      end do

      ! hub for SOLVE_Y(n)-dependencies
      solvename(1:len_short) = ' '
      solvename(1:10) = 'SOLVE_Y(n)'
      hubname(1:len_short) = ' '
      hubname(1:14) = 'HUB_SOLVE_Y(n)'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          x_max_ord = int((real(maxord)-1)/2+0.6)
        else
          op_parent = 'L'
          x_max_ord = int((real(maxord)-2)/2+0.6)
        end if
        solvename(7:7) = op_parent
        hubname(11:11) = op_parent
        do ord = 0,x_max_ord
          write(solvename(9:9),'(i1)') ord
          solvename(11:len_short) = ' '
          write(hubname(13:13),'(i1)') ord
          call add_target(trim(hubname),ttype_gen,.false.,tgt_info)
          allocate(ifreq(ord))
          do digit = 1,ord
            ifreq(digit) = digit
          end do
          comb_ok = .true.
          do while (comb_ok)
            do digit = 1,ord
              write(solvename(10+digit:10+digit),'(i1)') ifreq(digit)
            end do
            call set_dependency(trim(hubname),trim(solvename),tgt_info)
            comb_ok = next_comb(ifreq,ord,maxord)
          end do
          deallocate(ifreq)
        end do
      end do

      return
      end
