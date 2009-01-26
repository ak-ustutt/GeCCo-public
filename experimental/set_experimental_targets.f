*----------------------------------------------------------------------*
      subroutine set_experimental_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     set targets for response theory calculations
*
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'
      include 'mdef_me_list.h'
      include 'multd2h.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      character(*), intent(in) ::
     &     env_type

      integer, parameter ::
     &     len_short = 32, len_long = 256, maxsym = 8,
     &     maximum_order = 8, ntest = 100

      integer ::
     &     isim, ncat, nint, ndef, occ_def(ngastp,2,60),
     &     msc, sym, r12op, ansatz, maxexc, t1ext_mode,
     &     ord, op_par, len_op_exp, side, x_max_ord,
     &     freq_idx, digit, ilabels, ord2, op_par2, x_max_ord2,
     &     pertdir, freq_idxnew, ncnt, icnt, pos, maxexc_cur,icnt2

      character(len_target_name) ::
     &     mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(2)

      character(len_short) ::
     &     opname, melname, defmelname, formname, formname2,
     &     opname2, melname2, defmelname2, hubname, optname,
     &     solvename, evalname, lagname

      character(len_long) ::
     &     opexp, pert

      character ::
     &     op_name*4, op_parent*1, lagf_name*12, op_exp*50,
     &     pert_ord*1, lag_name*11, res_name*9, leftright*1,
     &     res_lag_name*15,resl_lag_name*16, resr_lag_name*16,
     &     lresp_name*8

      character(len=8) ::
     &     approx

      logical ::
     &     setr12, r12fix, set_zero, eval_dipmom(3)

      real(8) ::
     &     freqsum

      real(8), allocatable ::
     &     freq(:,:)

      integer, allocatable ::
     &     ifreq(:), isym(:), redun(:), ifreqnew(:), maxord(:),
     &     subset(:)

      logical, allocatable ::
     &     evaluate(:)

      type(me_list_array), allocatable ::
     &     mel_dummy(:)

      logical, external ::
     &     next_comb

      integer, external ::
     &     pert_sym, idx_target

*----------------------------------------------------------------------*
*     process keywords and allocate arrays:
*----------------------------------------------------------------------*

      ! skip this section if not requested
      ncnt = is_keyword_set('calculate.experimental')
      if (ncnt.eq.0) return

      ! set r12 targets
      setr12 = is_keyword_set('method.R12').gt.0
      if (setr12) then
        call get_argument_value('method.R12','ansatz',ival=ansatz)
        call get_argument_value('method.R12','approx',str=approx)
        call get_argument_value('method.R12','fixed',lval=r12fix)
        call get_argument_value('method.R12','r12op',ival=r12op)
        call get_argument_value('method.R12','T1ext',ival=t1ext_mode)
        if (r12op.le.1.and.r12fix) then
          call set_r12f_general_targets(tgt_info,orb_info,env_type)
        else
          call quit(1,'set_experimental_targets','wrong r12 method')
        end if
      else
        t1ext_mode = 0
      end if

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting experimental targets ...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

      allocate(maxord(ncnt),freq(ncnt,maximum_order),evaluate(ncnt))
      do icnt = 1,ncnt
        call get_argument_value('calculate.experimental','order',
     &       keycount=icnt,ival=maxord(icnt))
      end do
      allocate(redun(ncnt*maxval(maxord)),isym(ncnt*maxval(maxord)))
      redun = 0
      isym = 0
      if (maxval(maxord).gt.maximum_order) then
        write(pert_ord,'(i1)') maximum_order
        call quit(1,'set_experimental_targets',
     &        'maxord must not exceed '//pert_ord)
      end if
      pert(1:len_long) = ' '
      do icnt = 1,ncnt
        pos = (icnt-1)*maxval(maxord) + 1
        call get_argument_value('calculate.experimental','pert',
     &       keycount=icnt,str=pert(pos:len_long))

        ! duplicate values for pert if necessary and not specified
        do freq_idx = 2,maxord(icnt)
          pos = (icnt-1)*maxval(maxord) + freq_idx
          if (pert(pos:pos).eq.' ')
     &        pert(pos:pos) = pert(pos-1:pos-1)
        end do

        ! check if pert input is ok
        pos = (icnt-1)*maxval(maxord)
        do freq_idx = pos+1,pos+maxord(icnt)
          if (pert(freq_idx:freq_idx).ne.'X' .and.
     &          pert(freq_idx:freq_idx).ne.'Y' .and.
     &          pert(freq_idx:freq_idx).ne.'Z')
     &          call quit(1,'set_experimental_targets',
     &          'pert must contain X,Y,Z')

          ! determine irreps of perturbation operators
          isym(freq_idx) = pert_sym(pert(freq_idx:freq_idx),orb_info)
        end do
      end do

      ! determine redundancies
      if (maxval(maxord).gt.0) then

        do icnt = 1,ncnt

          ! get user defined frequency array, sum of frequencies is zero
          freq(icnt,1:maximum_order) = 0d0
          call get_argument_value('calculate.experimental','freq',
     &         keycount=icnt,xarr=freq(icnt,1:maximum_order))
          freq(icnt,maxord(icnt):maximum_order) = 0d0
          freq(icnt,maxord(icnt)) = -sum(freq(icnt,:))

          ! determine redundancies
          do digit = 1,maxord(icnt)
            pos = (icnt-1)*maxval(maxord) + digit
            do freq_idx = 1,pos-1
              icnt2 = (freq_idx-1)/maxval(maxord)+1
              if (pert(pos:pos).eq.pert(freq_idx:freq_idx) .and.
     &           abs(freq(icnt,digit) - freq(icnt2,freq_idx-(icnt2-1) *
     &           maxval(maxord))).lt.1d-12) redun(pos) = freq_idx
            end do
            if (redun(pos).eq.0) redun(pos) = pos
          end do
        end do
      end if

      evaluate = .true.

      ! maximum excitation
      maxexc = -1
      do icnt = 1,ncnt
        call get_argument_value('calculate.experimental','maxexc',
     &       keycount=icnt,ival=maxexc_cur)
        if (icnt.eq.1) then
          maxexc = maxexc_cur
        else if (maxexc.ne.maxexc_cur) then
          call quit(1,'set_experimental_targets',
     &                'maxexc must be the same for all requested props')
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'keywords processed:'
        write(luout,*) 'ncnt: ',ncnt
        write(luout,*) 'maxord: ',maxord
        write(luout,*) 'pert: ',pert(1:ncnt*maxval(maxord))
        write(luout,*) 'freq: ',freq
        write(luout,*) 'redun: ',redun
      end if

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! define op_ham as H(0)
      call ord_parameters(-1,parameters,0,3,-1)
      call set_rule(op_ham,ttype_op,SET_ORDER,op_ham,
     &              1,1,parameters,1,tgt_info)

      if (setr12) then
        ! set order of R12, V, B, Bh, X, R12-INT, C-INT, Xh to zero
        call ord_parameters(-1,parameters,0,3,-1)
        call set_rule(op_r12,ttype_op,SET_ORDER,op_r12,
     &                1,1,parameters,1,tgt_info)
        call set_rule(op_v_inter,ttype_op,SET_ORDER,op_v_inter,
     &                1,1,parameters,1,tgt_info)
        call set_rule(op_b_inter,ttype_op,SET_ORDER,op_b_inter,
     &                1,1,parameters,1,tgt_info)
        call set_rule(op_bh_inter,ttype_op,SET_ORDER,op_bh_inter,
     &                1,1,parameters,1,tgt_info)
        call set_rule(op_xh_inter,ttype_op,SET_ORDER,op_xh_inter,
     &                1,1,parameters,1,tgt_info)
        call set_rule(op_x_inter,ttype_op,SET_ORDER,op_x_inter,
     &                1,1,parameters,1,tgt_info)
        call set_rule(op_rint,ttype_op,SET_ORDER,op_rint,
     &                1,1,parameters,1,tgt_info)
        call set_rule(op_c_inter,ttype_op,SET_ORDER,op_c_inter,
     &                1,1,parameters,1,tgt_info)

      end if

      ! define V(1)
      call add_target('V(1)',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,1,1,setr12)
      call set_rule('V(1)',ttype_op,DEF_HAMILTONIAN,'V(1)',
     &              1,1,parameters,1,tgt_info)
      call ord_parameters(-1,parameters,1,3,-1)
      call set_rule('V(1)',ttype_op,SET_ORDER,'V(1)',
     &              1,1,parameters,1,tgt_info)

      ! define frequency components V(1)1, V(1)2, ...,V(1)maxord
      opname(1:len_short) = ' '
      opname(1:4) = 'V(1)'
      freq_idx = 0
      do while (next_comb(freq_idx,1,maxord,ncnt))
        call redundant_comb(freq_idx,freq_idxnew,redun,1,maxord,ncnt)
        write(opname(5:6),'(i2.2)') freq_idxnew
        if (idx_target(trim(opname),tgt_info).gt.0) cycle
        call add_target(trim(opname),ttype_op,.false.,tgt_info)
        call hop_parameters(-1,parameters,0,1,3,setr12)
        call set_rule(trim(opname),ttype_op,DEF_HAMILTONIAN,
     &                trim(opname),
     &                1,1,parameters,1,tgt_info)
        call ord_parameters(-1,parameters,1,3,freq_idxnew)
        call set_rule(trim(opname),ttype_op,SET_ORDER,trim(opname),
     &                1,1,parameters,1,tgt_info)
      end do

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
      if (t1ext_mode.eq.0) then
        call xop_parameters(-1,parameters,
     &                   .false.,1,maxexc,0,1)
        call set_rule(op_top,ttype_op,DEF_EXCITATION,
     &              op_top,1,1,
     &              parameters,1,tgt_info)
      else
        occ_def = 0
        ! 1 -- T1
        occ_def(IPART,1,1) = 1
        occ_def(IHOLE,2,1) = 1
        ! 2 -- T1'
        occ_def(IEXTR,1,2) = 1
        occ_def(IHOLE,2,2) = 1
        ndef = maxexc + 1
        ! 3,4,... -- T2,T3,...
        do nint = 2,maxexc
          occ_def(IPART,1,nint+1) = nint
          occ_def(IHOLE,2,nint+1) = nint
        end do
        call op_from_occ_parameters(-1,parameters,2,
     &                occ_def,ndef,1,(/.true.,.true./),ndef)
        call set_rule(op_top,ttype_op,DEF_OP_FROM_OCC,
     &                op_top,1,1,
     &                parameters,2,tgt_info)
      end if

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
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        else
          op_parent = 'L'
          op_name = 'L(x)'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
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
            allocate(ifreq(ord), ifreqnew(ord))
            ifreq = 0
            do while (next_comb(ifreq,ord,maxord,ncnt))
              call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
              do digit = 1,ord
                write(opname(3+2*digit:4+2*digit),'(i2.2)') 
     &                  ifreqnew(digit)
              end do
              if (idx_target(trim(opname),tgt_info).gt.0) cycle
              call add_target(trim(opname),ttype_op,.false.,tgt_info)
              call set_dependency(trim(opname),op_parent,tgt_info)
              call cloneop_parameters(-1,parameters,op_parent,.false.)
              call set_rule(trim(opname),ttype_op,CLONE_OP,
     &                      trim(opname),1,1,
     &                      parameters,1,tgt_info)
              call ord_parameters(-1,parameters,ord,op_par,ifreqnew)
              call set_rule(trim(opname),ttype_op,SET_ORDER,
     &                      trim(opname),
     &                      1,1,parameters,1,tgt_info)
            end do
            deallocate(ifreq, ifreqnew)
          end if
        end do
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'defining residuals'

      ! define residuals O(n)_L(0) and O(n)_T(0)
      ! following (2n+1) and (2n+2) rules regarding the maximum order
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          res_name = 'O(x)_L'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        else
          op_parent = 'L'
          res_name = 'O(x)_T'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
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
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        else
          op_parent = 'T'
          opname(6:6) = 'L'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
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
            ! is this one necessary (without frequency indices)???
            if (ord.gt.0) then
              call add_target(trim(opname),ttype_op,.false.,tgt_info)
              call set_dependency(trim(opname),op_parent,tgt_info)
              call cloneop_parameters(-1,parameters,op_parent,.false.)
              call set_rule(trim(opname),ttype_op,CLONE_OP,
     &                      trim(opname),1,1,
     &                      parameters,1,tgt_info)
            end if
            allocate(ifreq(ord),ifreqnew(ord))
            ifreq = 0
            set_zero = ord.eq.0
            do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
              set_zero = .false.
              call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
              do digit = 1,ord
                write(opname(5+2*digit:6+2*digit),'(i2.2)') 
     &                    ifreqnew(digit)
              end do
              if (idx_target(trim(opname),tgt_info).gt.0) cycle
              call add_target(trim(opname),ttype_op,.false.,tgt_info)
              call set_dependency(trim(opname),op_parent,tgt_info)
              call cloneop_parameters(-1,parameters,op_parent,.false.)
              call set_rule(trim(opname),ttype_op,CLONE_OP,
     &                      trim(opname),1,1,
     &                      parameters,1,tgt_info)
            end do
            deallocate(ifreq, ifreqnew)
          end do
        end do
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'defining lagrangians'

      ! define scalar response lagrangian LRESP(n) of order n
      lagname(1:len_short) = ' '
      lagname(1:8) = 'LRESP(n)'
      do ord=0,maxval(maxord)
        write(lagname(7:7),'(i1)') ord
        lagname(9:len_short) = ' '
        ! is this one necessary (without frequency indices)???
        if (ord.gt.0) then
          call add_target(trim(lagname),ttype_op,.false.,tgt_info)
          call set_dependency(trim(lagname),'LRESP',tgt_info)
          call cloneop_parameters(-1,parameters,'LRESP',.false.)
          call set_rule(trim(lagname),ttype_op,CLONE_OP,
     &                  trim(lagname),1,1,
     &                  parameters,1,tgt_info)
        end if
        allocate(ifreq(ord),ifreqnew(ord))
        ifreq = 0
        set_zero = ord.eq.0
        do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
          set_zero = .false.
          call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
          do digit = 1,ord
            write(lagname(7+2*digit:8+2*digit),'(i2.2)') ifreqnew(digit)
          end do
          if (idx_target(trim(lagname),tgt_info).gt.0) cycle
          call add_target(trim(lagname),ttype_op,.false.,tgt_info)
          call set_dependency(trim(lagname),'LRESP',tgt_info)
          call cloneop_parameters(-1,parameters,'LRESP',.false.)
          call set_rule(trim(lagname),ttype_op,CLONE_OP,
     &                  trim(lagname),1,1,
     &                  parameters,1,tgt_info)
        end do
        deallocate(ifreq, ifreqnew)
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'defining other operators'

      ! Diagonal
      call add_target(op_dia,ttype_op,.false.,tgt_info)
      call set_dependency(op_dia,'T(0)',tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        'T(0)',.false.)
      call set_rule(op_dia,ttype_op,CLONE_OP,
     &              op_dia,1,1,
     &              parameters,1,tgt_info)
c dbg
      print *,'m1'
c dbg

      ! T1 transformed Hamiltonian
      call add_target(op_hhat,ttype_op,.false.,tgt_info)
      if (.not.setr12) then
        call set_dependency(op_hhat,op_ham,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_ham,.false.)
        call set_rule(op_hhat,ttype_op,CLONE_OP,
     &                op_hhat,1,1,
     &                parameters,1,tgt_info)
      else
c dbg
      print *,'m1b'
c dbg
        call define_r12_hhat(tgt_info)
      end if
c dbg
      print *,'m2'
c dbg

      ! Hbar intermediate
      call add_target(op_hbar,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &                    .false.,1,2,0,1)
      call set_rule(op_hbar,ttype_op,DEF_CC_HBAR_OP,
     &              op_hbar,1,1,
     &              parameters,1,tgt_info)

      if (ntest.ge.100)
     &  write(luout,*) 'operators defined'

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      if (.not.setr12) then
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
     &       parameters,2,'response lagrange functional',0,'---')
        call set_rule('RESP_LAG',ttype_frm,DEF_EXP_FORMULA,
     &                labels,5,1,
     &                parameters,2,tgt_info)
      else
        ! define r12 response lagrangian
        call add_target('RESP_LAG',ttype_frm,.false.,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = 'RESP_LAG'
        labels(2) = 'LRESP'
        labels(3) = 'Hnew'
        labels(4) = op_r12
        labels(5) = op_r12
        labels(6) = 'L'
        labels(7) = 'T'
        ilabels = 7
        if (r12op.eq.1) then
          labels(8) = op_cexbar
          labels(9) = op_cex
          call set_dependency('RESP_LAG',op_cex,tgt_info)
          call set_dependency('RESP_LAG',op_cexbar,tgt_info)
          ilabels = 9
        end if
        call set_dependency('RESP_LAG','LRESP',tgt_info)
        call set_dependency('RESP_LAG','L',tgt_info)
        call set_dependency('RESP_LAG','Hnew',tgt_info)
        call set_dependency('RESP_LAG','T',tgt_info)
        call set_dependency('RESP_LAG',op_r12,tgt_info)
        call form_parameters(-1,
     &       parameters,2,'r12 response lagrange functional',
     &                    ansatz,'---')
        call set_rule('RESP_LAG',ttype_frm,DEF_CCR12_LAGRANGIAN,
     &                labels,ilabels,1,
     &                parameters,2,tgt_info)
        call form_parameters(-1,parameters,2,'stdout',1,'stdout')
        call set_rule('RESP_LAG',ttype_frm,PRINT_FORMULA,
     &                labels,2,1,parameters,2,tgt_info)
      end if

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
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        else
          op_parent = 'L'
          op_exp(1:len_op_exp) = 'L=L(0)'
          op_name = 'L(x)'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
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

      if (ntest.ge.100)
     &  write(luout,*) 'frequency expansion of perturbation operators'

      ! frequency expansion of V(1): V(1)=V(1)1+V(1)2+...

      formname(1:len_short) = ' '
      formname(1:6) = 'V(1)_F'
      do ord=0,maxval(maxord)
        allocate(ifreq(ord),ifreqnew(ord))
        ifreq = 0
        set_zero = ord.eq.0
        do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
          set_zero = .false.
          call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
          do digit = 1,ord
            write(formname(5+2*digit:6+2*digit),'(i2.2)')
     &                  ifreqnew(digit)
          end do
          if (idx_target(trim(formname),tgt_info).gt.0) cycle      
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(formname)
          call add_target(trim(formname),ttype_frm,.false.,tgt_info)
          opexp(1:len_long) = ' '
          opexp(1:5) = 'V(1)='
          len_op_exp = 5
          opname(1:len_short) = ' '
          opname(1:4) = 'V(1)'
          comb_loop: do digit = 1,ord
            write(opname(5:6),'(i2.2)') ifreqnew(digit)
            do pos = 1,len_op_exp-5
              if (opexp(pos:pos+5).eq.opname(1:6)) cycle comb_loop
            end do
            opexp(len_op_exp+1:len_op_exp+7) = trim(opname)//'+'
            len_op_exp = len_op_exp + 7
            call set_dependency(trim(formname),trim(opname),tgt_info)
          end do comb_loop
          opexp(len_op_exp:len_op_exp) = ' '
          len_op_exp = len_op_exp - 1
          if (ord.eq.0) then
            opexp(5:9) = '=V(1)'
            len_op_exp = 9
          end if
          call def_form_parameters(-1,
     &         parameters,2,opexp(1:len_op_exp),
     &         'freq exp of V(1)')
          call set_rule(trim(formname),ttype_frm,DEF_FORMULA,
     &                  labels,1,1,
     &                  parameters,2,tgt_info)
        end do
        deallocate(ifreq,ifreqnew)
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'frequency expansion of L and T operators'

      ! frequency expansion of X(n), n>0: X(1) = X(1)1+X(1)2+..., ...
      ! following the (2n+1) and (2n+2) rules regarding the maximum order
      formname(1:len_short) = ' '
      formname(1:6) = 'X(n)_F'
      do ord=1,maxval(maxord)
        do op_par = 1,2
          if (op_par.eq.1) then
            op_name = 'T(x)'
            x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
          else
            op_name = 'L(x)'
            x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
          end if
          do ord2=1,min(x_max_ord,ord)
            write(op_name(3:3),'(i1)') ord2
            formname(1:4) = op_name(1:4)
            allocate(ifreq(ord),ifreqnew(ord))
            ifreq = 0
            do while (next_comb(ifreq,ord,maxord,ncnt))
              call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
              do digit = 1,ord
                write(formname(5+2*digit:6+2*digit),'(i2.2)')
     &                      ifreqnew(digit)
              end do
              if (idx_target(trim(formname),tgt_info).gt.0) cycle
              labels(1:20)(1:len_target_name) = ' '
              labels(1) = trim(formname)
              call add_target(trim(formname),ttype_frm,.false.,tgt_info)
              call set_dependency(trim(formname),op_name,tgt_info)
              opexp(1:len_long) = ' '
              opexp(1:5) = op_name//'='
              len_op_exp = 5
              opname(1:len_short) = ' '
              opname(1:4) = op_name
              allocate(subset(ord2))
              subset = 0
              comb_loop2: do while (next_comb(subset,ord2,ord,1))
                do digit = 1,ord2
                  write(opname(3+2*digit:4+2*digit),'(i2.2)') 
     &                        ifreqnew(subset(digit))
                end do
                do pos = 1,len_op_exp-3-2*ord2
                  if (opexp(pos:pos+3+2*ord2).eq.opname(1:4+2*ord2))
     &                                  cycle comb_loop2
                end do
                opexp(len_op_exp+1:len_op_exp+5+2*ord2) = 
     &                              trim(opname)//'+'
                len_op_exp = len_op_exp + 5 + 2*ord2
                call set_dependency(trim(formname),trim(opname),
     &                              tgt_info)
              end do comb_loop2
              deallocate(subset)
              opexp(len_op_exp:len_op_exp) = ' '
              len_op_exp = len_op_exp - 1
              call def_form_parameters(-1,
     &             parameters,2,opexp(1:len_op_exp),
     &             'freq exp of '//op_name)
              call set_rule(trim(formname),ttype_frm,DEF_FORMULA,
     &                      labels,1,1,
     &                      parameters,2,tgt_info)
            end do
            deallocate(ifreq,ifreqnew)
          end do
        end do
      end do

      ! expand response lagrangian with Hnew=H+V(1)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'RESP_LAGF'
      labels(2) = 'RESP_LAG'
      labels(3) = 'H_FORM'
      call add_target('RESP_LAGF',ttype_frm,.false.,tgt_info)
      call set_dependency('RESP_LAGF','RESP_LAG',tgt_info)
      call set_dependency('RESP_LAGF','H_FORM',tgt_info)
      call set_dependency('RESP_LAGF','T_FORM',tgt_info)
      call set_dependency('RESP_LAGF','L_FORM',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'full response lagrangian',1,'---')
      call set_rule('RESP_LAGF',ttype_frm,EXPAND,
     &              labels,3,1,
     &              parameters,2,tgt_info)

      if (ntest.ge.100)
     &  write(luout,*) 'factoring out special intermediates'

      ! R12: factor out special intermediates
      if (setr12) then
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = 'RESP_LAGF' ! output formula (itself)
        labels(2) = 'RESP_LAGF' ! input formula
        labels(3) = form_r12_vint    ! the intermediates to be factored
        labels(4) = form_r12_vint//'^+'
        labels(5) = form_r12_bint
        labels(6) = form_r12_bhint
        labels(7) = form_r12_xint
        labels(8) = form_r12_xhint
        nint = 6
        call set_dependency('RESP_LAGF',form_r12_vint,tgt_info)
        call set_dependency('RESP_LAGF',form_r12_xint,tgt_info)
        call set_dependency('RESP_LAGF',form_r12_bint,tgt_info)
        call set_dependency('RESP_LAGF',form_r12_bhint,tgt_info)
        call set_dependency('RESP_LAGF',form_r12_xhint,tgt_info)
        if (ansatz.ne.1) then
          labels(9) = form_r12_cint
          labels(10) = trim(form_r12_cint)//'^+'
          call set_dependency('RESP_LAGF',form_r12_cint,tgt_info)
          nint = 8
        end if
        call form_parameters(-1,
     &       parameters,2,'r12 response lag., factored out',nint,'---')
        call set_rule('RESP_LAGF',ttype_frm,FACTOR_OUT,
     &                labels,nint+2,1,
     &                parameters,2,tgt_info)
c        call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c        call set_rule('RESP_LAGF',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

        ! replace r12 by the actual integrals
        if (ansatz.gt.1) then
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = 'RESP_LAGF'
          labels(2) = 'RESP_LAGF'
          labels(3) = op_r12
          labels(4) = op_rint
          call set_dependency('RESP_LAGF',op_rint,tgt_info)
          call form_parameters(-1,
     &         parameters,2,'complete r12 response lag.',1,'---')
          call set_rule('RESP_LAGF',ttype_frm,REPLACE,
     &                labels,4,1,
     &                parameters,2,tgt_info)
          call form_parameters(-1,parameters,2,'stdout',1,'stdout')
          call set_rule('RESP_LAGF',ttype_frm,PRINT_FORMULA,
     &                  labels,2,1,parameters,2,tgt_info)
        end if
      end if

      ! expand response lagrangian with T=T(0)+T(1)+...,L=...
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'RESP_LAGF'
      labels(2) = 'RESP_LAGF'
      labels(3) = 'T_FORM'
      labels(4) = 'L_FORM'
      call form_parameters(-1,
     &     parameters,2,'full response lagrangian',2,'---')
      call set_rule('RESP_LAGF',ttype_frm,EXPAND,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! expand RESP_LAG(n) with X(1)=X(1)1+X(1)2+..., X(2)=X(2)12+X(2)13+...
      ! and V(1)=V(1)1+V(1)2+...
      formname(1:len_short) = ' '
      formname(1:8) = 'F_LAG(n)'
      lag_name = 'RESP_LAG(n)'
      do ord = 1,maxval(maxord)
        write(formname(7:7),'(i1)') ord
        write(lag_name(10:10),'(i1)') ord
        allocate(ifreq(ord),ifreqnew(ord))
        ifreq = 0
        do while (next_comb(ifreq,ord,maxord,ncnt))
          call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
          formname2(1:len_short) = ' '
          formname2(1:6) = 'V(1)_F'
          do digit = 1,ord
            write(formname(7+2*digit:8+2*digit),'(i2.2)')
     &                       ifreqnew(digit)
            write(formname2(5+2*digit:6+2*digit),'(i2.2)')
     &                       ifreqnew(digit)
          end do
          if (idx_target(trim(formname),tgt_info).gt.0) cycle
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(formname)
          labels(2) = lag_name
          labels(3) = trim(formname2)
          ilabels = 3
          call add_target(trim(formname),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(formname),lag_name,tgt_info)
          call set_dependency(trim(formname),trim(formname2),tgt_info)
          do op_par = 1,2
            if (op_par.eq.1) then
              formname2(1:1) = 'T'
              x_max_ord = int((real(ord)-1)/2+0.6)
            else
              formname2(1:1) = 'L'
              x_max_ord = int((real(ord)-2)/2+0.6)
            end if
            do ord2=1,x_max_ord
              write(formname2(3:3),'(i1)') ord2
              ilabels = ilabels+1
              labels(ilabels) = trim(formname2)
              call set_dependency(trim(formname),trim(formname2),
     &                            tgt_info)
            end do
          end do
          call form_parameters(-1,
     &         parameters,2,'full response lagrangian with freqs',
     &         ilabels-2,'---')
          call set_rule(trim(formname),ttype_frm,EXPAND,
     &                  labels,ilabels,1,
     &                  parameters,2,tgt_info)
        end do
        deallocate(ifreq,ifreqnew)
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'expanding left and right residuals'

      ! expand left and right residuals RESS_LAG(n)_X
      ! with X(1)=X(1)1+X(1)2+..., X(2)=X(2)12+X(2)13+...
      ! and V(1)=V(1)1+V(1)2+...
      ! following (2n+1) and (2n+2) rules
      formname(1:len_short) = ' '
      formname(1:7) = 'RF(n)SX'
      resl_lag_name = 'RESS_LAG(n)_X'
      do op_par2 = 1,2
        if (op_par2.eq.1) then
          formname(7:7) = 'T'
          resl_lag_name(13:13) = 'T'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        else
          formname(7:7) = 'L'
          resl_lag_name(13:13) = 'L'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        end if
        do ord = op_par2-1,x_max_ord
          write(formname(4:4),'(i1)') ord
          write(resl_lag_name(10:10),'(i1)') ord
          do side = 1,2
            if (side.eq.1) then
              formname(6:6) = 'L'
              resl_lag_name(4:4) = 'L'
            else
              formname(6:6) = 'R'
              resl_lag_name(4:4) = 'R'
            end if
            allocate(ifreq(ord),ifreqnew(ord))
            ifreq = 0
            set_zero = ord.eq.0
            do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
              set_zero = .false.
              call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
              formname2(1:len_short) = ' '
              formname2(1:6) = 'V(1)_F'
              do digit = 1,ord
                write(formname(6+2*digit:7+2*digit),'(i2.2)')
     &                           ifreqnew(digit)
                write(formname2(5+2*digit:6+2*digit),'(i2.2)')
     &                           ifreqnew(digit)
              end do
              if (idx_target(trim(formname),tgt_info).gt.0) cycle
              labels(1:20)(1:len_target_name) = ' '
              labels(1) = trim(formname)
              labels(2) = resl_lag_name
              labels(3) = trim(formname2)
              ilabels = 3
              call add_target(trim(formname),ttype_frm,.false.,tgt_info)
              call set_dependency(trim(formname),resl_lag_name,tgt_info)
              call set_dependency(trim(formname),trim(formname2),
     &                            tgt_info)
              do op_par = 1,2
                if (op_par.eq.1) then
                  formname2(1:1) = 'T'
                  x_max_ord2 = int((real(maxval(maxord))-1)/2+0.6)
                else
                  formname2(1:1) = 'L'
                  x_max_ord2 = int((real(maxval(maxord))-2)/2+0.6)
                 end if
                do ord2=1,min(x_max_ord2,ord)
                  write(formname2(3:3),'(i1)') ord2
                  ilabels = ilabels+1
                  labels(ilabels) = trim(formname2)
                  call set_dependency(trim(formname),
     &                                trim(formname2),tgt_info)
                end do
              end do
              call form_parameters(-1,
     &             parameters,2,'left and right residuals with freqs',
     &             ilabels-2,'---')
              call set_rule(trim(formname),ttype_frm,EXPAND,
     &                      labels,ilabels,1,
     &                      parameters,2,tgt_info)
            end do
            deallocate(ifreq,ifreqnew)
          end do
        end do
      end do

      ! extract lagrangian of orders 0 to maxord
      lagf_name = 'RESP_LAGF(n)'
      do ord = 0,maxval(maxord)
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
      do ord = 0,maxval(maxord)
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

      if (ntest.ge.100)
     &  write(luout,*) 'extracting terms according to frequency pattern'

      ! extract terms of F_LAG(n) with correct freq. pattern
      formname(1:len_short) = ' '
      formname(1:6) = 'LAG(n)'
      formname2(1:len_short) = ' '
      formname2(1:8) = 'F_LAG(n)'
      lagname(1:len_short) = ' '
      lagname(1:8) = 'LRESP(n)'
      do ord = 1,maxval(maxord)
        write(formname(5:5),'(i1)') ord
        write(formname2(7:7),'(i1)') ord
        write(lagname(7:7),'(i1)') ord
        write(pert_ord,'(i1)') ord
        allocate(ifreq(ord),ifreqnew(ord))
        ifreq = 0
        do while (next_comb(ifreq,ord,maxord,ncnt))
          call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
          do digit = 1,ord
            write(formname(5+2*digit:6+2*digit),'(i2.2)') 
     &                       ifreqnew(digit)
            write(formname2(7+2*digit:8+2*digit),'(i2.2)')
     &                       ifreqnew(digit)
            write(lagname(7+2*digit:8+2*digit),'(i2.2)') ifreqnew(digit)
          end do
          if (idx_target(trim(formname),tgt_info).gt.0) cycle
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(formname)
          labels(2) = trim(formname2)
          labels(3) = trim(lagname)
          call add_target(trim(formname),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(formname),trim(formname2),tgt_info)
          call set_dependency(trim(formname),trim(lagname),tgt_info)
          call form_parameters2(-1,
     &         parameters,2,'lagrangian of pert order '//pert_ord//
     &         ' with freqs',
     &         ord,ifreqnew)
          call set_rule(trim(formname),ttype_frm,EXTRACT_FREQ,
     &                  labels,3,1,
     &                  parameters,2,tgt_info)
        end do
        deallocate(ifreq,ifreqnew)
      end do

      ! extract terms of left and right residuals RF(n)SX
      ! with correct freq. pattern
      ! following (2n+1) and (2n+2) rules
      formname(1:6) = 'R(n)SX'
      formname2(1:len_short) = ' '
      formname2(1:7) = 'RF(n)SX'
      opname(1:6) = 'O(n)SX'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        else
          op_parent = 'L'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        end if
        formname(6:6) = op_parent
        formname2(7:7) = op_parent
        opname(6:6) = op_parent
        do ord = op_par-1,x_max_ord
          write(formname(3:3),'(i1)') ord
          write(formname2(4:4),'(i1)') ord
          write(opname(3:3),'(i1)') ord
          write(pert_ord,'(i1)') ord
          do side = 1,2
            if (side.eq.1) then
              leftright = 'L'
            else
              leftright = 'R'
            end if
            formname(5:5) = leftright
            formname2(6:6) = leftright
            opname(5:5) = leftright
            allocate(ifreq(ord),ifreqnew(ord))
            ifreq = 0
            set_zero = ord.eq.0
            do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
              set_zero = .false.
              call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
              formname(7:len_short) = ' '
              opname(7:len_short) = ' '
              do digit = 1,ord
                write(formname(5+2*digit:6+2*digit),'(i2.2)')
     &                                          ifreqnew(digit)
                write(formname2(6+2*digit:7+2*digit),'(i2.2)')
     &                                          ifreqnew(digit)
                write(opname(5+2*digit:6+2*digit),'(i2.2)') 
     &                                          ifreqnew(digit)
              end do
              if (idx_target(trim(formname),tgt_info).gt.0) cycle
              labels(1:20)(1:len_target_name) = ' '
              labels(1) = trim(formname)
              labels(2) = trim(formname2)
              labels(3) = trim(opname)
              call add_target(trim(formname),ttype_frm,.false.,
     &                      tgt_info)
              call set_dependency(trim(formname),trim(formname2),
     &                      tgt_info)
              call set_dependency(trim(formname),trim(opname),
     &                      tgt_info)
              call form_parameters2(-1,
     &             parameters,2,'residual of pert order '//pert_ord//
     &             ' with freqs',
     &             ord,ifreqnew)
              call set_rule(trim(formname),ttype_frm,EXTRACT_FREQ,
     &                      labels,3,1,
     &                      parameters,2,tgt_info)
            end do
            deallocate(ifreq,ifreqnew)
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
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        else
          op_parent = 'L'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
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
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        else
          op_parent = 'L'
          op_name = 'T(n)'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
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

      if (ntest.ge.100)
     &  write(luout,*) 'formulae defined'

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
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        else
          op_parent = 'L'
          optname(1:8) = 'OPT_T(n)'
          defmelname2(1:11) = 'DEF_ME_T(n)'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
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
          allocate(ifreq(ord),ifreqnew(ord))
          ifreq = 0
          set_zero = ord.eq.0
          do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
            set_zero = .false.
            call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
            do digit = 1,ord
              write(optname(7+2*digit:8+2*digit),'(i2.2)') 
     &                             ifreqnew(digit)
              write(formname(5+2*digit:6+2*digit),'(i2.2)') 
     &                             ifreqnew(digit)
              write(formname2(5+2*digit:6+2*digit),'(i2.2)') 
     &                             ifreqnew(digit)
              write(defmelname(12+2*digit:13+2*digit),'(i2.2)')
     &                             ifreqnew(digit)
              write(defmelname2(10+2*digit:11+2*digit),'(i2.2)')
     &                             ifreqnew(digit)
            end do
            if (idx_target(trim(optname),tgt_info).gt.0) cycle
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
            call set_dependency(trim(optname),trim(formname2),
     &                    tgt_info)
            call set_dependency(trim(optname),mel_ham,tgt_info)
            if (ord.gt.0)
     &         call set_dependency(trim(optname),'HUB_DEFME_V(1)',
     &                             tgt_info)
            call set_dependency(trim(optname),trim(defmelname2),
     &                          tgt_info)
            defmelname(12:12) = 'L'
            call set_dependency(trim(optname),trim(defmelname),
     &                    tgt_info)
            defmelname(12:12) = 'R'
            call set_dependency(trim(optname),trim(defmelname),
     &                    tgt_info)
            call opt_parameters(-1,parameters,ncat,nint)
            call set_rule(trim(optname),ttype_frm,OPTIMIZE,
     &                    labels,ncat+nint+1,1,
     &                    parameters,1,tgt_info)
          end do
          deallocate(ifreq,ifreqnew)
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
      call get_argument_value('calculate.routes','simtraf',
     &                        ival=isim)
      if (isim.eq.1) then
        nint = 1
        call set_dependency('OPT_T(0)',form_cchhat,tgt_info)
        call set_dependency('OPT_T(0)',mel_hhatdef,tgt_info)
        labels(4) = form_cchhat
      else if (isim.eq.2) then
        nint = 1
        call set_dependency('OPT_T(0)',form_cchbar,tgt_info)
        call set_dependency('OPT_T(0)',meldef_hbar,tgt_info)
        labels(4) = form_cchbar
      end if
      call set_dependency('OPT_T(0)','RES_LAG(0)_L',tgt_info)
      call set_dependency('OPT_T(0)','RESP_LAG(0)',tgt_info)
      call set_dependency('OPT_T(0)','DEF_ME_LRESP(0)',tgt_info)
      call set_dependency('OPT_T(0)','DEF_ME_O(0)_L',tgt_info)
      call set_dependency('OPT_T(0)','DEF_ME_T(0)',tgt_info)
      call set_dependency('OPT_T(0)',mel_ham,tgt_info)
      ! r12 dependencies
      if (setr12) then
        call set_dependency('OPT_T(0)',mel_x_def,tgt_info)
        call set_dependency('OPT_T(0)',mel_bh_def,tgt_info)
        call set_dependency('OPT_T(0)',mel_b_def,tgt_info)
        call set_dependency('OPT_T(0)',mel_v_def,tgt_info)
        call set_dependency('OPT_T(0)',mel_c_def,tgt_info)
        call set_dependency('OPT_T(0)',mel_rint,tgt_info)
      end if 
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('OPT_T(0)',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! OPT_LRESP(n): for evaluation of LRESP(n)
      optname(1:len_short) =  ' '
      optname(1:12) = 'OPT_LRESP(n)'
      defmelname(1:len_short) = ' '
      defmelname(1:15) = 'DEF_ME_LRESP(n)'
      formname(1:len_short) = ' '
      formname(1:11) = 'RESP_LAG(0)'
      do ord = 0,maxval(maxord)
        write(optname(11:11),'(i1)') ord
        write(defmelname(14:14),'(i1)') ord
        allocate(ifreq(ord),ifreqnew(ord))
        ifreq = 0
        set_zero = ord.eq.0
        do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
          set_zero = .false.
          call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
          do digit = 1,ord
            write(defmelname(14+2*digit:15+2*digit),'(i2.2)')
     &                                   ifreqnew(digit)
            write(optname(11+2*digit:12+2*digit),'(i2.2)') 
     &                                   ifreqnew(digit)
            if (formname(1:3).eq.'LAG')
     &         write(formname(5+2*digit:6+2*digit),'(i2.2)') 
     &                                   ifreqnew(digit)
          end do
          if (idx_target(trim(optname),tgt_info).gt.0) cycle
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(optname)
          labels(2) = trim(formname)
          ncat = 1
          nint = 0
          call add_target(trim(optname),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(optname),trim(formname),tgt_info)
          call set_dependency(trim(optname),trim(defmelname),tgt_info)
          call opt_parameters(-1,parameters,ncat,nint)
          call set_rule(trim(optname),ttype_frm,OPTIMIZE,
     &                  labels,ncat+nint+1,1,
     &                  parameters,1,tgt_info)
        end do
        deallocate(ifreq,ifreqnew)
        formname(1:11) = 'LAG(n)     '
        write(formname(5:5),'(i1)') ord+1
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'opt. formulae defined'

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! for H(0) (op_ham), a ME-list is already defined and imported

      ! ME_V(1)1, ME_V(1)2, ...
      opname(1:len_short) = ' '
      opname(1:4) = 'V(1)'
      melname(1:len_short) = ' '
      melname(1:7) = 'ME_V(1)'
      defmelname(1:len_short) = ' '
      defmelname(1:11) = 'DEF_ME_V(1)'
      freq_idx = 0
      do while (next_comb(freq_idx,1,maxord,ncnt))
        call redundant_comb(freq_idx,freq_idxnew,redun,1,maxord,ncnt)
        write(opname(5:6),'(i2.2)') freq_idxnew
        write(melname(8:9),'(i2.2)') freq_idxnew
        write(defmelname(12:13),'(i2.2)') freq_idxnew
        if (idx_target(trim(defmelname),tgt_info).gt.0) cycle
        call add_target(trim(defmelname),ttype_opme,.false.,tgt_info)
        call set_dependency(trim(defmelname),trim(opname),tgt_info)
        ! (a) define
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = trim(melname)
        labels(2) = trim(opname)
        call me_list_parameters(-1,parameters,
     &       msc,0,isym(freq_idxnew),0,0,.false.)
        call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        ! (b) import
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = trim(melname)
        call import_parameters(-1,parameters,
     &                         pert(freq_idxnew:freq_idxnew)//
     &                         'DIPLEN','DALTON_SPECIAL')
        call set_rule(trim(defmelname),ttype_opme,IMPORT,
     &                labels,1,1,
     &                parameters,1,tgt_info)
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'ME_V lists defined'

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
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        else
          opname2(6:6) = 'L'
          opname(1:1) = 'T'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        end if
        do ord = 0,x_max_ord
          write(opname(3:3),'(i1)') ord
          opname(5:len_short) = ' '
          write(opname2(3:3),'(i1)') ord
          opname2(7:len_short) = ' '
          allocate(ifreq(ord),ifreqnew(ord))
          ifreq = 0
          set_zero = ord.eq.0
          do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
            set_zero = .false.
            call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
            sym = 1
            do digit = 1,ord
              write(opname(3+2*digit:4+2*digit),'(i2.2)') 
     &                                    ifreqnew(digit)
              sym = multd2h(sym,isym(ifreqnew(digit)))
            end do
            melname(4:len_short) = opname(1:len_short-3)
            defmelname(8:len_short) = opname(1:len_short-7)
            if (idx_target(trim(defmelname),tgt_info).gt.0) cycle
            call add_target(trim(defmelname),ttype_opme,.false.,
     &                      tgt_info)
            call set_dependency(trim(defmelname),trim(opname),
     &                          tgt_info)
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
            ! (b) ME_X(n)w
            labels(1:20)(1:len_target_name) = ' '
            labels(1) = trim(melname)
            labels(2) = trim(opname)
            call me_list_parameters(-1,parameters,
     &           msc,0,sym,0,0,.false.)
            call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                    labels,2,1,
     &                    parameters,1,tgt_info)
            call set_rule(trim(defmelname),ttype_opme,SET_FREQ,
     &                    labels,2,1,
     &                    parameters,1,tgt_info)
          end do
          deallocate(ifreq,ifreqnew)
          ! (c) ME_OS(n)_Xw (S=L,R; X=T,L; n=0,...,maxord)
          do side = 1,2
            if (side.eq.1) then
              leftright = 'L'
            else
              leftright = 'R'
            end if
            opname2(5:5) = leftright
            allocate(ifreq(ord),ifreqnew(ord))
            ifreq = 0
            set_zero = ord.eq.0
            do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
              set_zero = .false.
              call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
              sym = 1
              do digit = 1,ord
                write(opname2(5+2*digit:6+2*digit),'(i2.2)') 
     &                                       ifreqnew(digit)
                sym = multd2h(sym,isym(ifreqnew(digit)))
              end do
              melname2(4:len_short) = opname2(1:len_short-3)
              defmelname2(8:len_short) = opname2(1:len_short-7)
              if (idx_target(trim(defmelname2),tgt_info).gt.0) cycle
              if (ord.ge.op_par-1) then
                call add_target(trim(defmelname2),
     &                          ttype_opme,.false.,tgt_info)
                call set_dependency(trim(defmelname2),trim(opname2),
     &                              tgt_info)
                labels(1:20)(1:len_target_name) = ' '
                labels(1) = trim(melname2)
                labels(2) = trim(opname2)
                call me_list_parameters(-1,parameters,
     &               msc,0,sym,0,0,.false.)
                call set_rule(trim(defmelname2),ttype_opme,
     &                        DEF_ME_LIST,
     &                        labels,2,1,
     &                        parameters,1,tgt_info)
              end if
            end do
            deallocate(ifreq,ifreqnew)
          end do
        end do
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'ME_L, ME_T, ME_O-lists defined'

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
      defmelname(1:len_short) = ' '
      defmelname(1:15) = 'DEF_ME_LRESP(n)'
      melname(1:len_short) = ' '
      melname(1:11) = 'ME_LRESP(n)'
      lagname(1:len_short) = ' '
      lagname(1:8) = 'LRESP(n)'
      do ord = 0,maxval(maxord)
        write(defmelname(14:14),'(i1)') ord
        write(melname(10:10),'(i1)') ord
        write(lagname(7:7),'(i1)') ord
        allocate(ifreq(ord),ifreqnew(ord))
        ifreq = 0
        set_zero = ord.eq.0
        do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
          set_zero = .false.
          call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
          sym = 1
          do digit = 1,ord
            write(defmelname(14+2*digit:15+2*digit),'(i2.2)')
     &                                   ifreqnew(digit)
            write(melname(10+2*digit:11+2*digit),'(i2.2)') 
     &                                   ifreqnew(digit)
            write(lagname(7+2*digit:8+2*digit),'(i2.2)') ifreqnew(digit)
            sym = multd2h(sym,isym(ifreqnew(digit)))
          end do
          if (ord.gt.0) then
            icnt = (ifreq(1)-1)/maxval(maxord)+1
          else
            icnt = 1
          end if
          if (sym.ne.1.and.ord.eq.maxord(icnt)) then
            ! print zero as result and don't calculate any other property
            evaluate(icnt) = .false.
            allocate(mel_dummy(1))
            call print_result(ord,ifreqnew,
     &                        mel_dummy(1)%mel,.true.)
            deallocate(mel_dummy)
          end if
          if (idx_target(trim(defmelname),tgt_info).gt.0) cycle
          call add_target(trim(defmelname),ttype_opme,.false.,
     &                                                tgt_info)
          call set_dependency(trim(defmelname),trim(lagname),tgt_info)
          if (ord.gt.0) call set_dependency(trim(defmelname),
     &                       'HUB_DEFME_V(1)',tgt_info)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(melname)
          labels(2) = trim(lagname)
          call me_list_parameters(-1,parameters,
     &         msc,0,sym,0,0,.false.)
          call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                  labels,2,1,
     &                  parameters,1,tgt_info)
        end do
        deallocate(ifreq,ifreqnew)
      end do

      ! Preconditioner: one for each possible irrep
      do sym = 1,maxsym
        call me_list_label(mel_dia1,mel_dia,sym,0,0,0,.false.)
        call add_target(mel_dia1,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_dia1,mel_ham,tgt_info)
        call set_dependency(mel_dia1,op_dia,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_dia1
        labels(2) = op_dia
        call me_list_parameters(-1,parameters,
     &       0,0,sym,0,0,.false.)
        call set_rule(mel_dia1,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
        labels(1) = mel_dia1
        labels(2) = mel_ham
        call set_rule(mel_dia1,ttype_opme,PRECONDITIONER,
     &                labels,2,1,
     &                parameters,0,tgt_info)
      end do

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

      if (ntest.ge.100)
     &  write(luout,*) 'me-lists defined'

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*
  
      ! solve CC-equations 
      call add_target('SOLVE_T(0)',ttype_gen,.false.,tgt_info)
      if (setr12)
     &       call set_dependency('SOLVE_T(0)',eval_r12_inter,tgt_info)
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

      if (ntest.ge.100)
     &  write(luout,*) 'define solvers for linear equations'

      ! solve linear equations for X(n)w 
      ! (X=T,L, n=0,...,x_max_ord, not for T(0))
      do op_par = 1,2
        if (op_par.eq.1) then
          defmelname(1:11) = 'DEF_ME_L(n)'
          melname(1:7) = 'ME_L(n)'
          opname(1:6) = 'O(n)LT'
          optname(1:8) = 'OPT_L(n)'
          solvename(1:10) = 'SOLVE_L(n)'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        else
          defmelname(1:11) = 'DEF_ME_T(n)'
          melname(1:7) = 'ME_T(n)'
          opname(1:6) = 'O(n)LL'
          optname(1:8) = 'OPT_T(n)'
          solvename(1:10) = 'SOLVE_T(n)'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        end if
        do ord = op_par-1,x_max_ord
          write(solvename(9:9),'(i1)') ord
          write(defmelname(10:10),'(i1)') ord
          write(melname(6:6),'(i1)') ord
          write(opname(3:3),'(i1)') ord
          write(optname(7:7),'(i1)') ord
          allocate(ifreq(ord),ifreqnew(ord))
          ifreq = 0
          set_zero = ord.eq.0
          do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
            set_zero = .false.
            call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
            sym = 1
            solvename(11:len_short) = ' '
            defmelname(12:len_short) = ' '
            opname(7:len_short) = ' '
            melname(8:len_short) = ' '
            optname(9:len_short) = ' '
            do digit = 1,ord
              write(solvename(9+2*digit:10+2*digit),'(i2.2)')
     &                                   ifreqnew(digit)
              write(defmelname(10+2*digit:11+2*digit),'(i2.2)')
     &                                   ifreqnew(digit)
              write(optname(7+2*digit:8+2*digit),'(i2.2)') 
     &                                   ifreqnew(digit)
              write(opname(5+2*digit:6+2*digit),'(i2.2)') 
     &                                   ifreqnew(digit)
              write(melname(6+2*digit:7+2*digit),'(i2.2)') 
     &                                   ifreqnew(digit)
              sym = multd2h(sym,isym(ifreqnew(digit)))
            end do
            if (idx_target(trim(solvename),tgt_info).gt.0) cycle
            call add_target(trim(solvename),ttype_gen,.false.,
     &                          tgt_info)
            if (setr12)
     &              call set_dependency(trim(solvename),
     &              eval_r12_inter,tgt_info)
            call set_dependency(trim(solvename),trim(defmelname),
     &                          tgt_info)
            call set_dependency(trim(solvename),trim(optname),
     &                          tgt_info)
            call me_list_label(mel_dia1,mel_dia,sym,0,0,0,.false.)
            call solve_parameters(-1,parameters,2, 1,1,'DIA')
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
            ! SOLVE_Y(n)w
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
          end do
          deallocate(ifreq,ifreqnew)
        end do
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'define evaluaters'

      ! evaluate RESP_LAG(n)
      optname(1:len_short) = ' '
      optname(1:12) = 'OPT_LRESP(n)'
      evalname(1:len_short) = ' '
      evalname(1:11) = 'EVAL_LAG(n)'
      hubname(1:14) = 'HUB_SOLVE_X(n)'
      defmelname(1:len_short) = ' '
      defmelname(1:15) = 'DEF_ME_LRESP(n)'
      melname(1:len_short) = ' '
      melname(1:11) = 'ME_LRESP(n)'
      eval_dipmom = .false.
      do ord = 0,maxval(maxord)
        write(optname(11:11),'(i1)') ord
        write(evalname(10:10),'(i1)') ord
        write(defmelname(14:14),'(i1)') ord
        write(melname(10:10),'(i1)') ord
        allocate(ifreq(ord),ifreqnew(ord))
        ifreq = 0
        set_zero = ord.eq.0
        do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
          set_zero = .false.
          call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
          sym = 1
          freqsum = 0d0
          if (ord.gt.0) then
            icnt = (ifreq(1)-1)/maxval(maxord)+1
          else
            icnt = 1
          end if
          do digit = 1,ord
            pos = ifreq(digit)-(icnt-1)*maxval(maxord)
            write(evalname(10+2*digit:11+2*digit),'(i2.2)')
     &                                 ifreqnew(digit)
            write(defmelname(14+2*digit:15+2*digit),'(i2.2)')
     &                                 ifreqnew(digit)
            write(optname(11+2*digit:12+2*digit),'(i2.2)') 
     &                                 ifreqnew(digit)
            write(melname(10+2*digit:11+2*digit),'(i2.2)') 
     &                                 ifreqnew(digit)
            sym = multd2h(sym,isym(ifreqnew(digit)))
            freqsum = freqsum + freq(icnt,pos)
          end do
          ! only evaluate IF (requested prop not zero by symmetry (evaluate)
          !                   OR prop is energy (ord.eq.0))
          ! AND this prop not zero by symmetry (sym.eq.1)
          ! AND (sum of frequencies is zero (abs(freqsum).lt.1d-12)
          !      OR prop is dipole moment (ord.eq.1))
          if ((evaluate(icnt).or.(ord.eq.0.and..not.all(.not.evaluate)))
     &        .and.(sym.eq.1).and.
     &        (abs(freqsum).lt.1d-12.or.ord.eq.1)) then
            if (idx_target(trim(evalname),tgt_info).gt.0) cycle
            if (ord.eq.1) then
              ! evaluate dipole moment (ord.eq.1) only once
              if (pert(ifreqnew(1):ifreqnew(1)).eq.'X') then
                pertdir = 1
              else if (pert(ifreqnew(1):ifreqnew(1)).eq.'Y') then
                pertdir = 2
              else if (pert(ifreqnew(1):ifreqnew(1)).eq.'Z') then
                pertdir = 3
              else
                call quit(1,'set_experimental_targets',
     &                      'pert must be X, Y, or Z')
              end if
              if (eval_dipmom(pertdir)) cycle
              eval_dipmom(pertdir) = .true.
            end if
            call add_target(trim(evalname),ttype_gen,.true.,tgt_info)
            if (setr12)
     &           call set_dependency(trim(evalname),eval_r12_inter,
     &                                              tgt_info)
            ! first solve T(k), L(k) using (2n+1) and (2n+2) rules
            hubname(11:11) = 'T'
            x_max_ord = int((real(ord)-1)/2+0.6)
            write(hubname(13:13),'(i1)') x_max_ord
            call set_dependency(trim(evalname),hubname(1:14),tgt_info)
            if (maxord(icnt).gt.0) then
              hubname(11:11) = 'L'
              x_max_ord = int((real(ord)-2)/2+0.6)
              write(hubname(13:13),'(i1)') x_max_ord
              call set_dependency(trim(evalname),hubname(1:14),
     &                                           tgt_info)
            end if
            ! evaluate
            call set_dependency(trim(evalname),trim(optname),tgt_info)
            call set_dependency(trim(evalname),trim(defmelname),
     &                                         tgt_info)
            labels(1:20)(1:len_target_name) = ' '
            labels(1) = trim(optname)
            call set_rule(trim(evalname),ttype_opme,EVAL,
     &           labels,1,0,
     &           parameters,0,tgt_info)
            ! print result
            call ord_parameters(-1,parameters,ord,0,ifreqnew)
            call set_rule(trim(evalname),ttype_opme,PRINT_RES,
     &             trim(melname),1,1,
     &             parameters,1,tgt_info)
          end if
        end do
        deallocate(ifreq,ifreqnew)
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'define dependency hubs'

      ! hub for DEF_ME_Y(n)w-dependencies
      defmelname(1:len_short) = ' '
      defmelname(1:11) = 'DEF_ME_Y(n)'
      hubname(1:len_short) = ' '
      hubname(1:14) = 'HUB_DEFME_Y(n)'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        else
          op_parent = 'L'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        end if
        defmelname(8:8) = op_parent
        hubname(11:11) = op_parent
        do ord = 0,x_max_ord
          write(defmelname(10:10),'(i1)') ord
          defmelname(12:len_short) = ' '
          write(hubname(13:13),'(i1)') ord
          call add_target(trim(hubname),ttype_gen,.false.,tgt_info)
          allocate(ifreq(ord),ifreqnew(ord))
          ifreq = 0
          set_zero = ord.eq.0
          do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
            set_zero = .false.
            call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
            do digit = 1,ord
              write(defmelname(10+2*digit:11+2*digit),'(i2.2)')
     &                                     ifreqnew(digit)
            end do
            call set_dependency(trim(hubname),trim(defmelname),
     &                                     tgt_info)
          end do
          deallocate(ifreq,ifreqnew)
        end do
      end do   

      ! hub for DEF_ME_V(1)w-dependencies
      defmelname(1:len_short) = ' '
      defmelname(1:11) = 'DEF_ME_V(1)'
      call add_target('HUB_DEFME_V(1)',ttype_gen,.false.,tgt_info)
      freq_idx = 0
      do while (next_comb(freq_idx,1,maxord,ncnt))
        call redundant_comb(freq_idx,freq_idxnew,redun,1,maxord,ncnt)
        write(defmelname(12:13),'(i2.2)') freq_idxnew
        call set_dependency('HUB_DEFME_V(1)',trim(defmelname),
     &                                    tgt_info)
      end do

      ! hub for SOLVE_Y(n)w-dependencies
      solvename(1:len_short) = ' '
      solvename(1:10) = 'SOLVE_Y(n)'
      hubname(1:len_short) = ' '
      hubname(1:14) = 'HUB_SOLVE_Y(n)'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        else
          op_parent = 'L'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        end if
        solvename(7:7) = op_parent
        hubname(11:11) = op_parent
        do ord = 0,x_max_ord
          write(solvename(9:9),'(i1)') ord
          solvename(11:len_short) = ' '
          write(hubname(13:13),'(i1)') ord
          call add_target(trim(hubname),ttype_gen,.false.,tgt_info)
          allocate(ifreq(ord),ifreqnew(ord))
          ifreq = 0
          set_zero = ord.eq.0
          do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
            set_zero = .false.
            call redundant_comb(ifreq,ifreqnew,redun,ord,maxord,ncnt)
            do digit = 1,ord
              write(solvename(9+2*digit:10+2*digit),'(i2.2)')
     &                                     ifreqnew(digit)
            end do
            ! only solve if property (for which Y(n)w is required) is needed
            if (ord.gt.0) then
              icnt = (ifreq(1)-1)/maxval(maxord)+1
            else
              icnt = 1
            end if
            if (op_par.eq.1) then
              x_max_ord2 = int((real(maxord(icnt))-1)/2+0.6)
            else
              x_max_ord2 = int((real(maxord(icnt))-2)/2+0.6)
            end if
            if ((evaluate(icnt).or.ord.eq.0).and.ord.le.x_max_ord2)
     &        call set_dependency(trim(hubname),trim(solvename),
     &                                     tgt_info)
          end do
          deallocate(ifreq,ifreqnew)
        end do
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'phony targets processed'

*----------------------------------------------------------------------*
*     deallocate arrays
*----------------------------------------------------------------------*

      deallocate(maxord,freq,evaluate,redun,isym)

      return
      end
