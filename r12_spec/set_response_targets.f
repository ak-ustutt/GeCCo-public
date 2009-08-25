*----------------------------------------------------------------------*
      subroutine set_response_targets(tgt_info,orb_info,env_type)
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

      include 'def_pert_info.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      character(*), intent(in) ::
     &     env_type

      integer, parameter ::
     &     ntest = 100,
     &     len_short = 32, len_long = 256, maxsym = 8

      integer ::
     &     isim, ncat, nint, ndef, occ_def(ngastp,2,60),
     &     msc, sym, r12op, ansatz, maxexc, t1ext_mode,
     &     ord, op_par, len_op_exp, side, x_max_ord,
     &     digit, ilabels, ord2, op_par2, x_max_ord2,
     &     ncnt, icnt, pos,
     &     spec, trunc_type, ipop, npop, ncmp, restart

      character(len_target_name) ::
     &     mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)

      character(len_short) ::
     &     opname, melname, defmelname, formname, formname2,
     &     opname2, melname2, defmelname2, hubname, optname,
     &     solvename, evalname, lagname, formname3, formname4,
     &     defmelname3, opname3, opname4, defmelname4, melname3,
     &     defmelname5, formname5, optname2, melname4, evalname2

      character(len_long) ::
     &     opexp, method

      character ::
     &     op_name*4, op_parent*1, lagf_name*12, op_exp*50,
     &     pert_ord*1, lag_name*11, res_name*9, leftright*1,
     &     res_lag_name*15,resl_lag_name*16, resr_lag_name*16,
     &     lresp_name*8, approx_str*12

      character(len=8) ::
     &     approx

      character ::
     &     op_bprc*4 = 'BPRC',
     &     op_xprc*4 = 'XPRC',
     &     me_bprc*8 = 'BPRCLIST',
     &     me_xprc*8 = 'XPRCLIST',
     &     medef_bprc*12 = 'DEF-BPRCLIST',
     &     medef_xprc*12 = 'DEF-XPRCLIST'

      logical ::
     &     setr12, r12fix, set_zero, treat_bv, file_exists,
     &     truncate, half, userules, opt, screen

      real(8) ::
     &     freqsum

      type(pert_op_info) ::
     &     pop(3*maxpop)

      type(pert_component_info) ::
     &     cmp(maxcmp)

      integer, allocatable ::
     &     ifreq(:), ifreqnew(:), maxord(:),
     &     subset(:)

      logical, allocatable ::
     &     evaluate(:), eval_dipmom(:)

      type(me_list_array), allocatable ::
     &     mel_dummy(:)

      logical, external ::
     &     next_comb

      integer, external ::
     &     idx_target

*----------------------------------------------------------------------*
*     process keywords and allocate arrays:
*----------------------------------------------------------------------*

      ! skip this section if not requested
      ncnt = is_keyword_set('calculate.response')
      if (ncnt.eq.0) return

      ! set r12 targets
      setr12 = is_keyword_set('method.R12').gt.0
      if (setr12) then
        call get_argument_value('method.R12','ansatz',ival=ansatz)
        call get_argument_value('method.R12','approx',str=approx)
        call get_argument_value('method.R12','fixed',lval=r12fix)
        call get_argument_value('method.R12','r12op',ival=r12op)
        call get_argument_value('method.R12','T1ext',ival=t1ext_mode)
        call get_argument_value('method.R12','trunc',ival=trunc_type)
        call get_argument_value('method.R12','opt',lval=opt)
        call get_argument_value('method.R12','screen',lval=screen)
        truncate = trunc_type.ge.0
        if (is_keyword_set('method.truncate').gt.0) then
          truncate = is_keyword_set('method.truncate').gt.0
          if(truncate)
     &       call get_argument_value('method.truncate','trunc_type',
     &                                ival=trunc_type)
        end if
        call get_argument_value('calculate.response','BX_3C',
     &       keycount=1,lval=treat_bv)
        if (r12op.le.1.and.r12fix.and..not.(opt.and.r12op.ne.1)) then
          call set_r12f_general_targets(tgt_info,orb_info,env_type)
        else
          call quit(1,'set_response_targets','wrong r12 method')
        end if
      else
        t1ext_mode = 0
      end if

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting response targets ...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

      allocate(maxord(ncnt),evaluate(ncnt))

      call get_response_input(ncnt,maxord,ncmp,cmp,npop,pop,orb_info)

      allocate(eval_dipmom(npop))
      evaluate = .true.

      ! maximum excitation
      call get_argument_value('calculate.response','maxexc',
     &     keycount=1,ival=maxexc)
      ! use 2n+1 / 2n+2 rules?
      call get_argument_value('calculate.response','rules',
     &     keycount=1,lval=userules)
      if (.not.userules.and..not.all(abs(cmp(1:ncmp)%freq).lt.1d-12))
     &       call quit(1,'set_response_targets',
     &       'no-rules option only supported for static response')
      ! restart calculation? Requires amplitude mel files
      call get_argument_value('calculate.response','restart',
     &     keycount=1,ival=restart)
      melname(1:len_short) = ' '
      do op_par = 1,2
        if (op_par.eq.1) then
          melname(1:11) = 'ME_T(n)'
          x_max_ord = int((real(restart-1)-1)/2+0.6)
        else
          melname(1:11) = 'ME_L(n)'
          x_max_ord = int((real(restart-1)-2)/2+0.6)
        end if
        if (.not.userules.and.restart.gt.0) 
     &          call quit(1,'set_response_targets',
     &          'restart option currently only with rules=T')
        if (restart.lt.op_par) cycle
        do ord = 0,x_max_ord
          write(melname(6:6),'(i1)') ord
          melname(8:len_short) = ' '
          allocate(ifreq(ord),ifreqnew(ord))
          ifreq = 0
          set_zero = ord.eq.0
          do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
            set_zero = .false.
            call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                          maxord,ncnt)
            do digit = 1,ord
              write(melname(6+2*digit:7+2*digit),'(i2.2)')
     &                             ifreqnew(digit)
            end do
            inquire(file=trim(melname)//'_list.da',exist=file_exists)
            if (file_exists) then
              write(luout,*) 'Restart calculation using ',
     &                       trim(melname)//'_list.da'
            else
              call quit(1,'set_response_targets','Restart requires '//
     &                    'file '//trim(melname)//'_list.da')
            end if
          end do
          deallocate(ifreq,ifreqnew)
        end do
      end do

      if (ntest.ge.100) then
        write(luout,*) 'keywords processed:'
        write(luout,*) 'maxexc      : ',maxexc
        write(luout,*) 'userules    : ',userules
        write(luout,*) 'restart     : ',restart
        write(luout,*) 'BX_3C       : ',treat_bv
        write(luout,*) 'trunc_type  : ',trunc_type
        write(luout,*) 'ncnt        : ',ncnt
        write(luout,*) 'maxord      : ',maxord
        write(luout,*) 'freq        : ',cmp(1:ncmp)%freq
        write(luout,*) 'redun       : ',cmp(1:ncmp)%redun
        write(luout,*) 'pop_idx     : ',cmp(1:ncmp)%pop_idx
        write(luout,*) 'pop%name    : ',pop(1:npop)%name
        write(luout,*) 'pop%comp    : ',pop(1:npop)%comp
        write(luout,*) 'pop%int_name: ',pop(1:npop)%int_name
        write(luout,*) 'pop%isym    : ',pop(1:npop)%isym
      end if

      if (iprlvl.gt.0) then
        if (.not.userules)
     &        write(luout,*) 'No 2n+1 and 2n+2 rules will be used'
        method(1:6) = 'CCSDTQ'
        do ord = 5,maxexc
          write(method(ord+2:ord+2),'(i1)') ord
        end do
        pos = maxexc + 3
        method(pos:len_long) = ' '
        if (setr12.and.trunc_type.ne.2) then
          if (r12op.eq.0) then
            method(pos:pos+2) = '-SP'
            pos = pos+3
          else
            method(pos:pos+3) = '-XSP'
            pos = pos+4
            if (opt) then
              method(pos:pos+2) = 'opt'
              pos = pos+3
            end if
          end if
          if (screen) then
            method(pos:pos+7) = '(screen)'
            pos = pos+8
          end if
          if (treat_bv.and.maxval(maxord).gt.0) then
            write(luout,*) 'BX intermediate will be evaluated '
     &                     //'in approximation C'
          else if (maxval(maxord).gt.0) then
            write(luout,*) 'BX intermediate will be evaluated without '
     &                         //'special treatment'
          end if
        end if
        if (t1ext_mode.eq.1) method(pos:pos+3) = '-T1x'
        write(luout,*) 'Method to be used: '//trim(method)
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
        if (opt) then
          ! set operator block version of T12 and T12BAR to 2
          call ord_parameters(-1,parameters,1,-1,2)
          call set_rule(op_cex,ttype_op,SET_ORDER,op_cex,
     &                  1,1,parameters,1,tgt_info)
          call set_rule(op_cexbar,ttype_op,SET_ORDER,op_cexbar,
     &                  1,1,parameters,1,tgt_info)
        end if
        ! BV
        call add_target('BV',ttype_op,.false.,tgt_info)
        call set_dependency('BV',op_b_inter,tgt_info)
        call cloneop_parameters(-1,parameters,
     &       op_b_inter,.false.)
        call set_rule('BV',ttype_op,CLONE_OP,
     &                'BV',1,1,
     &                parameters,1,tgt_info)
        call ord_parameters(-1,parameters,1,3,-1)
        call set_rule('BV',ttype_op,SET_ORDER,'BV',
     &                1,1,parameters,1,tgt_info)
        ! CV
        call add_target('CV',ttype_op,.false.,tgt_info)
        call set_dependency('CV',op_c_inter,tgt_info)
        call cloneop_parameters(-1,parameters,
     &       op_c_inter,.false.)
        call set_rule('CV',ttype_op,CLONE_OP,
     &                'CV',1,1,
     &                parameters,1,tgt_info)
        call ord_parameters(-1,parameters,1,3,-1)
        call set_rule('CV',ttype_op,SET_ORDER,'CV',
     &                1,1,parameters,1,tgt_info)

        ! define direction components BVX, BVY, BVZ, CVX, CVY, CVZ
        ! and RDGB(V)X, RDGB(V)Y, RDGB(V)Z,
        ! and RBRV(V)X, RBRV(V)Y, RBRV(V)Z
        opname(1:len_short) = ' '
        opname(1:2) = 'BV'
        opname2(1:len_short) = ' '
        opname2(1:7) = 'RDGB(V)'
        opname3(1:len_short) = ' '
        opname3(1:7) = 'RBRV(V)'
        opname4(1:len_short) = ' '
        opname4(1:2) = 'CV'
        occ_def = 0
        ndef = 2
        ! 1
        occ_def(IHOLE,1,1) = 1
        occ_def(IPART,1,1) = 1
        occ_def(IHOLE,2,1) = 2
        ! 2
        occ_def(IHOLE,1,2) = 1
        occ_def(IEXTR,1,2) = 1
        occ_def(IHOLE,2,2) = 2
        if (r12op.eq.1) then
          ndef = 4
          ! 3
          occ_def(IHOLE,1,3) = 1
          occ_def(IPART,1,3) = 1
          occ_def(IHOLE,2,3) = 1
          occ_def(IPART,2,3) = 1
          ! 4
          occ_def(IHOLE,1,4) = 1
          occ_def(IEXTR,1,4) = 1
          occ_def(IHOLE,2,4) = 1
          occ_def(IPART,2,4) = 1
        end if
        do ipop = 1,npop
          opname(2:3) = pop(ipop)%name//pop(ipop)%comp
          opname2(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          opname3(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          opname4(2:3) = pop(ipop)%name//pop(ipop)%comp
          ! BVV and CVV:
          call add_target(trim(opname),ttype_op,.false.,tgt_info)
          call set_dependency(trim(opname),op_b_inter,tgt_info)
          call cloneop_parameters(-1,parameters,
     &         op_b_inter,.false.)
          call set_rule(trim(opname),ttype_op,CLONE_OP,
     &                  trim(opname),1,1,
     &                  parameters,1,tgt_info)
          call add_target(trim(opname4),ttype_op,.false.,tgt_info)
          call set_dependency(trim(opname4),op_c_inter,tgt_info)
          call cloneop_parameters(-1,parameters,
     &         op_c_inter,.false.)
          call set_rule(trim(opname4),ttype_op,CLONE_OP,
     &                  trim(opname4),1,1,
     &                  parameters,1,tgt_info)
          call ord_parameters(-1,parameters,1,ipop+3,0)
          call set_rule(trim(opname),ttype_op,SET_ORDER,trim(opname),
     &                  1,1,parameters,1,tgt_info)
          call set_rule(trim(opname4),ttype_op,SET_ORDER,trim(opname4),
     &                  1,1,parameters,1,tgt_info)
          ! RDGB(V)V:
          call add_target(trim(opname2),ttype_op,.false.,tgt_info)
          call set_dependency(trim(opname2),op_rint,tgt_info)
          call cloneop_parameters(-1,parameters,
     &         op_rint,.false.)
          call set_rule(trim(opname2),ttype_op,CLONE_OP,
     &                  trim(opname2),1,1,
     &                  parameters,1,tgt_info)
          ! RBRV(V)V:
          call add_target(trim(opname3),ttype_op,.false.,tgt_info)
          call op_from_occ_parameters(-1,parameters,2,
     &         occ_def,ndef,1,(/.false.,.true./),ndef)
          call set_rule(trim(opname3),ttype_op,DEF_OP_FROM_OCC,
     &                  trim(opname3),1,1,
     &                  parameters,2,tgt_info)
        end do

        ! define V(1)Vh and V(1)Vf (for BV intermediate, frozen core)
        if (orb_info%ngas.eq.4.and.treat_bv) then
          half = .true.
          opname(1:len_short) = ' '
          opname(1:6) = 'V(1)_h'
          do side = 1,2
            if (.not.half) opname(6:6) = 'f'
            do ipop = 1,npop
              opname(1:5) = pop(ipop)%name//'(1)'//pop(ipop)%comp
              call add_target(trim(opname),ttype_op,.false.,tgt_info)
              call genop_parameters(-1,parameters,3,
     &                        1,1,0,
     &                       -1,1,3,(/half,.false./),
     &                      (/0,2, 0,2, 0,0, 0,2/),
     &                      (/0,1, 0,1, 0,0, 0,1,
     &                        0,1, 0,1, 0,0, 0,1/),4,
     &                      (/0,0, 0,1, 0,1, 0,1,
     &                        0,1, 0,1, 0,1, 0,1,
     &                        0,0, 0,0, 0,0, 0,0,
     &                        0,0, 0,0, 0,0, 0,0/),4)
            call set_rule(trim(opname),ttype_op,DEF_GENERAL_OPERATOR,
     &                    trim(opname),1,1,
     &                    parameters,3,tgt_info)
            end do
            half = .false.
          end do
        else if (orb_info%ngas.ne.3.and.treat_bv) then
          call quit(1,'set_response_targets','special evaluation '//
     &              'BV intermediate requires ngas.eq.3/4')
        end if
      end if

      ! define V(1)
      call add_target('V(1)',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,1,1,setr12)
      call set_rule('V(1)',ttype_op,DEF_HAMILTONIAN,'V(1)',
     &              1,1,parameters,1,tgt_info)
      call ord_parameters(-1,parameters,1,3,-1)
      call set_rule('V(1)',ttype_op,SET_ORDER,'V(1)',
     &              1,1,parameters,1,tgt_info)

      ! define direction components V(1)X, V(1)Y, V(1)Z
      opname(1:len_short) = ' '
      opname(1:4) = 'V(1)'
      do ipop = 1,npop
        opname(1:5) = pop(ipop)%name//'(1)'//pop(ipop)%comp
        call add_target(trim(opname),ttype_op,.false.,tgt_info)
        call hop_parameters(-1,parameters,0,1,3,setr12)
        call set_rule(trim(opname),ttype_op,DEF_HAMILTONIAN,
     &                trim(opname),
     &                1,1,parameters,1,tgt_info)
        call ord_parameters(-1,parameters,1,ipop+3,0)
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
      if (t1ext_mode.eq.0.and..not.opt) then
        call xop_parameters(-1,parameters,
     &                   .false.,1,maxexc,0,1)
        call set_rule(op_top,ttype_op,DEF_EXCITATION,
     &              op_top,1,1,
     &              parameters,1,tgt_info)
      else
        occ_def = 0
        ! T1
        occ_def(IPART,1,1) = 1
        occ_def(IHOLE,2,1) = 1
        ndef = 2
        if (t1ext_mode.gt.0) then
          ! T1'
          occ_def(IEXTR,1,ndef) = 1
          occ_def(IHOLE,2,ndef) = 1
          ndef = ndef + 1
        end if
        if (opt) then
          ! another T1
          occ_def(IPART,1,ndef) = 1
          occ_def(IHOLE,2,ndef) = 1
          ndef = ndef + 1
        end if
        ! T2,T3,...
        do nint = 2,maxexc
          occ_def(IPART,1,nint+ndef-2) = nint
          occ_def(IHOLE,2,nint+ndef-2) = nint
        end do
        ndef = ndef + maxexc - 2
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
        if (.not.userules) x_max_ord = maxval(maxord)
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
              call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                            maxord,ncnt)
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
        if (.not.userules) x_max_ord = maxval(maxord)
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
      opname2(1:len_short) = ' '
      opname2(1:5) = 'SX(n)'
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
        if (.not.userules) x_max_ord = maxval(maxord)
        opname2(2:2) = op_parent
        do side = 1,2
          if (side.eq.1) then
            leftright = 'L'
          else
            leftright = 'R'
          end if
          opname(5:5) = leftright
          do ord=op_par-1,x_max_ord
            write(opname(3:3),'(i1)') ord
            write(opname2(4:4),'(i1)') ord
            opname(7:len_short) = ' '
            opname2(6:len_short) = ' '
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
              call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                            maxord,ncnt)
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
              if (setr12.and.side.eq.1.and.
     &            (ord.gt.0.or.opt.and.op_par.eq.1)) then
                ! SX(n)w: X(n)w vector times Metric
                opname2(6:len_short-1) = opname(7:len_short)
                call add_target(trim(opname2),ttype_op,.false.,tgt_info)
                call set_dependency(trim(opname2),op_parent,tgt_info)
                call cloneop_parameters(-1,parameters,
     &                                  op_parent,.false.)
                call set_rule(trim(opname2),ttype_op,CLONE_OP,
     &                        trim(opname2),1,1,
     &                        parameters,1,tgt_info)
              end if
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
          call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                        maxord,ncnt)
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

      ! Diagonal Preconditioner
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'L'
        else
          op_parent = 'T'
        end if
        call add_target(op_dia//'_'//op_parent,ttype_op,.false.,
     &                  tgt_info)
        call set_dependency(op_dia//'_'//op_parent,'T(0)',tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          'T(0)',op_par.eq.1)
        call set_rule(op_dia//'_'//op_parent,ttype_op,CLONE_OP,
     &                op_dia//'_'//op_parent,1,1,
     &                parameters,1,tgt_info)
      end do

      if (opt) then
        ! B-Matrix for Preconditioner
        call add_target(op_bprc,ttype_op,.false.,tgt_info)
        occ_def=0
        ndef = 1
        occ_def(IPART,1,1) = 1
        occ_def(IPART,2,1) = 1
        call op_from_occ_parameters(-1,parameters,2,
     &       occ_def,ndef,1,(/.true.,.true./),ndef)
        call set_rule(op_bprc,ttype_op,DEF_OP_FROM_OCC,
     &                op_bprc,1,1,
     &                parameters,2,tgt_info)
        ! X-Matrix for Preconditioner
        call add_target(op_xprc,ttype_op,.false.,tgt_info)
        occ_def=0
        ndef = 1
        occ_def(IPART,1,1) = 1
        occ_def(IPART,2,1) = 1
        call op_from_occ_parameters(-1,parameters,2,
     &       occ_def,ndef,1,(/.true.,.true./),ndef)
        call set_rule(op_xprc,ttype_op,DEF_OP_FROM_OCC,
     &                op_xprc,1,1,
     &                parameters,2,tgt_info)
      end if

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
        call define_r12_hhat(tgt_info,screen)
      end if

      ! Hbar intermediate
      call add_target(op_hbar,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &                    .false.,1,2,0,1)
      call set_rule(op_hbar,ttype_op,DEF_CC_HBAR_OP,
     &              op_hbar,1,1,
     &              parameters,1,tgt_info)

      if (setr12) then
        ! prototype formula for metric:
        call add_target(op_ccr12s0,ttype_op,.false.,tgt_info)
        call set_rule(op_ccr12s0,ttype_op,DEF_SCALAR,
     &                op_ccr12s0,1,1,
     &                parameters,0,tgt_info)

        ! S_X: L and T vectors times Metric
        opname(1:len_short) = ' '
        opname(1:3) = 'S_X'
        do op_par = 1,2
          if (op_par.eq.1) then
            op_parent = 'T'
            opname(3:3) = 'L'
          else
            op_parent = 'L'
            opname(3:3) = 'T'
          end if
          call add_target(trim(opname),ttype_op,.false.,tgt_info)
          call set_dependency(trim(opname),op_parent,tgt_info)
          call cloneop_parameters(-1,parameters,
     &                            op_parent,.false.)
          call set_rule(trim(opname),ttype_op,CLONE_OP,
     &                  trim(opname),1,1,
     &                  parameters,1,tgt_info)
        end do
      end if

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
        labels(3) = 'Hnew'
        labels(4) = 'L'
        labels(5) = 'T'
        call add_target('RESP_LAG',ttype_frm,.false.,tgt_info)
        call set_dependency('RESP_LAG','LRESP',tgt_info)
        call set_dependency('RESP_LAG','Hnew',tgt_info)
        call set_dependency('RESP_LAG','L',tgt_info)
        call set_dependency('RESP_LAG','T',tgt_info)
        call set_rule('RESP_LAG',ttype_frm,DEF_CC_LAGRANGIAN,
     &                labels,5,1,
     &                'conventional CC lagrangian',1,tgt_info)
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
      end if
      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
      call set_rule('RESP_LAG',ttype_frm,PRINT_FORMULA,
     &              labels,2,1,parameters,2,tgt_info)

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
        if (.not.userules) x_max_ord = maxval(maxord)
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

      ! frequency expansion of V(1): V(1)=V(1)X+V(1)Y+V(1)Z
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'V(1)_F'
      call add_target('V(1)_F',ttype_frm,.false.,tgt_info)
      opname(1:len_short) = ' '
      opname(1:4) = 'V(1)'
      opexp(1:len_long) = ' '
      opexp(1:4) = 'V(1)'
      len_op_exp = 4
      do ipop = 1,npop
        opname(1:5) = pop(ipop)%name//'(1)'//pop(ipop)%comp
        call set_dependency('V(1)_F',trim(opname),tgt_info)
        opexp(len_op_exp+1:len_op_exp+6) = '+'//opname(1:5)
        len_op_exp = len_op_exp + 6
      end do
      opexp(5:5) = '='
      call def_form_parameters(-1,
     &     parameters,2,trim(opexp),
     &     'expansion of V(1)')
      call set_rule('V(1)_F',ttype_frm,DEF_FORMULA,
     &              labels,1,1,
     &              parameters,2,tgt_info)

      if (setr12) then
        ! frequency expansion of BV: BV=BVX+BVY+BVZ
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = 'BV_F'
        call add_target('BV_F',ttype_frm,.false.,tgt_info)
        opname(1:len_short) = ' '
        opname(1:2) = 'BV'
        opexp(1:len_long) = ' '
        opexp(1:2) = 'BV'
        len_op_exp = 2
        do ipop = 1,npop
          opname(2:3) = pop(ipop)%name//pop(ipop)%comp
          call set_dependency('BV_F',trim(opname),tgt_info)
          opexp(len_op_exp+1:len_op_exp+4) = '+'//opname(1:3)
          len_op_exp = len_op_exp + 4
        end do
        opexp(3:3) = '='
        call def_form_parameters(-1,
     &       parameters,2,trim(opexp),
     &       'expansion of BV')
        call set_rule('BV_F',ttype_frm,DEF_FORMULA,
     &                labels,1,1,
     &                parameters,2,tgt_info)
        ! frequency expansion of CV: CV=CVX+CVY+CVZ
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = 'CV_F'
        call add_target('CV_F',ttype_frm,.false.,tgt_info)
        opname(1:len_short) = ' '
        opname(1:2) = 'CV'
        opexp(1:len_long) = ' '
        opexp(1:2) = 'CV'
        len_op_exp = 2
        do ipop = 1,npop
          opname(2:3) = pop(ipop)%name//pop(ipop)%comp
          call set_dependency('CV_F',trim(opname),tgt_info)
          opexp(len_op_exp+1:len_op_exp+4) = '+'//opname(1:3)
          len_op_exp = len_op_exp + 4
        end do
        opexp(3:3) = '='
        call def_form_parameters(-1,
     &       parameters,2,trim(opexp),
     &       'expansion of CV')
        call set_rule('CV_F',ttype_frm,DEF_FORMULA,
     &                labels,1,1,
     &                parameters,2,tgt_info)
      end if

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
          if (.not.userules) x_max_ord = maxval(maxord)
          do ord2=1,min(x_max_ord,ord)
            write(op_name(3:3),'(i1)') ord2
            formname(1:4) = op_name(1:4)
            allocate(ifreq(ord),ifreqnew(ord))
            ifreq = 0
            do while (next_comb(ifreq,ord,maxord,ncnt))
              call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                            maxord,ncnt)
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

      ! R12: define BV-intermediate and factor out special intermediates
      if (setr12) then
        ! define BV- and CV-intermediates
        do op_par = 1,2
          if (op_par.eq.1) then
            op_parent = 'B'
          else
            op_parent = 'C'
          end if
          formname(1:len_short) = ' '
          formname(1:9) = op_parent//'V_formal'
          opname(1:len_short) = ' '
          opname(1:4) = 'V(1)'
          opname2(1:len_short) = ' '
          opname2(1:2) = op_parent//'V'
          do ipop = 1,npop+1
            call add_target(trim(formname),ttype_frm,.false.,tgt_info)
            call set_dependency(trim(formname),trim(opname),tgt_info)
            call set_dependency(trim(formname),trim(opname2),tgt_info)
            call set_dependency(trim(formname),op_r12,tgt_info)
            labels(1) = trim(formname)
            labels(2) = trim(opname2)
            labels(3) = op_r12
            labels(4) = trim(opname)
            call form_parameters(-1,
     &           parameters,2,trim(formname),0,op_parent//'p')
            call set_rule(trim(formname),ttype_frm,DEF_R12INTM_FORMAL,
     &                    labels,4,1,
     &                    parameters,2,tgt_info)
            ! replace r12 by the actual integrals
            if (ipop.gt.1) then
              labels(1:20)(1:len_target_name) = ' '
              labels(1) = trim(formname)
              labels(2) = trim(formname)
              labels(3) = op_r12
              labels(4) = op_rint
              labels(5) = op_r12//'^+'
              labels(6) = op_rint//'^+'
              call set_dependency(trim(formname),op_rint,tgt_info)
              call form_parameters(-1,
     &             parameters,2,op_parent//'V-intermediate formula',2,
     &             '---')
              call set_rule(trim(formname),ttype_frm,REPLACE,
     &                    labels,6,1,
     &                    parameters,2,tgt_info)
            end if
            call form_parameters(-1,
     &           parameters,2,'stdout',0,'---')
            call set_rule(trim(formname),ttype_frm,PRINT_FORMULA,
     &                    labels,1,0,
     &                    parameters,2,tgt_info)
            formname(2:9) = pop(ipop)%name//pop(ipop)%comp//'_form '
            opname(1:5) = pop(ipop)%name//'(1)'//pop(ipop)%comp
            opname2(2:3) = pop(ipop)%name//pop(ipop)%comp
          end do
        end do

        ! R12: factor out special intermediates
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = 'RESP_LAGF' ! output formula (itself)
        labels(2) = 'RESP_LAGF' ! input formula
        labels(3) = form_r12_vint    ! the intermediates to be factored
        labels(4) = form_r12_vint//'^+'
        labels(5) = form_r12_bint
        labels(6) = form_r12_bhint
        labels(7) = form_r12_xint
        labels(8) = form_r12_xhint
        labels(9) = 'BV_formal'
        nint = 7
        call set_dependency('RESP_LAGF',form_r12_vint,tgt_info)
        call set_dependency('RESP_LAGF',form_r12_xint,tgt_info)
        call set_dependency('RESP_LAGF',form_r12_bint,tgt_info)
        call set_dependency('RESP_LAGF',form_r12_bhint,tgt_info)
        call set_dependency('RESP_LAGF',form_r12_xhint,tgt_info)
        call set_dependency('RESP_LAGF','BV_formal',tgt_info)
        if (ansatz.ne.1) then
          labels(10) = form_r12_cint
          labels(11) = trim(form_r12_cint)//'^+'
          labels(12) = 'CV_formal'
          labels(13) = 'CV_formal^+'
          call set_dependency('RESP_LAGF',form_r12_cint,tgt_info)
          call set_dependency('RESP_LAGF','CV_formal',tgt_info)
          nint = 11
        else
          call quit(1,'set_my_response_targets',
     &              'Currently not adapted for ansatz 1')
        end if
        call form_parameters(-1,
     &       parameters,2,'r12 response lag., factored out',nint,'---')
        call set_rule('RESP_LAGF',ttype_frm,FACTOR_OUT,
     &                labels,nint+2,1,
     &                parameters,2,tgt_info)

        ! screen all remaining R12 terms
        if (screen) then
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = 'RESP_LAGF'
          labels(2) = 'RESP_LAGF'
          labels(3) = 'LRESP'
          labels(4) = op_r12
          labels(5) = op_r12//'^+'
          call set_rule('RESP_LAGF',ttype_frm,INVARIANT,
     &                labels,5,1,
     &                title_ccr12lg0,1,tgt_info)
        end if

        ! replace r12 by the actual integrals
        if (ansatz.gt.1) then
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = 'RESP_LAGF'
          labels(2) = 'RESP_LAGF'
          labels(3) = op_r12
          labels(4) = op_rint
          labels(5) = op_r12//'^+'
          labels(6) = op_rint//'^+'
          call set_dependency('RESP_LAGF',op_rint,tgt_info)
          call form_parameters(-1,
     &         parameters,2,'complete r12 response lag.',2,'---')
          call set_rule('RESP_LAGF',ttype_frm,REPLACE,
     &                labels,6,1,
     &                parameters,2,tgt_info)
          call form_parameters(-1,parameters,2,'stdout',1,'stdout')
          call set_rule('RESP_LAGF',ttype_frm,PRINT_FORMULA,
     &                  labels,2,1,parameters,2,tgt_info)
c          call set_rule('RESP_LAGF',ttype_frm,TEX_FORMULA,
c     &                  labels,5,1,
c     &                  'ccr12_lag0.tex',1,tgt_info)
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
      ! and V(1)=V(1)X+V(1)Y+V(1)Z / BV=BVX+BVY+BVZ / CV=CVX+CVY+CVZ
      formname(1:len_short) = ' '
      formname(1:8) = 'F_LAG(n)'
      lag_name = 'RESP_LAG(n)'
      do ord = 1,maxval(maxord)
        write(formname(7:7),'(i1)') ord
        write(lag_name(10:10),'(i1)') ord
        allocate(ifreq(ord),ifreqnew(ord))
        ifreq = 0
        do while (next_comb(ifreq,ord,maxord,ncnt))
          call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                        maxord,ncnt)
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
          labels(3) = formname2(1:6)
          ilabels = 3
          call add_target(trim(formname),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(formname),lag_name,tgt_info)
          call set_dependency(trim(formname),formname2(1:6),tgt_info)
          if (setr12) then
            labels(4) = 'BV_F'
            labels(5) = 'CV_F'
            labels(6) = 'CV_F^+'
            ilabels = 6
            call set_dependency(trim(formname),'BV_F',tgt_info)
            call set_dependency(trim(formname),'CV_F',tgt_info)
          end if
          do op_par = 1,2
            if (op_par.eq.1) then
              formname2(1:1) = 'T'
              x_max_ord = int((real(ord)-1)/2+0.6)
            else
              formname2(1:1) = 'L'
              x_max_ord = int((real(ord)-2)/2+0.6)
            end if
            if (.not.userules) x_max_ord = ord
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
          call form_parameters(-1,parameters,2,'stdout',1,'stdout')
          call set_rule(trim(formname),ttype_frm,PRINT_FORMULA,
     &                  labels,2,1,parameters,2,tgt_info)
        end do
        deallocate(ifreq,ifreqnew)
      end do

      if (ntest.ge.100)
     &  write(luout,*) 'expanding left and right residuals'

      ! expand left and right residuals RESS_LAG(n)_X
      ! with X(1)=X(1)1+X(1)2+..., X(2)=X(2)12+X(2)13+...
      ! and V(1)=V(1)X+V(1)Y+V(1)Z / BV=BVX+BVY+BVZ
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
        if (.not.userules) x_max_ord = maxval(maxord)
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
              call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                            maxord,ncnt)
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
              labels(3) = formname2(1:6)
              ilabels = 3
              call add_target(trim(formname),ttype_frm,.false.,tgt_info)
              call set_dependency(trim(formname),resl_lag_name,tgt_info)
              call set_dependency(trim(formname),formname2(1:6),
     &                            tgt_info)
              if (setr12) then
                labels(4) = 'BV_F'
                labels(5) = 'CV_F'
                labels(6) = 'CV_F^+'
                ilabels = 6
                call set_dependency(trim(formname),'BV_F',tgt_info)
                call set_dependency(trim(formname),'CV_F',tgt_info)
              end if
              do op_par = 1,2
                if (op_par.eq.1) then
                  formname2(1:1) = 'T'
                  x_max_ord2 = int((real(maxval(maxord))-1)/2+0.6)
                else
                  formname2(1:1) = 'L'
                  x_max_ord2 = int((real(maxval(maxord))-2)/2+0.6)
                end if
                if (.not.userules) x_max_ord2 = maxval(maxord)
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

      if (.not.userules) then
        ! extract terms obeying (2n+1) and (2n+2) rules for NLEQ-solver
        lagf_name = 'RESP_LAGF(0)'
        lag_name = 'CCEQ_LAG(0)'
        lresp_name = 'LRESP(0)'
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = lag_name
        labels(2) = lagf_name
        labels(3) = lresp_name
        call add_target(lag_name,ttype_frm,.false.,tgt_info)
        call set_dependency(lag_name,lagf_name,tgt_info)
        call set_dependency(lag_name,lresp_name,tgt_info)
        call form_parameters(-1,
     &       parameters,2,'lagrangian of pert order 0 for CC-equations',
     &       -1,'---')
        call set_rule(lag_name,ttype_frm,EXTRACT_ORDER,
     &                labels,3,1,
     &                parameters,2,tgt_info)
      end if

      ! extract terms obeying (2n+1) and (2n+2) rules if userules
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
        nint = ord
        if (userules) nint = -1
        call form_parameters(-1,
     &       parameters,2,'lagrangian of pert order '//pert_ord,
     &       nint,'---')
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
          call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                        maxord,ncnt)
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
     &         parameters,3,'lagrangian of pert order '//pert_ord//
     &         ' with freqs',
     &         ord,ifreqnew,ncmp,cmp(1:ncmp)%pop_idx)
          call set_rule(trim(formname),ttype_frm,EXTRACT_FREQ,
     &                  labels,3,1,
     &                  parameters,3,tgt_info)
          if (.not.userules) then
            call form_parameters(-1,
     &           parameters,2,'stdout',1,'---')
            call set_rule(trim(formname),ttype_frm,CLASS_FORMULA,
     &                    labels,1,0,
     &                    parameters,2,tgt_info)
          end if
        end do
        deallocate(ifreq,ifreqnew)
      end do

      ! extract terms of left and right residuals RF(n)SX
      ! with correct freq. pattern
      ! following (2n+1) and (2n+2) rules
      formname(1:len_short) = ' '
      formname(1:6) = 'R(n)SX'
      formname2(1:len_short) = ' '
      formname2(1:7) = 'RF(n)SX'
      opname(1:len_short) = ' '
      opname(1:6) = 'O(n)SX'
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'T'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
        else
          op_parent = 'L'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
        end if
        if (.not.userules) x_max_ord = maxval(maxord)
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
              call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                            maxord,ncnt)
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
     &             parameters,3,'residual of pert order '//pert_ord//
     &             ' with freqs',
     &             ord,ifreqnew,ncmp,cmp(1:ncmp)%pop_idx)
              call set_rule(trim(formname),ttype_frm,EXTRACT_FREQ,
     &                      labels,3,1,
     &                      parameters,3,tgt_info)
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
        if (.not.userules) x_max_ord = maxval(maxord)
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
        if (.not.userules) x_max_ord = maxval(maxord)
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
      if (t1ext_mode.gt.0) then
        ! define a Hhat including T1' for linear equations
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_cchhat//'_LEQ'
        labels(2) = op_hhat
        labels(3) = op_ham
        labels(4) = 'T(0)'
        call add_target(form_cchhat//'_LEQ',ttype_frm,.false.,tgt_info)
        call set_dependency(form_cchhat//'_LEQ',op_hhat,tgt_info)
        call set_dependency(form_cchhat//'_LEQ',op_ham,tgt_info)
        call set_dependency(form_cchhat//'_LEQ','T(0)',tgt_info)
        call form_parameters(-1,
     &             parameters,2,title_cchhat,2,'---')
        call set_rule(form_cchhat//'_LEQ',ttype_frm,DEF_HHAT,
     &                labels,4,1,
     &                parameters,2,tgt_info)
      end if

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

      if (setr12) then
        ! formal Metric (in 'complete' basis)
        ! (a) Define
        call add_target(form_ccr12_s0,ttype_frm,.false.,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12_s0
        labels(2) = op_ccr12s0
        labels(3) = op_r12
        labels(4) = 'L'
        labels(5) = 'T'
        ilabels = 5
        if (r12op.eq.1) then
          call set_dependency(form_ccr12_s0,op_cex,tgt_info)
          call set_dependency(form_ccr12_s0,op_cexbar,tgt_info)
          labels(ilabels+1) = op_cexbar
          labels(ilabels+2) = op_cex
          ilabels = ilabels+2
        end if
        call set_dependency(form_ccr12_s0,op_ccr12s0,tgt_info)
        call set_dependency(form_ccr12_s0,op_r12,tgt_info)
        call set_dependency(form_ccr12_s0,'L',tgt_info)
        call set_dependency(form_ccr12_s0,'T',tgt_info)
        call form_parameters(-1,
     &       parameters,2,title_ccr12lg0,ansatz,'---')
        call set_rule(form_ccr12_s0,ttype_frm,DEF_CCR12_METRIC,
     &                labels,ilabels,1,
     &                parameters,2,tgt_info)
        ! (b) Factor out the R12 intermediates 
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12_s0 ! output formula (itself)
        labels(2) = form_ccr12_s0 ! input formula
        labels(3) = form_r12_xint
        nint = 1
        call set_dependency(form_ccr12_s0,form_r12_xint,tgt_info)
        call form_parameters(-1,
     &       parameters,2,title_ccr12lg0,nint,'---')
        call set_rule(form_ccr12_s0,ttype_frm,FACTOR_OUT,
     &                labels,nint+2,1,
     &                parameters,2,tgt_info)
        call form_parameters(-1,parameters,2,'stdout',1,'stdout')
        call set_rule(form_ccr12_s0,ttype_frm,PRINT_FORMULA,
     &                labels,2,1,parameters,2,tgt_info)
        ! (c) SF_X: derivation wrt L or T:
        opname(1:len_short) = ' '
        opname(1:3) = 'S_X'
        formname(1:len_short) = ' '
        formname(1:4) = 'SF_X'
        do op_par = 1,2
          if (op_par.eq.1) then
            op_parent = 'T'
          else
            op_parent = 'L'
          end if
          opname(3:3) = op_parent
          formname(4:4) = op_parent
          call add_target(trim(formname),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(formname),form_ccr12_s0,tgt_info)
          call set_dependency(trim(formname),trim(opname),tgt_info)
          call set_dependency(trim(formname),op_parent,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = trim(formname)
          labels(2) = form_ccr12_s0
          labels(3) = trim(opname)
          labels(4) = op_parent
          labels(5) = ' '
          call form_parameters(-1,
     &           parameters,2,
     &           'derivation of metric wrt '//op_parent,1,'-')
          call set_rule(trim(formname),ttype_frm,DERIVATIVE,
     &                  labels,5,1,
     &                  parameters,2,tgt_info)
          call form_parameters(-1,parameters,2,'stdout',1,'stdout')
          call set_rule(trim(formname),ttype_frm,PRINT_FORMULA,
     &                  labels,2,1,parameters,2,tgt_info)
        end do
        ! (d) derivation wrt T or L and replace with T(n)w or L(n)w
        formname(1:len_short) = ' '
        formname(1:6) = 'SFX(n)'
        formname2(1:len_short) = ' '
        formname2(1:4) = 'SF_Y'
        opname(1:len_short) = ' '
        opname(1:4) = 'X(n)'
        opname2(1:len_short) = ' '
        opname2(1:5) = 'SX(n)'
        do op_par = 1,2
          if (op_par.eq.1) then
            op_parent = 'L'
            formname2(4:4) = 'T'
            x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
          else
            op_parent = 'T'
            formname2(4:4) = 'L'
            x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
          end if
          if (.not.userules) x_max_ord = maxval(maxord)
          opname(1:1) = op_parent
          do ord = op_par-1,x_max_ord
            write(opname(3:3),'(i1)') ord
            allocate(ifreq(ord),ifreqnew(ord))
            ifreq = 0
            set_zero = ord.eq.0
            do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
              set_zero = .false.
              call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                            maxord,ncnt)
              formname(7:len_short) = ' '
              opname(5:len_short) = ' '
              opname2(6:len_short) = ' '
              do digit = 1,ord
                write(opname(3+2*digit:4+2*digit),'(i2.2)') 
     &                                          ifreqnew(digit)
              end do
              formname(3:len_short) = opname(1:len_short-2)
              opname2(2:len_short) = opname(1:len_short-1)
              if (idx_target(trim(formname),tgt_info).gt.0) cycle
              if (ord.eq.0.and..not.(opt.and.op_par.eq.1)) cycle
              labels(1:20)(1:len_target_name) = ' '
              labels(1) = trim(formname)
              labels(2) = trim(formname2)
              labels(3) = trim(opname2)
              labels(4) = op_parent
              labels(5) = trim(opname)
              call add_target(trim(formname),ttype_frm,.false.,
     &                      tgt_info)
              call set_dependency(trim(formname),trim(formname2),
     &                      tgt_info)
              call set_dependency(trim(formname),trim(opname),
     &                      tgt_info)
              call set_dependency(trim(formname),trim(opname2),
     &                      tgt_info)
              call set_dependency(trim(formname),op_parent,
     &                      tgt_info)
              call form_parameters(-1,
     &             parameters,2,
     &             'vector to solve times Metric',1,'-')
              call set_rule(trim(formname),ttype_frm,DERIVATIVE,
     &                    labels,5,1,
     &                    parameters,2,tgt_info)
              call form_parameters(-1,parameters,2,'stdout',1,'stdout')
              call set_rule(trim(formname),ttype_frm,PRINT_FORMULA,
     &                      labels,2,1,parameters,2,tgt_info)
            end do
            deallocate(ifreq,ifreqnew)
          end do
        end do
      end if

      if (setr12.and.treat_bv) then
        ! RDGB(V)V_F, RBRV(V)V_F, and BV_CABS_V:
        formname(1:len_short) = ' '
        formname(1:10) = 'RDGB(V)V_F'
        formname2(1:len_short) = ' '
        formname2(1:10) = 'RBRV(V)V_F'
        formname3(1:len_short) = ' '
        formname3(1:8) = 'BV_CABS_'
        opname(1:len_short) = ' '
        opname(1:7) = 'RDGB(V)'
        opname2(1:len_short) = ' '
        opname2(1:7) = 'RBRV(V)'
        opname3(1:len_short) = ' '
        opname3(1:4) = 'V(1)'
        opname4(1:len_short) = ' '
        opname4(1:2) = 'BV'
        do ipop = 1,npop
          formname(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          formname2(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          formname3(2:9) = pop(ipop)%name//'_CABS_'//pop(ipop)%comp
          opname(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          opname2(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          opname3(1:5) = pop(ipop)%name//'(1)'//pop(ipop)%comp
          opname4(2:3) = pop(ipop)%name//pop(ipop)%comp
          if (orb_info%ngas.eq.4) opname3(6:6) = 'f'
          ! RDGB(V)V_F:
          call add_target(trim(formname),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(formname),trim(opname3),tgt_info)
          call set_dependency(trim(formname),op_rint,tgt_info)
          call set_dependency(trim(formname),trim(opname),tgt_info)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(formname)
          labels(2) = trim(opname)
          labels(3) = op_rint
          labels(4) = '-'
          labels(5) = trim(opname3)
          call form_parameters(-1,
     &         parameters,2,'RBAR+(V)',3,'RD')
          call set_rule(trim(formname),ttype_frm,DEF_R12INTM_CABS,
     &                  labels,5,1,
     &                  parameters,2,tgt_info)
          call form_parameters(-1,
     &         parameters,2,'stdout',0,'---')
          call set_rule(trim(formname),ttype_frm,PRINT_FORMULA,
     &                  labels,1,0,
     &                  parameters,2,tgt_info)
          ! RBRV(V)V_F:
          call add_target(trim(formname2),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(formname2),trim(opname3),tgt_info)
          call set_dependency(trim(formname2),op_rint,tgt_info)
          call set_dependency(trim(formname2),trim(opname2),tgt_info)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(formname2)
          labels(2) = trim(opname2)
          labels(3) = op_rint
          labels(4) = '-'
          labels(5) = trim(opname3)
          call form_parameters(-1,
     &         parameters,2,'RBRV(V)',3,'RV')
          call set_rule(trim(formname2),ttype_frm,DEF_R12INTM_CABS,
     &                  labels,5,1,
     &                  parameters,2,tgt_info)
          call form_parameters(-1,
     &         parameters,2,'stdout',0,'---')
          call set_rule(trim(formname2),ttype_frm,PRINT_FORMULA,
     &                  labels,1,0,
     &                  parameters,2,tgt_info)
          ! BV_CABS_V:
          if (orb_info%ngas.eq.4) opname3(6:6) = 'h'
          call add_target(trim(formname3),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(formname3),trim(opname3),tgt_info)
          call set_dependency(trim(formname3),trim(opname4),tgt_info)
          call set_dependency(trim(formname3),op_rint,tgt_info)
          call set_dependency(trim(formname3),trim(opname),tgt_info)
          call set_dependency(trim(formname3),trim(opname2),tgt_info)
          call set_dependency(trim(formname3),'FF-X',tgt_info)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(formname3)
          labels(2) = trim(opname4)
          labels(3) = op_rint
          labels(4) = trim(opname)
          labels(5) = 'FF-X'
          labels(6) = trim(opname3)
          labels(7) = trim(opname2)
          approx_str(1:12) = ' '
          approx_str(1:1) = 'C'
          approx_str(12:12) = 'S'
          call form_parameters(-1,
     &         parameters,2,'BV CABS',3,'BV'//approx_str)
          call set_rule(trim(formname3),ttype_frm,DEF_R12INTM_CABS,
     &                  labels,7,1,
     &                  parameters,2,tgt_info)
          call form_parameters(-1,
     &         parameters,2,'stdout',0,'---')
          call set_rule(trim(formname3),ttype_frm,PRINT_FORMULA,
     &                  labels,1,0,
     &                  parameters,2,tgt_info)
        end do
      end if

      if (ntest.ge.100)
     &  write(luout,*) 'formulae defined'

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! OPT_Y(n)w: optimized formula for determination of Y(n)w
      ! depends on DEF_ME_H(0), DEF_ME_V(1)V (if n>0),
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
      formname3(1:len_short) = ' '
      formname3(1:6) = 'SFX(n)'
      defmelname3(1:len_short) = ' '
      defmelname3(1:12) = 'DEF_ME_SX(n)'
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
        if (.not.userules) x_max_ord = maxval(maxord)
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
            call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                          maxord,ncnt)
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
            formname3(3:len_short-2) =  optname(5:len_short)
            defmelname3(9:len_short) = optname(5:len_short-4)
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
            if (setr12.and.(ord.gt.0.or.opt.and.op_par.eq.1)) then
              ! add vector times Metric term:
              ncat = ncat + 1
              labels(ncat+1) = trim(formname3)
              call set_dependency(trim(optname),trim(formname3),
     &                            tgt_info)
              call set_dependency(trim(optname),trim(defmelname3),
     &                            tgt_info)
            end if
            if (isim.eq.1) then
              nint = 1
              if (t1ext_mode.eq.0) then
                call set_dependency(trim(optname),form_cchhat,tgt_info)
                labels(ncat+nint+1) = form_cchhat
              else
                call set_dependency(trim(optname),form_cchhat//'_LEQ',
     &                              tgt_info)
                labels(ncat+nint+1) = form_cchhat//'_LEQ'
              end if
              call set_dependency(trim(optname),mel_hhatdef,tgt_info)
            else if (isim.eq.2) then
              nint = 1
              call set_dependency(trim(optname),form_cchbar,tgt_info)
              call set_dependency(trim(optname),meldef_hbar,tgt_info)
              labels(ncat+nint+1) = form_cchbar
            end if
            call set_dependency(trim(optname),trim(formname),tgt_info)
            call set_dependency(trim(optname),trim(formname2),
     &                    tgt_info)
            call set_dependency(trim(optname),mel_ham,tgt_info)
            if (ord.gt.0)
     &         call set_dependency(trim(optname),'HUB_DEFME_V(1)',
     &                             tgt_info)
            if (ord.gt.0.and.setr12) then
              call set_dependency(trim(optname),'HUB_DEFME_BV',
     &                             tgt_info)
              call set_dependency(trim(optname),'HUB_DEFME_CV',
     &                             tgt_info)
            end if
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
      if (userules) then
        labels(3) = 'RESP_LAG(0)'
      else
        labels(3) = 'CCEQ_LAG(0)'
      end if
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
      if (userules) then
        call set_dependency('OPT_T(0)','RESP_LAG(0)',tgt_info)
        call set_dependency('OPT_T(0)','DEF_ME_LRESP(0)',tgt_info)
      else
        call set_dependency('OPT_T(0)','CCEQ_LAG(0)',tgt_info)
        call set_dependency('OPT_T(0)','DEF_ME_LCCEQ(0)',tgt_info)
      end if
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
      hubname(1:len_short) = ' '
      hubname(1:14) = 'HUB_DEFME_Y(n)'
      do ord = 0,maxval(maxord)
        write(optname(11:11),'(i1)') ord
        write(defmelname(14:14),'(i1)') ord
        allocate(ifreq(ord),ifreqnew(ord))
        ifreq = 0
        set_zero = ord.eq.0
        do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
          set_zero = .false.
          call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                        maxord,ncnt)
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
          do op_par = 1,2
            if (op_par.eq.1) then
              hubname(11:11) = 'T'
              x_max_ord2 = int((real(maxval(maxord))-1)/2+0.6)
            else
              hubname(11:11) = 'L'
              x_max_ord2 = int((real(maxval(maxord))-2)/2+0.6)
            end if
            if (.not.userules) x_max_ord2 = maxval(maxord)
            do ord2=1,min(x_max_ord2,ord)
              write(hubname(13:13),'(i1)') ord2
              call set_dependency(trim(optname),
     &                            trim(hubname),tgt_info)
            end do
          end do
          call opt_parameters(-1,parameters,ncat,nint)
          call set_rule(trim(optname),ttype_frm,OPTIMIZE,
     &                  labels,ncat+nint+1,1,
     &                  parameters,1,tgt_info)
        end do
        deallocate(ifreq,ifreqnew)
        formname(1:11) = 'LAG(n)     '
        write(formname(5:5),'(i1)') ord+1
      end do

      if (setr12) then
        ! optimize formulae for BV- and CV-intermediates
        optname(1:len_short) = ' '
        optname(1:6) = 'OPT_BV'
        optname2(1:len_short) = ' '
        optname2(1:6) = 'OPT_CV'
        defmelname(1:len_short) = ' '
        defmelname(1:9) = 'DEF_ME_BV'
        defmelname2(1:len_short) = ' '
        defmelname2(1:11) = 'DEF_ME_V(1)'
        formname(1:len_short) = ' '
        formname(1:8) = 'BV__form'
        formname2(1:len_short) = ' '
        formname2(1:10) = 'RDGB(V)V_F'
        formname3(1:len_short) = ' '
        formname3(1:10) = 'RBRV(V)V_F'
        formname4(1:len_short) = ' '
        formname4(1:8) = 'BV_CABS_'
        formname5(1:len_short) = ' '
        formname5(1:8) = 'CV__form'
        defmelname3(1:len_short) = ' '
        defmelname3(1:14) = 'DEF_ME_RDGB(V)'
        defmelname4(1:len_short) = ' '
        defmelname4(1:14) = 'DEF_ME_RBRV(V)'
        defmelname5(1:len_short) = ' '
        defmelname5(1:9) = 'DEF_ME_CV'
        do ipop = 1,npop
          optname(6:7) = pop(ipop)%name//pop(ipop)%comp
          optname2(6:7) = pop(ipop)%name//pop(ipop)%comp
          defmelname(9:10) = pop(ipop)%name//pop(ipop)%comp
          defmelname2(8:12) = pop(ipop)%name//'(1)'//pop(ipop)%comp
          if (orb_info%ngas.eq.4) defmelname2(13:13) = 'h'
          defmelname3(13:15) = pop(ipop)%name//')'//pop(ipop)%comp
          defmelname4(13:15) = pop(ipop)%name//')'//pop(ipop)%comp
          defmelname5(9:10) = pop(ipop)%name//pop(ipop)%comp
          formname(2:3) = pop(ipop)%name//pop(ipop)%comp
          formname2(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          formname3(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          formname4(2:9) = pop(ipop)%name//'_CABS_'//pop(ipop)%comp
          formname5(2:3) = pop(ipop)%name//pop(ipop)%comp
          labels(1:20)(1:len_target_name) = ' '
          ! OPT_BV
          labels(1) = trim(optname)
          call add_target(trim(optname),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(optname),trim(defmelname),tgt_info)
          call set_dependency(trim(optname),mel_rint,tgt_info)
          if (treat_bv) then
            labels(2) = trim(formname2)
            labels(3) = trim(formname3)
            labels(4) = trim(formname4)
            ncat = 3
            call set_dependency(trim(optname),trim(defmelname2),
     &                          tgt_info)
            if (orb_info%ngas.eq.4) then
              defmelname2(13:13) = 'f'
              call set_dependency(trim(optname),trim(defmelname2),
     &                            tgt_info)
            end if
            call set_dependency(trim(optname),trim(defmelname3),
     &                          tgt_info)
            call set_dependency(trim(optname),trim(defmelname4),
     &                          tgt_info)
            call set_dependency(trim(optname),'FF-X-INT',tgt_info)
            call set_dependency(trim(optname),trim(formname2),tgt_info)
            call set_dependency(trim(optname),trim(formname3),tgt_info)
            call set_dependency(trim(optname),trim(formname4),tgt_info)
          else
            labels(2) = trim(formname)
            ncat = 1
            call set_dependency(trim(optname),defmelname2(1:12),
     &                          tgt_info)
            call set_dependency(trim(optname),trim(formname),tgt_info)
          end if
          call opt_parameters(-1,parameters,ncat,0)
          call set_rule(trim(optname),ttype_frm,OPTIMIZE,
     &                  labels,ncat+1,1,
     &                  parameters,1,tgt_info)
          ! OPT_CV
          labels(1) = trim(optname2)
          call add_target(trim(optname2),ttype_frm,.false.,tgt_info)
          call set_dependency(trim(optname2),trim(defmelname5),tgt_info)
          call set_dependency(trim(optname2),mel_rint,tgt_info)
          labels(2) = trim(formname5)
          ncat = 1
          call set_dependency(trim(optname2),defmelname2(1:12),
     &                        tgt_info)
          call set_dependency(trim(optname2),trim(formname5),tgt_info)
          call opt_parameters(-1,parameters,ncat,0)
          call set_rule(trim(optname2),ttype_frm,OPTIMIZE,
     &                  labels,ncat+1,1,
     &                  parameters,1,tgt_info)
        end do
      end if

      if (ntest.ge.100)
     &  write(luout,*) 'opt. formulae defined'

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! for H(0) (op_ham), a ME-list is already defined and imported

      ! ME_V(1)X, ME_V(1)Y, ME_V(1)Z
      opname(1:len_short) = ' '
      opname(1:4) = 'V(1)'
      melname(1:len_short) = ' '
      melname(1:7) = 'ME_V(1)'
      defmelname(1:len_short) = ' '
      defmelname(1:11) = 'DEF_ME_V(1)'
      do ipop = 1,npop
        opname(1:5) = pop(ipop)%name//'(1)'//pop(ipop)%comp
        melname(4:8) = pop(ipop)%name//'(1)'//pop(ipop)%comp
        defmelname(8:12) = pop(ipop)%name//'(1)'//pop(ipop)%comp
        opname(1:1) = pop(ipop)%name
        melname(4:4) = pop(ipop)%name
        defmelname(8:8) = pop(ipop)%name
        call add_target(trim(defmelname),ttype_opme,.false.,tgt_info)
        call set_dependency(trim(defmelname),trim(opname),tgt_info)
        ! (a) define
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = trim(melname)
        labels(2) = trim(opname)
        call me_list_parameters(-1,parameters,msc,0,
     &       pop(ipop)%isym,0,0,.false.)
        call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        ! (b) import
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = trim(melname)
        call import_parameters(-1,parameters,
     &                         pop(ipop)%int_name,'DALTON_SPECIAL')
        call set_rule(trim(defmelname),ttype_opme,IMPORT,
     &                labels,1,1,
     &                parameters,1,tgt_info)
      end do

      if (setr12) then
        if (orb_info%ngas.eq.4.and.treat_bv) then
          ! ME_V(1)Vh, ME_V(1)Vf
          opname(1:len_short) = ' '
          opname(1:6) = 'V(1)_h'
          melname(1:len_short) = ' '
          melname(1:9) = 'ME_V(1)_h'
          defmelname(1:len_short) = ' '
          defmelname(1:13) = 'DEF_ME_V(1)_h'
          half = .true.
          do side = 1,2
            if (.not.half) then
              opname(6:6) = 'f'
              melname(9:9) = 'f'
              defmelname(13:13) = 'f'
            end if
            do ipop = 1,npop
              opname(1:5) = pop(ipop)%name//'(1)'//pop(ipop)%comp
              melname(4:8) = pop(ipop)%name//'(1)'//pop(ipop)%comp
              defmelname(8:12) = pop(ipop)%name//'(1)'//pop(ipop)%comp
              call add_target(trim(defmelname),ttype_opme,.false.,
     &                        tgt_info)
              call set_dependency(trim(defmelname),trim(opname),
     &                            tgt_info)
              ! (a) define
              labels(1:20)(1:len_target_name) = ' '
              labels(1) = trim(melname)
              labels(2) = trim(opname)
              call me_list_parameters(-1,parameters,msc,0,
     &            pop(ipop)%isym,0,0,.false.)
              call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                      labels,2,1,
     &                      parameters,1,tgt_info)
              ! (b) import
              labels(1:10)(1:len_target_name) = ' '
              labels(1) = trim(melname)
              call import_parameters(-1,parameters,
     &                            pop(ipop)%int_name,'DALTON_SPECIAL')
              call set_rule(trim(defmelname),ttype_opme,IMPORT,
     &                      labels,1,1,
     &                      parameters,1,tgt_info)
            end do
            half = .false.
          end do
        end if

        ! ME_BVX, ME_BVY, ME_BVZ; ME_RDGB(V)V and ME_BRV(V)V
        opname(1:len_short) = ' '
        opname(1:2) = 'BV'
        melname(1:len_short) = ' '
        melname(1:5) = 'ME_BV'
        defmelname(1:len_short) = ' '
        defmelname(1:9) = 'DEF_ME_BV'
        opname2(1:len_short) = ' '
        opname2(1:7) = 'RDGB(V)'
        melname2(1:len_short) = ' '
        melname2(1:10) = 'ME_RDGB(V)'
        defmelname2(1:len_short) = ' '
        defmelname2(1:14) = 'DEF_ME_RDGB(V)'
        opname3(1:len_short) = ' '
        opname3(1:7) = 'RBRV(V)'
        melname3(1:len_short) = ' '
        melname3(1:10) = 'ME_RBRV(V)'
        defmelname3(1:len_short) = ' '
        defmelname3(1:14) = 'DEF_ME_RBRV(V)'
        opname4(1:len_short) = ' '
        opname4(1:2) = 'CV'
        melname4(1:len_short) = ' '
        melname4(1:5) = 'ME_CV'
        defmelname4(1:len_short) = ' '
        defmelname4(1:9) = 'DEF_ME_CV'
        do ipop = 1,npop
          opname(2:3) = pop(ipop)%name//pop(ipop)%comp
          melname(5:6) = pop(ipop)%name//pop(ipop)%comp
          defmelname(9:10) = pop(ipop)%name//pop(ipop)%comp
          opname2(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          melname2(9:11) = pop(ipop)%name//')'//pop(ipop)%comp
          defmelname2(13:15) = pop(ipop)%name//')'//pop(ipop)%comp
          opname3(6:8) = pop(ipop)%name//')'//pop(ipop)%comp
          melname3(9:11) = pop(ipop)%name//')'//pop(ipop)%comp
          defmelname3(13:15) = pop(ipop)%name//')'//pop(ipop)%comp
          opname4(2:3) = pop(ipop)%name//pop(ipop)%comp
          melname4(5:6) = pop(ipop)%name//pop(ipop)%comp
          defmelname4(9:10) = pop(ipop)%name//pop(ipop)%comp
          ! DEF_ME_BVV:
          call add_target(trim(defmelname),ttype_opme,.false.,tgt_info)
          call set_dependency(trim(defmelname),trim(opname),tgt_info)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(melname)
          labels(2) = trim(opname)
          call me_list_parameters(-1,parameters,msc,0,
     &         pop(ipop)%isym,0,0,.false.)
          call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                  labels,2,1,
     &                  parameters,1,tgt_info)
          ! DEF_ME_CVV:
          call add_target(trim(defmelname4),ttype_opme,.false.,tgt_info)
          call set_dependency(trim(defmelname4),trim(opname4),tgt_info)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(melname4)
          labels(2) = trim(opname4)
          call me_list_parameters(-1,parameters,msc,0,
     &         pop(ipop)%isym,0,0,.false.)
          call set_rule(trim(defmelname4),ttype_opme,DEF_ME_LIST,
     &                  labels,2,1,
     &                  parameters,1,tgt_info)
          ! DEF_ME_RDGB(V)V:
          call add_target(trim(defmelname2),ttype_opme,.false.,tgt_info)
          call set_dependency(trim(defmelname2),trim(opname2),tgt_info)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(melname2)
          labels(2) = trim(opname2)
          call me_list_parameters(-1,parameters,msc,0,
     &         pop(ipop)%isym,0,0,.false.)
          call set_rule(trim(defmelname2),ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          ! DEF_ME_RBRV(V)V:
          call add_target(trim(defmelname3),ttype_opme,.false.,tgt_info)
          call set_dependency(trim(defmelname3),trim(opname3),tgt_info)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = trim(melname3)
          labels(2) = trim(opname3)
          call me_list_parameters(-1,parameters,msc,0,
     &         pop(ipop)%isym,0,0,.false.)
          call set_rule(trim(defmelname3),ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
        end do
      end if

      if (ntest.ge.100)
     &  write(luout,*) 'ME_V and ME_BV lists defined'

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
        if (.not.userules) x_max_ord = maxval(maxord)
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
            call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                          maxord,ncnt)
            sym = 1
            freqsum = 0d0
            do digit = 1,ord
              write(opname(3+2*digit:4+2*digit),'(i2.2)') 
     &                                    ifreqnew(digit)
              sym = multd2h(sym,pop(cmp(ifreqnew(digit))%pop_idx)%isym)
              freqsum = freqsum + cmp(ifreqnew(digit))%freq
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

            call freq_parameters(-1,parameters,freqsum)
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
              call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                            maxord,ncnt)
              sym = 1
              do digit = 1,ord
                write(opname2(5+2*digit:6+2*digit),'(i2.2)') 
     &                                       ifreqnew(digit)
                sym = multd2h(sym,
     &                        pop(cmp(ifreqnew(digit))%pop_idx)%isym)
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
          call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                        maxord,ncnt)
          sym = 1
          do digit = 1,ord
            write(defmelname(14+2*digit:15+2*digit),'(i2.2)')
     &                                   ifreqnew(digit)
            write(melname(10+2*digit:11+2*digit),'(i2.2)') 
     &                                   ifreqnew(digit)
            write(lagname(7+2*digit:8+2*digit),'(i2.2)') ifreqnew(digit)
            sym = multd2h(sym,pop(cmp(ifreqnew(digit))%pop_idx)%isym)
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
            call print_result(ord,ifreqnew,mel_dummy(1),.true.,orb_info)
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

      if (.not.userules) then
        ! ME_LCCEQ(n)
        defmelname(1:len_short) = ' '
        defmelname(1:15) = 'DEF_ME_LCCEQ(0)'
        melname(1:len_short) = ' '
        melname(1:11) = 'ME_LCCEQ(0)'
        lagname(1:len_short) = ' '
        lagname(1:8) = 'LRESP(0)'
        call add_target(trim(defmelname),ttype_opme,.false.,
     &                                              tgt_info)
        call set_dependency(trim(defmelname),trim(lagname),tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = trim(melname)
        labels(2) = trim(lagname)
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
      end if

      ! Diagonal Preconditioner: one for each possible irrep
      do sym = 1,maxsym
        do op_par = 1,2
          if (op_par.eq.1) then
            op_parent = 'L'
          else
            op_parent = 'T'
          end if
          call me_list_label(mel_dia1,mel_dia,sym,0,0,0,.false.)
          call add_target(trim(mel_dia1)//op_parent,ttype_opme,.false.,
     &                    tgt_info)
          call set_dependency(trim(mel_dia1)//op_parent,mel_ham,
     &                        tgt_info)
          call set_dependency(trim(mel_dia1)//op_parent,
     &                        op_dia//'_'//op_parent,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = trim(mel_dia1)//op_parent
          labels(2) = op_dia//'_'//op_parent
          call me_list_parameters(-1,parameters,
     &         0,0,sym,0,0,.false.)
          call set_rule(trim(mel_dia1)//op_parent,ttype_opme,
     &                  DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          labels(1) = trim(mel_dia1)//op_parent
          labels(2) = mel_ham
          call set_rule(trim(mel_dia1)//op_parent,ttype_opme,
     &                  PRECONDITIONER,
     &                  labels,2,1,
     &                  parameters,0,tgt_info)
        end do
      end do

      if (opt) then
        ! B-Matrix for Preconditioner (total symmetric)
        call add_target(me_bprc,ttype_opme,.false.,tgt_info)
        call set_dependency(me_bprc,op_bprc,tgt_info)
        call set_dependency(me_bprc,eval_r12_inter,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = me_bprc
        labels(2) = op_bprc
        call me_list_parameters(-1,parameters,
     &         msc,0,1,0,0,.false.)
        call set_rule(me_bprc,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = me_bprc
        labels(2) = mel_b_inter
        labels(3) = mel_bh_inter
        call add_parameters(-1,parameters,
     &       2,(/1d0,1d0/),2)
        call set_rule(me_bprc,ttype_opme,ADD,
     &                labels,3,1,
     &                parameters,1,tgt_info)
        ! X-Matrix for Preconditioner (total symmetric)
        call add_target(me_xprc,ttype_opme,.false.,tgt_info)
        call set_dependency(me_xprc,op_xprc,tgt_info)
        call set_dependency(me_bprc,eval_r12_inter,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = me_xprc
        labels(2) = op_xprc
        call me_list_parameters(-1,parameters,
     &         msc,0,1,0,0,.false.)
        call set_rule(me_xprc,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = me_xprc
        labels(2) = mel_x_inter
        call add_parameters(-1,parameters,
     &       1,(/1d0/),2)
        call set_rule(me_xprc,ttype_opme,ADD,
     &                labels,2,1,
     &                parameters,1,tgt_info)
      end if

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

      if (setr12) then
        ! ME_SX(n)w: ME for vector times Metric
        defmelname(1:len_short) = ' '
        defmelname(1:12) = 'DEF_ME_SX(n)'
        melname(1:len_short) = ' '
        melname(1:8) = 'ME_SX(n)'
        opname(1:len_short) = ' '
        opname(1:5) = 'SX(n)'
        do op_par = 1,2
          if (op_par.eq.1) then
            opname(2:2) = 'L'
            x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
          else
            opname(2:2) = 'T'
            x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
          end if
          if (.not.userules) x_max_ord = maxval(maxord)
          do ord = op_par-1,x_max_ord
            write(opname(4:4),'(i1)') ord
            opname(6:len_short) = ' '
            allocate(ifreq(ord),ifreqnew(ord))
            ifreq = 0
            set_zero = ord.eq.0.and.opt
            do while (next_comb(ifreq,ord,maxord,ncnt).or.set_zero)
              set_zero = .false.
              call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                            maxord,ncnt)
              sym = 1
              do digit = 1,ord
                write(opname(4+2*digit:5+2*digit),'(i2.2)') 
     &                                      ifreqnew(digit)
                sym = multd2h(sym,
     &                        pop(cmp(ifreqnew(digit))%pop_idx)%isym)
              end do
              melname(4:len_short) = opname(1:len_short-3)
              defmelname(8:len_short) = opname(1:len_short-7)
              if (idx_target(trim(defmelname),tgt_info).gt.0) cycle
              call add_target(trim(defmelname),ttype_opme,.false.,
     &                        tgt_info)
              call set_dependency(trim(defmelname),trim(opname),
     &                            tgt_info)
              labels(1:20)(1:len_target_name) = ' '
              labels(1) = trim(melname)
              labels(2) = trim(opname)
              call me_list_parameters(-1,parameters,
     &             msc,0,sym,0,0,.false.)
              call set_rule(trim(defmelname),ttype_opme,DEF_ME_LIST,
     &                      labels,2,1,
     &                      parameters,1,tgt_info)
            end do
            deallocate(ifreq,ifreqnew)
          end do
        end do
      end if

      if (ntest.ge.100)
     &  write(luout,*) 'me-lists defined'

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      if (setr12) then
        ! evaluate BV-Intermediate
        optname(1:len_short) = ' '
        optname(1:6) = 'OPT_BV'
        evalname(1:len_short) = ' '
        evalname(1:7) = 'EVAL_BV'
        do ipop = 1,npop
          optname(6:7) = pop(ipop)%name//pop(ipop)%comp
          evalname(7:8) = pop(ipop)%name//pop(ipop)%comp
          call add_target(trim(evalname),ttype_gen,.false.,tgt_info)
          call set_dependency(trim(evalname),trim(optname),tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = trim(optname)
          call set_rule(trim(evalname),ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)
        end do
        ! evaluate CV-Intermediate
        optname(1:len_short) = ' '
        optname(1:6) = 'OPT_CV'
        evalname(1:len_short) = ' '
        evalname(1:7) = 'EVAL_CV'
        do ipop = 1,npop
          optname(6:7) = pop(ipop)%name//pop(ipop)%comp
          evalname(7:8) = pop(ipop)%name//pop(ipop)%comp
          call add_target(trim(evalname),ttype_gen,.false.,tgt_info)
          call set_dependency(trim(evalname),trim(optname),tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = trim(optname)
          call set_rule(trim(evalname),ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)
        end do
      end if

      ! solve CC-equations 
      call add_target('SOLVE_T(0)',ttype_gen,.false.,tgt_info)
      if (setr12)
     &       call set_dependency('SOLVE_T(0)',eval_r12_inter,tgt_info)
      call set_dependency('SOLVE_T(0)','OPT_T(0)',tgt_info)
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
      call set_dependency('SOLVE_T(0)',trim(mel_dia1)//'T',tgt_info)
      call solve_parameters(-1,parameters,2, 1,1,'DIA')
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_T(0)'
      labels(2) = 'ME_O(0)_L'
      labels(3) = trim(mel_dia1)//'T'
      if (userules) then
        labels(4) = 'ME_LRESP(0)'
      else
        labels(4) = 'ME_LCCEQ(0)'
      end if
      labels(5) = 'OPT_T(0)'
      ilabels = 5
      if (opt) then
        call set_dependency('SOLVE_T(0)',me_bprc,tgt_info)
        call set_dependency('SOLVE_T(0)',me_xprc,tgt_info)
        call solve_parameters(-1,parameters,2,1,1,'MIX')
        labels(6) = me_bprc
        labels(7)= me_xprc
        labels(8)= mel_ham
        ilabels = 8
      end if
      if (restart.eq.0) 
     &       call set_rule('SOLVE_T(0)',ttype_opme,SOLVENLEQ,
     &       labels,ilabels,2,
     &       parameters,2,tgt_info)

      if (ntest.ge.100)
     &  write(luout,*) 'define solvers for linear equations'

      ! solve linear equations for X(n)w 
      ! (X=T,L, n=0,...,x_max_ord, not for T(0))
      do op_par = 1,2
        if (op_par.eq.1) then
          op_parent = 'L'
          defmelname(1:11) = 'DEF_ME_L(n)'
          melname(1:7) = 'ME_L(n)'
          opname(1:6) = 'O(n)LT'
          optname(1:8) = 'OPT_L(n)'
          solvename(1:10) = 'SOLVE_L(n)'
          x_max_ord = int((real(maxval(maxord))-2)/2+0.6)
          x_max_ord2 = int((real(restart-1)-2)/2+0.6)
        else
          op_parent = 'T'
          defmelname(1:11) = 'DEF_ME_T(n)'
          melname(1:7) = 'ME_T(n)'
          opname(1:6) = 'O(n)LL'
          optname(1:8) = 'OPT_T(n)'
          solvename(1:10) = 'SOLVE_T(n)'
          x_max_ord = int((real(maxval(maxord))-1)/2+0.6)
          x_max_ord2 = int((real(restart-1)-1)/2+0.6)
        end if
        if (.not.userules) x_max_ord = maxval(maxord)
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
            call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                          maxord,ncnt)
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
              sym = multd2h(sym,pop(cmp(ifreqnew(digit))%pop_idx)%isym)
            end do
            if (idx_target(trim(solvename),tgt_info).gt.0) cycle
            call add_target(trim(solvename),ttype_gen,.false.,
     &                          tgt_info)
            if (setr12)
     &              call set_dependency(trim(solvename),
     &              eval_r12_inter,tgt_info)
            if (setr12.and.ord.gt.0) then
              call set_dependency(trim(solvename),
     &              'HUB_EVAL_BV',tgt_info)
              call set_dependency(trim(solvename),
     &              'HUB_EVAL_CV',tgt_info)
            end if
            call set_dependency(trim(solvename),trim(defmelname),
     &                          tgt_info)
            call set_dependency(trim(solvename),trim(optname),
     &                          tgt_info)
            call me_list_label(mel_dia1,mel_dia,sym,0,0,0,.false.)
            call solve_parameters(-1,parameters,2, 1,1,'DIA')
            call set_dependency(trim(solvename),
     &                          trim(mel_dia1)//op_parent,tgt_info)
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
            if (setr12.and.(ord.gt.0.or.opt.and.op_par.eq.1)) then  ! Metric
              opname2(1:len_short-2) = 'S'//melname(4:len_short)
            else                           ! identity
              opname2(1:len_short-3) = melname(4:len_short)
            end if
            labels(1:20)(1:len_target_name) = ' '
            labels(1) = trim(melname)
            labels(2) = trim(mel_dia1)//op_parent
            labels(3) = trim(opname)
            labels(4) = trim(opname2)
            opname(5:5) = 'R'
            labels(5) = trim(opname)
            opname(5:5) = 'L'
            labels(6) = trim(optname)
            ilabels = 6
            if (opt) then
              call set_dependency(trim(solvename),me_bprc,tgt_info)
              call set_dependency(trim(solvename),me_xprc,tgt_info)
              call solve_parameters(-1,parameters,2,1,1,'MIX')
              labels(7) = me_bprc
              labels(8)= me_xprc
              labels(9)= mel_ham
              ilabels = 9
            end if
            if (ord.gt.x_max_ord2.or.op_par+restart.le.2)
     &            call set_rule(trim(solvename),ttype_opme,SOLVELEQ,
     &            labels,ilabels,1,
     &            parameters,2,tgt_info)
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
          call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                        maxord,ncnt)
          sym = 1
          freqsum = 0d0
          if (ord.gt.0) then
            icnt = (ifreq(1)-1)/maxval(maxord)+1
          else
            icnt = 1
          end if
          do digit = 1,ord
            write(evalname(10+2*digit:11+2*digit),'(i2.2)')
     &                                 ifreqnew(digit)
            write(defmelname(14+2*digit:15+2*digit),'(i2.2)')
     &                                 ifreqnew(digit)
            write(optname(11+2*digit:12+2*digit),'(i2.2)') 
     &                                 ifreqnew(digit)
            write(melname(10+2*digit:11+2*digit),'(i2.2)') 
     &                                 ifreqnew(digit)
            sym = multd2h(sym,pop(cmp(ifreqnew(digit))%pop_idx)%isym)
            freqsum = freqsum + cmp(ifreqnew(digit))%freq
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
              if (eval_dipmom(cmp(ifreqnew(1))%pop_idx)) cycle
              eval_dipmom(cmp(ifreqnew(1))%pop_idx) = .true.
            end if
            call add_target(trim(evalname),ttype_gen,.true.,tgt_info)
            if (setr12)
     &           call set_dependency(trim(evalname),eval_r12_inter,
     &                                              tgt_info)
            if (setr12.and.ord.gt.0) then
              call set_dependency(trim(evalname),
     &              'HUB_EVAL_BV',tgt_info)
              call set_dependency(trim(evalname),
     &              'HUB_EVAL_CV',tgt_info)
            end if
            ! first solve T(k), L(k) using (2n+1) and (2n+2) rules
            hubname(11:11) = 'T'
            x_max_ord = int((real(ord)-1)/2+0.6)
            if (.not.userules) x_max_ord = ord
            write(hubname(13:13),'(i1)') x_max_ord
            call set_dependency(trim(evalname),hubname(1:14),tgt_info)
            if (maxord(icnt).gt.0.or..not.userules) then
              hubname(11:11) = 'L'
              x_max_ord = int((real(ord)-2)/2+0.6)
              if (.not.userules) x_max_ord = ord
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
            if (restart.gt.ord) cycle
            call set_rule(trim(evalname),ttype_opme,EVAL,
     &           labels,1,0,
     &           parameters,0,tgt_info)
            ! print result
            call ord_parameters(-1,parameters,ord,0,ifreq)
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
        if (.not.userules) x_max_ord = maxval(maxord)
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
            call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                          maxord,ncnt)
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
      do ipop = 1,npop
        defmelname(8:12) = pop(ipop)%name//'(1)'//pop(ipop)%comp
        call set_dependency('HUB_DEFME_V(1)',trim(defmelname),
     &                                    tgt_info)
      end do

      if (setr12) then
        ! hub for DEF_ME_BV/CV-dependencies and EVAL_BV/CV-dependencies
        defmelname(1:len_short) = ' '
        defmelname(1:9) = 'DEF_ME_BV'
        defmelname2(1:len_short) = ' '
        defmelname2(1:9) = 'DEF_ME_CV'
        evalname(1:len_short) = ' '
        evalname(1:7) = 'EVAL_BV'
        evalname2(1:len_short) = ' '
        evalname2(1:7) = 'EVAL_CV'
        call add_target('HUB_DEFME_BV',ttype_gen,.false.,tgt_info)
        call add_target('HUB_DEFME_CV',ttype_gen,.false.,tgt_info)
        call add_target('HUB_EVAL_BV',ttype_gen,.false.,tgt_info)
        call add_target('HUB_EVAL_CV',ttype_gen,.false.,tgt_info)
        do ipop = 1,npop
          defmelname(9:10) = pop(ipop)%name//pop(ipop)%comp
          defmelname2(9:10) = pop(ipop)%name//pop(ipop)%comp
          evalname(7:8) = pop(ipop)%name//pop(ipop)%comp
          evalname2(7:8) = pop(ipop)%name//pop(ipop)%comp
          call set_dependency('HUB_DEFME_CV',trim(defmelname2),
     &                                      tgt_info)
          call set_dependency('HUB_EVAL_CV',trim(evalname2),
     &                                      tgt_info)
          if (pop(ipop)%isym.ne.1.and.
     &        r12op.eq.0) cycle
          call set_dependency('HUB_DEFME_BV',trim(defmelname),
     &                                      tgt_info)
          call set_dependency('HUB_EVAL_BV',trim(evalname),
     &                                      tgt_info)
        end do
      end if

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
        if (.not.userules) x_max_ord = maxval(maxord)
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
            call redundant_comb(ifreq,ifreqnew,cmp(1:ncmp)%redun,ord,
     &                          maxord,ncnt)
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
            if (.not.userules) x_max_ord2 = maxval(maxord)
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

      deallocate(maxord,evaluate,eval_dipmom)

      return
      end
