*----------------------------------------------------------------------*
      subroutine set_ic_mrcc_response_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed in MRCC excited state calculations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'multd2h.h'
      include 'def_orbinf.h'
      include 'opdim.h'
      include 'ifc_input.h'
      include 'ifc_targets.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     min_rank, max_rank,
     &     isim, ncat, nint, icnt, ncnt,
     &     s2, isym, ms, msc, lr_opt, sym_arr(8),
     &     s2_r, ms_r, isym_r, msc_r,
     &     s2_0, ms_0, isym_0, msc_0 
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     me_label_l, me_label_r, me_label_lr,
     &     labels(20)
      character(len=256),dimension(19) ::
     &     string   
      character(len=8) :: method='       '   
     

      ! skip this section if not requested
      ncnt = is_keyword_set('calculate.excitation')
      if (ncnt.eq.0) return

      if (iprlvl.gt.0)
     &    write(lulog,*) 'setting targets for icMRCC excited states ...'

      ! S2, MS and MS combi of reference state
      s2_0  = orb_info%imult
      isym_0 = orb_info%lsym
      ms_0 = orb_info%ims
      if (ms_0.eq.0.and.mod(orb_info%imult-1,4).eq.0) then
        msc_0 = 1
      else if (ms_0.eq.0.and.mod(orb_info%imult+1,4).eq.0) then
        msc_0 = -1
      else
        msc_0 = 0
      end if
 
      write(lulog,*) 'IRREP, S2, MS of reference state: ',isym_0, s2_0, 
     &               ms_0
      write(lulog,*) '  factor for spin-combinations: ',msc_0


        call get_argument_value('method.MRCC.excite','method',
     &       str=method)
  
        print*, trim(method)
        select case(trim(method))
        case('LR')
          write(lulog,*) 'doing ic_mrcc response in LR framework'
          lr_opt = 1
        case('EOM')
          write(lulog,*) 'doing ic_mrcc response in EOM framework'
          lr_opt = 2
        case default
          call quit(0,'set_ic_mrcc_response_target',
     &    'invalid method: '//trim(method))
        end select


      call add_target3([character(len=80) ::
     &      'target RSPNS_OP(                                  ',
     &      'depend (OMG,D,L,DENS,C0,A_C0,T)                   ',
     &      'CLONE_OPERATOR(label=R_prime_q,template=T)        ',
     &      'CLONE_OPERATOR(label=R_q,template=T)              ',
     &      'CLONE_OPERATOR(label=R_mu,template=C0)            ',
     &      'CLONE_OPERATOR(label=R_mu_x,adjoint=T,template=C0)',
     &      'CLONE_OPERATOR(label=AR_rspns_q,template=OMG)     ',
     &      'CLONE_OPERATOR(label=AR_rspns_mu,template=A_C0)   ',
     &      'CLONE_OPERATOR(label=SR_rspns_q,template=OMG)     ',
     &      'CLONE_OPERATOR(label=SR_rspns_mu,template=R_mu)   ',
     &      'DEF_SCALAR(label=den12)                           ',
     &      ')                                                 '],
     &      tgt_info)

      if(lr_opt.eq.1) then
      call add_target3([character(len=80) ::
     &      'target FORM_AR_RSPNS_Q(                              ',
     &      'depend (RSPNS_OP,F_OMG)                              ',
     &      'DERIVATIVE(label_res=F_AR_rspns_q,label_in=F_OMG,    ',
     &      'op_res=AR_rspns_q,op_deriv=(T,C0,C0^+),              ',
     &      'op_mult=(R_q,R_mu,R_mu_x))                           ',
     &      'REPLACE(label_res=F_AR_rspns_q,label_in=F_AR_rspns_q,',
     &      'op_list=(R_mu_x,R_mu^+)                              ',
     &      ')                                                    ',
     &      'PRINT_FORMULA(label=F_AR_rspns_q)                    ',
     &      ')                                                    '],
     &      tgt_info)
   
       else
      call add_target3([character(len=80) ::
     &  '    target FORM_AR_RSPNS_Q(                                 ',
     &  '    depend (RSPNS_OP,H0)                                    ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=T,op_res=den12,',
     &  '                 operators=(C0^+,L,H,R_mu),idx_sv=(1,2,3,4))',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '             operators=(C0^+,L,H,T,R_mu),idx_sv=(1,2,3,4,5),',
     &  '                      connect=(3,4))                        ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '             operators=(C0^+,L,T,H,R_mu),idx_sv=(1,2,3,4,5),',
     &  '                      connect=(3,4),fac=-1.0d0)             ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '         operators=(C0^+,L,H,T,T,R_mu),idx_sv=(1,2,3,4,5,6),',
     &  '                      fix_vtx=t,fac=0.5d0)                  ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '         operators=(C0^+,L,T,H,T,R_mu),idx_sv=(1,2,3,4,5,6),',
     &  '                      fix_vtx=t,fac=-1.0d0)                 ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '         operators=(C0^+,L,T,T,H,R_mu),idx_sv=(1,2,3,4,5,6),',
     &  '                      fix_vtx=t,fac=0.5d0)                  ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '         operators=(C0^+,L,H,R_q,C0),idx_sv=(1,2,3,4,5),    ',
     &  '                      connect=(3,4))                        ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '         operators=(C0^+,L,R_q,H,C0),idx_sv=(1,2,3,4,5),    ',
     &  '                      connect=(3,4),fac=-1.0d0)             ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '         operators=(C0^+,L,H,T,R_q,C0),idx_sv=(1,2,3,4,5,6))',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '         operators=(C0^+,L,R_q,H,T,C0),idx_sv=(1,2,3,4,5,6),',
     &  '                      fac=-1.0d0)                           ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '         operators=(C0^+,L,T,H,R_q,C0),idx_sv=(1,2,3,4,5,6),',
     &  '                      fac=-1.0d0)                           ',
     &  '    EXPAND_OP_PRODUCT(label=F_AR_rspns_t,new=F,op_res=den12,',
     &  '         operators=(C0^+,L,R_q,T,H,C0),idx_sv=(1,2,3,4,5,6))',
     &  'SELECT_SPECIAL(label_res=F_AR_rspns_t,label_in=F_AR_rspns_t,',
     &  '          type=nonzero,mode=sum)                            ',
     &  '    DERIVATIVE(label_res=F_AR_rspns_q,label_in=F_AR_rspns_t,',
     &  '      op_res=AR_rspns_q,op_deriv=L)                         ',
     &  '    PRINT_FORMULA(label=F_AR_rspns_q)                       ',
     &  '    )                                                       '
     &     ],tgt_info)

       endif

      ! trying fortran2003 style:
      !call add_target3((/
      call add_target3([character(len=80) ::
     &      'target F_prePPrint(',
     &      'depend (FORM_AR_RSPNS_Q,DEF_ME_INT_PP,H_PP)',
     &      'CLONE_OPERATOR(label=INT_PPr,template=INT_PP)',
     &      'REPLACE(label_res=F_prePPrint,label_in=F_AR_rspns_q,',
     &      'op_list=(H,H_PP))',
     &      'INVARIANT(label_res=F_prePPrint,label_in=F_prePPrint,',
     &      'op_res=AR_rspns_q,operators=H)',
     &      ')'],tgt_info)
C     !&      ')'/),tgt_info)

      !call add_target3((/
      call add_target3([character(len=80) ::
     &      'target F_PPrint(',
     &      'depend F_prePPrint',
     &      'DERIVATIVE(label_res=F_PPrint,label_in=F_prePPrint,',
     &      'op_res=INT_PPr,op_deriv=H_PP)',
     &      'PRINT_FORMULA(label=F_PPrint)',
     &      ')'],tgt_info)
C     !&      ')'/),tgt_info)
 
 
        if(lr_opt.eq.1) then
      !call add_target3((/
      call add_target3([character(len=80) ::
     &      'target FORM_AR_RSPNS_MU(',
     &      'depend(RSPNS_OP,F_OMG_C0,"E(MR)")',
     &      'DERIVATIVE(label_res=F_AR_rspns_mu,label_in=F_OMG_C0,',
     &      'op_res=AR_rspns_mu,op_deriv=(T,C0),',
     &      'op_mult=(R_q,R_mu))',
     &      'EXPAND_OP_PRODUCT(label=F_AR_rspns_mu,new=F,',
     &      'op_res=AR_rspns_mu,operators=(AR_rspns_mu,"E(MR)",',
     &      'R_mu,AR_rspns_mu),idx_sv=(1,2,3,1),fac=-1d0)',
     &      'PRINT_FORMULA(label=F_AR_rspns_mu)',
     &      ')'],tgt_info)
C     !&      ')'/),tgt_info)
        else
      !call add_target3((/
      call add_target3([character(len=80) ::
     &      'target FORM_AR_RSPNS_MU(',
     &      'depend(RSPNS_OP,F_E_C0,"E(MR)")',
     &      'DERIVATIVE(label_res=F_AR_rspns_c,label_in=F_E_C0,',
     &      'op_res=den12,op_deriv=C0,',
     &      'op_mult=R_mu)',
     &      'EXPAND_OP_PRODUCT(label=F_AR_rspns_c,new=F,op_res=den12,',
     &                 'operators=(C0^+,H,R_q,C0),idx_sv=(1,2,3,4),',
     &                 'connect=(2,3))',
     &      'EXPAND_OP_PRODUCT(label=F_AR_rspns_c,new=F,op_res=den12,',
     &              'operators=(C0^+,H,T,R_q,C0),idx_sv=(1,2,3,4,5))',
     & 'SELECT_SPECIAL(label_res=F_AR_rspns_c,label_in=F_AR_rspns_c,',
     &              'type=nonzero,mode=sum)',
     &      'DERIVATIVE(label_res=F_AR_rspns_mu,label_in=F_AR_rspns_c,',
     &                 'op_res=AR_rspns_mu,op_deriv=C0^+)',
     &      'EXPAND_OP_PRODUCT(label=F_AR_rspns_mu,new=F,',
     &      'op_res=AR_rspns_mu,operators=(AR_rspns_mu,"E(MR)",',
     &      'R_mu,AR_rspns_mu),idx_sv=(1,2,3,1),fac=-1d0)',
     &      'PRINT_FORMULA(label=F_AR_rspns_mu)',
     &      ')'],tgt_info)
C     !&      ')'/),tgt_info)
        endif 

        if(lr_opt.eq.1)then
      !call add_target3((/
      call add_target3([character(len=80) ::
     &      'target RSPNS_FORM(',
     &      'depend (Dtr,F_T,Ttr,RSPNS_OP,C0)',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=T,op_res=den12,',
     &      'operators=(den12,C0^+,L,R_q,C0,den12),',
     &      'idx_sv =(1,2,3,4,5,1)',
     &      ')',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=F,op_res=den12,',
     &         'operators=(den12,C0^+,L,R_q,T,C0,den12),',
     &         'idx_sv =(1,2,3,4,5,6,1),fac=0.5d0)', 
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=F,op_res=den12,',
     &         'operators=(den12,C0^+,L,T,R_q,C0,den12),',
     &         'idx_sv =(1,2,3,4,5,6,1),fac=-0.5d0)',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=F,op_res=den12,',
     &         'operators=(den12,C0^+,L,R_q,T,T,C0,den12),fix_vtx=t,',
     &         'idx_sv =(1,2,3,4,5,6,7,1),fac=0.1666666666666666d0)',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=F,op_res=den12,',
     &         'operators=(den12,C0^+,L,T,R_q,T,C0,den12),fix_vtx=t,',
     &         'idx_sv =(1,2,3,4,5,6,7,1),fac=-0.3333333333333333d0)',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=F,op_res=den12,',
     &         'operators=(den12,C0^+,L,T,T,R_q,C0,den12),fix_vtx=t,',
     &         'idx_sv =(1,2,3,4,5,6,7,1),fac=0.1666666666666666d0)',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=F,op_res=den12,',
     &         'operators=(den12,C0^+,L,R_q,T,T,T,C0,den12),fix_vtx=t,',
     &         'idx_sv =(1,2,3,4,5,6,7,8,1),fac=.04166666666666666d0)',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=F,op_res=den12,',
     &         'operators=(den12,C0^+,L,T,R_q,T,T,C0,den12),fix_vtx=t,',
     &         'idx_sv =(1,2,3,4,5,6,7,8,1),fac=-0.125d0)',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=F,op_res=den12,',
     &         'operators=(den12,C0^+,L,T,T,R_q,T,C0,den12),fix_vtx=t,',
     &         'idx_sv =(1,2,3,4,5,6,7,8,1),fac=0.125d0)',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=F,op_res=den12,',
     &         'operators=(den12,C0^+,L,T,T,T,R_q,C0,den12),fix_vtx=t,',
     &         'idx_sv =(1,2,3,4,5,6,7,8,1),fac=-.04166666666666666d0)',
     &      'SELECT_SPECIAL(label_res=F_den12,label_in=F_den12,',
     &              'type=nonzero,mode=sum)',
     &      'DERIVATIVE(label_res=F_SR_rspns_q,label_in=F_den12,',
     &      'op_res=SR_rspns_q,op_deriv=L)',
     &      'PRINT_FORMULA(label=F_SR_rspns_q)',
     &      'EXPAND_OP_PRODUCT(label=F_SR_rspns_mu,new=T,',
     &      'op_res=SR_rspns_mu,operators=(SR_rspns_mu,R_mu,',
     &      'SR_rspns_mu),idx_sv=(1,2,1))',
     &      'PRINT_FORMULA(label=F_SR_rspns_mu)',
     &      'INVARIANT(label_res=F_R_q,label_in=F_T,op_res=R_q,',
     &      'operators=D)',
     &      'REPLACE(label_res=F_R_q,label_in=F_R_q,',
     &      'op_list=(Ttr,R_prime_q)',
     &      ')',
     &      'PRINT_FORMULA(label=F_R_q)',
     &      ')'],tgt_info)
C     !&      ')'/),tgt_info)
 
         else 

      !call add_target3((/
      call add_target3([character(len=80) ::
     &      'target RSPNS_FORM(',
     &      'depend (Dtr,F_T,Ttr,RSPNS_OP,C0)',
     &      'EXPAND_OP_PRODUCT(label=F_den12,new=T,op_res=den12,',
     &      'operators=(den12,C0^+,L,R_q,C0,den12),',
     &      'idx_sv =(1,2,3,4,5,1)',
     &      ')',
     &      'DERIVATIVE(label_res=F_SR_rspns_q,label_in=F_den12,',
     &      'op_res=SR_rspns_q,op_deriv=L)',
     &      'PRINT_FORMULA(label=F_SR_rspns_q)',
     &      'EXPAND_OP_PRODUCT(label=F_SR_rspns_mu,new=T,',
     &      'op_res=SR_rspns_mu,operators=(SR_rspns_mu,R_mu,',
     &      'SR_rspns_mu),idx_sv=(1,2,1))',
     &      'PRINT_FORMULA(label=F_SR_rspns_mu)',
     &      'INVARIANT(label_res=F_R_q,label_in=F_T,op_res=R_q,',
     &      'operators=D)',
     &      'REPLACE(label_res=F_R_q,label_in=F_R_q,',
     &      'op_list=(Ttr,R_prime_q)',
     &      ')',
     &      'PRINT_FORMULA(label=F_R_q)',
     &      ')'],tgt_info)
C     &      ')'/),tgt_info)
          endif
              
      do icnt = 1, ncnt 
        call get_argument_value('calculate.excitation','sym',
     &       keycount=icnt,
     &       iarr=sym_arr)
        call get_argument_value('calculate.excitation','mult',
     &       keycount=icnt,
     &       ival=s2)
        ms = 0
        if (mod(s2,2).eq.0) ms = 1
        if (ms.eq.0.and.mod(s2,4).eq.1) msc = +1
        if (ms.eq.0.and.mod(s2,4).eq.3) msc = -1
        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle
          ! msc/isym is the final state MS/IRREP
          ! for the response we have to take care of the reference state symmetry
          isym_r = multd2h(isym,isym_0)
          if (s2.eq.s2_0) then
            s2_r = 1
          else if (abs(s2-s2_0).eq.2) then
            s2_r = 3
          else
            call quit(0,'set_ic_mrcc_response',
     &                  'cannot handle this S2 difference')
          end if
          ms_r = ms - ms_0
          msc_r = 0
          if (ms_r.eq.0.and.s2_r.eq.1) msc_r = +1
          if (ms_r.eq.0.and.s2_r.eq.3) msc_r = -1

c dbg
          write(lulog,*) 'isym, msc    ',isym, msc
          write(lulog,*) 'isym_r, msc_r',isym_r, msc_r
c dbg
          ! we will label most objects with isym and msc
          ! only the setup of the response lists needs the actual msc_r, isym_r

          call me_list_label(me_label_l,'ME_R_q',isym,0,0,msc,.false.)
          write(string(1),'("DEF_ME_LIST(list=",a,",operator=R_q,'//
     &   'irrep=",I1,",2ms=0,ab_sym=",I2,",min_rec=1,max_rec=",I3,")")')
     &     trim(me_label_l),isym_r,msc_r,sym_arr(isym)
          call me_list_label(me_label,'ME_INT_PPr',isym,0,0,msc,.false.)
          write(string(2),'("DEF_ME_LIST(list=",a,",operator=INT_PPr,'//
     &    'irrep=",I1,",2ms=0,ab_sym=",I2,")")')
     &     trim(me_label),isym_r,msc_r
          call me_list_label(me_label_r,'ME_R_mu',isym,0,0,msc,.false.)
          write(string(3),'("DEF_ME_LIST(list=",a,",operator=R_mu,'//
     &   'irrep=",I1,",2ms=0,ab_sym=",I2,",min_rec=1,max_rec=",I3,")")')
     &     trim(me_label_r),isym,msc,sym_arr(isym)
        call me_list_label(me_label_lr,'ME_R_prime_q',isym,0,0,msc,
     &      .false.)
        write(string(4),'("DEF_ME_LIST(list=",a,",operator=R_prime_q,'//
     &   'irrep=",I1,",2ms=0,ab_sym=",I2,",min_rec=1,max_rec=",I3,")")')
     &     trim(me_label_lr),isym_r,msc_r,sym_arr(isym)
       call me_list_label(me_label,'ME_AR_rspns_q',isym,0,0,msc,.false.)
          write(string(5),'("DEF_ME_LIST(list=",a,",'//
     &    'operator=AR_rspns_q,irrep=",I1,",2ms=0,ab_sym=",I2,")")')
     &     trim(me_label),isym_r,msc_r
      call me_list_label(me_label,'ME_AR_rspns_mu',isym,0,0,msc,.false.)
          write(string(6),'("DEF_ME_LIST(list=",a,",'//
     &    'operator=AR_rspns_mu,irrep=",I1,",2ms=0,ab_sym=",I2,")")')
     &     trim(me_label),isym,msc
       call me_list_label(me_label,'ME_SR_rspns_q',isym,0,0,msc,.false.)
          write(string(7),'("DEF_ME_LIST(list=",a,",'//
     &    'operator=SR_rspns_q,irrep=",I1,",2ms=0,ab_sym=",I2,")")')
     &     trim(me_label),isym_r,msc_r
      call me_list_label(me_label,'ME_SR_rspns_mu',isym,0,0,msc,.false.)
          write(string(8),'("DEF_ME_LIST(list=",a,",'//
     &    'operator=SR_rspns_mu,irrep=",I1,",2ms=0,ab_sym=",I2,")")')
     &     trim(me_label),isym,msc
        call me_list_label(medef_label,'ME_DIAG_q',isym,0,0,msc,.false.)
          write(string(9),'("DEF_ME_LIST(list=",a,",operator=DIA_T,'//
     &    'irrep=",I1,",2ms=0)")')
     &     trim(medef_label),isym_r
         call me_list_label(dia_label,'ME_DIAG_mu',isym,0,0,msc,.false.)
          write(string(10),'("DEF_ME_LIST(list=",a,",operator=DIA_C0,'//
     &    'irrep=",I1,",2ms=0)")')
     &     trim(dia_label),isym
          call me_list_label(mel_dia1,'ME_MINEN',isym,0,0,msc,.false.)
          write(string(11),'("DEF_ME_LIST(list=",a,",'//
     &    'operator=""E(MR)"",irrep=1,2ms=0)")')
     &     trim(mel_dia1)

          write(string(12),'("LIST_RSPNS_",I1,"_",I2.2)') isym,msc+1
          call add_target3([character(len=80) ::
     &        'target ', trim(string(12)), '(',
     &        'depend (RSPNS_OP,DEF_ME_C0,DEF_ME_Dtrdag,H0,DEF_ME_T,',
     &        'DIA_T,DIA_C0,"DEF_ME_E(MR)",F_prePPrint)',
     &        trim(string(1)),trim(string(2)),trim(string(3)),
     &        trim(string(4)),trim(string(5)),trim(string(6)),
     &        trim(string(7)),trim(string(8)),trim(string(9)),
     &        trim(string(10)),string(11),
     &        'ASSIGN_ME2OP(list="ME_E(MR)",operator="E(MR)")',
     &        ')'],tgt_info)     
C     !&        ')'/),tgt_info)     

 
        write(string(12),'("depend LIST_RSPNS_",I1,"_",I2.2)')isym,msc+1
         write(string(13),'("DIAG_CAL_q_",I1,"_",I2.2)') isym,msc+1
         write(string(17),'("EXTRACT_DIAG(list_res=",a,
     &              ",list_in=ME_A,mode=extend)")') trim(medef_label) 
         write(string(18),'("PRECONDITIONER(list_prc=",a,
     &              ",list_inp=ME_FREF)")') trim(medef_label)

         call add_target3([character(len=80) ::
     &        'target ', trim(string(13)), '(',
     &        'depend (EVAL_FREF,FOPT_Atr)',
     &        trim(string(12)),
     &        trim(string(18)),
*    &        'PRECONDITIONER(list_prc=ME_DIAG_q,list_inp=ME_FREF)',
     &        'ASSIGN_ME2OP(list=ME_Dtr,operator=Dtr)',
     &        'EVALUATE(form=FOPT_Atr)',
     &        string(17),            
*    &        'EXTRACT_DIAG(list_res=',medef_label,',list_in=ME_A,',
*    &        'mode=extend)',
     &        ')'],tgt_info)
          
         write(string(13),'("DIAG_CAL_mu_",I1,"_",I2.2)') isym,msc+1
         write(string(17),'("EXTRACT_DIAG(list_res=",a,
     &        ",list_in=",a,",mode=ext_act)")') 
     &        trim(dia_label),trim(mel_dia1) 
         write(string(18),'("PRECONDITIONER(list_prc=",a,
     &                     ",list_inp=H0,mode=dia-H)")') trim(dia_label)
         write(string(19),'("SCALE_COPY(list_res=",a,
     &            ",list_inp=""ME_E(MR)"",fac=-1d0)")') trim(mel_dia1)
         call add_target3([character(len=80) ::
     &         'target',trim(string(13)),'(',
     &         trim(string(12)),
     &         trim(string(18)),
     &         trim(string(19)),string(17),                      ! trim(string(17)) FIXME 
     &         ')'],tgt_info)


         write(string(13),'("RSPNS_OPT_",I1,"_",I2.2)') isym,msc+1
         write(string(14),'("OPTIMIZE(label_opt=RSPNS_OPT_",I1,"_",I2.2,
     &   ",labels_in=(F_AR_rspns_q,F_AR_rspns_mu,F_SR_rspns_q,",
     &    "F_SR_rspns_mu,F_R_q),interm=F_PPrint)")') isym,msc+1

           call add_target3([character(len=80) ::
     &     'target',trim(string(13)),'(',
     &     'depend (RSPNS_FORM,FORM_AR_RSPNS_Q,',
     &     'FORM_AR_RSPNS_MU,"DEF_ME_E(MR)",F_PPrint)',
     &     string(12),string(14),
C    &    'OPTIMIZE(label_opt=',trim(string(14)),'labels_in=(F_AR_rspns_q,',
C    &     'F_AR_rspns_mu,F_SR_rspns_q,F_SR_rspns_mu,F_R_q),',
C    &     'interm=F_PPrint)',
     &     ')'],tgt_info)

       write(string(12),'("mode=""TRF DIA"",n_roots=",I1,")")')
     &                   sym_arr(isym)
       write(string(13),'("MY_TARGET_",I1,"_",I2.2)') isym,msc+1
       write(string(14),'("depend RSPNS_OPT_",I1,"_",I2.2)') isym,msc+1
       write(string(15),'("depend DIAG_CAL_q_",I1,"_",I2.2)') isym,msc+1
       write(string(16),'("depend DIAG_CAL_mu_",I1,"_",I2.2)')isym,msc+1
       write(string(17),'("list_opt=(",a,",",a,"),")')
     &                    trim(me_label_l),trim(me_label_r)
       write(string(18),'("list_prc=(",a,",",a,"),")')
     &                   trim( medef_label),trim(dia_label)
       write(string(19),'("form=RSPNS_OPT_",I1,"_",I2.2,
     &          ",list_spc=(",a,",ME_Dtr,ME_Dtrdag),")')
     &                     isym,msc+1,trim(me_label_lr)

       call add_target3([character(len=80) ::
     &     'target', trim(string(13)), '(',
     &     'required',
     &     trim(string(14)),trim(string(15)),trim(string(16)),
     &     'SOLVE_EVP(',
     &     trim(string(17)),
     &     string(18),
     &     'op_mvp=(AR_rspns_q,AR_rspns_mu),',
     &     'op_svp=(SR_rspns_q,SR_rspns_mu),',
     &     string(19),
     &     string(12),
     &     ')'],tgt_info)


        end do
      end do
       
      end subroutine
