*----------------------------------------------------------------------*
      subroutine set_ic_mrcc_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets for internally contracted MRCC
*
*     matthias 2010
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'ifc_targets.h'
      include 'def_orbinf.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'

      integer, parameter ::
     &     ntest = 100

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ndef, occ_def(ngastp,2,60),
     &     icnt,
     &     msc, maxexc, ip, ih, ivv, iv, ivv2, jvv,
     &     minh, maxh,
     &     minp, maxp, maxv, maxvv, minexc, maxcom, maxcom_en,
     &     n_t_cls, i_cls,
     &     n_tred_cls, len_form, optref, idef, ciroot,
     &     version(60), ivers, stndT(2,60), stndD(2,60), nsupT, nsupD
      logical ::
     &     pure_vv, update_prc, skip, preopt, singrm, first, cheap_prc
      character(len_target_name) ::
     &     dia_label, dia_label2,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character ::
     &     op_ht*3, f_ht*5, op_ht0to*6, f_ht0to*8, form_str*50,
     &     def_ht*10

      ! first set targets for CASSCF or uncontracted CI wave function
      ! (if not done already)
      if (.not.is_keyword_set('method.MR').gt.0)
     &      call quit(1,'set_ic_mrcc_targets',
     &      'MRCC requires MR wave function')

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting multireference targets #3...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

      ! get minimum and maximum numbers of excitations, holes, particles,
      ! valence-valence excitations
      call get_argument_value('method.MR','minh',
     &     ival=minh)
      call get_argument_value('method.MR','maxh',
     &     ival=maxh)
      call get_argument_value('method.MR','minp',
     &     ival=minp)
      call get_argument_value('method.MR','maxp',
     &     ival=maxp)
      call get_argument_value('method.MR','maxv',
     &     ival=maxv)
      call get_argument_value('method.MR','maxvv',
     &     ival=maxvv)
      call get_argument_value('method.MR','minexc',
     &     ival=minexc)
      call get_argument_value('method.MR','maxexc',
     &     ival=maxexc)
      if (maxh.lt.0) maxh = maxexc
      if (maxp.lt.0) maxp = maxexc
      if (maxv.lt.0) maxv = 2*maxexc
      if (maxvv.lt.0) maxvv = maxexc
      call get_argument_value('method.MR','pure_vv',
     &     lval=pure_vv)
      call get_argument_value('method.MR','ciroot',
     &     ival=ciroot)
      call get_argument_value('method.MR','cheap_prc',
     &     lval=cheap_prc)
      call get_argument_value('calculate.solve.non_linear','optref',
     &     ival=optref)
      call get_argument_value('calculate.solve.non_linear','update_prc',
     &     lval=update_prc)
      call get_argument_value('calculate.solve.non_linear','preopt',
     &     lval=preopt)
      call get_argument_value('calculate.solve.non_linear','singrm',
     &     lval=singrm)
      call get_argument_value('method.MRCC','maxcom_res',
     &     ival=maxcom)
      call get_argument_value('method.MRCC','maxcom_en',
     &     ival=maxcom_en)

      if (ntest.ge.100) then
        write(luout,*) 'maxcom_en  = ', maxcom_en
        write(luout,*) 'maxcom_res = ', maxcom
        write(luout,*) 'preopt     = ', preopt
      end if
      
*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define particle conserving deexcitation operator L (for icMRCC)
      call add_target('L',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
c dbg
c            if (ih.gt.ip) cycle
            if (singrm.and.ivv.eq.0.and.(ih+ip.eq.1.or.ih.eq.ip)) cycle
c dbgend
cmh         no blocks which can be modeled by a contraction of two operators
c            if (ip.ge.1.and.ih.ge.1) cycle
cmh end
c            ! only hole-particle singles
c            if (max(ip,ih)+ivv.eq.1.and.(ip.ne.1.or.ih.ne.1)) cycle
            ndef = ndef + 1
            occ_def(IHOLE,1,ndef) = ih
            occ_def(IPART,2,ndef) = ip
            occ_def(IVALE,2,ndef) = max(ih-ip,0) + ivv
            occ_def(IVALE,1,ndef) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      if (singrm) then ! define conventional blocks separately
        do ip = minp, maxp ! ih=ip
          if (ip.lt.minh.or.ip.gt.maxh) cycle
          if (ip.lt.minexc) cycle
          ndef = ndef + 1 
          occ_def(IHOLE,1,ndef) = ip
          occ_def(IPART,2,ndef) = ip
        end do
      end if
      n_t_cls = ndef
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('L',ttype_op,DEF_OP_FROM_OCC,
     &              'L',1,1,
     &              parameters,2,tgt_info)

      ! define particle conserving excitation operator T (for icMRCC)
      call add_target('T',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      nsupT=0
      do ip = minp, maxp
        do ih = minh, maxh
          first = .true.
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
c dbg
c            if (ih.gt.ip) cycle
            if (singrm.and.ivv.eq.0.and.(ih+ip.eq.1.or.ih.eq.ip)) cycle
c dbgend
cmh         no blocks which can be modeled by a contraction of two operators
c            if (ip.ge.1.and.ih.ge.1) cycle
cmh end
c            ! only hole-particle singles
c            if (max(ip,ih)+ivv.eq.1.and.(ip.ne.1.or.ih.ne.1)) cycle
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef) = max(ip-ih,0) + ivv
            if (first) then
              nsupT = nsupT + 1
              stndT(1,nsupT) = ndef
              first = .false.
            end if
            stndT(2,nsupT) = ndef
          end do
        end do
      end do
      if (singrm) then ! define conventional blocks separately
        do ip = minp, maxp ! ih=ip
          if (ip.lt.minh.or.ip.gt.maxh) cycle
          if (ip.lt.minexc) cycle
          ndef = ndef + 1
          occ_def(IHOLE,2,ndef) = ip
          occ_def(IPART,1,ndef) = ip
          nsupT = nsupT + 1
          stndT(1,nsupT) = ndef
          stndT(2,nsupT) = ndef
        end do
      end if
      n_t_cls = ndef
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('T',ttype_op,DEF_OP_FROM_OCC,
     &              'T',1,1,
     &              parameters,2,tgt_info)
cmh      ! T excitation operator
cmh      call add_target2('T',.false.,tgt_info)
      call set_dependency('T','L',tgt_info)
cmh      call set_rule2('T',CLONE_OP,tgt_info)
cmh      call set_arg('T',CLONE_OP,'LABEL',1,tgt_info,
cmh     &     val_label=(/'T'/))
cmh      call set_arg('T',CLONE_OP,'TEMPLATE',1,tgt_info,
cmh     &     val_label=(/'L'/))
cmh      call set_arg('T',CLONE_OP,'ADJOINT',1,tgt_info,
cmh     &     val_log=(/.true./))
c      call set_rule2('T',SET_ORDER,tgt_info)
c      call set_arg('T',SET_ORDER,'LABEL',1,tgt_info,
c     &     val_label=(/'T'/))
c      call set_arg('T',SET_ORDER,'SPECIES',1,tgt_info,
c     &             val_int=(/1/))

      ! we need to set "super-block" information of metric
      ! since metric is defined in set_ic_mrci_targets,
      ! we just run over the same loops here:
      ndef = 0
      nsupD = 0
      do ip = minp, maxp
        do ih = minh, maxh
          first = .true.
          do ivv = min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv),0,-1
           do jvv =min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv),0,-1
            if (ip.eq.ih.and.ip.eq.maxexc) cycle
            if (.not.pure_vv.and.ip.eq.0.and.ih.eq.0.and.
     &          ivv+jvv.gt.0) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv.or.
     &          abs(ih-ip)+2*jvv.gt.maxv) cycle
            if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
            if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
            if (max(ip,ih)+ivv.lt.minexc.or.
     &          max(ip,ih)+jvv.lt.minexc) cycle
c dbg
            if (singrm.and.max(ih,ip).eq.1.and.min(ivv,jvv).eq.0) cycle
c dbgend
c            ! skip if block already exists
c            skip = .false.
c            do idef = 1, ndef
c             skip = skip.or.(occ_def(IVALE,1,idef*3).eq.max(ip-ih,0)+ivv
c     &               .and.occ_def(IVALE,1,idef*3-1).eq.max(ih-ip,0)+ivv)
c            end do
c            if (skip) cycle
            ndef = ndef + 1
            if (first) then
              nsupD = nsupD + 1
              stndD(1,nsupD) = ndef
              first = .false.
            end if
            stndD(2,nsupD) = ndef
           end do
          end do
        end do
      end do
      if (nsupD.gt.nsupT) call quit(1,'set_ic_mrcc_targets',
     &       'more super-blocks for metric than for T?')

      ! TT intermediate
      call add_target('TT',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
            if (ip.lt.1.or.ih.lt.1) cycle
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('TT',ttype_op,DEF_OP_FROM_OCC,
     &              'TT',1,1,
     &              parameters,2,tgt_info)

      ! transformed T
      call add_target2('Ttr',.false.,tgt_info)
      call set_dependency('Ttr','T',tgt_info)
      call set_rule2('Ttr',CLONE_OP,tgt_info)
      call set_arg('Ttr',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Ttr'/))
      call set_arg('Ttr',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'T'/))

      ! transformed L
      call add_target2('Ltr',.false.,tgt_info)
      call set_dependency('Ltr','L',tgt_info)
      call set_rule2('Ltr',CLONE_OP,tgt_info)
      call set_arg('Ltr',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Ltr'/))
      call set_arg('Ltr',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'L'/))

      ! subset of L for non-redundant valence-only metric
      call add_target('Lred',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (ip.eq.ih.and.ip.eq.maxexc) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
c dbg
c            if (ih.gt.ip) cycle
            if (singrm.and.ivv.eq.0.and.max(ih,ip).eq.1) cycle
c dbgend
c            ! exclude blocks with one or more hole-part. excitations
c            if (min(ih,ip).ge.1) cycle
            ! same valence structure already exists?
            ivers = 1
            do idef = 1, ndef
              if (occ_def(IVALE,1,idef).eq.max(ip-ih,0)+ivv
     &                 .and.occ_def(IVALE,2,idef).eq.max(ih-ip,0)+ivv)
     &           ivers = ivers + 1
            end do
c            if (skip) cycle
            ndef = ndef + 1
            occ_def(IHOLE,1,ndef) = ih
            occ_def(IPART,2,ndef) = ip
            occ_def(IVALE,2,ndef) = max(ih-ip,0) + ivv
            occ_def(IVALE,1,ndef) = max(ip-ih,0) + ivv
            ! distinguish ops with same valence part by blk_version
            version(ndef) = ivers
          end do
        end do
      end do
      n_tred_cls = ndef
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('Lred',ttype_op,DEF_OP_FROM_OCC,
     &              'Lred',1,1,
     &              parameters,2,tgt_info)
      call set_rule2('Lred',SET_ORDER,tgt_info)
      call set_arg('Lred',SET_ORDER,'LABEL',1,tgt_info,
     &     val_label=(/'Lred'/))
      call set_arg('Lred',SET_ORDER,'SPECIES',1,tgt_info,
     &             val_int=(/1/))
      call set_rule2('Lred',SET_ORDER,tgt_info)
      call set_arg('Lred',SET_ORDER,'LABEL',1,tgt_info,
     &     val_label=(/'Lred'/))
      call set_arg('Lred',SET_ORDER,'SPECIES',1,tgt_info,
     &             val_int=(/-1/))
      call set_arg('Lred',SET_ORDER,'ORDER',1,tgt_info,
     &             val_int=(/ndef/))
      call set_arg('Lred',SET_ORDER,'IDX_FREQ',ndef,tgt_info,
     &             val_int=version(1:ndef))

      ! subset of T operator
      call add_target2('Tred',.false.,tgt_info)
      call set_dependency('Tred','Lred',tgt_info)
      call set_rule2('Tred',CLONE_OP,tgt_info)
      call set_arg('Tred',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Tred'/))
      call set_arg('Tred',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'Lred'/))
      call set_arg('Tred',CLONE_OP,'ADJOINT',1,tgt_info,
     &     val_log=(/.true./))

      ! subset of Residual
      call add_target('OMGred',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (ip.eq.ih.and.ip.eq.maxexc) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
c dbg
c            if (ih.gt.ip) cycle
            if (singrm.and.ivv.eq.0.and.max(ih,ip).eq.1) cycle
c dbgend
            ! same valence structure already exists?
            ivers = 1
            do idef = 1, ndef
              if (occ_def(IVALE,2,idef*2-1).eq.max(ip-ih,0)+ivv
     &            .and.occ_def(IVALE,1,idef*2).eq.max(ih-ip,0)+ivv)
     &           ivers = ivers + 1
            end do
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef*2) = ih
            occ_def(IPART,1,ndef*2) = ip
            occ_def(IVALE,1,ndef*2) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef*2-1) = max(ip-ih,0) + ivv
            ! distinguish ops with same valence part by blk_version
            version(ndef) = ivers
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('OMGred',ttype_op,DEF_OP_FROM_OCC,
     &              'OMGred',1,1,
     &              parameters,2,tgt_info)
      call set_rule2('OMGred',SET_ORDER,tgt_info)
      call set_arg('OMGred',SET_ORDER,'LABEL',1,tgt_info,
     &     val_label=(/'OMGred'/))
      call set_arg('OMGred',SET_ORDER,'SPECIES',1,tgt_info,
     &             val_int=(/-1/))
      call set_arg('OMGred',SET_ORDER,'ORDER',1,tgt_info,
     &             val_int=(/ndef/))
      call set_arg('OMGred',SET_ORDER,'IDX_FREQ',ndef,tgt_info,
     e             val_int=version(1:ndef))

      ! define Residual
      call add_target('OMG',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
c dbg
c            if (ih.gt.ip) cycle
            if (singrm.and.ivv.eq.0.and.(ih+ip.eq.1.or.ih.eq.ip)) cycle
c dbgend
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef*2) = ih
            occ_def(IPART,1,ndef*2) = ip
            occ_def(IVALE,1,ndef*2) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef*2-1) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      if (singrm) then ! define conventional blocks separately
        do ip = minp, maxp ! ih=ip
          if (ip.lt.minh.or.ip.gt.maxh) cycle
          if (ip.lt.minexc) cycle
          ndef = ndef + 1
          occ_def(IHOLE,2,ndef*2) = ip
          occ_def(IPART,1,ndef*2) = ip
        end do
      end if
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('OMG',ttype_op,DEF_OP_FROM_OCC,
     &              'OMG',1,1,
     &              parameters,2,tgt_info)

      ! define transformed Residual (for preconditioner)
      call add_target('OMGtr',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
            if (ip.eq.ih.and.ip.eq.maxexc) cycle ! no conv. blocks (for PREC)
c dbg
c            if (ih.gt.ip) cycle
            if (singrm.and.ivv.eq.0.and.max(ih,ip).eq.1) cycle
c dbgend
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef*2) = ih
            occ_def(IPART,1,ndef*2) = ip
            occ_def(IVALE,1,ndef*2) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef*2-1) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('OMGtr',ttype_op,DEF_OP_FROM_OCC,
     &              'OMGtr',1,1,
     &              parameters,2,tgt_info)

      ! define Jacobian (diagonal blocks only)
      call add_target('A(CC)',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
            if (ip.eq.ih.and.ip.eq.maxexc) cycle ! no conv. blocks (for PREC)
c dbg
c            if (ih.gt.ip) cycle
            if (singrm.and.ivv.eq.0.and.max(ih,ip).eq.1) cycle
c dbgend
            ndef = ndef + 1
            occ_def(IHOLE,1,3*ndef-1) = ih
            occ_def(IHOLE,2,3*ndef-1) = ih
            occ_def(IPART,1,3*ndef-1) = ip
            occ_def(IPART,2,3*ndef-1) = ip
            occ_def(IVALE,1,3*ndef-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,3*ndef-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,3*ndef-2) = max(ip-ih,0) + ivv
            occ_def(IVALE,1,3*ndef) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
      call set_rule('A(CC)',ttype_op,DEF_OP_FROM_OCC,
     &              'A(CC)',1,1,
     &              parameters,2,tgt_info)

      ! Metric operator (for testing, also non-diagonal blocks)
      call add_target('S',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
           do ivv2 = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*max(ivv,ivv2).gt.maxv) cycle
            if (max(ip,ih).eq.0.and.
     &          (min(ivv,ivv2).eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+min(ivv,ivv2).lt.minexc) cycle
            if (abs(ih-ip).eq.0.and.min(ivv,ivv2).eq.0) cycle ! no conv. blocks
            ndef = ndef + 1
            occ_def(IHOLE,1,3*ndef-1) = ih
            occ_def(IHOLE,2,3*ndef-1) = ih
            occ_def(IPART,1,3*ndef-1) = ip
            occ_def(IPART,2,3*ndef-1) = ip
            occ_def(IVALE,1,3*ndef-1) = max(ih-ip,0) + ivv2
            occ_def(IVALE,2,3*ndef-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,3*ndef-2) = max(ip-ih,0) + ivv2
            occ_def(IVALE,1,3*ndef) = max(ip-ih,0) + ivv
           end do
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
      call set_rule('S',ttype_op,DEF_OP_FROM_OCC,
     &              'S',1,1,
     &              parameters,2,tgt_info)

      ! Diagonal Preconditioner
      call add_target(op_dia//'_T',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(op_dia//'_T','T',tgt_info)
      call cloneop_parameters(-1,parameters,'T',.false.)
      call set_rule(op_dia//'_T',ttype_op,CLONE_OP,op_dia//'_T',1,1,
     &              parameters,1,tgt_info)

      ! Diagonal Preconditioner for L
      call add_target(op_dia//'_L',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(op_dia//'_L','L',tgt_info)
      call cloneop_parameters(-1,parameters,'L',.false.)
      call set_rule(op_dia//'_L',ttype_op,CLONE_OP,op_dia//'_L',1,1,
     &              parameters,1,tgt_info)

      ! C0 dagger operator (for testing)
      call add_target2('C0dag',.false.,tgt_info)
      call set_dependency('C0dag','C0',tgt_info)
      call set_rule2('C0dag',CLONE_OP,tgt_info)
      call set_arg('C0dag',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'C0dag'/))
      call set_arg('C0dag',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'C0'/))
      call set_arg('C0dag',CLONE_OP,'ADJOINT',1,tgt_info,
     &     val_log=(/.true./))

      ! HT commutator: blocks up to rank 2
      op_ht = 'HT '
      do icnt = 1, 6
        write(op_ht(3:3),'(i1)') icnt
        call add_target2(op_ht,.false.,tgt_info)
        call set_dependency(op_ht,op_ham,tgt_info)
        call set_rule2(op_ht,CLONE_OP,tgt_info)
        call set_arg(op_ht,CLONE_OP,'LABEL',1,tgt_info,
     &       val_label=(/op_ht/))
        call set_arg(op_ht,CLONE_OP,'TEMPLATE',1,tgt_info,
     &       val_label=(/op_ham/))
      end do

      ! sums of HT commutators
      op_ht0to = 'HT0to '
      do icnt = 1, 4
        write(op_ht0to(6:6),'(i1)') icnt
        call add_target2(op_ht0to,.false.,tgt_info)
        call set_dependency(op_ht0to,op_ham,tgt_info)
        call set_rule2(op_ht0to,CLONE_OP,tgt_info)
        call set_arg(op_ht0to,CLONE_OP,'LABEL',1,tgt_info,
     &       val_label=(/op_ht0to/))
        call set_arg(op_ht0to,CLONE_OP,'TEMPLATE',1,tgt_info,
     &       val_label=(/op_ham/))
      end do

      ! define scalar unit operator
      call add_target('1scal',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 1
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('1scal',ttype_op,DEF_OP_FROM_OCC,
     &              '1scal',1,1,
     &              parameters,2,tgt_info)
      call opt_parameters(-1,parameters,+1,0)
      call set_rule('1scal',ttype_op,SET_HERMIT,
     &              '1scal',1,1,
     &              parameters,1,tgt_info)
c dbg
cmh    ****************************************************************
cmh    Want to activate calculation of Heff?
cmh    1.) comment in everything related to Heff and 1v in this routine
cmh    2.) Extend definition of 1v in set_ic_mrci_targets up to nactel
cmh    3.) new case in topo_remove_vtxs:
cmh        xlines_new = xlines_new + xlines_scr_u 
cmh    ****************************************************************
c      ! define effective hamiltonian
c      call add_target('Heff',ttype_op,.false.,tgt_info)
c      occ_def = 0
c      ndef = 1
c      occ_def(IVALE,1,ndef) = orb_info%nactel
c      occ_def(IVALE,2,ndef) = orb_info%nactel
c      call op_from_occ_parameters(-1,parameters,2,
c     &              occ_def,ndef,1,(/0,0/),ndef)
c      call set_rule('Heff',ttype_op,DEF_OP_FROM_OCC,
c     &              'Heff',1,1,
c     &              parameters,2,tgt_info)
c dbgend
*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      ! multireference CC lagrangian
      ! a) set up
      call add_target2('F_MRCC_LAG',.false.,tgt_info)
      call set_dependency('F_MRCC_LAG','NORM',tgt_info)
      call set_dependency('F_MRCC_LAG','H',tgt_info)
      call set_dependency('F_MRCC_LAG','C0',tgt_info)
      call set_dependency('F_MRCC_LAG','T',tgt_info)
      call set_dependency('F_MRCC_LAG','L',tgt_info)
      call set_rule2('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,tgt_info)
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &     tgt_info,val_label=(/'L','H','T','C0'/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &     tgt_info,val_int=(/maxcom/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &     tgt_info,val_int=(/maxcom_en/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &     val_str='---')
        call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'TITLE',1,
     &     tgt_info,val_str='ic-MRCC Lagrangian')
      if (optref.eq.-1.or.optref.eq.-2) then
        call set_dependency('F_MRCC_LAG','E(MR)',tgt_info)
        call set_rule2('F_MRCC_LAG',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'LABEL',1,
     &       tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'OP_RES',1,
     &       tgt_info,
     &       val_label=(/'NORM'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'E(MR)','C0^+','C0'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'IDX_SV',3,
     &       tgt_info,
     &       val_int=(/2,3,4/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &       val_rl8=(/-1d0/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
      end if
      if (.false..and.maxh.gt.0) then
        ! replace T-T double contractions by TT intermediate
        call set_dependency('F_MRCC_LAG','F_TT',tgt_info)
        call set_rule2('F_MRCC_LAG',FACTOR_OUT,tgt_info)
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_TT'/))
      end if
      ! factor out HT intermediates
c      call set_dependency('F_MRCC_LAG','F_HT1',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT2',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT0to1',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT0to2',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT3',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT4',tgt_info)
c      call set_rule2('F_MRCC_LAG',FACTOR_OUT,tgt_info)
c      call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_MRCC_LAG'/))
c      call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_MRCC_LAG'/))
c      call set_arg('F_MRCC_LAG',FACTOR_OUT,'INTERM',4,tgt_info,
c     &     val_label=(/'F_HT2','F_HT1','F_HT0to2','F_HT0to1'/))
      call set_rule2('F_MRCC_LAG',SELECT_SPECIAL,tgt_info)
      call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &     val_label=(/'H','T'/))
      call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC2')
      call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str='CHECK_FAC')
      call set_rule2('F_MRCC_LAG',PRINT_FORMULA,tgt_info)
      call set_arg('F_MRCC_LAG',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))

      ! Residual
      call add_target2('F_OMG',.false.,tgt_info)
      call set_dependency('F_OMG','F_MRCC_LAG',tgt_info)
      call set_dependency('F_OMG','OMG',tgt_info)
      call set_rule2('F_OMG',DERIVATIVE,tgt_info)
      call set_arg('F_OMG',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG'/))
      call set_arg('F_OMG',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_OMG',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      call set_arg('F_OMG',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'L'/))
      call set_rule2('F_OMG',SELECT_SPECIAL,tgt_info)
      call set_arg('F_OMG',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG'/))
      call set_arg('F_OMG',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_OMG'/))
      call set_arg('F_OMG',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &     val_label=(/'H','T'/))
      call set_arg('F_OMG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC2')
c dbg
c      call set_rule2('F_OMG',PRINT_FORMULA,tgt_info)
c      call set_arg('F_OMG',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_OMG'/))
c dbgend

      ! Lagrangian without Lambda...
      call add_target2('F_E_C0',.false.,tgt_info)
      call set_dependency('F_E_C0','F_MRCC_LAG',tgt_info)
      call set_dependency('F_E_C0','L',tgt_info)
      call set_rule2('F_E_C0',INVARIANT,tgt_info)
      call set_arg('F_E_C0',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E_C0'/))
      call set_arg('F_E_C0',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_E_C0',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_E_C0',INVARIANT,'OPERATORS',1,tgt_info,
     &     val_label=(/'L'/))
      call set_arg('F_E_C0',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='MRCC energy expression for C0 equations')
c dbg
c      ! insert unit operators to allow for double differentiation
c      call set_dependency('F_E_C0','1v',tgt_info)
c      call set_rule2('F_E_C0',INSERT,tgt_info)
c      call set_arg('F_E_C0',INSERT,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_E_C0'/))
c      call set_arg('F_E_C0',INSERT,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_E_C0'/))
c      call set_arg('F_E_C0',INSERT,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_E_C0',INSERT,'OP_INS',1,tgt_info,
c     &     val_label=(/'1v'/))
c      call set_arg('F_E_C0',INSERT,'OP_INCL',2,tgt_info,
c     &     val_label=(/'C0^+','C0'/))
c dbgend

      ! Residual for C0
      call add_target2('F_OMG_C0',.false.,tgt_info)
      call set_dependency('F_OMG_C0','F_E_C0',tgt_info)
      call set_dependency('F_OMG_C0','A_C0',tgt_info)
      call set_rule2('F_OMG_C0',DERIVATIVE,tgt_info)
      call set_arg('F_OMG_C0',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG_C0'/))
      call set_arg('F_OMG_C0',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E_C0'/))
      call set_arg('F_OMG_C0',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'A_C0'/))
      call set_arg('F_OMG_C0',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'C0^+'/))
c      call set_rule2('F_OMG_C0',INVARIANT,tgt_info)
c      call set_arg('F_OMG_C0',INVARIANT,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_OMG_C0'/))
c      call set_arg('F_OMG_C0',INVARIANT,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_OMG_C0'/))
c      call set_arg('F_OMG_C0',INVARIANT,'OP_RES',1,tgt_info,
c     &     val_label=(/'A_C0'/))
c      call set_arg('F_OMG_C0',INVARIANT,'OPERATORS',1,tgt_info,
c     &     val_label=(/'L'/))
c      call set_arg('F_OMG_C0',INVARIANT,'TITLE',1,tgt_info,
c     &     val_str='Residual for Reference function')
      call set_rule2('F_OMG_C0',SELECT_SPECIAL,tgt_info)
      call set_arg('F_OMG_C0',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG_C0'/))
      call set_arg('F_OMG_C0',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_OMG_C0'/))
      call set_arg('F_OMG_C0',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &     val_label=(/'H','T'/))
      call set_arg('F_OMG_C0',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC2')
c dbg
c      call set_rule2('F_OMG_C0',PRINT_FORMULA,tgt_info)
c      call set_arg('F_OMG_C0',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_OMG_C0'/))
c dbgend

      ! Energy
      call add_target2('F_MRCC_E',.false.,tgt_info)
      call set_dependency('F_MRCC_E','F_MRCC_LAG',tgt_info)
      call set_dependency('F_MRCC_E','E(MR)',tgt_info)
      call set_dependency('F_MRCC_E','L',tgt_info)
      call set_rule2('F_MRCC_E',INVARIANT,tgt_info)
      call set_arg('F_MRCC_E',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_E'/))
      call set_arg('F_MRCC_E',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_MRCC_E',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'E(MR)'/))
      call set_arg('F_MRCC_E',INVARIANT,'OPERATORS',2,tgt_info,
     &     val_label=(/'L','E(MR)'/))
      call set_arg('F_MRCC_E',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='MRCC energy expression')
      call set_rule2('F_MRCC_E',SELECT_SPECIAL,tgt_info)
      call set_arg('F_MRCC_E',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_E'/))
      call set_arg('F_MRCC_E',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_E'/))
      call set_arg('F_MRCC_E',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &     val_label=(/'H','T'/))
      call set_arg('F_MRCC_E',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC2')
      call set_rule2('F_MRCC_E',PRINT_FORMULA,tgt_info)
      call set_arg('F_MRCC_E',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_E'/))

      ! multireference CC norm
      ! a) set up using reduced densities
      call add_target2('F_MRCC_NORM',.false.,tgt_info)
      call set_dependency('F_MRCC_NORM','NORM',tgt_info)
      call set_dependency('F_MRCC_NORM','Tred',tgt_info)
      call set_dependency('F_MRCC_NORM','Lred',tgt_info)
      call set_dependency('F_MRCC_NORM','DENS',tgt_info)
      call set_rule2('F_MRCC_NORM',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'OPERATORS',2,
     &     tgt_info,
     &     val_label=(/'Lred','Tred'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'IDX_SV',2,tgt_info,
     &     val_int=(/2,3/))
      call set_rule2('F_MRCC_NORM',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'DENS','Lred','Tred','DENS'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/2,3,4,2/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &     val_int=(/1,4/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
     &     val_int=(/orb_info%nactel,-1,-1,orb_info%nactel/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! b) insert unit operators to allow for differentiation
      call set_dependency('F_MRCC_NORM','1v',tgt_info)
      call set_rule2('F_MRCC_NORM',INSERT,tgt_info)
      call set_arg('F_MRCC_NORM',INSERT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',INSERT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',INSERT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_MRCC_NORM',INSERT,'OP_INS',1,tgt_info,
     &     val_label=(/'1v'/))
      call set_arg('F_MRCC_NORM',INSERT,'OP_INCL',2,tgt_info,
     &     val_label=(/'Lred','Tred'/))
      call set_dependency('F_MRCC_NORM','1scal',tgt_info)
      call set_rule2('F_MRCC_NORM',INSERT,tgt_info)
      call set_arg('F_MRCC_NORM',INSERT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',INSERT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',INSERT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_MRCC_NORM',INSERT,'OP_INS',1,tgt_info,
     &     val_label=(/'1scal'/))
      call set_arg('F_MRCC_NORM',INSERT,'OP_INCL',2,tgt_info,
     &     val_label=(/'Lred','Tred'/))
      ! c) replace 1v by 1 (was used because we only needed valence blocks)
      call set_dependency('F_MRCC_NORM','1',tgt_info)
      call set_rule2('F_MRCC_NORM',REPLACE,tgt_info)
      call set_arg('F_MRCC_NORM',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'1v','1'/))
c dbg
c      call set_rule2('F_MRCC_NORM',PRINT_FORMULA,tgt_info)
c      call set_arg('F_MRCC_NORM',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_MRCC_NORM'/))
c dbgend

      ! transformed excitation operator
      call add_target2('F_T',.false.,tgt_info)
      call set_dependency('F_T','T',tgt_info)
      call set_dependency('F_T','Dtr',tgt_info)
      call set_dependency('F_T','Ttr',tgt_info)
      do i_cls = 1, nsupD
        call set_rule2('F_T',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_T',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'T','Dtr','Ttr','Dtr','T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/1,2,3,2,1/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MIN',5,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndD(1,i_cls),stndT(1,i_cls),
     &                 stndD(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MAX',5,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndD(2,i_cls),stndT(2,i_cls),
     &                 stndD(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/4/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'AVOID',8,tgt_info,
     &       val_int=(/3,5,2,4,1,4,2,5/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
      do i_cls = nsupD+1, nsupT
        call set_rule2('F_T',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_T',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'T','Ttr','T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &       val_int=(/1,2,1/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndT(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndT(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
      end do
      ! no active open lines from Ttr
      call set_rule2('F_T',SELECT_LINE,tgt_info)
      call set_arg('F_T',SELECT_LINE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_T'/))
      call set_arg('F_T',SELECT_LINE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_T'/))
      call set_arg('F_T',SELECT_LINE,'OP_RES',1,tgt_info,
     &     val_label=(/'T'/))
      call set_arg('F_T',SELECT_LINE,'OP_INCL',1,tgt_info,
     &     val_label=(/'Ttr'/))
      call set_arg('F_T',SELECT_LINE,'IGAST',1,tgt_info,
     &     val_int=(/3/))
      call set_arg('F_T',SELECT_LINE,'MODE',1,tgt_info,
     &     val_str='no_ext')
c dbg
c      call set_rule2('F_T',PRINT_FORMULA,tgt_info)
c      call set_arg('F_T',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_T'/))
c dbgend

      ! transformed deexcitation operator
      call add_target2('F_L',.false.,tgt_info)
      call set_dependency('F_L','L',tgt_info)
      call set_dependency('F_L','Dtr',tgt_info)
      call set_dependency('F_L','Ltr',tgt_info)
      do i_cls = 1, nsupD
        call set_rule2('F_L',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_L',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'L','Dtr^+','Ltr','Dtr^+','L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/1,2,3,2,1/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MIN',5,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndD(1,i_cls),stndT(1,i_cls),
     &                 stndD(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MAX',5,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndD(2,i_cls),stndT(2,i_cls),
     &                 stndD(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/4/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'AVOID',8,tgt_info,
     &       val_int=(/1,3,2,4,1,4,2,5/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
      do i_cls = nsupD+1, nsupT
        call set_rule2('F_L',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_L',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'L','Ltr','L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &       val_int=(/1,2,1/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndT(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndT(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
      end do
      ! no active open lines from Ttr
      call set_rule2('F_L',SELECT_LINE,tgt_info)
      call set_arg('F_L',SELECT_LINE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_L'/))
      call set_arg('F_L',SELECT_LINE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_L'/))
      call set_arg('F_L',SELECT_LINE,'OP_RES',1,tgt_info,
     &     val_label=(/'L'/))
      call set_arg('F_L',SELECT_LINE,'OP_INCL',1,tgt_info,
     &     val_label=(/'Ltr'/))
      call set_arg('F_L',SELECT_LINE,'IGAST',1,tgt_info,
     &     val_int=(/3/))
      call set_arg('F_L',SELECT_LINE,'MODE',1,tgt_info,
     &     val_str='no_ext')
c dbg
c      call set_rule2('F_L',PRINT_FORMULA,tgt_info)
c      call set_arg('F_L',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_L'/))
c dbgend

      ! transformed multireference CC energy
      ! a) set up
      call add_target2('F_E(MRCC)tr',.false.,tgt_info)
      call set_dependency('F_E(MRCC)tr','E(MR)',tgt_info)
      call set_dependency('F_E(MRCC)tr','H',tgt_info)
c      call set_dependency('F_E(MRCC)tr','FREF',tgt_info)
      call set_dependency('F_E(MRCC)tr','C0',tgt_info)
      call set_dependency('F_E(MRCC)tr','T',tgt_info)
      call set_dependency('F_E(MRCC)tr','L',tgt_info)
      call set_rule2('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,tgt_info)
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'LABEL',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'OP_RES',1,
     &     tgt_info,val_label=(/'E(MR)'/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &     tgt_info,val_label=(/'L','H','T','C0'/))
c     &     tgt_info,val_label=(/'L','FREF','T','C0'/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &     tgt_info,val_int=(/1/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &     tgt_info,val_int=(/0/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &     val_str='NOSCAL')
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'TITLE',1,tgt_info,
     &     val_str='Precursor for linearized Jacobian')
      ! f) insert 1 (particle/hole space) for later differentiation
      call set_dependency('F_E(MRCC)tr','1ph',tgt_info)
      call set_rule2('F_E(MRCC)tr',INSERT,tgt_info)
      call set_arg('F_E(MRCC)tr',INSERT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',INSERT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',INSERT,'OP_RES',1,tgt_info,
     &     val_label=(/'E(MR)'/))
      call set_arg('F_E(MRCC)tr',INSERT,'OP_INS',1,tgt_info,
     &     val_label=(/'1ph'/))
      call set_arg('F_E(MRCC)tr',INSERT,'OP_INCL',2,tgt_info,
     &     val_label=(/'L','T'/))
      ! replace 1ph by 1
      call set_dependency('F_E(MRCC)tr','1',tgt_info)
      call set_rule2('F_E(MRCC)tr',REPLACE,tgt_info)
      call set_arg('F_E(MRCC)tr',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'1ph','1'/))
      ! g) expand T and L
      call set_dependency('F_E(MRCC)tr','F_T',tgt_info)
      call set_dependency('F_E(MRCC)tr','F_L',tgt_info)
      call set_rule2('F_E(MRCC)tr',EXPAND,tgt_info)
      call set_arg('F_E(MRCC)tr',EXPAND,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',EXPAND,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',EXPAND,'INTERM',2,tgt_info,
     &     val_label=(/'F_T','F_L'/))
      ! h) let only diagonal blocks of Jacobian survive
      call set_rule2('F_E(MRCC)tr',SELECT_SPECIAL,tgt_info)
      call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &     val_label=(/'Ltr','Ttr'/))
      call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='SAME')
c dbg
      call set_rule2('F_E(MRCC)tr',PRINT_FORMULA,tgt_info)
      call set_arg('F_E(MRCC)tr',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
c dbgend

      ! (transformed) Jacobian times vector
      call add_target2('F_A_Ttr',.false.,tgt_info)
      call set_dependency('F_A_Ttr','F_E(MRCC)tr',tgt_info)
      call set_dependency('F_A_Ttr','OMGtr',tgt_info)
      call set_rule2('F_A_Ttr',DERIVATIVE,tgt_info)
      call set_arg('F_A_Ttr',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_A_Ttr'/))
      call set_arg('F_A_Ttr',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_A_Ttr',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMGtr'/))
      call set_arg('F_A_Ttr',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'Ltr'/))
c dbg
c      call set_rule2('F_A_Ttr',PRINT_FORMULA,tgt_info)
c      call set_arg('F_A_Ttr',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_A_Ttr'/))
c dbgend

      ! (transformed) Jacobian
      call add_target2('F_Atr',.false.,tgt_info)
      call set_dependency('F_Atr','F_A_Ttr',tgt_info)
      call set_dependency('F_Atr','A(CC)',tgt_info)
      call set_rule2('F_Atr',DERIVATIVE,tgt_info)
      call set_arg('F_Atr',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_Atr'/))
      call set_arg('F_Atr',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_A_Ttr'/))
      call set_arg('F_Atr',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'A(CC)'/))
      call set_arg('F_Atr',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'Ttr'/))
c dbg
c      call set_rule2('F_Atr',KEEP_TERMS,tgt_info)
c      call set_arg('F_Atr',KEEP_TERMS,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_Atr'/))
c      call set_arg('F_Atr',KEEP_TERMS,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_Atr'/))
c      call set_arg('F_Atr',KEEP_TERMS,'TERMS',2,tgt_info,
c     &     val_int=(/1,7/))
c
      call set_rule2('F_Atr',PRINT_FORMULA,tgt_info)
      call set_arg('F_Atr',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_Atr'/))
c dbgend

      ! Metric times amplitudes
      call add_target2('F_MRCC_SC',.false.,tgt_info)
      call set_dependency('F_MRCC_SC','F_MRCC_NORM',tgt_info)
      call set_dependency('F_MRCC_SC','OMGred',tgt_info)
      call set_rule2('F_MRCC_SC',DERIVATIVE,tgt_info)
      call set_arg('F_MRCC_SC',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_SC'/))
      call set_arg('F_MRCC_SC',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_SC',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMGred'/))
      call set_arg('F_MRCC_SC',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'Lred'/))
c dbg
c      call set_rule2('F_MRCC_SC',PRINT_FORMULA,tgt_info)
c      call set_arg('F_MRCC_SC',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_MRCC_SC'/))
c dbgend

      ! Metric (valence-only part)
      call add_target2('F_MRCC_D',.false.,tgt_info)
      call set_dependency('F_MRCC_D','F_MRCC_SC',tgt_info)
      call set_dependency('F_MRCC_D','D',tgt_info)
      call set_rule2('F_MRCC_D',DERIVATIVE,tgt_info)
      call set_arg('F_MRCC_D',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_D'/))
      call set_arg('F_MRCC_D',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_SC'/))
      call set_arg('F_MRCC_D',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'D'/))
      call set_arg('F_MRCC_D',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'Tred'/))
c dbg
      call set_rule2('F_MRCC_D',PRINT_FORMULA,tgt_info)
      call set_arg('F_MRCC_D',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_D'/))
c dbgend

      ! HT intermediate
      f_ht = 'F_HT '
      op_ht = 'HT '
      do icnt = 1, 2
        write(f_ht(5:5),'(i1)') icnt
        write(op_ht(3:3),'(i1)') icnt
        call add_target2(f_ht,.false.,tgt_info)
        call set_dependency(f_ht,op_ht,tgt_info)
        call set_dependency(f_ht,op_ham,tgt_info)
        call set_dependency(f_ht,'T',tgt_info)
        call set_rule2(f_ht,DEF_MRCC_LAGRANGIAN,tgt_info)
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'LABEL',1,tgt_info,
     &       val_label=(/f_ht/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'OP_RES',1,tgt_info,
     &       val_label=(/op_ht/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'OPERATORS',4,tgt_info,
     &       val_label=(/'T  ','H  ','T  ',op_ht/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,tgt_info,
     &       val_int=(/icnt/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,tgt_info,
     &       val_int=(/icnt/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &       val_str='HBAR')
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'TITLE',1,tgt_info,
     &       val_str='Commutator intermediate')
        if (icnt.gt.1.and.maxh.gt.0) then
          ! replace T-T double contractions by TT intermediate
          call set_dependency(f_ht,'F_TT',tgt_info)
          call set_rule2(f_ht,FACTOR_OUT,tgt_info)
          call set_arg(f_ht,FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &         val_label=(/f_ht/))
          call set_arg(f_ht,FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &         val_label=(/f_ht/))
          call set_arg(f_ht,FACTOR_OUT,'INTERM',1,tgt_info,
     &         val_label=(/'F_TT'/))
        end if
c dbg
c        if (icnt.eq.2) then
c        ! factor out lower HT intermediates
c        call set_dependency(f_ht,'F_HT1',tgt_info)
c        call set_rule2(f_ht,FACTOR_OUT,tgt_info)
c        call set_arg(f_ht,FACTOR_OUT,'LABEL_RES',1,tgt_info,
c     &       val_label=(/f_ht/))
c        call set_arg(f_ht,FACTOR_OUT,'LABEL_IN',1,tgt_info,
c     &       val_label=(/f_ht/))
c        call set_arg(f_ht,FACTOR_OUT,'INTERM',1,tgt_info,
c     &       val_label=(/'F_HT1'/))
c        end if
c dbgend
        call set_rule2(f_ht,PRINT_FORMULA,tgt_info)
        call set_arg(f_ht,PRINT_FORMULA,'LABEL',1,tgt_info,
     &       val_label=(/f_ht/))
      end do

      ! sums of HT intermediates
      f_ht0to = 'F_HT0to '
      op_ht0to = 'HT0to '
      form_str = 'HT0to =H'
      len_form = 8
      do icnt = 1, 4
        write(f_ht0to(8:8),'(i1)') icnt
        write(op_ht0to(6:6),'(i1)') icnt
        call add_target2(f_ht0to,.false.,tgt_info)
        call set_dependency(f_ht0to,op_ht0to,tgt_info)
        call set_dependency(f_ht0to,op_ham,tgt_info)
        op_ht = 'HT '
        do ih = 1, icnt
          write(op_ht(3:3),'(i1)') ih
          call set_dependency(f_ht0to,op_ht,tgt_info)
        end do
        write(form_str(6:6),'(i1)') icnt
        write(form_str(len_form+1:len_form+4),'(a1,a3)') '+',op_ht
        len_form = len_form + 4
        call set_rule2(f_ht0to,DEF_FORMULA,tgt_info)
        call set_arg(f_ht0to,DEF_FORMULA,'LABEL',1,tgt_info,
     &       val_label=(/f_ht0to/))
        call set_arg(f_ht0to,DEF_FORMULA,'FORMULA',1,tgt_info,
     &       val_str=form_str(1:len_form))
        call set_rule2(f_ht0to,PRINT_FORMULA,tgt_info)
        call set_arg(f_ht0to,PRINT_FORMULA,'LABEL',1,tgt_info,
     &       val_label=(/f_ht0to/))
      end do

      call add_target2('F_TT',.false.,tgt_info)
      call set_dependency('F_TT','TT',tgt_info)
      call set_dependency('F_TT','T',tgt_info)
      call set_rule2('F_TT',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_TT',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_TT'/))
      call set_arg('F_TT',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'TT'/))
      call set_arg('F_TT',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'TT','T','T','TT'/))
      call set_arg('F_TT',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,3,1/))
      call set_rule2('F_TT',PRINT_FORMULA,tgt_info)
      call set_arg('F_TT',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_TT'/))

c dbg
c      ! transformed multireference CC energy
c      ! a) set up
c      call add_target2('F_NORMtr',.false.,tgt_info)
c      call set_dependency('F_NORMtr','NORM',tgt_info)
c      call set_dependency('F_NORMtr','T',tgt_info)
c      call set_dependency('F_NORMtr','L',tgt_info)
c      call set_dependency('F_NORMtr','DENS',tgt_info)
c      call set_rule2('F_NORMtr',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'OPERATORS',2,
c     &     tgt_info,
c     &     val_label=(/'L','T'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'IDX_SV',2,tgt_info,
c     &     val_int=(/2,3/))
c      call set_rule2('F_NORMtr',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'OPERATORS',4,
c     &     tgt_info,
c     &     val_label=(/'DENS','L','T','DENS'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
c     &     val_int=(/2,3,4,2/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
c     &     val_int=(/1,4/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
c     &     val_int=(/orb_info%nactel,-1,-1,orb_info%nactel/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      ! f) insert 1 (particle/hole space) for later differentiation
c      call set_dependency('F_NORMtr','1ph',tgt_info)
c      call set_rule2('F_NORMtr',INSERT,tgt_info)
c      call set_arg('F_NORMtr',INSERT,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',INSERT,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',INSERT,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_NORMtr',INSERT,'OP_INS',1,tgt_info,
c     &     val_label=(/'1ph'/))
c      call set_arg('F_NORMtr',INSERT,'OP_INCL',2,tgt_info,
c     &     val_label=(/'L','T'/))
c      ! replace 1ph by 1
c      call set_dependency('F_NORMtr','1',tgt_info)
c      call set_rule2('F_NORMtr',REPLACE,tgt_info)
c      call set_arg('F_NORMtr',REPLACE,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',REPLACE,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',REPLACE,'OP_LIST',2,tgt_info,
c     &     val_label=(/'1ph','1'/))
c      ! g) expand T and L
c      call set_dependency('F_NORMtr','F_T',tgt_info)
c      call set_dependency('F_NORMtr','F_L',tgt_info)
c      call set_rule2('F_NORMtr',EXPAND,tgt_info)
c      call set_arg('F_NORMtr',EXPAND,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',EXPAND,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',EXPAND,'INTERM',2,tgt_info,
c     &     val_label=(/'F_T','F_L'/))
cc dbg
c      call set_rule2('F_NORMtr',PRINT_FORMULA,tgt_info)
c      call set_arg('F_NORMtr',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
cc dbgend
c
c      ! transformed Metric times vector
c      call add_target2('F_S_Ttr',.false.,tgt_info)
c      call set_dependency('F_S_Ttr','F_NORMtr',tgt_info)
c      call set_dependency('F_S_Ttr','OMG',tgt_info)
c      call set_rule2('F_S_Ttr',DERIVATIVE,tgt_info)
c      call set_arg('F_S_Ttr',DERIVATIVE,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_S_Ttr'/))
c      call set_arg('F_S_Ttr',DERIVATIVE,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_S_Ttr',DERIVATIVE,'OP_RES',1,tgt_info,
c     &     val_label=(/'OMG'/))
c      call set_arg('F_S_Ttr',DERIVATIVE,'OP_DERIV',1,tgt_info,
c     &     val_label=(/'Ltr'/))
c      call set_rule2('F_S_Ttr',PRINT_FORMULA,tgt_info)
c      call set_arg('F_S_Ttr',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_S_Ttr'/))
c
c      ! transformed Metric (should be unity)
c      call add_target2('F_Str',.false.,tgt_info)
c      call set_dependency('F_Str','F_S_Ttr',tgt_info)
c      call set_dependency('F_Str','A(CC)',tgt_info)
c      call set_rule2('F_Str',DERIVATIVE,tgt_info)
c      call set_arg('F_Str',DERIVATIVE,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_Str'/))
c      call set_arg('F_Str',DERIVATIVE,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_S_Ttr'/))
c      call set_arg('F_Str',DERIVATIVE,'OP_RES',1,tgt_info,
c     &     val_label=(/'A(CC)'/))
c      call set_arg('F_Str',DERIVATIVE,'OP_DERIV',1,tgt_info,
c     &     val_label=(/'Ttr'/))
c      call set_rule2('F_Str',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Str',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Str'/))
c dbgend

c dbg
c      ! effective Hamiltonian
c      call add_target2('F_Heff',.false.,tgt_info)
c      call set_dependency('F_Heff','F_OMG_C0',tgt_info)
c      call set_dependency('F_Heff','Heff',tgt_info)
c      call set_rule2('F_Heff',DERIVATIVE,tgt_info)
c      call set_arg('F_Heff',DERIVATIVE,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_Heff'/))
c      call set_arg('F_Heff',DERIVATIVE,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_OMG_C0'/))
c      call set_arg('F_Heff',DERIVATIVE,'OP_RES',1,tgt_info,
c     &     val_label=(/'Heff'/))
c      call set_arg('F_Heff',DERIVATIVE,'OP_DERIV',1,tgt_info,
c     &     val_label=(/'C0'/))
c      call set_rule2('F_Heff',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Heff',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Heff'/))
c dbgend
*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! density matrix
      call add_target2('FOPT_MRCC_D',.false.,tgt_info)
      call set_dependency('FOPT_MRCC_D','F_DENS0',tgt_info)
      call set_dependency('FOPT_MRCC_D','F_MRCC_D',tgt_info)
      call set_dependency('FOPT_MRCC_D','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_MRCC_D','DEF_ME_1scal',tgt_info)
      call set_dependency('FOPT_MRCC_D','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_MRCC_D','DEF_ME_D',tgt_info)
      call set_dependency('FOPT_MRCC_D','DEF_ME_DENS',tgt_info)
      call set_rule2('FOPT_MRCC_D',OPTIMIZE,tgt_info)
      call set_arg('FOPT_MRCC_D',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_MRCC_D'/))
      call set_arg('FOPT_MRCC_D',OPTIMIZE,'LABELS_IN',2,tgt_info,
     &             val_label=(/'F_DENS0','F_MRCC_D'/))

      ! transformation
      call add_target2('FOPT_T',.false.,tgt_info)
      call set_dependency('FOPT_T','F_T',tgt_info)
      call set_dependency('FOPT_T','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_T','DEF_ME_Ttr',tgt_info)
      call set_dependency('FOPT_T','DEF_ME_Dtr',tgt_info)
      call set_rule2('FOPT_T',OPTIMIZE,tgt_info)
      call set_arg('FOPT_T',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_T'/))
      call set_arg('FOPT_T',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_T'/))

      ! transformed Hessian
      call add_target2('FOPT_Atr',.false.,tgt_info)
      call set_dependency('FOPT_Atr','F_Atr',tgt_info)
      call set_dependency('FOPT_Atr',mel_ham,tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_A(CC)',tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_Dtr',tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_C0',tgt_info)
      call set_rule2('FOPT_Atr',ASSIGN_ME2OP,tgt_info)
      call set_arg('FOPT_Atr',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_Dtr'/))
      call set_arg('FOPT_Atr',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'Dtr'/))
      call set_rule2('FOPT_Atr',OPTIMIZE,tgt_info)
      call set_arg('FOPT_Atr',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_Atr'/))
      call set_arg('FOPT_Atr',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_Atr'/))

      ! transformed Hessian times vector
      call add_target2('FOPT_A_Ttr',.false.,tgt_info)
      call set_dependency('FOPT_A_Ttr','F_A_Ttr',tgt_info)
      call set_dependency('FOPT_A_Ttr',mel_ham,tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_OMGtr',tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_Ttr',tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_Dtr',tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_C0',tgt_info)
      call set_rule2('FOPT_A_Ttr',ASSIGN_ME2OP,tgt_info)
      call set_arg('FOPT_A_Ttr',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_Dtr'/))
      call set_arg('FOPT_A_Ttr',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'Dtr'/))
      call set_rule2('FOPT_A_Ttr',OPTIMIZE,tgt_info)
      call set_arg('FOPT_A_Ttr',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_A_Ttr'/))
      call set_arg('FOPT_A_Ttr',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_A_Ttr'/))

      ! Residual
      call add_target2('FOPT_OMG',.false.,tgt_info)
      call set_dependency('FOPT_OMG','F_OMG',tgt_info)
      call set_dependency('FOPT_OMG','F_MRCC_E',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_T',tgt_info)
      if (.false..and.maxh.gt.0)
     &    call set_dependency('FOPT_OMG','DEF_ME_TT',tgt_info)
c      call set_dependency('FOPT_OMG','DEF_ME_HT1',tgt_info)
c      call set_dependency('FOPT_OMG','DEF_ME_HT2',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_OMG',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_E(MR)',tgt_info)
      call set_dependency('FOPT_OMG',mel_ham,tgt_info)
      if (optref.eq.-1.or.optref.eq.-2) then
        call set_dependency('FOPT_OMG','F_OMG_C0',tgt_info)
        call set_dependency('FOPT_OMG','DEF_ME_A_C0',tgt_info)
c dbg
c      call set_dependency('FOPT_OMG','DEF_ME_1v',tgt_info)
c dbgend
        call set_rule2('FOPT_OMG',ASSIGN_ME2OP,tgt_info)
        call set_arg('FOPT_OMG',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_A_C0'/))
        call set_arg('FOPT_OMG',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &             val_label=(/'A_C0'/))
      end if
      call set_rule2('FOPT_OMG',OPTIMIZE,tgt_info)
      call set_arg('FOPT_OMG',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_OMG'/))
      if (optref.eq.-1.or.optref.eq.-2) then
        call set_arg('FOPT_OMG',OPTIMIZE,'LABELS_IN',3,tgt_info,
     &             val_label=(/'F_MRCC_E','F_OMG','F_OMG_C0'/))
      else
        call set_arg('FOPT_OMG',OPTIMIZE,'LABELS_IN',2,tgt_info,
     &             val_label=(/'F_MRCC_E','F_OMG'/))
      end if

      ! Residual for C0
      call add_target2('FOPT_OMG_C0',.false.,tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_T',tgt_info)
c dbg
c      call set_dependency('FOPT_OMG_C0','DEF_ME_1v',tgt_info)
c dbgend
      if (.false..and.maxh.gt.0)
     &    call set_dependency('FOPT_OMG_C0','DEF_ME_TT',tgt_info)
c      call set_dependency('FOPT_OMG_C0','DEF_ME_HT1',tgt_info)
c      call set_dependency('FOPT_OMG_C0','DEF_ME_HT2',tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_E(MR)',tgt_info)
      call set_dependency('FOPT_OMG_C0',mel_ham,tgt_info)
      call set_dependency('FOPT_OMG_C0','F_OMG_C0',tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_A_C0',tgt_info)
      call set_rule2('FOPT_OMG_C0',ASSIGN_ME2OP,tgt_info)
      call set_arg('FOPT_OMG_C0',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_A_C0'/))
      call set_arg('FOPT_OMG_C0',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'A_C0'/))
      call set_rule2('FOPT_OMG_C0',OPTIMIZE,tgt_info)
      call set_arg('FOPT_OMG_C0',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_OMG_C0'/))
      call set_arg('FOPT_OMG_C0',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_OMG_C0'/))

c dbg
c      ! transformed metric
c      call add_target2('FOPT_Str',.false.,tgt_info)
c      call set_dependency('FOPT_Str','F_Str',tgt_info)
c      call set_dependency('FOPT_Str','DEF_ME_A(CC)',tgt_info)
c      call set_dependency('FOPT_Str','DEF_ME_1',tgt_info)
c      call set_dependency('FOPT_Str','DEF_ME_Dtr',tgt_info)
c      call set_rule2('FOPT_Str',OPTIMIZE,tgt_info)
c      call set_arg('FOPT_Str',OPTIMIZE,'LABEL_OPT',1,tgt_info,
c     &             val_label=(/'FOPT_Str'/))
c      call set_arg('FOPT_Str',OPTIMIZE,'LABELS_IN',1,tgt_info,
c     &             val_label=(/'F_Str'/))
c dbgend

c dbg
c      ! effective Hamiltonian
c      call add_target2('FOPT_Heff',.false.,tgt_info)
c      call set_dependency('FOPT_Heff','F_Heff',tgt_info)
c      call set_dependency('FOPT_Heff','DEF_ME_1v',tgt_info)
c      call set_dependency('FOPT_Heff','DEF_ME_Heff',tgt_info)
c      call set_dependency('FOPT_Heff',mel_ham,tgt_info)
c      call set_dependency('FOPT_Heff','DEF_ME_T',tgt_info)
c      call set_rule2('FOPT_Heff',OPTIMIZE,tgt_info)
c      call set_arg('FOPT_Heff',OPTIMIZE,'LABEL_OPT',1,tgt_info,
c     &             val_label=(/'FOPT_Heff'/))
c      call set_arg('FOPT_Heff',OPTIMIZE,'LABELS_IN',1,tgt_info,
c     &             val_label=(/'F_Heff'/))
c dbgend
*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! ME for Hessian
      call add_target2('DEF_ME_A(CC)',.false.,tgt_info)
      call set_dependency('DEF_ME_A(CC)','A(CC)',tgt_info)
      call set_rule2('DEF_ME_A(CC)',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_A(CC)'/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'A(CC)'/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'DIAG_IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'DIAG_MS',1,tgt_info,
     &     val_int=(/0/))

      ! Diagonal Preconditioner
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call add_target(trim(dia_label),ttype_opme,.false.,tgt_info)
      call set_dependency(trim(dia_label),'EVAL_FREF',tgt_info)
      call set_dependency(trim(dia_label),
     &                    op_dia//'_'//'T',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = trim(dia_label)
      labels(2) = op_dia//'_'//'T'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule(trim(dia_label),ttype_opme,
     &              DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
      ! use effective Fock op. (needed only for pure inactive exc.)
      labels(1) = trim(dia_label)
      labels(2) = 'ME_FREF'
      call set_rule(trim(dia_label),ttype_opme,
     &              PRECONDITIONER,
     &              labels,2,1,
     &              parameters,1,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Preconditioner (a):',0,'LIST')
c      call set_rule(trim(dia_label),ttype_opme,PRINT_MEL,
c     &     trim(dia_label),1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! ME for T
      call add_target2('DEF_ME_T',.false.,tgt_info)
      call set_dependency('DEF_ME_T','T',tgt_info)
      call set_rule2('DEF_ME_T',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_T',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_T'/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'T'/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for Ttr
      call add_target2('DEF_ME_Ttr',.false.,tgt_info)
      call set_dependency('DEF_ME_Ttr','Ttr',tgt_info)
      call set_rule2('DEF_ME_Ttr',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Ttr'/))
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'Ttr'/))
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))

      ! ME-List for HT intermediates
      op_ht = 'HT '
      op_ht0to = 'ME_HT '
      def_ht = 'DEF_ME_HT '
      do icnt = 1, 4
        write(op_ht(3:3),'(i1)') icnt
        write(op_ht0to(6:6),'(i1)') icnt
        write(def_ht(10:10),'(i1)') icnt
        call add_target2(def_ht,.false.,tgt_info)
        call set_dependency(def_ht,op_ht,tgt_info)
        call set_rule2(def_ht,DEF_ME_LIST,tgt_info)
        call set_arg(def_ht,DEF_ME_LIST,'LIST',1,tgt_info,
     &               val_label=(/op_ht0to/))
        call set_arg(def_ht,DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &               val_label=(/op_ht/))
        call set_arg(def_ht,DEF_ME_LIST,'MS',1,tgt_info,
     &               val_int=(/0/))
        call set_arg(def_ht,DEF_ME_LIST,'IRREP',1,tgt_info,
     &               val_int=(/1/))
      end do

      ! ME for TT intermediate
      call add_target2('DEF_ME_TT',.false.,tgt_info)
      call set_dependency('DEF_ME_TT','TT',tgt_info)
      call set_rule2('DEF_ME_TT',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_TT',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_TT'/))
      call set_arg('DEF_ME_TT',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'TT'/))
      call set_arg('DEF_ME_TT',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_TT',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))

      ! ME for Residual
      call add_target2('DEF_ME_OMG',.false.,tgt_info)
      call set_dependency('DEF_ME_OMG','OMG',tgt_info)
      call set_rule2('DEF_ME_OMG',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_OMG'/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'OMG'/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for transformed Residual
      call add_target2('DEF_ME_OMGtr',.false.,tgt_info)
      call set_dependency('DEF_ME_OMGtr','OMGtr',tgt_info)
      call set_rule2('DEF_ME_OMGtr',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_OMGtr',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_OMGtr'/))
      call set_arg('DEF_ME_OMGtr',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'OMGtr'/))
      call set_arg('DEF_ME_OMGtr',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_OMGtr',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))

      ! reordered projector matrix ME_Dproj (to eliminate lin. dep.)
      call add_target('DEF_ME_Dproj',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_Dproj','Dtr',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_Dproj'
      labels(2) = 'Dtr'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_Dproj',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call set_dependency('DEF_ME_Dproj','DEF_ME_Dinv',tgt_info)
      labels(1) = 'ME_Dproj'
      labels(2) = 'ME_D'
      call form_parameters(-1,parameters,2,
     &     '---',13,'---')
      call set_rule('DEF_ME_Dproj',ttype_opme,
     &              REORDER_MEL,
     &              labels,2,1,
     &              parameters,2,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Reordered projector matrix :',0,'LIST')
c      call set_rule('DEF_ME_Dproj',ttype_opme,PRINT_MEL,
c     &     'ME_Dproj',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! ME_1scal
      call add_target2('DEF_ME_1scal',.false.,tgt_info)
      call set_dependency('DEF_ME_1scal','1scal',tgt_info)
      call set_rule2('DEF_ME_1scal',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_1scal'/))
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'1scal'/))
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &     val_int=(/1/))
      call dens_parameters(-1,parameters,0,0,0)
      call set_rule('DEF_ME_1scal',ttype_opme,UNITY,
     &     'ME_1scal',1,1,
     &     parameters,1,tgt_info)

c dbg
c      ! ME_Heff
c      call add_target2('DEF_ME_Heff',.false.,tgt_info)
c      call set_dependency('DEF_ME_Heff','Heff',tgt_info)
c      call set_rule2('DEF_ME_Heff',DEF_ME_LIST,tgt_info)
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'LIST',1,tgt_info,
c     &     val_label=(/'ME_Heff'/))
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'OPERATOR',1,tgt_info,
c     &     val_label=(/'Heff'/))
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'IRREP',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'MS',1,tgt_info,
c     &     val_int=(/0/))
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'DIAG_IRREP',1,tgt_info,
c     &     val_int=(/orb_info%lsym/))
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'DIAG_MS',1,tgt_info,
c     &     val_int=(/mod(orb_info%imult-1,2)/))
c
c      ! ME_1v
c      call add_target2('DEF_ME_1v',.false.,tgt_info)
c      call set_dependency('DEF_ME_1v','1v',tgt_info)
c      call set_rule2('DEF_ME_1v',DEF_ME_LIST,tgt_info)
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'LIST',1,tgt_info,
c     &     val_label=(/'ME_1v'/))
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'OPERATOR',1,tgt_info,
c     &     val_label=(/'1v'/))
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'IRREP',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'MS',1,tgt_info,
c     &     val_int=(/0/))
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
c     &     val_int=(/1/))
c      call dens_parameters(-1,parameters,0,0,0)
c      call set_rule('DEF_ME_1v',ttype_opme,UNITY,
c     &     'ME_1v',1,1,
c     &     parameters,1,tgt_info)
c dbgend
*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! evaluate valence-only metric
      call add_target2('EVAL_MRCC_D',.false.,tgt_info)
      call set_dependency('EVAL_MRCC_D','FOPT_MRCC_D',tgt_info)
      call set_dependency('EVAL_MRCC_D','EVAL_REF_S(S+1)',tgt_info)
c      call set_dependency('EVAL_MRCC_D','EVAL_DENS0',tgt_info)
      call set_rule2('EVAL_MRCC_D',EVAL,tgt_info)
      call set_arg('EVAL_MRCC_D',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_MRCC_D'/))
c dbg
c      call set_rule2('EVAL_MRCC_D',PRINT_MEL,tgt_info)
c      call set_arg('EVAL_MRCC_D',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_DENS'/))
c      call set_rule2('EVAL_MRCC_D',PRINT_MEL,tgt_info)
c      call set_arg('EVAL_MRCC_D',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_D'/))
c dbgend

      ! Evaluate diagonal elements of Jacobian
      call add_target('EVAL_Atr',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_Atr','FOPT_Atr',tgt_info)
c      call set_dependency('EVAL_Atr','EVAL_FREF',tgt_info)
      call set_rule('EVAL_Atr',ttype_opme,EVAL,
     &     'FOPT_Atr',1,0,
     &     parameters,0,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'transformed Jacobian :',0,'LIST')
c      call set_rule('EVAL_Atr',ttype_opme,PRINT_MEL,
c     &     'ME_A(CC)',1,0,
c     &     parameters,2,tgt_info)
c dbgend
      ! put diagonal elements to preconditioner
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      labels(1) = trim(dia_label)
      labels(2) = 'ME_A(CC)'
      call set_dependency('EVAL_Atr',trim(dia_label),tgt_info)
      call set_rule('EVAL_Atr',ttype_opme,
     &              EXTRACT_DIAG,
     &              labels,2,1,
     &              parameters,0,tgt_info)
c dbg
      call form_parameters(-1,parameters,2,
     &     'Preconditioner (b) :',0,'LIST')
      call set_rule('EVAL_Atr',ttype_opme,PRINT_MEL,
     &     trim(dia_label),1,0,
     &     parameters,2,tgt_info)
c dbgend

      ! Evaluate approximation of diagonal elements of Jacobian
      call add_target('EVAL_A_Ttr',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_A_Ttr','FOPT_A_Ttr',tgt_info)
      ! (a) set all elements of right hand vector to one
      call set_rule2('EVAL_A_Ttr',SET_MEL,tgt_info)
      call set_arg('EVAL_A_Ttr',SET_MEL,'LIST',1,tgt_info,
     &             val_label=(/'ME_Ttr'/))
      call set_arg('EVAL_A_Ttr',SET_MEL,'VAL_LIST',2,tgt_info,
     &             val_rl8=(/1d0,1d0/))
      call set_arg('EVAL_A_Ttr',SET_MEL,'IDX_LIST',2,tgt_info,
     &             val_int=(/-1,-2/))
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Ttr set to :',0,'LIST')
c      call set_rule('EVAL_A_Ttr',ttype_opme,PRINT_MEL,
c     &     'ME_Ttr',1,0,
c     &     parameters,2,tgt_info)
c dbgend
      ! (b) evaluate
      call set_rule('EVAL_A_Ttr',ttype_opme,EVAL,
     &     'FOPT_A_Ttr',1,0,
     &     parameters,0,tgt_info)
      ! (c) scale back and copy onto preconditioner list
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call set_dependency('EVAL_A_Ttr',trim(dia_label),tgt_info)
      call set_rule2('EVAL_A_Ttr',SCALE_COPY,tgt_info)
      call set_arg('EVAL_A_Ttr',SCALE_COPY,'LIST_RES',1,tgt_info,
     &             val_label=(/trim(dia_label)/))
      call set_arg('EVAL_A_Ttr',SCALE_COPY,'LIST_INP',1,tgt_info,
     &             val_label=(/'ME_OMGtr'/))
      call set_arg('EVAL_A_Ttr',SCALE_COPY,'FAC',2,tgt_info,
     &             val_rl8=(/1d0,1d0/))
      ! (d) free up Ttr list
      call set_rule('EVAL_A_Ttr',ttype_opme,RES_ME_LIST,
     &     'ME_Ttr',1,0,
     &     parameters,0,tgt_info)
c dbg
      call form_parameters(-1,parameters,2,
     &     'Preconditioner (b) :',0,'LIST')
      call set_rule('EVAL_A_Ttr',ttype_opme,PRINT_MEL,
     &     trim(dia_label),1,0,
     &     parameters,2,tgt_info)
c dbgend

      ! Solve MR coupled cluster equations
      call add_target2('SOLVE_MRCC',.true.,tgt_info)
      call set_dependency('SOLVE_MRCC','EVAL_REF_S(S+1)',tgt_info)
      call set_dependency('SOLVE_MRCC','FOPT_OMG',tgt_info)
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call set_dependency('SOLVE_MRCC',trim(dia_label),tgt_info)
      if (cheap_prc) then
        call set_dependency('SOLVE_MRCC','EVAL_A_Ttr',tgt_info)
      else
        call set_dependency('SOLVE_MRCC','EVAL_Atr',tgt_info)
      end if
      call set_dependency('SOLVE_MRCC','EVAL_MRCC_D',tgt_info)
      call set_dependency('SOLVE_MRCC','DEF_ME_Dtrdag',tgt_info)
      call set_dependency('SOLVE_MRCC','FOPT_T',tgt_info)
      if (optref.ne.0) then
        call me_list_label(dia_label2,mel_dia,orb_info%lsym,
     &                     0,0,0,.false.)
        dia_label2 = trim(dia_label2)//'C0'
        call set_dependency('SOLVE_MRCC',trim(dia_label2),tgt_info)
        if (optref.ne.-1.and.optref.ne.-2) 
     &     call set_dependency('SOLVE_MRCC','FOPT_OMG_C0',tgt_info)
        call set_dependency('SOLVE_MRCC','DEF_ME_Dproj',tgt_info)
      end if
      do icnt = 1, max(1,optref)
      call set_rule2('SOLVE_MRCC',SOLVENLEQ,tgt_info)
      if (optref.lt.0) then
        if (preopt) then
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',1,tgt_info,
     &         val_label=(/'ME_T'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
     &         val_str='TRF')
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',1,tgt_info,
     &         val_label=(/'ME_OMG'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',1,tgt_info,
     &         val_label=(/trim(dia_label)/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_E',1,tgt_info,
     &       val_label=(/'ME_E(MR)'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',3,tgt_info,
     &       val_label=(/'ME_Ttr','ME_Dtr','ME_Dtrdag'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',1,tgt_info,
     &       val_label=(/'FOPT_T'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM',1,tgt_info,
     &         val_label=(/'FOPT_OMG'/))
          call set_rule2('SOLVE_MRCC',SOLVENLEQ,tgt_info)
        end if
      end if
      if (optref.eq.-1.or.optref.eq.-2) then
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',2,tgt_info,
     &       val_label=(/'ME_T','ME_C0'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
     &       val_str='TRF/NRM')
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',2,tgt_info,
     &       val_label=(/'ME_OMG','ME_A_C0'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',2,tgt_info,
     &       val_label=(/trim(dia_label),trim(dia_label2)/))
      else
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',1,tgt_info,
     &       val_label=(/'ME_T'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
     &       val_str='TRF')
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',1,tgt_info,
     &       val_label=(/'ME_OMG'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',1,tgt_info,
     &       val_label=(/trim(dia_label)/))
      end if
      call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_E',1,tgt_info,
     &     val_label=(/'ME_E(MR)'/))
      if (optref.ne.0.and.update_prc) then
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',7,tgt_info,
     &     val_label=(/'ME_Ttr','ME_Dtr','ME_Dtrdag','ME_Dproj',
     &                 'ME_D','ME_Dinv',
     &                 'ME_A(CC)'/))
      else if (optref.ne.0.and..not.update_prc) then
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',6,tgt_info,
     &     val_label=(/'ME_Ttr','ME_Dtr','ME_Dtrdag','ME_Dproj',
     &                 'ME_D','ME_Dinv'/))
      else
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',3,tgt_info,
     &     val_label=(/'ME_Ttr','ME_Dtr','ME_Dtrdag'/))
      end if
      if (optref.ne.0) then
        if (update_prc) then
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',3,tgt_info,
     &         val_label=(/'FOPT_T','FOPT_MRCC_D','FOPT_Atr'/))
        else
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',2,tgt_info,
     &         val_label=(/'FOPT_T','FOPT_MRCC_D'/))
        end if
      else
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',1,tgt_info,
     &       val_label=(/'FOPT_T'/))
      end if
      call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM',1,tgt_info,
     &     val_label=(/'FOPT_OMG'/))
c dbg
c      ! print effective Hamiltonian
c      ! Evaluate transformed metric
c      call set_dependency('SOLVE_MRCC','FOPT_Heff',tgt_info)
c      call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
c     &     'FOPT_Heff',1,0,
c     &     parameters,0,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'effective Hamiltonian :',0,'LIST')
c      call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
c     &     'ME_Heff',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      if (optref.gt.0.and.icnt.ne.optref) then !not in last iteration
        call set_rule2('SOLVE_MRCC',SOLVEEVP,tgt_info)
        call set_arg('SOLVE_MRCC',SOLVEEVP,'LIST_OPT',1,tgt_info,
     &       val_label=(/'ME_C0'/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'MODE',1,tgt_info,
     &       val_str='DIA')
        call set_arg('SOLVE_MRCC',SOLVEEVP,'N_ROOTS',1,tgt_info,
     &       val_int=(/ciroot/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'OP_MVP',1,tgt_info,
     &       val_label=(/'A_C0'/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'LIST_PRC',1,tgt_info,
     &       val_label=(/trim(dia_label2)/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'OP_SVP',1,tgt_info,
     &     val_label=(/'C0'/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'FORM',1,tgt_info,
     &       val_label=(/'FOPT_OMG_C0'/))
c dbg
        call form_parameters(-1,parameters,2,
     &       'CI coefficients :',0,'LIST')
        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
     &       'ME_C0',1,0,
     &       parameters,2,tgt_info)
        call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
     &       'FOPT_REF_S(S+1)',1,0,
     &       parameters,0,tgt_info)
c dbgend
      end if
      end do
c dbg
      if (optref.ne.0) then
        call form_parameters(-1,parameters,2,
     &       'final CI coefficients :',0,'LIST')
        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
     &       'ME_C0',1,0,
     &       parameters,2,tgt_info)
        call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
     &       'FOPT_REF_S(S+1)',1,0,
     &       parameters,0,tgt_info)
      end if
c dbgend

c dbg
c        call form_parameters(-1,parameters,2,
c     &       'final T amplitudes :',0,'LIST')
c        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
c     &       'ME_T',1,0,
c     &       parameters,2,tgt_info)
c dbgend

c dbg
c      ! Evaluate transformed metric
c      call add_target('EVAL_Str',ttype_gen,.false.,tgt_info)
c      call set_dependency('EVAL_Str','FOPT_Str',tgt_info)
c      call set_dependency('EVAL_Str','SOLVE_REF',tgt_info)
c      call set_rule('EVAL_Str',ttype_opme,RES_ME_LIST,
c     &     'ME_A(CC)',1,0,
c     &     parameters,0,tgt_info)
c      call set_rule('EVAL_Str',ttype_opme,EVAL,
c     &     'FOPT_Str',1,0,
c     &     parameters,0,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'transformed metric :',0,'LIST')
c      call set_rule('EVAL_Str',ttype_opme,PRINT_MEL,
c     &     'ME_A(CC)',1,0,
c     &     parameters,2,tgt_info)
c dbgend

c dbg
c      call add_target2('PRINT_Heff',.false.,tgt_info)
c      call set_dependency('PRINT_Heff','FOPT_Heff',tgt_info)
c      call set_dependency('PRINT_Heff','SOLVE_REF',tgt_info)
c      call set_rule2('PRINT_Heff',PRINT_MEL,tgt_info)
c      call set_arg('PRINT_Heff',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_T'/))
c      call set_rule('PRINT_Heff',ttype_opme,EVAL,
c     &     'FOPT_Heff',1,0,
c     &     parameters,0,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'effective Hamiltonian :',0,'LIST')
c      call set_rule('PRINT_Heff',ttype_opme,PRINT_MEL,
c     &     'ME_Heff',1,0,
c     &     parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'bare Hamiltonian :',0,'LIST')
c      call set_rule('PRINT_Heff',ttype_opme,PRINT_MEL,
c     &     mel_ham,1,0,
c     &     parameters,2,tgt_info)
c dbgend

c dbg
c      ! Evaluate Residual
c      call add_target('EVAL_OMG',ttype_gen,.false.,tgt_info)
c      call set_dependency('EVAL_OMG','FOPT_OMG',tgt_info)
c      call set_dependency('EVAL_OMG','SOLVE_REF',tgt_info)
cc      call set_dependency('EVAL_OMG','PRINT_Heff',tgt_info)
c      call set_rule('EVAL_OMG',ttype_opme,RES_ME_LIST,
c     &     'ME_OMG',1,0,
c     &     parameters,0,tgt_info)
c      call set_rule('EVAL_OMG',ttype_opme,EVAL,
c     &     'FOPT_OMG',1,0,
c     &     parameters,0,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'Residual :',0,'LIST')
c      call set_rule('EVAL_OMG',ttype_opme,PRINT_MEL,
c     &     'ME_OMG',1,0,
c     &     parameters,2,tgt_info)
c dbgend

c dbg
c      ! Evaluate projected T
c      call add_target('EVAL_Tproj',ttype_gen,.false.,tgt_info)
c      call set_dependency('EVAL_Tproj','FOPT_T',tgt_info)
c      call set_dependency('EVAL_Tproj','SOLVE_REF',tgt_info)
c      call set_dependency('EVAL_Tproj','DEF_ME_Dproj',tgt_info)
c      call set_rule2('EVAL_Tproj',ASSIGN_ME2OP,tgt_info)
c      call set_arg('EVAL_Tproj',ASSIGN_ME2OP,'LIST',1,tgt_info,
c     &           val_label=(/'ME_T'/))
c      call set_arg('EVAL_Tproj',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
c     &           val_label=(/'Ttr'/))
c      call set_rule2('EVAL_Tproj',ASSIGN_ME2OP,tgt_info)
c      call set_arg('EVAL_Tproj',ASSIGN_ME2OP,'LIST',1,tgt_info,
c     &           val_label=(/'ME_Ttr'/))
c      call set_arg('EVAL_Tproj',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
c     &           val_label=(/'T'/))
c      call set_rule2('EVAL_Tproj',ASSIGN_ME2OP,tgt_info)
c      call set_arg('EVAL_Tproj',ASSIGN_ME2OP,'LIST',1,tgt_info,
c     &           val_label=(/'ME_Dproj'/))
c      call set_arg('EVAL_Tproj',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
c     &           val_label=(/'Dtr'/))
c      call set_rule('EVAL_Tproj',ttype_opme,EVAL,
c     &     'FOPT_T',1,0,
c     &     parameters,0,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'projected T :',0,'LIST')
c      call set_rule('EVAL_Tproj',ttype_opme,PRINT_MEL,
c     &     'ME_Ttr',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      return
      end
