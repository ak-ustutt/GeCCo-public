*----------------------------------------------------------------------*
      subroutine frm_sched1(xret,flist,depend_info,idxselect,nselect,
     &         op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     schedule the evaluation of a formula
*     
*     version 1 -- complete replacement of version 0
*----------------------------------------------------------------------*

      implicit none

      include 'routes.h'

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'par_opnames_gen.h'
      include 'def_reorder_info.h'
      include 'def_dependency_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      character, parameter ::
     &     name_scr0*6 = 'cntscr'

      character ::
     &     name_scr*13

      real(8), intent(inout) ::
     &     xret(*)
      integer, intent(in) ::
     &     nselect, idxselect(nselect)
      type(formula_item), intent(in), target ::
     &     flist
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info
      type(dependency_info) ::
     &     depend_info

      logical ::
     &     update, reo_op1op2, reo_other, possible, skip, new,
     &     tra_op1, tra_op2, tra_op1op2, set_reo, make_contr_red, self
      integer ::
     &     idxopres, idxres, nres, type_xret, type_xret_cur,
     &     idxme_res, idxmel, absym12,
     &     n_occ_cls, maxvtx, maxarc, maxfac, nblk_res,
     &     nfact, idxop1op2, iblkop1op2, iops, iblkres, ifree,
     &     ninter, idx, nsym, ngas, nexc, ndis, iprint, iterm, len,
     &     idoffop1, idoffop2, idoffop1op2, ivtx_new, nvtx,
     &     iarc, idum, idxop_intm
      real(8) ::
     &     fac, facc, xnrm, bc_sign, xret_last
      character ::
     &     title*256, opscrnam*8

      type(reorder_info) ::
     &     reo_info
      type(me_list_array), pointer ::
     &     mel_arr(:)
      type(operator_array), pointer ::
     &     ops(:)
      integer ::
     &     mstop(2), igamtop(2), idxop(2), iblkop(2),
     &     mstop1op2, igamtop1op2, njoined,
     &     njoined_op(2), njoined_op1op2,
     &     njoined_cnt, njoined_res, idxinp
      real(8), pointer ::
     &     xret_blk(:), xret_pnt(:)
      real(8), target ::
     &     xret_scr(1)
      integer, pointer ::
     &     op2list(:),
     &     occ_vtx(:,:,:), irestr_vtx(:,:,:,:,:), info_vtx(:,:),
     &     occ_vtx_red(:,:,:), irestr_vtx_red(:,:,:,:,:),
     &                                        info_vtx_red(:,:),
     &     merge_op1(:), merge_op2(:), merge_op1op2(:), merge_op2op1(:),
     &     iocc_op1(:,:,:), iocc_op2(:,:,:),
     &     irst_op1(:,:,:,:,:),
     &     irst_op2(:,:,:,:,:),
     &     irst_res(:,:,:,:,:),
     &     iocc_ex1(:,:,:), iocc_ex2(:,:,:), iocc_cnt(:,:,:),
     &     iocc_op1op2(:,:,:),
     &     irst_op1op2(:,:,:,:,:),
     &     iocc_op1op2tmp(:,:,:),
     &     irst_op1op2tmp(:,:,:,:,:)
      type(operator), pointer ::
     &     opscr(:), optmp
      type(me_list), pointer ::
     &     melscr(:), meltmp
      type(operator), pointer ::
     &     opres
      type(me_list), pointer ::
     &     me_op1, me_op2, me_op1op2, me_res, me_op1op2tmp
      type(filinf), pointer ::
     &     ffop1, ffop2, ffop1op2, ffres
      type(formula_item), pointer ::
     &     cur_form
      type(contraction) ::
     &     cur_contr, cur_contr_red

      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      integer, external ::
     &     idxlist, idx_oplist2
      logical, external ::
     &     me_list_uptodate
      real(8), external ::
     &     xnormop

      ifree = mem_setmark('frm_sched0')
      
      iprint = max(ntest,iprlvl)
      if (ntest.gt.0.or.iprint.ge.5) then
        write(luout,*) '============================='
        write(luout,*) '= Entered formula scheduler ='
        write(luout,*) '============================='
      end if

      nsym = orb_info%nsym
      ngas = orb_info%ngas

      ops => op_info%op_arr
      mel_arr => op_info%mel_arr
      op2list => op_info%op2list

      cur_form => flist

      skip = .false.
      nres  = 0
      idxres = 0
      iterm = 0
      nullify(xret_blk)
      xret_last = 0d0

      call init_contr(cur_contr)
      call init_contr(cur_contr_red)

      ! loop over entries
      term_loop: do 

        if (nres.gt.0) cur_form => cur_form%next

        ! skip to next entry, if requested
        if (skip) then
          do while(cur_form%command.eq.command_add_contribution)
            cur_form => cur_form%next
            if (.not.associated(cur_form))
     &           call quit(1,'frm_sched1',
     &           'formula list is not terminated properly!')
          end do
        end if
      
        if (.not.associated(cur_form))
     &       call quit(1,'frm_sched1',
     &       'formula list is not terminated properly!')

        if (nres.eq.0 .and.
     &      cur_form%command.ne.command_set_target_init)
     &     call quit(1,'frm_sched1','first command must define target')

        select case(cur_form%command)
        case(command_end_of_formula,command_set_target_init)

          ! for previous target: assemble xret value
          if (idxres.gt.0.and..not.skip) then

            ! in case of splin-flip symmetry exploitation:
            ! symmetrize final result here (if not scalar)
            if (use_tr.and.me_res%absym.ne.0.and.type_xret.eq.1) then
              call sym_ab_list(0.5d0,me_res,me_res,
     &             xret_blk,.true.,
     &             op_info,str_info,strmap_info,orb_info)
            end if

            if (type_xret.eq.1) then
              xret(idxres) = sqrt(sum(xret_blk(1:nblk_res)))
            else if (type_xret.eq.2) then
              xret(idxres) = xret_blk(1)
            end if

            call atim_csw(cpu,sys,wall)
            call prtim(luout,'time for target',
     &           cpu-cpu0,sys-sys0,wall-wall0)

          end if

          call atim_csw(cpu0,sys0,wall0)

          if (cur_form%command.eq.command_end_of_formula)
     &         exit term_loop

          xret_last = 0d0

          ! initialize result
          nres = nres+1
          idxres = nres

          ! requested?
          skip = nselect.gt.0.and.
     &           idxlist(idxres,idxselect,nselect,1).le.0
          ! check dependency
          skip = skip.or.me_list_uptodate(idxres,depend_info,op_info)

          idxopres = cur_form%target      ! op index of result
          if (idxopres.eq.0)
     &         call quit(1,'frm_sched1','idxopres==0 is obsolete!')

          idxme_res = op2list(idxopres)  ! list index of result
          if (iprint.ge.10.and.skip) then
            write(luout,*) 'Skipping target: ',
     &           trim(mel_arr(idxme_res)%mel%label)
          else if (iprint.ge.10) then
            write(luout,*) 'New target: ',
     &           trim(mel_arr(idxme_res)%mel%label)
          end if

          if (skip) cycle term_loop

c          ffres => op_info%opfil_arr(idxopres)%fhand
          me_res => op_info%mel_arr(idxme_res)%mel
          opres  => me_res%op
          ffres  => me_res%fhand
          nblk_res = opres%n_occ_cls
          type_xret = 2
          if (me_res%len_op.gt.1) type_xret = 1
          if (ffres%unit.le.0)
     &           call file_open(ffres)
          call zeroop(me_res)
 
          if (associated(xret_blk))
     &           ifree = mem_dealloc('xret_blk')
          ifree = mem_alloc_real(xret_blk,nblk_res,'xret_blk')
          xret_blk(1:nblk_res) = 0d0
 
          ! mark current term as updated:
          call touch_file_rec(ffres)

          cycle term_loop
c        case(command_set_target_update)
        case(command_symmetrise)
          
          call symmetrise(1d0,me_res,me_res,
     &         xret_blk,
     &         op_info,str_info,orb_info)

          cycle term_loop

        case(command_add_contribution)
        case default
          write(luout,*) 'command = ',cur_form%command
          call quit(1,'frm_sched','command not defined/implemented')
        end select

        ! get *copy* of contr as we are going to modify that
        call copy_contr(cur_form%contr,cur_contr)

        iterm = iterm+1
        if (iprint.ge.20)
     &     write(luout,*) '   term #',iterm

        if (ntest.ge.50)
     &       call prt_contr2(luout,cur_contr,op_info)
         
        ! process info
        fac = cur_contr%fac
        nfact = cur_contr%nfac
        
        if (idxopres.ne.cur_contr%idx_res)
     &       call quit(1,'frm_sched1','inconsistency in result index!')
        njoined_res = opres%njoined

        if (ntest.ge.100) write(luout,*) 'nfact, fac: ',nfact,fac

        idxmel = op2list(cur_contr%idx_res)
        me_op1op2 => mel_arr(idxmel)%mel
        iblkres = (cur_contr%iblk_res-1)/njoined_res + 1
        if (me_op1op2%len_op_occ(iblkres).eq.0)
     &       cycle term_loop
        ! check here for other zero blocks as well ...

        if (nfact.eq.0) then

c??          iblkres = (cur_contr%iblk_res-1)/njoined_res + 1
          iblkres = cur_contr%iblk_res
          idxop(1) = cur_contr%vertex(1)%idx_op
          iblkop(1) = cur_contr%vertex(1)%iblk_op
          tra_op1 = cur_contr%vertex(1)%dagger
          tra_op1op2 = cur_contr%dagger

          ! special: unit operator
          if (ops(idxop(1))%op%name.eq.op_unity) then
            call add_unity(fac,1,me_res,iblkres,orb_info,str_info)
          else
            idxmel = op2list(idxop(1))
            if (mel_arr(idxmel)%mel%fhand%unit.le.0)
     &           call file_open(mel_arr(idxmel)%mel%fhand)
c fix:
            njoined = mel_arr(idxmel)%mel%op%njoined
            iblkop(1) = (iblkop(1)-1)/njoined + 1
c fix:
            if (     tra_op1.and..not.tra_op1op2.or.
     &          .not.tra_op1.and.     tra_op1op2) then
              call add_opblk_transp(xret_blk(iblkres),type_xret,fac,
     &             mel_arr(idxmel)%mel,me_res,tra_op1,tra_op1op2,
     &             iblkop(1),iblkres,
     &             op_info,str_info,orb_info,.false.)
            else
              call add_opblk(xret_blk(iblkres),type_xret,fac,
     &             mel_arr(idxmel)%mel,me_res,
     &             iblkop(1),iblkres,orb_info,.false.)
            end if

            if (type_xret.eq.2.and.iprint.ge.3)
     &           write(luout,'(1x,"term # ",i5,":",2(x,g19.10))')
     &           iterm, xret_blk(1)-xret_last, xret_blk(1)
            if (type_xret.eq.2) xret_last = xret_blk(1)

            if (ntest.ge.50) then
              write(luout,*) 'xret after term ',iterm
              write(luout,'(x,4g19.10)') xret_blk(1:nblk_res)
            end if

          end if

          cycle term_loop
        end if

        nvtx = cur_form%contr%nvtx

        ! allocate arrays for occupations and restrictions
        allocate(
     &       iocc_op1(ngastp,2,nvtx), iocc_op2(ngastp,2,nvtx),
     &       irst_op1(2,ngas,2,2,nvtx),
     &       irst_op2(2,ngas,2,2,nvtx),
     &       irst_res(2,ngas,2,2,nvtx),
     &       iocc_ex1(ngastp,2,nvtx),
     &       iocc_ex2(ngastp,2,nvtx), iocc_cnt(ngastp,2,2*nvtx),
     &       iocc_op1op2(ngastp,2,nvtx),
     &       irst_op1op2(2,orb_info%ngas,2,2,nvtx),
     &       iocc_op1op2tmp(ngastp,2,nvtx),
     &       irst_op1op2tmp(2,orb_info%ngas,2,2,nvtx))        

        ! preliminary fix to set irst_res (needed in get_bc_info)
        call set_restr_prel(irst_res,
     &       cur_contr,op_info,orb_info%ihpvgas,orb_info%ngas)
        
        ! allocate arrays for intermediates
        allocate(
     &       occ_vtx(ngastp,2,nvtx+njoined_res),
     &       irestr_vtx(2,orb_info%ngas,2,2,nvtx+njoined_res),
     &       info_vtx(2,nvtx+njoined_res),
     &       merge_op1(10*nvtx*nvtx), ! a bit too large, I guess ...
     &       merge_op2(10*nvtx*nvtx),
     &       merge_op1op2(10*nvtx*nvtx),
     &       merge_op2op1(10*nvtx*nvtx))
        allocate(
     &       occ_vtx_red(ngastp,2,nvtx+njoined_res),
     &       irestr_vtx_red(2,orb_info%ngas,2,2,nvtx+njoined_res),
     &       info_vtx_red(2,nvtx+njoined_res))
        if (nfact.gt.1) then
          allocate(opscr(nfact-1),optmp,melscr(nfact-1),meltmp)
          do idx = 1, nfact-1
            melscr(idx)%fhand => null()
          end do
        end if

        ! reset intermediate counter
        ninter = 0

        call occvtx4contr(0,occ_vtx,cur_contr,op_info)
        call vtxinf4contr(irestr_vtx,info_vtx,
     &                            cur_contr,op_info,ngas)

        fac = cur_contr%fac

        ! add 0-contractions, if necessary
        call check_disconnected(cur_contr)

        ! loop over binary contractions
        bin_loop: do idx = 1, nfact

          ! use factor in last contraction
          facc = 1d0
          if (idx.eq.nfact) facc=fac

          if (iprint.ge.20) write(luout,*) '    contr #',idx
          iarc = cur_contr%inffac(5,idx)
          ninter = ninter + 1

          ! set up info for binary contraction
c          new = .false.!cur_contr%nvtx.ge.4
          make_contr_red = cur_contr%nsupvtx.gt.2 .or.
     &                    (cur_contr%nsupvtx.eq.2 .and.
     &                     cur_contr%narc.gt.1)
          set_reo = make_contr_red
          if (set_reo) then
            ! reset reo_info
            call init_reo_info(reo_info)
          end if

          call get_bc_info3(bc_sign,possible,
     &         idxop,iblkop,
     &         iocc_ex1,iocc_ex2,iocc_cnt,
     &         iocc_op1,iocc_op2,iocc_op1op2,
     &         irst_op1,irst_op2,irst_op1op2,
     &         tra_op1,tra_op2,tra_op1op2,
     &         mstop,mstop1op2,
     &         igamtop,igamtop1op2,
     &         njoined_op, njoined_op1op2, njoined_cnt,
     &         merge_op1,merge_op2,merge_op1op2, merge_op2op1,
     &         cur_contr,occ_vtx,irestr_vtx,info_vtx,
     &         make_contr_red,
     &         cur_contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &         set_reo,reo_info, reo_info,
     &         iarc,.false.,-ninter,
     &         irst_res,njoined_res,orb_info,op_info)
          if (.not.possible) then
            call prt_contr3(luout,cur_contr,-1)
            write(luout,*) 'get_bc_info did not raise "possible"-flag!'
            call quit(1,'frm_sched1','could not continue ...')
          end if

          ! set up reduced contraction after 
          ! current binary contraction
          if (idx.ne.nfact) then
            ! process reordering info
            call get_reo_info(reo_op1op2,reo_other,
     &           iocc_op1op2,iocc_op1op2tmp,
     &           irst_op1op2,irst_op1op2tmp,
     &           njoined_op1op2,
     &           reo_info,str_info,orb_info)
            if (reo_other) then
              call quit(1,'frm_sched1',
     &             'reordering of operator other than OP1OP2 requested')
            end if

            call copy_contr(cur_contr_red,cur_contr)
            occ_vtx    = occ_vtx_red
            irestr_vtx = irestr_vtx_red
            info_vtx   = info_vtx_red
            ! add 0-contractions, if necessary
            call check_disconnected(cur_contr)

          else
            reo_op1op2 = .false.
            reo_other = .false.
            iocc_op1op2tmp = iocc_op1op2
            irst_op1op2tmp = irst_op1op2
          end if

          ! set up operator 1 and 2
          self = .false.
          do iops = 1, 2
            if (idxop(iops).gt.0) then
              idxmel = op2list(idxop(iops))
              ! primary operator or long-term intermediate
              if (iops.eq.1) me_op1 => mel_arr(idxmel)%mel
              if (iops.eq.2) me_op2 => mel_arr(idxmel)%mel
            else if (idxop(iops).lt.0) then
              ! intermediate for current contraction only
              if (iops.eq.1) me_op1 => melscr(-idxop(iops))
              if (iops.eq.2) me_op2 => melscr(-idxop(iops))
            else if (iops.eq.2) then
              self = .true.
            else
              call quit(1,'frm_sched1','inconsistent idxop occurred!')
            end if
          end do

          if (me_op1%fhand%unit.le.0)
     &             call file_open(me_op1%fhand)
          if (.not.self.and.me_op2%fhand%unit.le.0)
     &             call file_open(me_op2%fhand)

          ! set up result
          if (idx.eq.nfact) then
            ! last operation: store on result array
            update = .true.
            idxop1op2 = cur_contr%idx_res
            iblkop1op2 = cur_contr%iblk_res
            idxmel = op2list(idxop1op2)
            me_op1op2 => mel_arr(idxmel)%mel
            xret_pnt => xret_blk(iblkop1op2:iblkop1op2)
            type_xret_cur = type_xret
          else
            ! new intermediate
            update = .false.
            iblkop1op2 = 1

            ! set up pseudo-operator for current intermediate
            ! refers to reordered operator (if this matters)
            write(opscrnam,'("INT",i3.3)') ninter
            call set_ps_op(opscr(ninter),opscrnam,
     &           iocc_op1op2,irst_op1op2,njoined_op1op2,orb_info)
            melscr(ninter)%op => opscr(ninter)

            ! Fix Ms?
            melscr(ninter)%fix_vertex_ms = .false.

            absym12 = me_op1%absym*me_op2%absym 
            call set_ps_list(melscr(ninter),opscrnam,
     &           absym12,0,mstop1op2,igamtop1op2,0,
     &           str_info,strmap_info,orb_info)
            call file_open(melscr(ninter)%fhand)

            me_op1op2 => melscr(ninter)
            xret_pnt => xret_scr
            type_xret_cur = 0

          end if

          ! if the result is reordered as a final step of the
          ! contraction, we need a pseudo-operator for the
          ! raw contraction result:
          if (reo_op1op2) then
            write(opscrnam,'("INT",i3.3,"RW")') ninter
            call set_ps_op(optmp,opscrnam,
     &           iocc_op1op2tmp,irst_op1op2tmp,njoined_op1op2,
     &           orb_info)
            meltmp%op => optmp

            ! Should Ms be fixed?
            meltmp%fix_vertex_ms = me_op1op2%fix_vertex_ms

            call set_ps_list(meltmp,opscrnam,
     &           0,0,mstop1op2,igamtop1op2,0,
     &           str_info,strmap_info,orb_info)
            me_op1op2tmp => meltmp
c            me_op1op2tmp%fix_vertex_ms = me_op1op2%fix_vertex_ms
          else
            me_op1op2tmp => me_op1op2
            me_op1op2tmp%fix_vertex_ms = me_op1op2%fix_vertex_ms
          end if

          ! translate records to offset in file:
          ! (makes life easier in case we once decide to use
          ! one scratch file only: no changes to contr_op1op2 necessary)
          ffop1 => me_op1%fhand
          if (.not.self) ffop2 => me_op2%fhand
          ffop1op2 => me_op1op2%fhand
          idoffop1 = ffop1%length_of_record*(ffop1%current_record-1)
          if (.not.self)
     &        idoffop2 = ffop2%length_of_record*(ffop2%current_record-1)
          idoffop1op2 = ffop1op2%length_of_record*
     &                                   (ffop1op2%current_record-1)

          if (ntest.ge.100)
     &         write(luout,*) 'calling contraction kernel'
          ! do the contraction
          call contr_op1op2(facc,bc_sign,
     &       update,self,xret_pnt,type_xret_cur,
     &       me_op1,me_op2,me_op1op2, me_op1op2tmp,
     &       tra_op1, tra_op2, tra_op1op2,
     &       iblkop(1),iblkop(2),iblkop1op2,iblkop1op2,
     &       idoffop1,idoffop2,idoffop1op2,
     &       iocc_ex1,iocc_ex2,iocc_cnt,
     &       idum, idum, idum, 0,
     &       iocc_op1, iocc_op2, iocc_op1op2, iocc_op1op2tmp,
     &       irst_op1,irst_op2,irst_op1op2, irst_op1op2tmp,
     &       merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &       njoined_op(1), njoined_op(2),njoined_op1op2, njoined_cnt,
     &       mstop(1),mstop(2),mstop1op2,
     &       igamtop(1),igamtop(2),igamtop1op2,
     &       reo_info,
     &       str_info,strmap_info,orb_info)
          if (ntest.ge.100)
     &         write(luout,*) 'returned from contraction kernel'

          if (reo_op1op2) then
            call dealloc_me_list(meltmp)
            call dealloc_operator(optmp)
            call dealloc_reo_info(reo_info)
          end if

        end do bin_loop

        if (type_xret.eq.2.and.iprint.ge.3)
     &       write(luout,'(1x,"term # ",i5,":",2(x,g19.10))')
     &       iterm, xret_blk(1)-xret_last, xret_blk(1)
        if (type_xret.eq.2) xret_last = xret_blk(1)

        if (ntest.ge.50) then
          write(luout,*) 'xret after term ',iterm
          write(luout,'(x,4f19.10)') xret_blk(1:nblk_res)
        end if

        deallocate(
     &       occ_vtx,irestr_vtx,info_vtx,
     &       occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &       merge_op1,merge_op2,merge_op1op2,merge_op2op1)
        deallocate(
     &       iocc_op1, iocc_op2,
     &       irst_op1, irst_op2,
     &       irst_res,
     &       iocc_ex1,iocc_ex2, iocc_cnt,
     &       iocc_op1op2,irst_op1op2)        

        if (nfact.gt.1) then
          ! get rid of intermediate definitions and files
          do idx = 1, ninter-1
c            call file_close_delete(ffscr(idx))
            call dealloc_me_list(melscr(idx))
            call dealloc_operator(opscr(idx))
          end do
          deallocate(opscr,optmp,melscr,meltmp)
        end if

      end do term_loop

      call dealloc_contr(cur_contr)
      call dealloc_contr(cur_contr_red)

      ifree = mem_flushmark()

      if (ntest.ge.100)
     &     write(luout,*) 'returning from frm_sched1'

      return
      end
