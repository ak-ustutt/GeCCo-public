*----------------------------------------------------------------------*
      subroutine frm_sched1(xret,flist,depend_info,idxselect,nselect,
     &         op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     schedule the evaluation of a formula
*     
*     version 1 -- complete replacement of version 0
*----------------------------------------------------------------------*

      implicit none

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
     &     update, reo_op1op2, reo_other, possible, skip,
     &     tra_op1, tra_op2, tra_op1op2
      integer ::
     &     lufrm, idxopres, idxres, nres, type_xret, type_xret_cur,
     &     idxme_res, idxmel,
     &     n_occ_cls, maxvtx, maxarc, maxfac, nblk_res,
     &     nfact, idxop1op2, iblkop1op2, iops, iblkres, ifree,
     &     ninter, idx, nsym, ngas, nexc, ndis, iprint, iterm, len,
     &     idoffop1, idoffop2, idoffop1op2, ivtx_new, nvtx,
     &     iarc, idum, idxop_intm
      real(8) ::
     &     fac, facc, xnrm, bc_sign
      character ::
     &     title*256, opscrnam*8

      type(reorder_info) ::
     &     reo_info
      type(me_list_array), pointer ::
     &     mel_arr(:)
c      type(file_array), pointer ::
c     &     ffops(:)
      type(operator_array), pointer ::
     &     ops(:)
      integer ::
     &     mstop(2), igamtop(2), idxop(2), iblkop(2),
     &     mstop1op2, igamtop1op2, njoined,
     &     njoined_op(2), njoined_op1op2,
     &     njoined_cnt, njoined_res
      real(8), pointer ::
     &     xret_blk(:), xret_pnt(:)
      real(8), target ::
     &     xret_scr(1)
      integer, pointer ::
     &     op2list(:),
     &     occ_vtx(:,:,:), irestr_vtx(:,:,:,:,:), info_vtx(:,:),
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
c      type(filinf), pointer ::
c     &     ffscr(:)
      type(operator), pointer ::
     &     opres
c      type(operator), pointer ::
c     &     op1, op2, op1op2, opres, op1op2tmp
      type(me_list), pointer ::
     &     me_op1, me_op2, me_op1op2, me_res, me_op1op2tmp
      type(filinf), pointer ::
     &     ffop1, ffop2, ffop1op2, ffres
      type(formula_item), pointer ::
     &     cur_form
      type(contraction) ::
     &     cur_contr

      integer, external ::
     &     idxlist
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

c      ! open files
c      call file_open(fffrm)

c      if (iprint.ge.10) then
c        write(luout,*) 'formula file: ',
c     &       fffrm%name(1:len_trim(fffrm%name))
c      end if

c      lufrm = fffrm%unit
c      rewind lufrm

c      read(lufrm) len,title(1:len)

c      if (iprint.ge.5)
c     &     write(luout,*) 'Evaluating:',title(1:len)

      cur_form => flist
c      allocate(cur_form%contr)
c      call init_contr(cur_form%contr)

      skip = .false.
      nres  = 0
      idxres = 0
      iterm = 0
      nullify(xret_blk)

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
     &     call quit(1,'frm_sched','first command must define target')

        select case(cur_form%command)
        case(command_end_of_formula)
          ! get xret value for final target
          if (type_xret.eq.1) then
            xret(idxres) = sqrt(sum(xret_blk(1:nblk_res)))
          else if (type_xret.eq.2) then
            xret(idxres) = xret_blk(1)
          end if

          exit term_loop

        case(command_set_target_init)

          ! for previous target: assemble xret value
          if (idxres.gt.0.and..not.skip) then
            if (type_xret.eq.1) then
              xret(idxres) = sqrt(sum(xret_blk(1:nblk_res)))
            else if (type_xret.eq.2) then
              xret(idxres) = xret_blk(1)
            end if
          end if

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
          
          call symmetrise(1d0,me_res,me_res,xret_blk,op_info,orb_info)

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
        if (me_op1op2%len_op_occ(cur_contr%iblk_res).eq.0)
     &       cycle term_loop
        ! check here for other zero blocks as well ...

        if (nfact.eq.0) then

          iblkres = (cur_contr%iblk_res-1)/njoined_res + 1
          idxop(1) = cur_contr%vertex(1)%idx_op
          iblkop(1) = cur_contr%vertex(1)%iblk_op
          tra_op1 = cur_contr%vertex(1)%dagger
          tra_op1op2 = cur_contr%dagger

          ! special: unit operator
          if (ops(idxop(1))%op%name.eq.op_unity) then
            call add_unity(fac,me_res,iblkres,orb_info)
          else
            idxmel = op2list(idxop(1))
            if (mel_arr(idxmel)%mel%fhand%unit.le.0)
     &           call file_open(mel_arr(idxmel)%mel%fhand)
c fix:
            njoined = mel_arr(idxmel)%mel%op%njoined
            iblkop(1) = (iblkop(1)-1)/njoined + 1
c fix:
            if (tra_op1.xor.tra_op1op2) then
              call add_opblk_transp(xret_blk(iblkres),fac,
     &             mel_arr(idxmel)%mel,me_res,tra_op1,tra_op1op2,
     &             iblkop(1),iblkres,
     &             op_info,str_info,orb_info)
            else
              call add_opblk(xret_blk(iblkres),fac,
     &             mel_arr(idxmel)%mel,me_res,
     &             iblkop(1),iblkres,orb_info)
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
     &       iocc_ex2(ngastp,2,nvtx), iocc_cnt(ngastp,2,nvtx),
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
c dbg
c        print *,'info from vtxinf4contr:'
c        print *,'result:'
cc        do idum = 1, njoined_res
c          call wrt_rstr(6,irestr_vtx(1,1,1,1,idum),ngas)
c        end do
c        print *,'op-vertices:'
c        do idum = njoined_res+1, njoined_res+nvtx
c          call wrt_rstr(6,irestr_vtx(1,1,1,1,idum),ngas)
c        end do
c dbg

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
          call get_bc_info2(bc_sign,
     &         idxop,iblkop,
     &         iocc_ex1,iocc_ex2,iocc_cnt,
     &         iocc_op1,iocc_op2,iocc_op1op2,
     &         irst_op1,irst_op2,irst_op1op2,
     &         tra_op1,tra_op2,tra_op1op2,
     &         mstop,mstop1op2,
     &         igamtop,igamtop1op2,
     &         njoined_op, njoined_op1op2, njoined_cnt,
     &         merge_op1,merge_op2,merge_op1op2, merge_op2op1,
     &         cur_contr,njoined_res,
     &                        occ_vtx,irestr_vtx,info_vtx,iarc,
     &         irst_res,orb_info)

          ! set up reduced contraction after 
          ! current binary contraction
          if (idx.ne.nfact) then
            ivtx_new = cur_contr%inffac(3,idx)
            idxop_intm = -ninter

            ! reset reo_info
            call init_reo_info(reo_info)
            
            call reduce_contr(cur_contr,occ_vtx,
     &           possible,
     &           iarc,idxop_intm,ivtx_new,
     &           njoined_res,
     &           .false.,idum,idum,
     &           .true.,irestr_vtx,info_vtx,irst_res,
     &           .true.,reo_info,orb_info)
            if (.not.possible)
     &           call quit(1,'frm_sched1',
     &           'inconsistency: reduce_contr is in difficulties')

            ! add 0-contractions, if necessary
            call check_disconnected(cur_contr)
            ! process reordering info
            call get_reo_info(reo_op1op2,reo_other,
     &           iocc_op1op2,iocc_op1op2tmp,
     &           irst_op1op2,irst_op1op2tmp,
     &           njoined_op1op2,
     &           cur_contr,reo_info,str_info,orb_info)
          else
            reo_op1op2 = .false.
            reo_other = .false.
            iocc_op1op2tmp = iocc_op1op2
            irst_op1op2tmp = irst_op1op2
          end if

          ! set up operator 1 and 2
          do iops = 1, 2
            if (idxop(iops).gt.0) then
              idxmel = op2list(idxop(iops))
              ! primary operator or long-term intermediate
c              if (iops.eq.1) ffop1 => ffops(idxop(iops))%fhand
c              if (iops.eq.2) ffop2 => ffops(idxop(iops))%fhand
              if (iops.eq.1) me_op1 => mel_arr(idxmel)%mel
              if (iops.eq.2) me_op2 => mel_arr(idxmel)%mel
            else
              ! intermediate for current contraction only
c              if (iops.eq.1) ffop1 => ffscr(-idxop(iops))
c              if (iops.eq.2) ffop2 => ffscr(-idxop(iops))
              if (iops.eq.1) me_op1 => melscr(-idxop(iops))
              if (iops.eq.2) me_op2 => melscr(-idxop(iops))
            end if
          end do

          if (me_op1%fhand%unit.le.0)
     &             call file_open(me_op1%fhand)
          if (me_op2%fhand%unit.le.0)
     &             call file_open(me_op2%fhand)

          ! set up result
          if (idx.eq.nfact) then
            ! last operation: store on result array
            update = .true.
            idxop1op2 = cur_contr%idx_res
            iblkop1op2 = cur_contr%iblk_res
c dbg
c            print *,'SET TO ',iblkop1op2,' IN RES BRANCH:'
c dbg
            idxmel = op2list(idxop1op2)
            me_op1op2 => mel_arr(idxmel)%mel
            xret_pnt => xret_blk(iblkop1op2:iblkop1op2)
            type_xret_cur = type_xret
          else
            ! new intermediate
            update = .false.
c            write(name_scr,'(a,i3.3,".da ")') name_scr0,ninter
c            call file_init(ffscr(ninter),name_scr,1,lblk_da)
c            call file_open(ffscr(ninter))
c            if (ntest.ge.100)
c     &           write(luout,*) 'new intermediate file: ',
c     &           ffscr(ninter)%name(1:len_trim(ffscr(ninter)%name))

c            ffop1op2 => ffscr(ninter)
            iblkop1op2 = 1

            ! set up pseudo-operator for current intermediate
            ! refers to reordered operator (if this matters)
            write(opscrnam,'("INT",i3.3)') ninter
            call set_ps_op(opscr(ninter),opscrnam,
     &           iocc_op1op2,irst_op1op2,njoined_op1op2,orb_info)
            melscr(ninter)%op => opscr(ninter)
            call set_ps_list(melscr(ninter),opscrnam,
     &           0,0,mstop1op2,igamtop1op2,0,
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
            call set_ps_list(meltmp,opscrnam,
     &           0,0,mstop1op2,igamtop1op2,0,
     &           str_info,strmap_info,orb_info)
            me_op1op2tmp => meltmp
          else
            me_op1op2tmp => me_op1op2
          end if

          ! translate records to offset in file:
          ! (makes life easier in case we once decide to use
          ! one scratch file only: no changes to contr_op1op2 necessary)
          ffop1 => me_op1%fhand
          ffop2 => me_op2%fhand
          ffop1op2 => me_op1op2%fhand
          idoffop1 = ffop1%length_of_record*(ffop1%current_record-1)
          idoffop2 = ffop2%length_of_record*(ffop2%current_record-1)
          idoffop1op2 = ffop1op2%length_of_record*
     &                                   (ffop1op2%current_record-1)

          if (ntest.ge.100)
     &         write(luout,*) 'calling contraction kernel'
          ! do the contraction
          call contr_op1op2(facc,bc_sign,
     &       update,xret_pnt,type_xret_cur,
     &       me_op1,me_op2,me_op1op2, me_op1op2tmp,
     &       tra_op1, tra_op2, tra_op1op2,
     &       iblkop(1),iblkop(2),iblkop1op2,iblkop1op2,
     &       idoffop1,idoffop2,idoffop1op2,
     &       iocc_ex1,iocc_ex2,iocc_cnt,
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

c dbg
c          if(idxopres.eq.8)then
c            write(luout,*)'xret', xret_pnt(1)
c          endif
c dbg

          if (reo_op1op2) then
            call dealloc_me_list(meltmp)
            call dealloc_operator(optmp)
            call dealloc_reo_info(reo_info)
          end if

        end do bin_loop

        deallocate(
     &       occ_vtx,irestr_vtx,info_vtx,
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

      ifree = mem_flushmark()

      if (ntest.ge.100)
     &     write(luout,*) 'returning from frm_sched1'

      return
      end
