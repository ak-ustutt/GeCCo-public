*----------------------------------------------------------------------*
      subroutine frm_sched1(xret,fffrm,
     &         op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     schedule the evaluation of a formula file
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
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      character, parameter ::
     &     name_scr0*6 = 'cntscr'

      character ::
     &     name_scr*13

      real(8), intent(inout) ::
     &     xret(*)
      type(filinf), intent(inout) ::
     &     fffrm
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info

      type(formula_item) ::
     &     cur_form

      logical ::
     &     update, reo_op1op2, reo_other, possible
      integer ::
     &     lufrm, idxopres, idxres, nres, type_xret, type_xret_cur,
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
      type(file_array), pointer ::
     &     ffops(:)
      type(operator_array), pointer ::
     &     ops(:)
      integer ::
     &     mstop(2), igamtop(2), idxop(2), iblkop(2),
     &     mstop1op2, igamtop1op2,
     &     njoined_op(2), njoined_op1op2,
     &     njoined_cnt, njoined_res
      real(8), pointer ::
     &     xret_blk(:), xret_pnt(:)
      real(8), target ::
     &     xret_scr(1)
      integer, pointer ::
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
      type(filinf), pointer ::
     &     ffscr(:)
      type(operator), pointer ::
     &     op1, op2, op1op2, opres, op1op2tmp
      type(filinf), pointer ::
     &     ffop1, ffop2, ffop1op2, ffres
      
      logical, external ::
     &     rd_formula
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

      ffops => op_info%opfil_arr
      ops => op_info%op_arr

      ! open files
      call file_open(fffrm)

      if (iprint.ge.10) then
        write(luout,*) 'formula file: ',
     &       fffrm%name(1:len_trim(fffrm%name))
      end if

      lufrm = fffrm%unit
      rewind lufrm

      read(lufrm) len,title(1:len)

      if (iprint.ge.5)
     &     write(luout,*) 'Evaluating:',title(1:len)

      
      allocate(cur_form%contr)
      call init_contr(cur_form%contr)

      nres  = 0
      idxres = 0
      iterm = 0
      nullify(xret_blk)

      ! loop over entries
      term_loop: do while(rd_formula(fffrm,cur_form))

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
          if (idxres.gt.0) then
            if (type_xret.eq.1) then
              xret(idxres) = sqrt(sum(xret_blk(1:nblk_res)))
            else if (type_xret.eq.2) then
              xret(idxres) = xret_blk(1)
            end if
          end if

          ! initialize result
          nres = nres+1
          idxres = nres
          idxopres = cur_form%target
          if (iprint.ge.10) then
            write(luout,*) 'New target: ',idxopres
          end if
          if (idxopres.gt.0) then
            ffres => op_info%opfil_arr(idxopres)%fhand
            opres => op_info%op_arr(idxopres)%op
            nblk_res = opres%n_occ_cls
            type_xret = 2
            if (opres%len_op.gt.1) type_xret = 1
            if (ffres%unit.le.0)
     &           call file_open(ffres)
            call zeroop(ffres,opres)
          else
            nblk_res = 1
            type_xret = 2
          end if
 
          if (associated(xret_blk))
     &           ifree = mem_dealloc('xret_blk')
          ifree = mem_alloc_real(xret_blk,nblk_res,'xret_blk')
          xret_blk(1:nblk_res) = 0d0
 
c          xret(idxres) = 0d0
          cycle term_loop
c        case(command_set_target_update)
        case(command_add_contribution)
        case default
          write(luout,*) 'command = ',cur_form%command
          call quit(1,'frm_sched','command not defined/implemented')
        end select

        iterm = iterm+1
        if (iprint.ge.20)
     &     write(luout,*) '   term #',iterm

        if (ntest.ge.50)
     &       call prt_contr2(luout,cur_form%contr,op_info)
         
        ! process info
        fac = cur_form%contr%fac
        nfact = cur_form%contr%nfac
        
        idxopres = cur_form%contr%idx_res 
        njoined_res = ops(idxopres)%op%njoined

        if (ntest.ge.100) write(luout,*) 'nfact, fac: ',nfact,fac

        if (nfact.eq.0) then

          iblkres = cur_form%contr%iblk_res
          idxop(1) = cur_form%contr%vertex(1)%idx_op
          iblkop(1) = cur_form%contr%vertex(1)%iblk_op

          ! special: unit operator
          if (ops(idxop(1))%op%name.eq.op_unity) then
            call add_unity(fac,ffres,opres,iblkres,orb_info)
          else
            if (ffops(idxop(1))%fhand%unit.le.0)
     &         call file_open(ffops(idxop(1))%fhand)

            call add_opblk(fac,ffops(idxop(1))%fhand,ffres,
     &         ops(idxop(1))%op,opres,
     &         iblkop(1),iblkres,orb_info)
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
     &       cur_form%contr,op_info,orb_info%ihpvgas,orb_info%ngas)
        
        ! allocate arrays for intermediates
        allocate(
     &       occ_vtx(ngastp,2,nvtx+njoined_res),
     &       irestr_vtx(2,orb_info%ngas,2,2,nvtx+njoined_res),
     &       info_vtx(2,nvtx+njoined_res),
     &       merge_op1(nvtx*nvtx*10), ! a bit too large, I guess ...
     &       merge_op2(nvtx*nvtx*10),
     &       merge_op1op2(nvtx*nvtx*10),
     &       merge_op2op1(nvtx*nvtx*10))
        if (nfact.gt.1)
     &       allocate(opscr(nfact-1),optmp,ffscr(nfact-1))

        ! reset intermediate counter
        ninter = 0

        call occvtx4contr(0,occ_vtx,cur_form%contr,op_info)
        call vtxinf4contr(irestr_vtx,info_vtx,
     &                            cur_form%contr,op_info,ngas)

        fac = cur_form%contr%fac

        ! add 0-contractions, if necessary
        call check_disconnected(cur_form%contr)

        ! loop over binary contractions
        bin_loop: do idx = 1, nfact

          ! use factor in last contraction
          facc = 1d0
          if (idx.eq.nfact) facc=fac

          if (iprint.ge.20) write(luout,*) '    contr #',idx
          iarc = cur_form%contr%inffac(5,idx)

          ninter = ninter + 1

          ! set up info for binary contraction
          call get_bc_info2(bc_sign,
     &         idxop,iblkop,
     &         iocc_ex1,iocc_ex2,iocc_cnt,
     &         iocc_op1,iocc_op2,iocc_op1op2,
     &         irst_op1,irst_op2,irst_op1op2,
     &         mstop,mstop1op2,
     &         igamtop,igamtop1op2,
     &         njoined_op, njoined_op1op2, njoined_cnt,
     &         merge_op1,merge_op2,merge_op1op2, merge_op2op1,
     &         cur_form%contr,njoined_res,
     &                        occ_vtx,irestr_vtx,info_vtx,iarc,
     &         irst_res,orb_info%ihpvgas,ngas)

          ! set up reduced contraction after 
          ! current binary contraction
          if (idx.ne.nfact) then
            ivtx_new = cur_form%contr%inffac(3,idx)
            idxop_intm = -ninter

            ! reset reo_info
            call init_reo_info(reo_info)
            
            call reduce_contr(cur_form%contr,occ_vtx,
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
            call check_disconnected(cur_form%contr)
            ! process reordering info
            call get_reo_info(reo_op1op2,reo_other,
     &           iocc_op1op2,iocc_op1op2tmp,
     &           irst_op1op2,irst_op1op2tmp,
     &           njoined_op1op2,
     &           cur_form%contr,reo_info,str_info,orb_info)
          else
            reo_op1op2 = .false.
            reo_other = .false.
            iocc_op1op2tmp = iocc_op1op2
            irst_op1op2tmp = irst_op1op2
          end if

          ! set up operator 1 and 2
          do iops = 1, 2
            if (idxop(iops).gt.0) then
              ! primary operator or long-term intermediate
              if (iops.eq.1) ffop1 => ffops(idxop(iops))%fhand
              if (iops.eq.2) ffop2 => ffops(idxop(iops))%fhand
              if (iops.eq.1) op1 => ops(idxop(iops))%op
              if (iops.eq.2) op2 => ops(idxop(iops))%op
            else
              ! intermediate for current contraction only
              if (iops.eq.1) ffop1 => ffscr(-idxop(iops))
              if (iops.eq.2) ffop2 => ffscr(-idxop(iops))
              if (iops.eq.1) op1 => opscr(-idxop(iops))
              if (iops.eq.2) op2 => opscr(-idxop(iops))
            end if
          end do

          if (ffop1%unit.le.0)
     &             call file_open(ffop1)
          if (ffop2%unit.le.0)
     &             call file_open(ffop2)

          ! set up result
          if (idx.eq.nfact) then
            ! last operation: store on result array
            update = .true.
            idxop1op2 = cur_form%contr%idx_res
            iblkop1op2 = cur_form%contr%iblk_res
            ffop1op2 => ffres
            op1op2 => ops(idxop1op2)%op
            xret_pnt => xret_blk(iblkop1op2:iblkop1op2)
            type_xret_cur = type_xret
          else
            ! new intermediate
            update = .false.
            write(name_scr,'(a,i3.3,".da ")') name_scr0,ninter
            call file_init(ffscr(ninter),name_scr,1,lblk_da)
            call file_open(ffscr(ninter))
            if (ntest.ge.100)
     &           write(luout,*) 'new intermediate file: ',
     &           ffscr(ninter)%name(1:len_trim(ffscr(ninter)%name))

            ffop1op2 => ffscr(ninter)
            iblkop1op2 = 1

            ! allocate further stuff in operator structure
c            opscr(ninter)%n_occ_cls = 1
c            opscr(ninter)%njoined = 1
c            call init_operator(0,opscr(ninter),orb_info)
            ! set up pseudo-operator for current intermediate
            ! refers to reordered operator (if this matters)
            write(opscrnam,'("INT",i3.3)') ninter
            call set_ps_op(opscr(ninter),opscrnam,
     &           iocc_op1op2,irst_op1op2,njoined_op1op2,
     &           mstop1op2,igamtop1op2,
     &           orb_info,str_info)
            ! allocate sub-arrays
            call init_operator(1,opscr(ninter),orb_info)
            ! set up dimensions (pass 1)
            call set_op_dim2(1,opscr(ninter),str_info,nsym)
            ! some more sub-arrays
            call init_operator(2,opscr(ninter),orb_info)
            ! set up dimensions (pass 2)
            call set_op_dim2(2,opscr(ninter),str_info,nsym)
            op1op2 => opscr(ninter)
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
     &           mstop1op2,igamtop1op2,
     &           orb_info,str_info)
            ! and the allocation sequence
            call init_operator(1,optmp,orb_info)
            call set_op_dim2(1,optmp,str_info,nsym)
            call init_operator(2,optmp,orb_info)
            call set_op_dim2(2,optmp,str_info,nsym)
            op1op2tmp => optmp
          else
            op1op2tmp => op1op2
          end if

          ! translate records to offset in file:
          ! (makes life easier in case we once decide to use
          ! one scratch file only: no changes to contr_op1op2 necessary)
          idoffop1 = ffop1%length_of_record*(ffop1%current_record-1)
          idoffop2 = ffop2%length_of_record*(ffop2%current_record-1)
          idoffop1op2 = ffop1op2%length_of_record*
     &                                   (ffop1op2%current_record-1)

          if (ntest.ge.100)
     &         write(luout,*) 'calling contraction kernel'
          ! do the contraction
          call contr_op1op2(facc,bc_sign,ffop1,ffop2,
     &       update,ffop1op2,xret_pnt,type_xret_cur,
     &       op1,op2,op1op2, op1op2tmp,
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

          if (reo_op1op2) then
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
            call file_close_delete(ffscr(idx))
            call dealloc_operator(opscr(idx))
          end do
          deallocate(opscr,optmp,ffscr)
        end if

      end do term_loop

      call dealloc_contr(cur_form%contr)
      deallocate(cur_form%contr)
 
      call file_close_keep(fffrm)

      ifree = mem_flushmark()

      if (ntest.ge.100)
     &     write(luout,*) 'returning from frm_sched1'

      return
      end
