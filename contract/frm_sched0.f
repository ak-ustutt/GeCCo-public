*----------------------------------------------------------------------*
      subroutine frm_sched0(xret,fffrm,
     &         op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     schedule the evaluation of a formula file
*     
*     version 0
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula.h'

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
      type(orbinf) ::
     &     orb_info

      type(formula) ::
     &     cur_form

      logical ::
     &     update
      integer ::
     &     lufrm, idxopres, idxres, nres,
     &     n_occ_cls, maxvtx, maxarc, maxfac,
     &     nfact, idxop1op2, iblkop1op2, iops, iblkres,
     &     ninter, idx, nsym, ngas, nexc, ndis, iprint, iterm, len
      real(8) ::
     &     fac, facc, xnrm
      character ::
     &     title*256, opscrnam*8

      type(file_array), pointer ::
     &     ffops(:)
      type(operator_array), pointer ::
     &     ops(:)
      integer ::
     &     iocc_op(ngastp,2,2), irst_op(2,orb_info%ngas,2,2,2),
     &     irst_res(2,orb_info%ngas,2,2),
     &     mstop(2), igamtop(2), idxop(2), iblkop(2),
     &     iocc_ext(ngastp,2,2), iocc_cnt(ngastp,2)
      integer, allocatable ::
     &     interm(:), iocc_op1op2(:,:,:), irst_op1op2(:,:,:,:,:),
     &     mstop1op2(:), igamtop1op2(:)
      type(operator), pointer ::
     &     opscr(:)
      type(filinf), pointer ::
     &     ffscr(:)
      type(operator), pointer ::
     &     op1, op2, op1op2, opres
      type(filinf), pointer ::
     &     ffop1, ffop2, ffop1op2, ffres

      logical, external ::
     &     rd_formula
      real(8), external ::
     &     xnormop

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

c      call file_init(ffinterm,name_interm,1,lblk_da)
c      call file_open(ffinterm)

      lufrm = fffrm%unit
      rewind lufrm

      read(lufrm) len,title(1:len)
c      read(lufrm) idxopres,n_occ_cls,maxvtx,maxarc,maxfac

      if (iprint.ge.5)
     &     write(luout,*) 'Evaluating:',title(1:len)
c      if (ntest.ge.10)
c     &     write(luout,'(x,a,5i6)')
c     &     'idxopres,n_occ_cls,maxvtx,maxarc,maxfac:',
c     &     idxopres,n_occ_cls,maxvtx,maxarc,maxfac

      
      allocate(cur_form%contr)
      cur_form%contr%mxvtx = 0
      cur_form%contr%mxarc = 0
      cur_form%contr%mxfac = 0

      nres  = 0
      idxres = 0
      iterm = 0
      
      ! loop over entries
      term_loop: do while(rd_formula(fffrm,cur_form))

        if (nres.eq.0 .and.
     &      cur_form%command.ne.command_set_target_init)
     &     call quit(1,'frm_sched','first command must define target')

        select case(cur_form%command)
        case(command_end_of_formula)
          exit term_loop
        case(command_set_target_init)
c quick hack:
          if (idxres.gt.0) then
            xret(idxres) = xnormop(ffres,ops(idxopres)%op)
          end if
c
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
            if (ffres%unit.le.0)
     &           call file_open(ffres)
            call zeroop(ffres,opres)
          end if
          xret(idxres) = 0d0
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
     &       call prt_contr2(luout,cur_form%contr,op_info%op_arr)
         
        ! process info
        fac = cur_form%contr%fac
        nfact = cur_form%contr%nfac
        
        if (ntest.ge.100) write(luout,*) 'nfact, fac: ',nfact,fac

        if (nfact.eq.0) then

          iblkres = cur_form%contr%iblk_res
          idxop(1) = cur_form%contr%vertex(1)%idx_op
          iblkop(1) = cur_form%contr%vertex(1)%iblk_op

          if (ffops(idxop(1))%fhand%unit.le.0)
     &         call file_open(ffops(idxop(1))%fhand)

          call add_opblk(fac,ffops(idxop(1))%fhand,ffres,
     &         ops(idxop(1))%op,opres,
     &         iblkop(1),iblkres,orb_info)

          cycle term_loop
        end if

        ! preliminary fix to set irst_res (needed in get_bc_info)
        call set_restr_prel(irst_res,
     &       cur_form%contr,op_info,orb_info%ihpvgas,orb_info%ngas)

        ! allocate arrays for intermediates
        allocate(interm(nfact),iocc_op1op2(ngastp,2,nfact),
     &       irst_op1op2(2,ngas,2,2,nfact),
     &       mstop1op2(nfact),igamtop1op2(nfact))
        if (nfact.gt.1)
     &       allocate(opscr(nfact-1),ffscr(nfact-1))

        ! reset intermediate counter
        ninter = 0

        ! loop over binary contractions
        bin_loop: do idx = 1, nfact

          ! use factor in last contraction
          facc = 1d0
          if (idx.eq.nfact) facc=fac

          if (iprint.ge.20) write(luout,*) '    contr #',idx

          ! set up info for binary contraction
          call get_bc_info(iocc_op,irst_op,
     &       mstop,igamtop,idxop,iblkop,
     &       iocc_ext,iocc_cnt,
     &       cur_form%contr%inffac(1,idx),cur_form%contr,
     &                                    op_info,irst_res,
     &       interm,ninter,iocc_op1op2,irst_op1op2,
     &       mstop1op2,igamtop1op2,
     &       orb_info%ihpvgas,ngas)

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
            allocate(opscr(ninter)%ihpvca_occ(ngastp,2,1),
     &               opscr(ninter)%ica_occ(2,1),
     &               opscr(ninter)%igasca_restr(2,ngas,2,2,1),
     &               opscr(ninter)%len_op_occ(1),
     &               opscr(ninter)%idx_graph(ngastp,2,1))
            ! set up pseudo-operator for current intermediate
            write(opscrnam,'("INT",i3.3)') ninter
            call set_ps_op(opscr(ninter),opscrnam,
     &           iocc_op1op2(1,1,ninter),irst_op1op2(1,1,1,1,ninter),
     &           mstop1op2(ninter),igamtop1op2(ninter),
     &           ngas,orb_info%ihpvgas,str_info)
            ! allocate sub-arrays
            allocate(opscr(ninter)%off_op_occ(1),
     &               opscr(ninter)%off_op_gmo(1),
     &               opscr(ninter)%len_op_gmo(1),
     &               opscr(ninter)%off_op_gmox(1),
     &               opscr(ninter)%len_op_gmox(1))
            nexc = min(opscr(ninter)%ica_occ(1,1),
     &                 opscr(ninter)%ica_occ(2,1))
            allocate(opscr(ninter)%len_op_gmo(1)%gam_ms(nsym,nexc+1),
     &               opscr(ninter)%off_op_gmo(1)%gam_ms(nsym,nexc+1))
            ! set up dimensions (pass 1)
            call set_op_dim(1,.false.,opscr(ninter),
     &           str_info,nsym)
            ! some more sub-arrays
            ndis = opscr(ninter)%off_op_gmox(1)%maxd
            allocate(
     &          opscr(ninter)%len_op_gmox(1)%d_gam_ms(ndis,nsym,nexc+1),
     &          opscr(ninter)%off_op_gmox(1)%d_gam_ms(ndis,nsym,nexc+1),
     &          opscr(ninter)%off_op_gmox(1)%did(ndis,nsym,nexc+1),
     &          opscr(ninter)%off_op_gmox(1)%ndis(nsym,nexc+1))
            ! set up dimensions (pass 2)
            call set_op_dim(2,.false.,opscr(ninter),
     &           str_info,nsym)
            op1op2 => opscr(ninter)
          end if

          if (ntest.ge.100)
     &         write(luout,*) 'calling contraction kernel'
          ! do the contraction
          call contr_op1op2(facc,ffop1,ffop2,
     &       update,ffop1op2,xret(idxres),
     &       op1,op2,op1op2,
     &       iblkop(1),iblkop(2),iblkop1op2,
     &       iocc_ext(1,1,1),iocc_ext(1,1,2),iocc_cnt,
     &       irst_op(1,1,1,1,1),irst_op(1,1,1,1,2),
     &                 irst_op1op2(1,1,1,1,ninter),
     &       mstop(1),mstop(2),mstop1op2(ninter),
     &       igamtop(1),igamtop(2),igamtop1op2(ninter),
     &       str_info,orb_info)
          if (ntest.ge.100)
     &         write(luout,*) 'returned from contraction kernel'

        end do bin_loop

        deallocate(interm,iocc_op1op2,irst_op1op2,
     &       mstop1op2,igamtop1op2)
        if (nfact.gt.1) then
          ! a rather awkward deallocation section follows
          ! looking for more elegant solutions ....
          do idx = 1, ninter-1
            call file_close_delete(ffscr(idx))
            deallocate(
     &          opscr(idx)%len_op_gmox(1)%d_gam_ms,
     &          opscr(idx)%off_op_gmox(1)%d_gam_ms,
     &          opscr(idx)%off_op_gmox(1)%did,
     &          opscr(idx)%off_op_gmox(1)%ndis)
            deallocate(
     &          opscr(idx)%len_op_gmo(1)%gam_ms,
     &          opscr(idx)%off_op_gmo(1)%gam_ms)
            deallocate(
     &          opscr(idx)%off_op_occ,
     &          opscr(idx)%off_op_gmo,
     &          opscr(idx)%len_op_gmo,
     &          opscr(idx)%off_op_gmox,
     &          opscr(idx)%len_op_gmox)
            deallocate(
     &          opscr(idx)%ihpvca_occ,
     &          opscr(idx)%ica_occ,
     &          opscr(idx)%igasca_restr,
     &          opscr(idx)%len_op_occ,
     &          opscr(idx)%idx_graph)
          end do
          deallocate(opscr,ffscr)
        end if

      end do term_loop

      call dealloc_contr(cur_form%contr)
      deallocate(cur_form%contr)

      if (idxopres.gt.0) then
        ! return norm of result
        xret(idxres) = xnormop(ffres,ops(idxopres)%op)
      end if
c      do iops = 1, op_info%nops
c        if (ffops(iops)%fhand%unit.gt.0)
c     &       call file_close_keep(ffops(iops)%fhand)
c      end do
      call file_close_keep(fffrm)

      if (ntest.ge.100)
     &     write(luout,*) 'returning from frm_sched0'

      return
      end
