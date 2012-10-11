*----------------------------------------------------------------------*
      subroutine frm_sched2(xret,flist,depend_info,idxselect,nselect,
     &         op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     schedule the evaluation of a formula
*     
*     version 2 -- for new optimized formulae
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
     &     update, add, reo, skip, call_sti_remover, measure
      integer ::
     &     idxopres, idxres, nres, type_xret,
     &     idxme_res, nblk_res, ifree,
     &     iterm, icmd, icmd_pfi, iprint
      real(8) ::
     &     xret_last

      type(me_list_array), pointer ::
     &     mel_arr(:)
      real(8), pointer ::
     &     xret_blk(:), xret_pnt(:)
      real(8), target ::
     &     xret_scr(1)
      integer, pointer ::
     &     op2list(:)
      type(operator), pointer ::
     &     opres
      type(me_list), pointer ::
     &     me_res
      type(filinf), pointer ::
     &     ffres
      type(formula_item), pointer ::
     &     cur_form

      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0,
     &     cpus,syss,walls,cpus0,syss0,walls0

      integer, external ::
     &     idxlist, idx_oplist2
      logical, external ::
     &     me_list_uptodate

      ifree = mem_setmark('frm_sched0')
      
      iprint = max(ntest,iprlvl)
      if (ntest.gt.0.or.iprint.ge.5) then
        write(luout,*) '============================='
        write(luout,*) '= Entered formula scheduler ='
        write(luout,*) '============================='
      end if

      mel_arr => op_info%mel_arr
      op2list => op_info%op2list

      cur_form => flist

      skip = .false.
      nres  = 0
      idxres = 0
      iterm = 0
      icmd = 0
      icmd_pfi = 0 ! only for ntest.ge.100
      nullify(xret_blk)
      xret_last = 0d0


      ! loop over entries
      main_loop: do 

        if (lustat.gt.0) call atim_csw(cpus0,syss0,walls0)
        measure = .false.

        if (nres.eq.0 .and.
     &      cur_form%command.ne.command_set_target_init)
     &     call quit(1,'frm_sched2','first command must define target')

        update = .false.
        call_sti_remover = .false.

        if (ntest.ge.100)
     &       call print_form_item(luout,icmd_pfi,cur_form,op_info)

        select case(cur_form%command)
        case(command_end_of_formula,command_set_target_init)

          ! set pointers again, as the targets may change
          mel_arr => op_info%mel_arr
          op2list => op_info%op2list

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
     &         exit main_loop

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

          idxme_res = op2list(idxopres)  ! list index of result
          if (iprint.ge.10.and.skip) then
            write(luout,*) 'Skipping target: ',
     &           trim(mel_arr(idxme_res)%mel%label)
          else if (iprint.ge.10) then
            write(luout,*) 'New target: ',
     &           trim(mel_arr(idxme_res)%mel%label)
          end if

          if (skip) then
            ! fast-forward to next [init target]
            do while(
     &           cur_form%next%command.ne.command_set_target_init.and.
     &           cur_form%next%command.ne.command_end_of_formula)
              cur_form => cur_form%next
            end do

          else

            me_res => op_info%mel_arr(idxme_res)%mel
            opres  => me_res%op
            ffres  => me_res%fhand
            nblk_res = opres%n_occ_cls
            type_xret = 2
            if (me_res%len_op.gt.1) type_xret = 1
            if (ffres%unit.le.0)
     &           call file_open(ffres)
            call zeroop(me_res)
 
            if (lustat.gt.0) write(lustat,*) 'NEW TARGET: ',
     &           trim(me_res%label)

            if (associated(xret_blk))
     &           ifree = mem_dealloc('xret_blk')
            ifree = mem_alloc_real(xret_blk,nblk_res,'xret_blk')
            xret_blk(1:nblk_res) = 0d0
 
            ! mark current term as updated:
            call touch_file_rec(ffres)

          end if

        case(command_new_intermediate)

          call fs_newintm_drv(cur_form,
     &         op_info,str_info,strmap_info,orb_info)

        case(command_del_intermediate)

          call fs_delintm_drv(cur_form,op_info)

        case(command_reorder,command_add_reo)

          measure = .true.
          !update = cur_form%command.eq.command_add_reo
          update = idx_oplist2(cur_form%reo%label_out,op_info)
     &                 .eq.idxopres
          if (cur_form%command.eq.command_add_reo) icmd = icmd+1
          call fs_reo_drv(xret_blk,type_xret,idxopres,me_res,
     &         cur_form,update,
     &         op_info,str_info,strmap_info,orb_info)

          call_sti_remover = .true.!cur_form%command.eq.command_reorder

        case(command_add_intm,command_cp_intm)

          measure = .true.
          icmd = icmd+1
          !update = cur_form%command.eq.command_add_intm
          update = idx_oplist2(cur_form%bcontr%label_res,op_info)
     &                 .eq.idxopres

          call fs_add_drv(xret_blk,type_xret,idxopres,me_res,
     &         cur_form,update,
     &         op_info,str_info,strmap_info,orb_info)

          call_sti_remover = .true.

        case(command_bc,command_add_bc,
     &       command_bc_reo,command_add_bc_reo)

          measure = .true.
          icmd = icmd+1
          add    = (cur_form%command.eq.command_add_bc.or.
     &              cur_form%command.eq.command_add_bc_reo)
          update = idx_oplist2(cur_form%bcontr%label_res,op_info)
     &                 .eq.idxopres
          reo    = (cur_form%command.eq.command_bc_reo.or.
     &              cur_form%command.eq.command_add_bc_reo)
          
          call fs_contr_drv(xret_blk,
     &                             type_xret,idxopres,me_res,
     &         cur_form,update,add,reo,
     &         op_info,str_info,strmap_info,orb_info)

          call_sti_remover = .true.

        case(command_symmetrise)
          
          measure = .true.
          call symmetrise(1d0,me_res,me_res,
     &         xret_blk,
     &         op_info,str_info,orb_info)

        case(command_add_contribution)
          call quit(1,'frm_sched2',
     &       'command not valid in this version of frm_sched!: [CONTR]')

        case default
          write(luout,*) 'command = ',cur_form%command
          call quit(1,'frm_sched','command not defined/implemented')
        end select

        ! remove short-term intermediates
        if (call_sti_remover)
     &       call fs_sti_remover(cur_form,op_info)

        if (update) then
          ! we count only terms that actually update a target
          iterm = iterm + 1
          if (type_xret.eq.2.and.iprint.ge.3)
     &         write(luout,'(1x,"term # ",i5,":",2(x,g19.10))')
     &         iterm, xret_blk(1)-xret_last, xret_blk(1)
          if (type_xret.eq.2) xret_last = xret_blk(1)

          if (ntest.ge.50) then
            write(luout,*) 'xret after term ',iterm
            write(luout,'(x,4f19.10)') xret_blk(1:nblk_res)
          end if
        end if

        if (lustat.gt.0.and.measure) then
          call atim_csw(cpus,syss,walls)
          write(lustat,'(i8,3f18.4)')
     &         icmd,cpus-cpus0,syss-syss0,walls-walls0
        end if

        if (.not.associated(cur_form%next))
     &       call quit(1,'frm_sched2',
     &       'unexpected end of formula list ([END] missing)')

        cur_form => cur_form%next

      end do main_loop


      ifree = mem_flushmark()

      if (ntest.ge.100)
     &     write(luout,*) 'returning from frm_sched2'

      return
      end
