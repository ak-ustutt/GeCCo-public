*----------------------------------------------------------------------*
      subroutine optc_macit(imacit,imicit,imicit_tot,
     &       task,iroute,nopt,
     &       me_opt,me_grd,me_dia,
     &       me_trv,me_h_trv,
     &       n_states,
     &       me_special,nspecial,
c     &       ffopt,ffgrd,ffdia,ffmet,
c     &       ff_trv,ff_h_trv,
     &       nincore,lenbuf,ffscr,
     &       xbuf1,xbuf2,xbuf3,
     &       fspc,nspcfrm,energy,xngrd,
     &       opti_info,opti_stat_ini,
     &       orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     driver for macro-iterations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
c      include 'def_filinf.h'
c      include 'mdef_me_list.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_file_array.h'
      include 'def_optimize_info.h'
      include 'def_optimize_status.h'
      include 'ifc_memman.h'
      include 'def_orbinf.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'ifc_adv.h'

      integer, intent(inout) ::
     &     task
      integer, intent(inout) ::
     &     imacit, imicit, imicit_tot
      integer, intent(in) ::
     &     iroute, nincore, lenbuf, nopt, nspecial, nspcfrm, n_states

      type(me_list_array), intent(inout) ::
     &     me_opt(nopt), me_grd(nopt), me_dia(nopt),
     &     me_trv(nopt), me_h_trv(nopt),
     &     me_special(nspecial)
      type(filinf), intent(in) ::
     &     ffscr

      type(formula_item), intent(in) ::
     &     fspc(nspcfrm)

      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout), target ::
     &     opti_stat_ini

      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*), energy(0:*), xngrd(*)

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info

* local
      type(file_array), pointer ::
     &     ffopt(:), ffgrd(:), ffdia(:)
      logical ::
     &     accept, shift, init, normalize
      integer ::
     &     irecr, irecv, klsmat,
     &     imet, idamp,
     &     ndel, iopt, lenscr, ifree,
     &     nspcfrm_eff, i_spc, i_state, i_state_energy,
     &     nopt_state, iopt_state, idx
      integer, allocatable ::
     &     ndim_save(:)
      real(8) ::
     &     xnrm
      real(8), pointer ::
     &     xscr(:), xscr2(:), vec(:)
      integer, pointer ::
     &     ivec(:)
      type(optimize_status), pointer ::
     &     opti_stat
      logical ::
     &     lzero
      real(8), external :: dnrm2

      integer, external ::
     &     idx_mel_list
      integer :: iii
      
      ! set file arrays for calls to "old" routines
      allocate(ffopt(nopt),ffgrd(nopt),ffdia(nopt))
      do iopt = 1, nopt
        ffopt(iopt)%fhand => me_opt(iopt)%mel%fhand
        ffgrd(iopt)%fhand => me_grd(iopt)%mel%fhand
        ffdia(iopt)%fhand => me_dia(iopt)%mel%fhand
      end do

      allocate(ndim_save(n_states))
      ! for using the special formulas specifically for each state
      nspcfrm_eff = nspcfrm/n_states
      nopt_state = nopt/n_states
* check last step

      accept = .true.
      ! adapt trust radius
      ! ....

* accept step?
      if (accept) then
* step accepted:

        if (iroute.eq.1.or.iroute.eq.2) then

          ! for DIIS or precond. ASSJ: make D^-1|gradient(n)>
          ! (for DIIS, damping is allowed)
          if (iroute.eq.1) then

           opti_stat => opti_stat_ini
           do i_state = 1,n_states
            ndim_save(i_state) = opti_stat%ndim_rsbsp
            opti_stat => opti_stat%next_state
           end do

           opti_stat => opti_stat_ini
           i_state = 0
           iopt_state = nopt_state
           do iopt = 1, opti_info%nopt
            if(iopt_state.eq.nopt_state)then
             i_state = i_state+1
             iopt_state = 1
            else
             iopt_state = iopt_state + 1
            end if

c dbg
c            write(lulog,*) "MACIT, first part:"
c            write(lulog,*) " imacit    : ",imacit
c            write(lulog,*) " iopt      : ",iopt
c            write(lulog,*) " iopt_state: ",iopt_state
c            write(lulog,*) " i_state   : ",i_state
c dbg end            

            if(n_states.gt.1)then
             i_state_energy = i_state
            else
             i_state_energy = 0
            end if

            init = iopt_state.eq.1
c dbg
c              print *,'xbuf2 (2)',xbuf2(1:opti_info%nwfpar(iopt))
c dbg

c dbg
c$$$              if(imacit.gt.1)then
c$$$              write(lulog,*) 'MACIT 1: dump of result for ',
c$$$     &             trim( me_opt(iopt)%mel%label)
c$$$              call wrt_mel_file(lulog,5,
c$$$     &             me_opt(iopt)%mel,
c$$$     &             1,me_opt(iopt)%mel%op%n_occ_cls,
c$$$     &             str_info,orb_info)
c$$$              write(lulog,*) 'MACIT 1: dump of result for ',
c$$$     &             trim( me_grd(iopt)%mel%label)
c$$$              call wrt_mel_file(lulog,5,
c$$$     &             me_grd(iopt)%mel,
c$$$     &             1,me_grd(iopt)%mel%op%n_occ_cls,
c$$$     &             str_info,orb_info)
c$$$              write(lulog,*) 'MACIT 1: dump of result for ',
c$$$     &             trim( me_dia(iopt)%mel%label)
c$$$              call wrt_mel_file(lulog,5,
c$$$     &             me_dia(iopt)%mel,
c$$$     &             1,me_dia(iopt)%mel%op%n_occ_cls,
c$$$     &             str_info,orb_info)
c$$$              idx = idx_mel_list('ME_DENS',op_info)
c$$$              write(lulog,*) 'MACIT 1: dump of result for ',
c$$$     &             trim(op_info%mel_arr(idx)%mel%label)
c$$$              call wrt_mel_file(lulog,5,
c$$$     &             op_info%mel_arr(idx)%mel,
c$$$     &             1,op_info%mel_arr(idx)%mel%op%n_occ_cls,
c$$$     &             str_info,orb_info)
c$$$              idx = idx_mel_list('ME_C0',op_info)
c$$$              write(lulog,*) 'MACIT 1: dump of result for ',
c$$$     &             trim(op_info%mel_arr(idx)%mel%label)
c$$$              call wrt_mel_file(lulog,5,
c$$$     &             op_info%mel_arr(idx)%mel,
c$$$     &             1,op_info%mel_arr(idx)%mel%op%n_occ_cls,
c$$$     &             str_info,orb_info)
c$$$              do iii=1,nspecial
c$$$               write(lulog,*) 'MACIT 1: dump of result for ',
c$$$     &              trim( me_special(iii)%mel%label)
c$$$               call wrt_mel_file(lulog,5,
c$$$     &              me_special(iii)%mel,
c$$$     &              1,me_special(iii)%mel%op%n_occ_cls,
c$$$     &              str_info,orb_info)
c$$$              end do
c$$$              end if
c dbg

              call optc_diis_sbsp_add(opti_stat%ndim_rsbsp,
     &             opti_stat%ndim_vsbsp,opti_stat%mxdim_sbsp,
     &             init,
     &             opti_stat%iord_vsbsp, opti_stat%iord_rsbsp,
     &             me_opt(iopt)%mel,me_grd(iopt)%mel,me_dia(iopt)%mel,
     &             me_u,use_u,
     &             me_special,nspecial,
     &             opti_stat_ini%ffrsbsp(iopt)%fhand,
     &             opti_stat_ini%ffvsbsp(iopt)%fhand,
     &             opti_info%typ_prc(iopt),
     &             nincore,opti_info%nwfpar(iopt),
     &             lenbuf,xbuf1,xbuf2,xbuf3,
     &             fspc((i_state-1)*nspcfrm_eff+1:i_state*nspcfrm_eff),
     &             nspcfrm_eff,
     &             energy(i_state_energy),xngrd,iopt,imacit,i_state,
     &             opti_info,
     &             orb_info,op_info,str_info,strmap_info)

              shift = ndim_save(i_state).eq.opti_stat%ndim_rsbsp.and.
     &             iopt_state.eq.1

              call optc_update_redsp1(opti_stat%sbspmat,
     &             opti_stat%ndim_rsbsp,opti_stat%mxdim_sbsp,
     &             shift,init,
     &             opti_stat%iord_rsbsp,
     &             opti_stat_ini%ffrsbsp(iopt)%fhand,
     &             nincore,opti_info%nwfpar(iopt),
     &             lenbuf,xbuf1,xbuf2,xbuf3)

              if(n_states.gt.1.and.iopt_state.eq.nopt_state) then
               do i_spc = 1,nspecial
                call mel_adv_state(me_special(i_spc)%mel,n_states)
               end do
               idx = idx_mel_list('ME_DENS',op_info) ! quick & dirty
               call mel_adv_state(op_info%mel_arr(idx)%mel,n_states)
               idx = idx_mel_list('ME_C0',op_info) ! quick & dirty -> one can pass another "me_spacial", just to advance the states...
               call mel_adv_state(op_info%mel_arr(idx)%mel,n_states)
               call op_adv_state(["C0"],1,n_states,op_info,.false.)
               if(i_state.lt.n_states.and.iopt_state.eq.nopt_state)
     &              opti_stat => opti_stat%next_state
              end if

            end do

          else
           if(n_states.gt.1)
     &          call quit(1,'optc_macit',
     &          'this else part is not adapt for n_states>1')

            if (opti_info%nopt.gt.1)
     &           call quit(1,'optc_macit','adapt for nopt.gt.1')
            ! add |vec(n)> , on xbuf1 if incore, xbuf2 unused
            ndim_save(i_state) = opti_stat%ndim_vsbsp
            call optc_sbsp_add(opti_stat%ndim_vsbsp,
     &           opti_stat%mxdim_sbsp,opti_stat%iord_vsbsp,
     &           ffopt(1)%fhand,1,ffdia(1)%fhand,.false.,
     &           opti_stat%ffvsbsp(1)%fhand,
     &           nincore,opti_info%nwfpar(1),lenbuf,xbuf1,xbuf2)

            if (nincore.ge.2)
     &         call vec_from_da(ffopt(1)%fhand,1,
     &                          xbuf1,opti_info%nwfpar(1))
c dbg
c            print *,'amp: ',opti_info%nwfpar(1),
c     6              dnrm2(opti_info%nwfpar(1),xbuf1,1)
c dbg
            klsmat = opti_stat%mxdim_sbsp**2 + 1
            ! if incore: xbuf1 remains unchanged
            shift = ndim_save(i_state).eq.opti_stat%ndim_vsbsp
            call optc_update_redsp1(opti_stat%sbspmat(klsmat),
     &           opti_stat%ndim_vsbsp,
     &           opti_stat%mxdim_sbsp,shift,.true.,
     &           opti_stat%iord_vsbsp,opti_stat%ffvsbsp(1)%fhand,
     &           nincore,opti_info%nwfpar(1),
     &           lenbuf,xbuf1,xbuf2,xbuf3)

            ! add |gradient(n)>
            ! project out linear dependencies if required

            opti_stat => opti_stat_ini
            i_state = 0
            iopt_state = nopt_state
            do iopt = 1, opti_info%nopt
             if(iopt_state.eq.nopt_state)then
              i_state = i_state+1
              iopt_state = 1
             else
              iopt_state = iopt_state + 1
             end if

             if (opti_info%typ_prc(iopt).eq.optinf_prc_traf.and.
     &          opti_info%optref.ne.0.and.nspecial.ge.4) then

Cc dbg
C          print *,'call to projector' !HIERHIER
Cc dbg
C                call optc_project(me_grd(iopt)%mel,
C     &               me_opt(iopt)%mel,me_dia(iopt)%mel,
C     &            me_special,nspecial,
C     &            opti_info%nwfpar(iopt),xbuf2,
C     &            fspc,nspcfrm,iopt,opti_info,
C     &            orb_info,op_info,str_info,strmap_info)

              end if
            end do
            ! add |gradient(n)>, on xbuf2 if incore, xbuf1 unused
            ndim_save(i_state) = opti_stat%ndim_rsbsp
            call optc_sbsp_add(opti_stat%ndim_rsbsp, ! yaa, attention: maybe advance or something is necessary: ffgrd
     &           opti_stat%mxdim_sbsp,opti_stat%iord_rsbsp,
     &           ffgrd(1)%fhand,1,ffdia(1)%fhand,.false.,
     &           opti_stat_ini%ffrsbsp(1)%fhand,
     &           nincore,opti_info%nwfpar(1),lenbuf,xbuf2,xbuf1)

            if (nincore.ge.2)
     &         call vec_from_da(ffgrd(1)%fhand,1,
     &                          xbuf2,opti_info%nwfpar(1))
            shift = ndim_save(i_state).eq.opti_stat%ndim_rsbsp
            call optc_update_redsp2(
     &           opti_stat%sbspmat,opti_stat%ndim_vsbsp,
     &              opti_stat%ndim_rsbsp,opti_stat%mxdim_sbsp,shift,
     &           opti_stat%iord_vsbsp,opti_stat%ffvsbsp(1)%fhand,
     &           opti_stat%iord_rsbsp,opti_stat_ini%ffrsbsp(1)%fhand,
     &           nincore,opti_info%nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)

          end if

        end if

* step not accepted:
      else
      ! DIIS:
        call quit(1,'optc_macit', 'reject step: not yet ready')
      ! restore |vec(n-1)> from subspace info
      
      ! modify D^-1|gradient(n-1)> for larger damping

      ! modify |vec(n-1)> + D^-1|gradient(n-1)> for larger damping

      ! ASSJ: signal that subspace is unchanged

      ! 2nd order: line-searching would be appropriate now ...
      end if

      if (iroute.eq.0) then
        ! only for testing ...
        do iopt = 1, opti_info%nopt
          call da_diavec(ffopt(iopt)%fhand,1,1d0,
     &               ffgrd(iopt)%fhand,1,-1d0,
     &               ffdia(iopt)%fhand,1,0d0,-1d0,
     &               opti_info%nwfpar(iopt),xbuf1,xbuf2,lenbuf)
        end do

      else if (iroute.eq.1) then
* do DIIS extrapolation for current subspace
        if (opti_stat_ini%ndim_rsbsp.ne.opti_stat_ini%ndim_vsbsp)
     &       call quit(1,'optc_macit','inconsistency (DIIS)')

        opti_stat => opti_stat_ini
        i_state = 0
        iopt_state = nopt_state
        do iopt = 1, opti_info%nopt
         if(iopt_state.eq.nopt_state)then
          i_state = i_state+1
          iopt_state = 1
         else
          iopt_state = iopt_state + 1
         end if

c dbg
c         write(lulog,*) "MACIT, second part:"
c         write(lulog,*) " imacit    : ",imacit
c         write(lulog,*) " iopt      : ",iopt
c         write(lulog,*) " iopt_state: ",iopt_state
c         write(lulog,*) " i_state   : ",i_state
c dbg end

         if(iopt_state.eq.1)then
          lenscr = (opti_stat%ndim_rsbsp+1)*(opti_stat%ndim_rsbsp+2)/2
          ifree = mem_alloc_real(xscr,lenscr,'DIIS_mat')
          ifree = mem_alloc_real(vec,opti_stat%ndim_rsbsp+1,'DIIS_vec')
          ifree = mem_alloc_int(ivec,opti_stat%ndim_rsbsp+1,'DIIS_piv')
          call optc_diis_extr(opti_stat%sbspmat,
     &         vec,opti_stat%ndim_rsbsp,ndel,
     &         xscr,ivec)
         end if

         normalize = opti_info%typ_prc(iopt).eq.optinf_prc_norm
         call optc_expand_vec(vec,opti_stat%ndim_vsbsp,xnrm,normalize,
     &        ffopt(iopt)%fhand,
     &        ffopt(iopt)%fhand%current_record,0d0,
     &        opti_stat_ini%ffvsbsp(iopt)%fhand,opti_stat%iord_vsbsp,
     &        nincore,opti_info%nwfpar(iopt),lenbuf,xbuf1,xbuf2)
         if (normalize) then
          call dscal(opti_info%nwfpar(iopt),1d0/xnrm,xbuf1,1)
          call vec_to_da(ffopt(iopt)%fhand,
     &         ffopt(iopt)%fhand%current_record,xbuf1,
     &         opti_info%nwfpar(iopt))
         end if

         call touch_file_rec(ffopt(iopt)%fhand)
         if(iopt_state.eq.nopt_state)then
          if(n_states.GT.1)then
           do i_spc = 1,nspecial
            call mel_adv_state(me_special(i_spc)%mel,n_states)
           end do
           idx = idx_mel_list('ME_DENS',op_info) ! quick & dirty
           call mel_adv_state(op_info%mel_arr(idx)%mel,n_states)
           idx = idx_mel_list('ME_C0',op_info) ! quick & dirty
           call mel_adv_state(op_info%mel_arr(idx)%mel,n_states)
           call op_adv_state(["C0"],1,n_states,op_info,.false.)
           if(i_state.lt.n_states)
     &          opti_stat => opti_stat%next_state
          end if
          ifree = mem_dealloc('DIIS_mat')
          ifree = mem_dealloc('DIIS_vec')
          ifree = mem_dealloc('DIIS_piv')
         end if

        end do

        ! project out linear dependencies if required
        opti_stat => opti_stat_ini
        i_state = 0
        iopt_state = nopt_state
        do iopt = 1, opti_info%nopt
         if(iopt_state.eq.nopt_state)then
          i_state = i_state+1
          iopt_state = 1
         else
          iopt_state = iopt_state + 1
         end if
         if ((opti_info%typ_prc(iopt).eq.optinf_prc_traf.or.
     &        opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc.or.
     &        opti_info%typ_prc(iopt).eq.optinf_prc_invH0 ).and.
     &        opti_info%optref.ne.0.and.nspecial.ge.5) then
              lzero =opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc
              call optc_project(me_opt(iopt)%mel,me_opt(iopt)%mel,
     &        me_dia(iopt)%mel,me_special,nspecial,
     &        opti_info%nwfpar(iopt),xbuf1,
     &        fspc((i_state-1)*nspcfrm_eff+1:i_state*nspcfrm_eff),
     &        nspcfrm_eff,iopt,imacit,i_state,lzero,
     &        opti_info,orb_info,op_info,str_info,strmap_info)
         end if
         if(n_states.gt.1.and.iopt_state.eq.nopt_state)then
          do i_spc = 1,nspecial
           call mel_adv_state(me_special(i_spc)%mel,n_states)
          end do
          idx = idx_mel_list('ME_DENS',op_info) ! quick & dirty
          call mel_adv_state(op_info%mel_arr(idx)%mel,n_states)
          idx = idx_mel_list('ME_C0',op_info) ! quick & dirty
          call mel_adv_state(op_info%mel_arr(idx)%mel,n_states)
          call op_adv_state(["C0"],1,n_states,op_info,.false.)
         end if
        end do

c dbg
c$$$        do iopt=1,opti_info%nopt
c$$$         if(imacit.gt.1)then
c$$$          write(lulog,*) 'MACIT final: dump of result for ',
c$$$     &         trim( me_opt(iopt)%mel%label)
c$$$          call wrt_mel_file(lulog,5,
c$$$     &         me_opt(iopt)%mel,
c$$$     &         1,me_opt(iopt)%mel%op%n_occ_cls,
c$$$     &         str_info,orb_info)
c$$$          write(lulog,*) 'MACIT final: dump of result for ',
c$$$     &         trim( me_grd(iopt)%mel%label)
c$$$          call wrt_mel_file(lulog,5,
c$$$     &         me_grd(iopt)%mel,
c$$$     &         1,me_grd(iopt)%mel%op%n_occ_cls,
c$$$     &         str_info,orb_info)
c$$$          write(lulog,*) 'MACIT final: dump of result for ',
c$$$     &         trim( me_dia(iopt)%mel%label)
c$$$          call wrt_mel_file(lulog,5,
c$$$     &         me_dia(iopt)%mel,
c$$$     &         1,me_dia(iopt)%mel%op%n_occ_cls,
c$$$     &         str_info,orb_info)
c$$$          idx = idx_mel_list('ME_DENS',op_info)
c$$$          write(lulog,*) 'MACIT final: dump of result for ',
c$$$     &         trim(op_info%mel_arr(idx)%mel%label)
c$$$          call wrt_mel_file(lulog,5,
c$$$     &         op_info%mel_arr(idx)%mel,
c$$$     &         1,op_info%mel_arr(idx)%mel%op%n_occ_cls,
c$$$     &         str_info,orb_info)
c$$$          idx = idx_mel_list('ME_C0',op_info)
c$$$          write(lulog,*) 'MACIT final: dump of result for ',
c$$$     &         trim(op_info%mel_arr(idx)%mel%label)
c$$$          call wrt_mel_file(lulog,5,
c$$$     &         op_info%mel_arr(idx)%mel,
c$$$     &         1,op_info%mel_arr(idx)%mel%op%n_occ_cls,
c$$$     &         str_info,orb_info)
c$$$          do iii=1,nspecial
c$$$           write(lulog,*) 'MACIT final: dump of result for ',
c$$$     &          trim( me_special(iii)%mel%label)
c$$$           call wrt_mel_file(lulog,5,
c$$$     &          me_special(iii)%mel,
c$$$     &          1,me_special(iii)%mel%op%n_occ_cls,
c$$$     &          str_info,orb_info)
c$$$          end do
c$$$         end if
c$$$        end do
c dbg

* do ASSJ step for current subspace
      else if (iroute.eq.2) then

       if(n_states.gt.1)
     &      call quit(1,'optc_macit','ASSJ: adapt for n_states>1')

        if (opti_stat%ndim_rsbsp.ne.opti_stat%ndim_vsbsp)
     &       call quit(1,'optc_macit','inconsistency (ASSJ)')

        lenscr = opti_stat%mxdim_sbsp**2
        ifree = mem_alloc_real(xscr,lenscr,'ASSJ_mat1')
        ifree = mem_alloc_real(xscr2,lenscr,'ASSJ_mat2')
        ifree = mem_alloc_real(vec,opti_stat%mxdim_sbsp,'DIIS_vec')

        imet = 0
        idamp = 1
        klsmat = opti_stat%mxdim_sbsp**2 + 1
        call optc_assj_step(imet,
     &       opti_stat%sbspmat,opti_stat%sbspmat(klsmat),
     &       vec,opti_stat%ndim_vsbsp,opti_stat%mxdim_sbsp,ndel,
     &       idamp,opti_stat%trrad,xscr,xscr2)

        
        if (opti_info%nopt.gt.1)
     &       call quit(1,'optc_macit','ASSJ: adapt for nopt>1')
        ! add internal step to vector
        if (opti_stat%ndim_vsbsp.gt.1) then
          ! to be optimized later ...
          if (nincore.ge.2)
     &         call vec_from_da(ffopt(1)%fhand,1,
     &                          xbuf1,opti_info%nwfpar(1))
          call optc_expand_vec(vec,opti_stat%ndim_vsbsp,xnrm,.false.,
     &         ffopt(1)%fhand,1,1d0,
     &       opti_stat%ffvsbsp(1)%fhand,opti_stat%iord_vsbsp,
     &       nincore,opti_info%nwfpar(1),lenbuf,xbuf1,xbuf2)
        end if
        ! add correction to gradient

        if (opti_stat%ndim_rsbsp.gt.1) then
          ! to be optimized later ...
          if (nincore.ge.2)
     &         call vec_from_da(ffgrd(1)%fhand,1,
     &                          xbuf1,opti_info%nwfpar(1))
          call optc_expand_vec(vec,opti_stat%ndim_rsbsp,xnrm,.false.,
     &         ffgrd(1)%fhand,1,1d0,
     &       opti_stat_ini%ffrsbsp(1)%fhand,opti_stat%iord_rsbsp,
     &       nincore,opti_info%nwfpar(1),lenbuf,xbuf1,xbuf2)
        end if

        ! make external step
        opti_stat => opti_stat_ini
        i_state = 0
        iopt_state = nopt_state
        do iopt = 1, opti_info%nopt
         if(iopt_state.eq.nopt_state)then
          i_state = i_state+1
          iopt_state = 1
         else
          iopt_state = iopt_state + 1
         end if

          if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then

            if (nincore.lt.2) 
     &              call quit(1,'optc_macit','need more memory')

            call optc_prc_traf(
     &              me_opt(iopt)%mel,me_grd(iopt)%mel,me_dia(iopt)%mel,
     &              me_special,nspecial,
     &              opti_info%nwfpar(iopt),xbuf1,xbuf2,
     &              fspc((i_state-1)*nspcfrm_eff+1:i_state*nspcfrm_eff),
     &              nspcfrm_eff,
     &              xngrd,iopt,imacit,i_state,opti_info,
     &              orb_info,op_info,str_info,strmap_info)

            if(n_states.gt.1.and.iopt_state.eq.nopt_state)then
             do i_spc = 1,nspecial
              call mel_adv_state(me_special(i_spc)%mel,n_states)
             end do
            end if
          else
            call optc_pert_step(ffopt(iopt)%fhand,
     &         ffgrd(iopt)%fhand,ffdia(iopt)%fhand,opti_stat%trrad,
     &         nincore,opti_info%nwfpar(iopt),lenbuf,xbuf1,xbuf2,xbuf3,
     &         ffscr)
          end if
        end do


      end if


      deallocate(ffopt,ffdia,ffgrd)

      return
      end

