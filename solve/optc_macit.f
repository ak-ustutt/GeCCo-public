*----------------------------------------------------------------------*
      subroutine optc_macit(imacit,imicit,imicit_tot,
     &       task,iroute,nopt,
     &       me_opt,me_grd,me_dia,
     &       me_trv,me_h_trv,
     &       me_special,nspecial,
c     &       ffopt,ffgrd,ffdia,ffmet,
c     &       ff_trv,ff_h_trv,
     &       nincore,lenbuf,ffscr,
     &       xbuf1,xbuf2,xbuf3,
     &       flist,depend,energy,
     &       opti_info,opti_stat,
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
      include 'def_dependency_info.h'

      integer, intent(inout) ::
     &     task
      integer, intent(inout) ::
     &     imacit, imicit, imicit_tot
      integer, intent(in) ::
     &     iroute, nincore, lenbuf, nopt, nspecial

      type(me_list_array), intent(in) ::
     &     me_opt(nopt), me_grd(nopt), me_dia(nopt),
     &     me_trv(nopt), me_h_trv(nopt),
     &     me_special(nspecial)
      type(filinf), intent(in) ::
     &     ffscr

      type(formula_item), intent(inout) ::
     &     flist
      type(dependency_info) ::
     &     depend

      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout) ::
     &     opti_stat

      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*), energy

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
     &     ndim_save, ndel, iopt, lenscr, ifree, nselect, idx
      real(8) ::
     &     xnrm
      real(8), pointer ::
     &     xscr(:), xscr2(:), vec(:), xret(:)
      integer, pointer ::
     &     ivec(:), idxselect(:)
      integer, external ::
     &     idx_mel_list

      ! set file arrays for calls to "old" routines
      allocate(ffopt(nopt),ffgrd(nopt),ffdia(nopt))
      do iopt = 1, nopt
        ffopt(iopt)%fhand => me_opt(iopt)%mel%fhand
        ffgrd(iopt)%fhand => me_grd(iopt)%mel%fhand
        ffdia(iopt)%fhand => me_dia(iopt)%mel%fhand
      end do

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

            ndim_save = opti_stat%ndim_rsbsp
            
            do iopt = 1, opti_info%nopt

              init = iopt.eq.1
c dbg
c              print *,'xbuf2 (2)',xbuf2(1:opti_info%nwfpar(iopt))
c dbg
              call optc_diis_sbsp_add(opti_stat%ndim_rsbsp,
     &             opti_stat%ndim_vsbsp,opti_stat%mxdim_sbsp,
     &             init,
     &             opti_stat%iord_vsbsp, opti_stat%iord_rsbsp,
     &             me_opt(iopt)%mel,me_grd(iopt)%mel,me_dia(iopt)%mel,
     &             me_special,nspecial,
     &             opti_stat%ffrsbsp(iopt)%fhand,
     &             opti_stat%ffvsbsp(iopt)%fhand,
     &             opti_info%typ_prc(iopt),
     &             nincore,opti_info%nwfpar(iopt),
     &             lenbuf,xbuf1,xbuf2,xbuf3,
     &             flist,depend,energy,iopt,opti_info,
     &             orb_info,op_info,str_info,strmap_info)

              shift = ndim_save.eq.opti_stat%ndim_rsbsp.and.iopt.eq.1
              call optc_update_redsp1(opti_stat%sbspmat,
     &             opti_stat%ndim_rsbsp,opti_stat%mxdim_sbsp,
     &             shift,init,
     &             opti_stat%iord_rsbsp,opti_stat%ffrsbsp(iopt)%fhand,
     &             nincore,opti_info%nwfpar(iopt),
     &             lenbuf,xbuf1,xbuf2,xbuf3)


            end do

          else
            if (opti_info%nopt.gt.1)
     &           call quit(1,'optc_macit','adapt for nopt.gt.1')
            ! add |vec(n)> , on xbuf1 if incore, xbuf2 unused
            ndim_save = opti_stat%ndim_vsbsp
            call optc_sbsp_add(opti_stat%ndim_vsbsp,
     &           opti_stat%mxdim_sbsp,opti_stat%iord_vsbsp,
     &           ffopt(1)%fhand,1,ffdia(1)%fhand,.false.,
     &           opti_stat%ffvsbsp(1)%fhand,
     &           nincore,opti_info%nwfpar(1),lenbuf,xbuf1,xbuf2)

            klsmat = opti_stat%mxdim_sbsp**2 + 1
            ! if incore: xbuf1 remains unchanged
            shift = ndim_save.eq.opti_stat%ndim_vsbsp
            call optc_update_redsp1(opti_stat%sbspmat(klsmat),
     &           opti_stat%ndim_vsbsp,
     &           opti_stat%mxdim_sbsp,shift,.true.,
     &           opti_stat%iord_vsbsp,opti_stat%ffvsbsp(1)%fhand,
     &           nincore,opti_info%nwfpar(1),
     &           lenbuf,xbuf1,xbuf2,xbuf3)

            ! add |gradient(n)>, on xbuf2 if incore, xbuf1 unused
            ndim_save = opti_stat%ndim_rsbsp
            call optc_sbsp_add(opti_stat%ndim_rsbsp,
     &           opti_stat%mxdim_sbsp,opti_stat%iord_rsbsp,
     &           ffgrd(1)%fhand,1,ffdia(1)%fhand,.false.,
     &           opti_stat%ffrsbsp(1)%fhand,
     &           nincore,opti_info%nwfpar(1),lenbuf,xbuf2,xbuf1)

            shift = ndim_save.eq.opti_stat%ndim_rsbsp
            call optc_update_redsp2(
     &           opti_stat%sbspmat,opti_stat%ndim_vsbsp,
     &              opti_stat%ndim_rsbsp,opti_stat%mxdim_sbsp,shift,
     &           opti_stat%iord_vsbsp,opti_stat%ffvsbsp(1)%fhand,
     &           opti_stat%iord_rsbsp,opti_stat%ffrsbsp(1)%fhand,
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
          call da_diavec(ffopt(iopt)%fhand,1,1,1d0,
     &               ffgrd(iopt)%fhand,1,1,-1d0,
     &               ffdia(iopt)%fhand,1,1,0d0,-1d0,
     &               opti_info%nwfpar(iopt),xbuf1,xbuf2,lenbuf)
        end do

      else if (iroute.eq.1) then
* do DIIS extrapolation for current subspace
        if (opti_stat%ndim_rsbsp.ne.opti_stat%ndim_vsbsp)
     &       call quit(1,'optc_macit','inconsistency (DIIS)')

        lenscr = (opti_stat%ndim_rsbsp+1)*(opti_stat%ndim_rsbsp+2)/2
        ifree = mem_alloc_real(xscr,lenscr,'DIIS_mat')
        ifree = mem_alloc_real(vec,opti_stat%ndim_rsbsp+1,'DIIS_vec')
        ifree = mem_alloc_int(ivec,opti_stat%ndim_rsbsp+1,'DIIS_piv')
        call optc_diis_extr(opti_stat%sbspmat,
     &       vec,opti_stat%ndim_rsbsp,ndel,
     &       xscr,ivec)

        do iopt = 1, opti_info%nopt
          normalize = opti_info%typ_prc(iopt).eq.optinf_prc_norm
          call optc_expand_vec(vec,opti_stat%ndim_vsbsp,xnrm,normalize,
     &         ffopt(iopt)%fhand,1,0d0,
     &         opti_stat%ffvsbsp(iopt)%fhand,opti_stat%iord_vsbsp,
     &         nincore,opti_info%nwfpar(iopt),lenbuf,xbuf1,xbuf2)
          if (normalize) then
            call dscal(opti_info%nwfpar(iopt),1d0/xnrm,xbuf1,1)
            call vec_to_da(ffopt(iopt)%fhand,1,xbuf1,
     &                     opti_info%nwfpar(iopt))
          end if
          call touch_file_rec(ffopt(iopt)%fhand)
        end do

        ! project out linear dependencies if required
        do iopt = 1, opti_info%nopt
          if (opti_info%typ_prc(iopt).eq.optinf_prc_traf.and.
     &        nspecial.ge.6) then

cmh   comment this in, if you wish to use optimization algorithm "1a"
c          ! update metric, trafo matrices and projector if not up to date
c          idx = idx_mel_list('ME_C0',op_info)  ! quick & dirty
c          if (nspecial.ge.6.and.
c     &        op_info%mel_arr(idx)%mel%fhand%last_mod(1).gt.
c     &        me_special(2)%mel%fhand%last_mod(1))
c     &        call update_metric(me_dia(iopt)%mel,me_special,nspecial,
c     &          flist,depend,orb_info,op_info,str_info,strmap_info)

          ! put new vector to special list for transformation
          call list_copy(me_opt(iopt)%mel,me_special(1)%mel)

          ! use projection matrix
          call assign_me_list(me_special(4)%mel%label,
     &                           me_special(2)%mel%op%name,op_info)

          ! calculate transformed vector
          allocate(xret(depend%ntargets),idxselect(depend%ntargets))
          nselect = 0
          call select_formula_target(idxselect,nselect,
     &                me_opt(iopt)%mel%label,depend,op_info)
          ! pretend that vector is not up to date
          ffopt(iopt)%fhand%last_mod(
     &           ffopt(iopt)%fhand%current_record) = -1
          call frm_sched(xret,flist,depend,idxselect,nselect,
     &                op_info,str_info,strmap_info,orb_info)
          deallocate(xret,idxselect)

          end if
        end do

* do ASSJ step for current subspace
      else if (iroute.eq.2) then

        if (opti_stat%ndim_rsbsp.ne.opti_stat%ndim_vsbsp)
     &       call quit(1,'optc_macit','inconsistency (ASSJ)')

        lenscr = opti_stat%mxdim_sbsp**2
        ifree = mem_alloc_real(xscr,lenscr,'ASSJ_mat1')
        ifree = mem_alloc_real(xscr2,lenscr,'ASSJ_mat2')
        ifree = mem_alloc_real(vec,opti_stat%mxdim_sbsp,'DIIS_vec')

        imet = 0
        idamp = 1
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
     &       opti_stat%ffrsbsp(1)%fhand,opti_stat%iord_rsbsp,
     &       nincore,opti_info%nwfpar(1),lenbuf,xbuf1,xbuf2)
        end if

        ! make external step
        call optc_pert_step(ffopt(1)%fhand,
     &       ffgrd(1)%fhand,ffdia(1)%fhand,opti_stat%trrad,
     &       nincore,opti_info%nwfpar(1),lenbuf,xbuf1,xbuf2,xbuf3,ffscr)

      end if

      deallocate(ffopt,ffdia,ffgrd)

      return
      end

