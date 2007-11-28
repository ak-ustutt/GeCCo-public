*----------------------------------------------------------------------*
      subroutine optc_macit(imacit,imicit,imicit_tot,
     &       task,iroute,
     &       ffopt,ffgrd,ffdia,
     &       ff_trv,ff_h_trv,
     &       nincore,lenbuf,ffscr,
     &       xbuf1,xbuf2,xbuf3,
     &       opti_info,opti_stat)
*----------------------------------------------------------------------*
*     driver for macro-iterations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_file_array.h'
      include 'def_optimize_info.h'
      include 'def_optimize_status.h'
      include 'ifc_memman.h'

      integer, intent(inout) ::
     &     task
      integer, intent(inout) ::
     &     imacit, imicit, imicit_tot
      integer, intent(in) ::
     &     iroute, nincore, lenbuf

      type(file_array), intent(in) ::
     &     ffopt(*), ffgrd(*), ffdia(*),
     &     ff_trv(*), ff_h_trv(*)
      type(filinf), intent(in) ::
     &     ffscr

      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout) ::
     &     opti_stat

      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

* local
      logical ::
     &     accept, shift, init
      integer ::
     &     irecr, irecv, klsmat,
     &     imet, idamp,
     &     ndim_save, ndel, iopt, lenscr, ifree
      real(8) ::
     &     xdum
      real(8), pointer ::
     &     xscr(:), xscr2(:), vec(:)
      integer, pointer ::
     &     ivec(:)

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

              call optc_diis_sbsp_add(opti_stat%ndim_rsbsp,
     &             opti_stat%ndim_vsbsp,opti_stat%mxdim_sbsp,
     &             init,
     &             opti_stat%iord_vsbsp, opti_stat%iord_rsbsp,
     &             ffopt(iopt)%fhand,ffgrd(iopt)%fhand,
     &                                      ffdia(iopt)%fhand,
     &             opti_stat%ffrsbsp(iopt)%fhand,
     &             opti_stat%ffvsbsp(iopt)%fhand,
     &             nincore,opti_info%nwfpar(iopt),lenbuf,xbuf1,xbuf2)

              shift = ndim_save.eq.opti_stat%ndim_rsbsp.and.iopt.eq.1
              call optc_update_redsp1(opti_stat%sbspmat,
     &             opti_stat%ndim_rsbsp,opti_stat%mxdim_sbsp,
     &             shift,init,
     &             opti_stat%iord_rsbsp,opti_stat%ffrsbsp(iopt)%fhand,
     &             nincore,opti_info%nwfpar(iopt),lenbuf,xbuf1,xbuf2)

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
     &           nincore,opti_info%nwfpar(1),lenbuf,xbuf1,xbuf2)

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
          call da_diavec(ffopt(iopt)%fhand,1,1d0,
     &               ffgrd(iopt)%fhand,1,-1d0,
     &               ffdia(iopt)%fhand,1,0d0,-1d0,
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
          call optc_expand_vec(vec,opti_stat%ndim_vsbsp,xdum,.false.,
     &         ffopt(iopt)%fhand,1,0d0,
     &         opti_stat%ffvsbsp(iopt)%fhand,opti_stat%iord_vsbsp,
     &         nincore,opti_info%nwfpar(iopt),lenbuf,xbuf1,xbuf2)
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
          call optc_expand_vec(vec,opti_stat%ndim_vsbsp,xdum,.false.,
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
          call optc_expand_vec(vec,opti_stat%ndim_rsbsp,xdum,.false.,
     &         ffgrd(1)%fhand,1,1d0,
     &       opti_stat%ffrsbsp(1)%fhand,opti_stat%iord_rsbsp,
     &       nincore,opti_info%nwfpar(1),lenbuf,xbuf1,xbuf2)
        end if

        ! make external step
        call optc_pert_step(ffopt(1)%fhand,
     &       ffgrd(1)%fhand,ffdia(1)%fhand,opti_stat%trrad,
     &       nincore,opti_info%nwfpar(1),lenbuf,xbuf1,xbuf2,xbuf3,ffscr)

      end if

      return
      end

