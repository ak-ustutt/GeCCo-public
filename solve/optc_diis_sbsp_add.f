*----------------------------------------------------------------------*
      subroutine optc_diis_sbsp_add(ndim_rsbsp,ndim_vsbsp,mxdim_sbsp,
     &     get_new_rec,
     &     iord_vsbsp,iord_rsbsp,
     &     me_amp,me_grd,me_dia,
     &     me_special,nspecial,
     &     ff_rsbsp,ff_vsbsp,
     &     typ_prc,
     &     nincore,nwfpar,
     &     lenbuf,xbuf1,xbuf2,xbuf3,
     &     fspc,nspcfrm,energy,xngrd,iopt,imacit,i_state,opti_info,
     &     orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
*     add (preconditioned) gradient and sum of current vector and
*      minus the prec. gradient to the respective subspace files
*
*     incore.ge.2: on output is the precond. gradient on xbuf1
*                      and vector - precond. gradient on xbuf2
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_optimize_info.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nspecial, iopt, nspcfrm, imacit,i_state
      logical, intent(in) ::
     &     get_new_rec
      type(me_list_array), intent(inout) ::
     &     me_special(nspecial)
      type(me_list), intent(in) ::
     &     me_amp,me_grd,me_dia
      type(filinf), intent(in) ::
     &     ff_rsbsp,ff_vsbsp
      integer, intent(inout) ::
     &     ndim_rsbsp,ndim_vsbsp,
     &     iord_vsbsp(mxdim_sbsp),iord_rsbsp(mxdim_sbsp)
      integer, intent(in) ::
     &     mxdim_sbsp,
     &     nincore, nwfpar,
     &     lenbuf, typ_prc
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*), energy, xngrd(*)

      type(formula_item), intent(in) ::
     &     fspc(nspcfrm)

      type(optimize_info), intent(in) ::
     &     opti_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(strmapinf),intent(in) ::
     &     strmap_info
      type(operator_info), intent(inout) ::
     &     op_info

      logical::
     &     lzero_flag
      integer ::
     &     irecr, irecv, inum, idx_inv
      character(len_opname) ::
     &     op_grd_name, op_trf_name, op_amp_name

      real(8) ::
     &     xdum

      type(filinf), pointer ::
     &     ffamp, ffgrd, ffdia

      integer, external ::
     &     ioptc_get_sbsp_rec, idx_oplist2, idx_mel_list
      real(8), external ::
     &     dnrm2

      ! pointers to file handle
      ffamp => me_amp%fhand
      ffgrd => me_grd%fhand
      ffdia => me_dia%fhand

      ! get record-number for new vector in subspace
      inum = ndim_rsbsp   ! get last entry
      if (get_new_rec) inum = 0 ! request new entry
      irecr = ioptc_get_sbsp_rec(inum,iord_rsbsp,ndim_rsbsp,mxdim_sbsp)
      inum = ndim_vsbsp   ! get last entry
      if (get_new_rec) inum = 0 ! request new entry
      irecv = ioptc_get_sbsp_rec(inum,iord_vsbsp,ndim_vsbsp,mxdim_sbsp)

      if (ntest.ge.100) then
        write(lulog,*) 'this is optc_diis_sbsp_add.f'
        write(lulog,*) 'added records: ',irecv, irecr
        write(lulog,*) 'nwfpar: ',nwfpar
        write(lulog,*) 'typ_prc',typ_prc
        write(lulog,*) 'nincore',nincore
        write(lulog,*) 'set T1 zero',lzero_flag
      end if

      if (nincore.ge.2) then

!        select case(typ_prc) ...
c        if(trim(ffdia%name).ne.'op_B_INV_elements.da')then
c dbg
c          print *,'prc name: ',trim(ffdia%name)
c          print *,'|g|,|d|:',dnrm2(nwfpar,xbuf1,1),dnrm2(nwfpar,xbuf2,1)
c          print *,'g: ', xbuf1(1:nwfpar)
c          print *,'d: ', xbuf2(1:nwfpar)
c          xbuf2(1:nwfpar) = 0d0
c dbg

        select case(typ_prc)
        ! prelim. w/o damping
        case(optinf_prc_file,optinf_prc_norm)
          call vec_from_da(ffgrd,1,xbuf1,nwfpar)
          if (ntest.ge.100) then
            write(lulog,*) 'aaaaaa gradient vector before:'
            write(lulog,*) xbuf1(1:nwfpar)
          end if

          call vec_from_da(ffdia,1,xbuf2,nwfpar)

          ! preconditioner for eigenvalue equation: substract energy
          if (typ_prc.eq.optinf_prc_norm)
     &       xbuf2(1:nwfpar) = xbuf2(1:nwfpar) - energy

          call diavc(xbuf1,xbuf1,1d0,xbuf2,0d0,nwfpar)

          if (ntest.ge.100) then
            write(lulog,*) 'gradient vector afterwards:'
            write(lulog,*) xbuf1(1:nwfpar)
          end if
          call vec_from_da(ffamp,ffamp%current_record,xbuf2,nwfpar)
        ! the preconditioning by solving a set of LEq's was done before,
        ! so only load the (preconditioned) gradient here
        case(optinf_prc_invH0)
          call vec_from_da(ffgrd,1,xbuf1,nwfpar)
          call vec_from_da(ffamp,ffamp%current_record,xbuf2,nwfpar)
c dbg
c          print *,'read gradient (already prec): ',dnrm2(nwfpar,xbuf1,1)
c          print *,'read amplitude: ',dnrm2(nwfpar,xbuf2,1)
c dbg

        case(optinf_prc_blocked)
          call vec_from_da(ffgrd,1,xbuf1,nwfpar)

          call optc_prc_special2(me_grd,me_special,nspecial,
     &                           me_amp%op%name,0d0,
     &                          nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                          orb_info,op_info,str_info,strmap_info)
          call vec_from_da(ffamp,1,xbuf2,nwfpar)
        case(optinf_prc_mixed)
          call vec_from_da(ffgrd,1,xbuf1,nwfpar)
          call vec_from_da(ffdia,1,xbuf2,nwfpar)
          call optc_prc_mixed(me_grd,me_special,nspecial,
     &                           me_amp%op%name,0d0,
     &                          nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                          orb_info,op_info,str_info,strmap_info)
          call vec_from_da(ffamp,1,xbuf2,nwfpar)

        case(optinf_prc_traf,optinf_prc_traf_spc)
           lzero_flag=.false.
           if (typ_prc .eq. optinf_prc_traf_spc) lzero_flag=.true.

           call optc_prc_traf(me_amp,me_grd,me_dia,
     &                       me_special,nspecial,
     &                       nwfpar,xbuf1,xbuf2,
     &                       fspc,nspcfrm,xngrd,iopt,imacit,i_state,
     &                       opti_info,
     &                       orb_info,op_info,str_info,strmap_info,
     &                       lzero_flag)

        case default
          call quit(1,'optc_diis_sbsp_add','unknown route')
        end select

        xbuf2(1:nwfpar) = xbuf2(1:nwfpar) - xbuf1(1:nwfpar)

        call vec_to_da(ff_rsbsp,irecr,xbuf1,nwfpar)
        call vec_to_da(ff_vsbsp,irecv,xbuf2,nwfpar)

      else

        call quit(1,'optc_diis_sbsp_add','nincore<3 route not debugged')

        if (typ_prc.eq.optinf_prc_blocked)
     &        call quit(1,'optc_diis_sbsp_add',
     &       '(2): blocked preconditioning for nincore==3, only')
        if (typ_prc.eq.optinf_prc_mixed)
     &        call quit(1,'optc_diis_sbsp_add',
     &       '(2b): blocked preconditioning for nincore==3, only')

        ! prelim. w/o damping
        ! add D^-1|gradient(n)> to subspace
        call da_diavec(ff_rsbsp,irecr,0d0,
     &       ffgrd,1,1d0,
     &       ffdia,1,0d0,-1d0,
     &       nwfpar,xbuf1,xbuf2,lenbuf)

        ! add |vec(n)> - D^-1|gradient(n)> to subspace
        call da_vecsum(ff_vsbsp,irecv,
     &       ff_rsbsp,irecr,-1d0,
     &       ffamp,1,1d0,
     &       nwfpar,xbuf1,xbuf2,lenbuf)

      end if
c dbg      write (lulog,*) "optc_diis_sbsp_add ended"
      return
      end
