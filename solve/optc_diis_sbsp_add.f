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

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nspecial
      logical, intent(in) ::
     &     get_new_rec
      type(me_list_array), intent(in) ::
     &     me_special(nspecial)
      type(me_list), intent(in) ::
     &     me_amp,me_grd,me_dia
      type(filinf), intent(in) ::
     &     ff_rsbsp,ff_vsbsp
      integer, intent(inout) ::
     &     ndim_rsbsp,ndim_vsbsp,iord_vsbsp,iord_rsbsp
      integer, intent(in) ::
     &     mxdim_sbsp,
     &     nincore, nwfpar, 
     &     lenbuf, typ_prc
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(strmapinf),intent(in) ::
     &     strmap_info
      type(operator_info), intent(inout) ::
     &     op_info


      integer ::
     &     irecr, irecv, inum, idx_inv, idx, len, iblk

      type(filinf), pointer ::
     &     ffamp, ffgrd, ffdia

      integer, external ::
     &     ioptc_get_sbsp_rec, idx_oplist2
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
        write(luout,*) 'added records: ',irecv, irecr
        write(luout,*) 'nwfpar: ',nwfpar
      end if

      if (nincore.ge.2) then

        call vec_from_da(ffgrd,1,xbuf1,nwfpar)

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
        case(optinf_prc_file)
c dbg
c          print *,'Raw grad:'
c          call wrt_mel_buf(luout,5,xbuf1,me_grd,1,2,
c     &     str_info,orb_info)
c dbg
          call vec_from_da(ffdia,1,xbuf2,nwfpar)
          xbuf1(1:nwfpar) = xbuf1(1:nwfpar)/xbuf2(1:nwfpar)
c dbg
c          print *,'Precond. grad:'
c          call wrt_mel_buf(luout,5,xbuf1,me_grd,1,2,
c     &         str_info,orb_info)
c dbg
        case(optinf_prc_blocked)
c dbg
c          print *,'Raw grad:'
c          call wrt_mel_buf(luout,5,xbuf1,me_grd,1,1,
c     &     str_info,orb_info)
c dbg
c          call optc_prc_special(me_grd,me_special,nspecial,
c     &                          nincore,xbuf1,xbuf2,xbuf3,lenbuf,
c     &                          orb_info,str_info)
c          xbuf1(1:nwfpar) = xbuf3(1:nwfpar)
          call optc_prc_special2(me_grd,me_special,nspecial,
     &                           me_amp%op%name,0d0,
     &                          nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                          orb_info,op_info,str_info,strmap_info)
c          call mem_check('after prc_special2')
c dbg
c          print *,'Precond. grad:'
c          call wrt_mel_buf(luout,5,xbuf1,me_grd,1,1,
c     &         str_info,orb_info)
c dbg
        end select
c dbg
c          print *,'|g/d|:' ,dnrm2(nwfpar,xbuf1,1)
c          print *,'g/d: ', xbuf1(1:nwfpar)
c dbg
        call vec_from_da(ffamp,1,xbuf2,nwfpar)
c dbg
c          print *,'t norm:',dnrm2(nwfpar,xbuf2,1)
c          print *,'t before: ', xbuf2(1:nwfpar)
c dbg
        xbuf2(1:nwfpar) = xbuf2(1:nwfpar) - xbuf1(1:nwfpar)
c dbg
c          print *,'t norm new:',dnrm2(nwfpar,xbuf2,1)
c          print *,'t before new: ', xbuf2(1:nwfpar)
c dbg

c        else
cc dbg
c          print *,'prc name: ',trim(ffdia%name)
c          print *,'g,d norm:',
c     &         dnrm2(nwfpar,xbuf1,1),dnrm2(nwfpar,xbuf2,1)
cc dbg
c
c          idx_inv = idx_oplist2(op_b_inv,op_info)
c          op => op_info%op_arr(idx_inv)%op
cc dbg
c          print *,'g norm:',dnrm2(nwfpar,xbuf1,1)
c          print *,'g: ',xbuf1(1:nwfpar)
cc          print *,'g: ',xbuf1(1:20)
cc          print *,'g: ',xbuf1(72:92)          
cc dbg
c          call op_vec_mult(1d0,ffdia,op,xbuf1,nwfpar,
c     &         op_info,orb_info)
cc dbg
c          print *,'g/d norm:',dnrm2(nwfpar,xbuf1,1)
c          print *,'g/d: ',xbuf1(1:nwfpar)
c          print *,'g/d: ',xbuf1(1:20)
c          print *,'g/d: ',xbuf1(72:92)
cc dbg
c          call vec_from_da(ffamp,1,xbuf2,nwfpar)
cc dbg
c          print *,'t norm:',dnrm2(nwfpar,xbuf2,1)
c          print *,'t after: ',xbuf2(1:nwfpar)
cc          print *,'t: ',xbuf2(72:92)
cc dbg
c          xbuf2(1:nwfpar) = xbuf2(1:nwfpar) - xbuf1(1:nwfpar)
c
c        endif

        call vec_to_da(ff_rsbsp,irecr,xbuf1,nwfpar)
        call vec_to_da(ff_vsbsp,irecv,xbuf2,nwfpar)

      else
        
        if (typ_prc.eq.optinf_prc_blocked)
     &        call quit(1,'optc_diis_sbsp_add',
     &       '(2): blocked preconditioning for nincore==3, only')

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

      return
      end
