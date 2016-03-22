*----------------------------------------------------------------------*
      subroutine optc_prc_traf(me_amp,me_grd,me_dia,me_special,nspecial,
     &     nwfpar,xbuf1,xbuf2,
     &     fspc,nspcfrm,xngrd,iopt,imacit,i_state,opti_info,
     &     orb_info,op_info,str_info,strmap_info,
     &     lzero_flag)
*----------------------------------------------------------------------*
*
*     preconditioning step when transformations are involved
*     fspc(1)         : transformation formula
*     me_special(1)   : transformed vector or vector to be transformed
*     me_special(2,3) : transformation matrices
*     higher fspc, me_special: for update of metric (if needed)
*     xbuf1 : contains new step on exit
*     xbuf2 : contains old (current) guess vector on exit
*
*     matthias, Nov. 2010
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
      include 'ifc_input.h'
      include 'mdef_target_info.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nspecial, iopt, nspcfrm, nwfpar, imacit, i_state
      type(me_list_array), intent(inout) ::
     &     me_special(nspecial)
      type(me_list), intent(in) ::
     &     me_amp,me_grd,me_dia
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xngrd(*)

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
      logical,intent(in) ::
     &     lzero_flag


      integer ::
     &     idx, nopt, prc_iter
      character(len_opname) ::
     &     op_grd_name, op_trf_name, op_amp_name

      real(8) ::
     &     xdum, prc_impfac

      type(filinf), pointer ::
     &     ffamp, ffgrd, ffdia

      integer, external ::
     &     idx_mel_list

      character(len_target_name) ::
     &     c_st
      character(len_target_name), external ::
     &     state_label
      integer :: iii
      real(8),external::
     & xnormop


      ! pointers to file handle
      ffamp => me_amp%fhand
      ffgrd => me_grd%fhand
      ffdia => me_dia%fhand

      nopt = opti_info%nopt

      ! Quick & dirty over quick & dirty:
      ! The list for the reference is ME_C0_<state> for opt_ref=-1,-2
      if(opti_info%optref.EQ.-1.or.opti_info%optref.EQ.-2)then
       c_st = state_label(i_state,.false.)
      else
       c_st = ""
      end if
      idx = idx_mel_list('ME_C0'//trim(c_st),op_info) ! quick & dirty
      ! update me lists for transformation matrices if required
      if (opti_info%optref.ne.0.and.
     &    op_info%mel_arr(idx)%mel%fhand%last_mod(
     &    op_info%mel_arr(idx)%mel%fhand%current_record).gt.
     &    me_special(2)%mel%fhand%last_mod( ! 1
     &    me_special(2)%mel%fhand%current_record)) then
        call update_metric(me_dia,me_special,nspecial,
     &      fspc,nspcfrm,orb_info,op_info,str_info,strmap_info,
     &      opti_info%update_prc.gt.0.and.
     &      mod(imacit,max(opti_info%update_prc,1)).eq.0)
      else if (opti_info%optref.ne.0.and.iprlvl.ge.10) then
        write(lulog,*) ' Metric is already up to date!'
      end if

c dbg
c      call vec_from_da(me_dia%fhand,1,xbuf1,nwfpar)
c      print *,'Preconditioner:'
c      call wrt_mel_buf(lulog,5,xbuf1,me_dia,1,
c     &     me_dia%op%n_occ_cls,
c     &     str_info,orb_info)
c dbgend

      ! Transform residual into orthogonal basis:

      op_grd_name = trim(me_grd%op%name)

!     whenever a zeroing of the t1 part is requested,
!     also transform this according to the second formula 
!     and use the fourth ME_special for the transformed vector
!     d.h.: when using 'TR0'
!     set the transformation formula for Tges to Fspec 2
!     and the list to me_special(4)
!     very ugly hack
      if(lzero_flag)then
         op_trf_name = trim(me_special(4)%mel%op%name)
      else
         op_trf_name = trim(me_special(1)%mel%op%name)
      end if
      op_amp_name = trim(me_amp%op%name)

c dbg
c      call vec_from_da(me_grd%fhand,1,xbuf1,nwfpar)
c        write(lulog,*) 'untransformed gradient vector:'
c        write(lulog,*) xbuf1(1:nwfpar)
c dbg

      ! apply sign-fix
c      write(lulog,*) 'Fixing sign of residual for iopt =',iopt
      call optc_fix_signs2(me_grd%fhand,1,
     &                    opti_info,iopt,
     &                    nwfpar,xbuf1)
      ! assign op. to be transformed with list of gradient
      call assign_me_list(me_grd%label,
     &                    trim(op_trf_name),op_info)

      if (lzero_flag)then
      ! assign op. with list containing special vector
         call assign_me_list(me_special(4)%mel%label,
     &                    trim(op_amp_name),op_info)
         else
      ! assign op. with list containing special vector
            call assign_me_list(me_special(1)%mel%label,
     &                    trim(op_amp_name),op_info)
         end if
      ! use daggered transformation matrix if requested
      if (nspecial.ge.3)
     &   call assign_me_list(me_special(3)%mel%label,
     &                       me_special(3)%mel%op%name,op_info)

      ! calculate transformed residual

      if (lzero_flag)then
         call evaluate2(fspc(2),.true.,.true.,
     &            op_info,str_info,strmap_info,orb_info,
     &            xngrd(iopt),.true.) !get transformed res. norm
         else
         call evaluate2(fspc(1),.true.,.true.,
     &            op_info,str_info,strmap_info,orb_info,
     &            xngrd(iopt),.true.) !get transformed res. norm
       end if
cdbg
!      if (lzero_flag)then
!       call print_list('before zeroing',me_special(4)%mel,"NORM",
!     &                  -1d0,0d0,
!     &                  orb_info,str_info)
!       end if
cdbg

!     set all single excitations to zero if requestes
!     after long deliberation, I decided to also include V,V so one can be sure to include
!     **all** singular excitations
!     originally i set the T1 part after preconditioning to zero. 
!     But we need the gradient of the zeroed omega andd since preconditioning doesn't mix the
!     T_parts, it is instead done here.
      if (lzero_flag) then
         if (ntest.ge.100) write (lulog,*) " setting O1 part to 0.0"
         call set_blks(me_special(4)%mel,"P,H|P,V|V,H|V,V",0d0)

      xngrd(iopt)=xnormop(me_special(4)%mel)
      endif

      write (lulog,*) "Norm of transformed Gradient ",xngrd(iopt)

      if (lzero_flag)then
         call vec_from_da(me_special(4)%mel%fhand,
     &     me_special(4)%mel%fhand%current_record,
     &     xbuf1,nwfpar)
      else
         call vec_from_da(me_special(1)%mel%fhand,
     &     me_special(1)%mel%fhand%current_record,
     &     xbuf1,nwfpar)
      end if
      call vec_from_da(ffdia,
     &     ffdia%current_record,
     &     xbuf2,nwfpar)

      if (ntest.ge.100) then
        write(lulog,*) 'transformed gradient vector:'
c        write(lulog,*) xbuf1(1:nwfpar)
        if (lzero_flag)then
        call wrt_mel_buf(lulog,5,xbuf1,me_special(4)%mel,1,
     &       me_special(4)%mel%op%n_occ_cls,
     &       str_info,orb_info)
        else
        call wrt_mel_buf(lulog,5,xbuf1,me_special(1)%mel,1,
     &       me_special(1)%mel%op%n_occ_cls,
     &       str_info,orb_info)
        end if
      end if

      ! preconditioning step (optional: with iterative improvement)
      ! CAUTION: relies on that Aoff is on last special list
      call get_argument_value('method.MR','prc_iter',ival=prc_iter)
      call get_argument_value('method.MR','prc_impfac',xval=prc_impfac)
      if (prc_iter.ge.1) then
        call assign_me_list(me_special(nspecial)%mel%label,
     &                      me_special(nspecial)%mel%op%name,op_info)
        write(lulog,'(a,i4,a)') 'Trying to improve preconditioner in',
     &                 prc_iter,' iteration(s)'
      end if
      if (lzero_flag)then
      call prc_iterimp_rec(prc_iter,prc_iter,prc_impfac,
     &                     xbuf1,xbuf2,nwfpar,iopt,ffdia,
     &                     me_special(4)%mel,me_special(nspecial)%mel,
     &                     me_grd,fspc(2),opti_info,
     &                     orb_info,op_info,str_info,strmap_info)
      else
         call prc_iterimp_rec(prc_iter,prc_iter,prc_impfac,
     &                     xbuf1,xbuf2,nwfpar,iopt,ffdia,
     &                     me_special(1)%mel,me_special(nspecial)%mel,
     &                     me_grd,fspc(1),opti_info,
     &                     orb_info,op_info,str_info,strmap_info)
      end if
      if (ntest.ge.100) then
        write(lulog,*) 'preconditioned gradient vector:'
c        write(lulog,*) xbuf1(1:nwfpar)
      if (lzero_flag)then
        call wrt_mel_buf(lulog,5,xbuf1,me_special(4)%mel,1,
     &       me_special(4)%mel%op%n_occ_cls,
     &       str_info,orb_info)
        else
        call wrt_mel_buf(lulog,5,xbuf1,me_special(1)%mel,1,
     &       me_special(1)%mel%op%n_occ_cls,
     &       str_info,orb_info)
        end if
      end if

      ! get current trial vector (list will be overwritten)
      ! attention !! If, for some reason, the amplitudes ME
      ! be stored in different records, the 1 has to be changed
      call vec_from_da(ffamp,1,xbuf2,nwfpar)


      ! Transform new vector into original basis

      ! assign op. with its proper list
      call assign_me_list(me_amp%label,
     &                    trim(op_amp_name),op_info)
      ! assign op. to be transformed with special list
      if (lzero_flag) then
         call assign_me_list(me_special(4)%mel%label,
     &                    trim(op_trf_name),op_info)
      else
         call assign_me_list(me_special(1)%mel%label,
     &                    trim(op_trf_name),op_info)
      end if
      ! use non-daggered transformation matrix if requested
      if (nspecial.ge.3)
     &   call assign_me_list(me_special(2)%mel%label,
     &                       me_special(2)%mel%op%name,op_info)
      ! assign gradient list to its proper op.
      call assign_me_list(me_grd%label,
     &                    trim(op_grd_name),op_info)

      ! calculate transformed vector
      if (lzero_flag) then
         call evaluate2(fspc(2),.true.,.true.,
     &            op_info,str_info,strmap_info,orb_info,xdum,.false.)
      else
         call evaluate2(fspc(1),.true.,.true.,
     &            op_info,str_info,strmap_info,orb_info,xdum,.false.)
      end if   

      call vec_from_da(ffamp,1,xbuf1,nwfpar)
c dbg
c      ! vector could be put back to list here:
c      call vec_to_da(ffamp,1,xbuf2,nwfpar)
c dbgend

cdbg - the following lines do not give any sensible output:
c      if (lzero_flag)then
c       call print_list('before zeroing',me_special(4)%mel,"NORM",
c     &                  -1d0,0d0,
c     &                  orb_info,str_info)
c      else 
c       call print_list('before zeroing',me_special(1)%mel,"NORM",
c     &                  -1d0,0d0,
c     &                  orb_info,str_info)
c      end if
cdbg
      if (ntest.ge.100) then
        write(lulog,*) 'gradient vector afterwards:'
      do iii=1,nwfpar
        write(lulog,*) xbuf1(iii)
      end do
c      call wrt_mel_buf(lulog,5,xbuf1,me_amp,1,
c     &     me_amp%op%n_occ_cls,
c     &     str_info,orb_info)
      end if


      return
      end
