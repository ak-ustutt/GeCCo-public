*----------------------------------------------------------------------*
      subroutine optc_prc_traf(me_amp,me_grd,me_dia,me_special,nspecial,
     &     nwfpar,xbuf1,xbuf2,
     &     fspc,nspcfrm,xngrd,iopt,opti_info,
     &     orb_info,op_info,str_info,strmap_info)
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

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nspecial, iopt, nspcfrm, nwfpar
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


      integer ::
     &     idx, nopt, nsec, isec, stsec, ndsec
      character(len_opname) ::
     &     op_grd_name, op_trf_name, op_amp_name

      integer, pointer ::
     &     nwfpsec(:), idstsec(:), nsec_arr(:)

      real(8), pointer ::
     &     signsec(:)
      real(8) ::
     &     xdum

      type(filinf), pointer ::
     &     ffamp, ffgrd, ffdia

      integer, external ::
     &     idx_mel_list

      ! pointers to file handle
      ffamp => me_amp%fhand
      ffgrd => me_grd%fhand
      ffdia => me_dia%fhand

      nopt = opti_info%nopt

      ! update me lists for transformation matrices if required
      idx = idx_mel_list('ME_C0',op_info)  ! quick & dirty
      
      if (opti_info%optref.ne.0.and.
     &    op_info%mel_arr(idx)%mel%fhand%last_mod(1).gt.
     &    me_special(2)%mel%fhand%last_mod(1)) then
        call update_metric(me_dia,me_special,nspecial,
     &      fspc,nspcfrm,orb_info,op_info,str_info,strmap_info,
     &      opti_info%update_prc,opti_info%singrm)
      else if (opti_info%optref.ne.0.and.iprlvl.ge.10) then
        write(luout,*) ' Metric is already up to date!'
      end if

c dbg
c      call vec_from_da(me_dia%fhand,1,xbuf1,nwfpar)
c      print *,'Preconditioner:'
c      call wrt_mel_buf(luout,5,xbuf1,me_dia,1,
c     &     me_dia%op%n_occ_cls,
c     &     str_info,orb_info)
c dbgend

      ! Transform residual into orthogonal basis:

      op_grd_name = trim(me_grd%op%name)
      op_trf_name = trim(me_special(1)%mel%op%name)
      op_amp_name = trim(me_amp%op%name)
      ! assign op. to be transformed with list of gradient
      call assign_me_list(me_grd%label,
     &                    trim(op_trf_name),op_info)
      ! assign op. with list containing special vector
      call assign_me_list(me_special(1)%mel%label,
     &                    trim(op_amp_name),op_info)
      ! use daggered transformation matrix if requested
      if (nspecial.ge.3)
     &   call assign_me_list(me_special(3)%mel%label,
     &                       me_special(3)%mel%op%name,op_info)

      ! calculate transformed residual
      call evaluate2(fspc(1),
     &            op_info,str_info,strmap_info,orb_info,
     &            xngrd(iopt),.true.) !get transformed res. norm

      call vec_from_da(me_special(1)%mel%fhand,1,xbuf1,nwfpar)

      if (ntest.ge.100) then
        write(luout,*) 'transformed gradient vector:'
        write(luout,*) xbuf1(1:nwfpar)
      end if

      write(luout,'(x,a,i1,a,x,g10.4)')
     &   'Norm of transformed residual for vector ',
     &   iopt,':',xngrd(iopt)

      call vec_from_da(ffdia,1,xbuf2,nwfpar)

      nsec_arr => opti_info%nsec(1:nopt)
      nsec = sum(nsec_arr)
      nwfpsec => opti_info%nwfpsec(1:nsec)
      idstsec => opti_info%idstsec(1:nsec)
      signsec => opti_info%signsec2(1:nsec)
      stsec = 1
      ndsec = 0
      if (iopt.gt.1) stsec = stsec + nsec_arr(iopt-1)
      ndsec = ndsec + nsec_arr(iopt)

      ! Divide by precond., account for sign changes if necessary

      do isec = stsec, ndsec
        call diavc(xbuf1(idstsec(isec)),xbuf1(idstsec(isec)),
     &             signsec(isec),xbuf2(idstsec(isec)),
     &             0d0,nwfpsec(isec))
      end do

      ! put new vector to special list for transformation
      call vec_to_da(me_special(1)%mel%fhand,1,xbuf1,nwfpar)
      ! get current trial vector (list will be overwritten)
      call vec_from_da(ffamp,1,xbuf2,nwfpar)

      ! Transform new vector into original basis

      ! assign op. with its proper list
      call assign_me_list(me_amp%label,
     &                    trim(op_amp_name),op_info)
      ! assign op. to be transformed with special list
      call assign_me_list(me_special(1)%mel%label,
     &                    trim(op_trf_name),op_info)
      ! use non-daggered transformation matrix if requested
      if (nspecial.ge.3)
     &   call assign_me_list(me_special(2)%mel%label,
     &                       me_special(2)%mel%op%name,op_info)
      ! assign gradient list to its proper op.
      call assign_me_list(me_grd%label,
     &                    trim(op_grd_name),op_info)

      ! calculate transformed vector
      call evaluate2(fspc(1),
     &            op_info,str_info,strmap_info,orb_info,xdum,.false.)

      call vec_from_da(ffamp,1,xbuf1,nwfpar)
c dbg
c      ! vector could be put back to list here:
c      call vec_to_da(ffamp,1,xbuf2,nwfpar)
c dbgend

      if (ntest.ge.100) then
        write(luout,*) 'gradient vector afterwards:'
        write(luout,*) xbuf1(1:nwfpar)
      end if


      return
      end
