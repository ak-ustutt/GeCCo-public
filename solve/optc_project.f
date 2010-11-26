*----------------------------------------------------------------------*
      subroutine optc_project(me_amp,me_dia,me_special,nspecial,
     &     nwfpar,xbuf1,fspc,nspcfrm,iopt,opti_info,
     &     orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
*     projects out linear dependencies in me_amp
*     fspc(1)         : transformation formula
*     me_special(1)   : transformed vector or vector to be transformed
*     me_special(2,3) : transformation matrices
*     me_special(4)   : projector
*     higher fspc, me_special: for update of metric (if needed)
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
     &     nspecial, nspcfrm, nwfpar, iopt
      type(me_list_array), intent(inout) ::
     &     me_special(nspecial)
      type(me_list), intent(in) ::
     &     me_amp,me_dia
      real(8), intent(inout) ::
     &     xbuf1(*)

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
     &     idx, nsec, isec, stsec, ndsec, ioff

      integer, pointer ::
     &     nwfpsec(:), idstsec(:), nsec_arr(:)

      real(8), pointer ::
     &     signsec(:)
      real(8) ::
     &     xdum

      type(filinf), pointer ::
     &     ffamp

      integer, external ::
     &     idx_mel_list

      ! pointers to file handle
      ffamp => me_amp%fhand

      if (opti_info%optref.eq.-2) then
        ! update metric, trafo matrices and projector if not up to date
        idx = idx_mel_list('ME_C0',op_info)  ! quick & dirty
        if (op_info%mel_arr(idx)%mel%fhand%last_mod(1).gt.
     &      me_special(2)%mel%fhand%last_mod(1))
     &      call update_metric(me_dia,me_special,nspecial,
     &        fspc,nspcfrm,orb_info,op_info,str_info,strmap_info,
     &        opti_info%update_prc,opti_info%singrm)
      end if

      ! put new vector to special list for transformation
      call list_copy(me_amp,me_special(1)%mel)

      ! use projection matrix
      call assign_me_list(me_special(4)%mel%label,
     &                       me_special(2)%mel%op%name,op_info)

      ! evaluate projected vector
      call evaluate2(fspc(1),
     &         op_info,str_info,strmap_info,orb_info,xdum,.false.)

      ! Since formally we get a transposed vector, we need to
      ! account for sign changes when reordering
      nsec = opti_info%nsec(iopt)
      if (nsec.gt.1) then
        ioff = sum(opti_info%nsec(1:iopt))-nsec
        nwfpsec => opti_info%nwfpsec(ioff+1:ioff+nsec)
        idstsec => opti_info%idstsec(ioff+1:ioff+nsec)
        signsec => opti_info%signsec2(ioff+1:ioff+nsec)
        call vec_from_da(ffamp,1,xbuf1,nwfpar)
        do isec = 1, nsec
          xbuf1(idstsec(isec):idstsec(isec)+nwfpsec(isec)-1) =
     &       signsec(isec)
     &       *xbuf1(idstsec(isec):idstsec(isec)+nwfpsec(isec)-1)
        end do
        call vec_to_da(ffamp,1,xbuf1,nwfpar)
      end if

      return
      end
