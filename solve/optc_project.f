*----------------------------------------------------------------------*
      subroutine optc_project(me_amp,me_trf,me_dia,me_special,nspecial,
     &     nwfpar,xbuf1,fspc,nspcfrm,iopt,imacit,opti_info,
     &     orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
*     projects out linear dependencies in me_amp
*     fspc(1)         : transformation formula
*     me_trf          : the list that is the result of fspc (if not 
*                       identical to me_amp)
*     me_special(1)   : transformed vector or vector to be transformed
*     me_special(2,3) : transformation matrices
*     me_special(4)   : projector
*     higher fspc, me_special: for update of metric (if needed)
*
*     after projection, the vector resides on the list me_amp again
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
     &     nspecial, nspcfrm, nwfpar, iopt, imacit
      type(me_list_array), intent(inout) ::
     &     me_special(nspecial)
      type(me_list), intent(in) ::
     &     me_amp,me_dia,me_trf
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
     &     idx, nsec, isec, stsec, ndsec, ioff,
     &     nblk, iblk, nj, iblkoff, jblk
      character(len_opname) ::
     &     op_spc_name, op_trf_name, op_amp_name

      integer, pointer ::
     &     nwfpsec(:), idstsec(:), nsec_arr(:), occ(:,:,:)

      real(8), pointer ::
     &     signsec(:)
      real(8) ::
     &     xdum

      integer, external ::
     &     idx_mel_list, iblk_occ

      real(8), external :: dnrm2

      if (opti_info%optref.eq.-2) then
        ! update metric, trafo matrices and projector if not up to date
        idx = idx_mel_list('ME_C0',op_info)  ! quick & dirty
        if (op_info%mel_arr(idx)%mel%fhand%last_mod(
     &      op_info%mel_arr(idx)%mel%fhand%current_record).gt.
     &      me_special(2)%mel%fhand%last_mod(1))
     &      call update_metric(me_dia,me_special,nspecial,
     &        fspc,nspcfrm,orb_info,op_info,str_info,strmap_info,
     &        opti_info%update_prc.gt.0.and.
     &        mod(imacit,max(opti_info%update_prc,1)).eq.0)
      end if

      call assign_me_list(me_special(1)%mel%label,
     &                      me_special(1)%mel%op%name,op_info)

      ! put new vector to special list for transformation
      call list_copy(me_amp,me_special(1)%mel,.false.)

      ! use projection matrix
      call assign_me_list(me_special(4)%mel%label,
     &                      me_special(2)%mel%op%name,op_info)

      ! evaluate projected vector
      call evaluate2(fspc(1),.true.,.true.,
     &         op_info,str_info,strmap_info,orb_info,xdum,.true.)


c      ! Since formally we get a transposed vector, we need to
c      ! account for sign changes when reordering
c      nsec = opti_info%nsec(iopt)
c      if (nsec.gt.1) then
c        ioff = sum(opti_info%nsec(1:iopt))-nsec
c        nwfpsec => opti_info%nwfpsec(ioff+1:ioff+nsec)
c        idstsec => opti_info%idstsec(ioff+1:ioff+nsec)
c        signsec => opti_info%signsec2(ioff+1:ioff+nsec)
c        call vec_from_da(ffamp,1,xbuf1,nwfpar)
c        do isec = 1, nsec
c          xbuf1(idstsec(isec):idstsec(isec)+nwfpsec(isec)-1) =
c     &       signsec(isec)
c     &       *xbuf1(idstsec(isec):idstsec(isec)+nwfpsec(isec)-1)
c        end do
c        call vec_to_da(ffamp,1,xbuf1,nwfpar)
c      end if

c      ! Now add "redundant" T components?
c      if (opti_info%update_prc.gt.0.and.
c     &    nspecial.ge.8.and.nspcfrm.ge.4.or.
c     &    opti_info%update_prc.le.0.and.nspecial.ge.7.and.nspcfrm.ge.3)
c     &   then
c        ! does T(3)red exist? If so, we need two steps
c        if (opti_info%update_prc.gt.0.and.nspecial.eq.9.or.
c     &      opti_info%update_prc.le.0.and.nspecial.eq.8) then
c          ! evaluate T(2)red (projector list is already assigned)
c          call evaluate2(fspc(nspcfrm-1),.true.,.false.,
c     &           op_info,str_info,strmap_info,orb_info,xdum,.false.)
c          ! Here we need to add this to T
c          nj = me_amp%op%njoined
c          nblk = me_special(nspecial-1)%mel%op%n_occ_cls
c          do iblk = 1, nblk
c           iblkoff = (iblk-1)*nj
c           occ => me_special(nspecial-1)%mel%op%ihpvca_occ(1:ngastp,1:2,
c     &                                          iblkoff+1:iblkoff+nj)
c           jblk = iblk_occ(occ,.false.,me_amp%op,
c     &                  me_special(nspecial-1)%mel%op%blk_version(iblk))
c           if (jblk.lt.1) call quit(1,'optc_project','blk not found(1)')
c           call add_opblk(xdum,0,1d0,me_special(nspecial-1)%mel,me_amp,
c     &                    iblk,jblk,orb_info,.false.)
c          end do
c        end if
c        ! Now evaluate T(3)red (or T(2)red if not done above)
c        call evaluate2(fspc(nspcfrm),.true.,.false.,
c     &         op_info,str_info,strmap_info,orb_info,xdum,.false.)
c        ! and add this to T as well
c        nblk = me_special(nspecial)%mel%op%n_occ_cls
c        do iblk = 1, nblk
c          iblkoff = (iblk-1)*nj
c          occ => me_special(nspecial)%mel%op%ihpvca_occ(1:ngastp,1:2,
c     &                                         iblkoff+1:iblkoff+nj)
c          jblk = iblk_occ(occ,.false.,me_amp%op,
c     &                    me_special(nspecial)%mel%op%blk_version(iblk))
c          if (jblk.lt.1) call quit(1,'optc_project','blk not found(2)')
c          call add_opblk(xdum,0,1d0,me_special(nspecial)%mel,me_amp,
c     &                   iblk,jblk,orb_info,.false.)
c        end do
c      end if

      return
      end
