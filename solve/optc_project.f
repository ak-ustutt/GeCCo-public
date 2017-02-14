*----------------------------------------------------------------------*
      subroutine optc_project(me_amp,me_trf,me_dia,
     &     me_u,use_u,
     &     me_special,nspecial,
     &     nwfpar,xbuf1,fspc,nspcfrm,iopt,imacit,i_state,lzero,
     &     opti_info,orb_info,op_info,str_info,strmap_info)
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
      include 'mdef_target_info.h'

      integer, parameter ::
     &     ntest = 100
      character(len=*),parameter::
     &     i_am = "optc_project"
      integer, intent(in) ::
     &     nspecial, nspcfrm, nwfpar, iopt, imacit, i_state,use_u
      type(me_list_array), intent(inout) ::
     &     me_special(nspecial),me_u(3)
      type(me_list), intent(in) ::
     &     me_amp,me_dia,me_trf
      real(8), intent(inout) ::
     &     xbuf1(*)
      logical, intent(in)::
     &     lzero

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
     &     idx, ioff,
     &     nblk, iblk, nj, iblkoff, jblk
      character(len_opname) ::
     &     op_spc_name, op_trf_name, op_amp_name

      integer, pointer ::
     &     occ(:,:,:)

      real(8) ::
     &     xdum

      integer, external ::
     &     idx_mel_list, iblk_occ

      real(8), external :: dnrm2
      type(operator),pointer::
     &     op_amp, op_orth, op_trf
      character(len_target_name) ::
     &     c_st
      character(len_target_name), external ::
     &     state_label
      if(ntest.ge.10)
     &     write(lulog,*) "Projecting now "
      print *, i_am," called for ",iopt
      if (opti_info%optref.eq.-2) then
        ! update metric, trafo matrices and projector if not up to date
        c_st = state_label(i_state,.false.)
        idx = idx_mel_list('ME_C0'//trim(c_st),op_info)  ! quick & dirty
        if (op_info%mel_arr(idx)%mel%fhand%last_mod(
     &      op_info%mel_arr(idx)%mel%fhand%current_record).gt.
     &      me_special(2)%mel%fhand%last_mod(
     &      me_special(2)%mel%fhand%current_record))
     &       call update_metric(me_dia,
     &       me_u,use_u,
     &       me_special,nspecial,
     &        fspc,nspcfrm,orb_info,op_info,str_info,strmap_info,
     &        opti_info%update_prc.gt.0.and.
     &        mod(imacit,max(opti_info%update_prc,1)).eq.0)
      end if
      if(lzero)then
         call print_list("unitary_matrix", me_u(1)%mel,
     &        "LIST",0d0,0d0,
     &        orb_info,str_info)
         call print_list("unitary_matrix reordered", me_u(2)%mel,
     &        "LIST",0d0,0d0,
     &        orb_info,str_info)
         call print_list("unitary_matrix reo_transposed", me_u(3)%mel,
     &        "LIST",0d0,0d0,
     &        orb_info,str_info)
      end if
      
      
      if (lzero)then            !also use_u is >0 and optref == -3
         op_amp => me_amp%op
         op_orth => me_special(7)%mel%op
         op_trf => me_special(2)%mel%op
         print *,"current record", me_u(2)%mel%fhand%current_record
         call print_list("untransformed", me_amp,
     &        "LIST",0d0,0d0,
     &        orb_info,str_info)
         call assign_me_list(me_special(7)%mel%label,
     &        op_amp%name,op_info)
         call assign_me_list(me_amp%label,
     &        op_orth%name,op_info)
         call assign_me_list(me_u(2)%mel%label,
     &        op_trf%name,op_info)
         call evaluate2(fspc(3),.true.,.true.,
     &        op_info,str_info,strmap_info,orb_info,xdum,.true.)
         call set_blks(me_special(7)%mel,"P,H|P,V|V,H|V,V",0d0)
c         call print_list("transformed vector", me_special(7)%mel,
c     &        "LIST",0d0,0d0,
c     &        orb_info,str_info)
         call assign_me_list(me_special(7)%mel%label,
     &        op_orth%name,op_info)
         call assign_me_list(me_amp%label,
     &        op_amp%name,op_info)
         call assign_me_list(me_u(3)%mel%label,
     &        op_trf%name,op_info)
         call evaluate2(fspc(3),.true.,.true.,
     &        op_info,str_info,strmap_info,orb_info,xdum,.true.)
         call print_list("backtranformed", me_amp,
     &        "LIST",0d0,0d0,
     &        orb_info,str_info)
      else
         call assign_me_list(me_special(1)%mel%label,
     &        me_special(1)%mel%op%name,op_info)
         
! put new vector to special list for transformation
         call list_copy(me_amp,me_special(1)%mel,.false.)
         
      ! use projection matrix
         call assign_me_list(me_special(4)%mel%label,
     &        me_special(2)%mel%op%name,op_info)

      ! evaluate projected vector
         call evaluate2(fspc(1),.true.,.true.,
     &        op_info,str_info,strmap_info,orb_info,xdum,.true.)
      end if

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
