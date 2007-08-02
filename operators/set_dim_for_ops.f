*----------------------------------------------------------------------*
      subroutine set_dim_for_ops(op_list,nops,str_info,orb_info)
*----------------------------------------------------------------------*
*
*     andreas, april 2007
*
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 100

      type(operator_list), intent(inout), target ::
     &     op_list
      integer, intent(in) ::
     &     nops
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(operator_list), pointer ::
     &     current

      integer ::
     &     iprint, ipass, iop, nexc, nsym, ndis, iocc_cls,
     &     len, i

      iprint = max(ntest,iprlvl)

      if (iprint.gt.3) then
        write(luout,*) 'Setting up operator dimensions'
      end if

      nsym = orb_info%nsym

      ! scan through operator list
      current => op_list
      do iop = 1, nops
        if (.not.associated(current%op))
     &       call quit(0,'set_dim_for_ops','buggy operator list (a)')
        if (current%op%formal) then
          if (iop.lt.nops) current => current%next
          cycle
        end if

        call init_operator(1,current%op,orb_info)
c        allocate(current%op%off_op_occ(current%op%n_occ_cls),
c     &           current%op%len_op_occ(current%op%n_occ_cls),
c     &           current%op%off_op_gmo(current%op%n_occ_cls),
c     &           current%op%len_op_gmo(current%op%n_occ_cls),
c     &           current%op%off_op_gmox(current%op%n_occ_cls),
c     &           current%op%len_op_gmox(current%op%n_occ_cls))
c        do iocc_cls = 1, current%op%n_occ_cls
c          nexc = min(current%op%ica_occ(1,iocc_cls),
c     &               current%op%ica_occ(2,iocc_cls))
c          allocate(current%op%len_op_gmo(iocc_cls)%gam_ms(nsym,nexc+1),
c     &             current%op%off_op_gmo(iocc_cls)%gam_ms(nsym,nexc+1))
c        end do
        ! set up length info for operator
        ipass = 1
        call set_op_dim(ipass,.false.,current%op,str_info,nsym)

        call init_operator(2,current%op,orb_info)
c        do iocc_cls = 1, current%op%n_occ_cls
c          nexc = min(current%op%ica_occ(1,iocc_cls),
c     &               current%op%ica_occ(2,iocc_cls))
c          ndis = current%op%off_op_gmox(iocc_cls)%maxd
c          allocate(current%op%len_op_gmox(iocc_cls)%
c     &                        d_gam_ms(ndis,nsym,nexc+1),
c     &             current%op%off_op_gmox(iocc_cls)%
c     &                        d_gam_ms(ndis,nsym,nexc+1),
c     &             current%op%off_op_gmox(iocc_cls)%
c     &                        did(ndis,nsym,nexc+1),
c     &             current%op%off_op_gmox(iocc_cls)%ndis(nsym,nexc+1))
c        end do

        ! extended length info for operator
        ipass = 2
        call set_op_dim(ipass,.false.,current%op,str_info,nsym)

        if (iprint.gt.3) then
          len = len_trim(current%op%name)
          write(luout,*) trim(current%op%name),(' ',i=1,len_opname-len),
     &         ' -- total length: ',current%op%len_op
        end if
        
        if (iop.lt.nops.and..not.associated(current%next))
     &       call quit(0,'set_dim_for_ops','buggy operator list (b)')
        if (iop.lt.nops) current => current%next
      end do

      end
