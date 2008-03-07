*----------------------------------------------------------------------*
      subroutine set_contrav_op(op_contra,name,
     &     op_co,orb_info)
*----------------------------------------------------------------------*
*     set the contravariant operator which compensates any open lines
*     of the given operator     
*     WARNING: not yet rigorously tested, but works so far for the call
*         in set_r12intm_formal()
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'stdunit.h'
      include 'ifc_baserout.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 100

      type(operator), intent(inout) ::
     &     op_contra
      character, intent(in) ::
     &     name*(*)
      type(operator), intent(in) ::
     &     op_co

      type(orbinf), intent(in), target ::
     &     orb_info

      logical ::
     &     add_before, add_after
      integer ::
     &     iblk, ioffblk, ijoin, idef, idx, idx_co, idx_contra,
     &     njoined_co, njoined_contra, nblk
      integer, pointer ::
     &     hpvxgas(:,:), ngas, nspin

      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,
     &     'set_contrav_op speaking')
      
      nspin => orb_info%nspin
      ngas => orb_info%ngas
      hpvxgas => orb_info%ihpvgas

      if (len_trim(name).gt.len_opname)
     &     call quit(1,'set_contrav_op',
     &     'name too long: "'//trim(name)//'"')

      op_contra%name = '        '
      op_contra%name = name
      
      njoined_co = op_co%njoined
      nblk    = op_co%n_occ_cls
      ! at least njoined_co-1:
      njoined_contra = njoined_co-1
      ! scan through occupations of input operator and get
      ! njoined for output operator
      idx = 0 
      add_before = .false.
      add_after  = .false.
      do iblk = 1, nblk
        ! EX from first?
        idx = idx+1
        add_before = add_before .or.
     &      iocc_nonzero(iocc_xdn(1,op_co%ihpvca_occ(1:ngastp,1:2,idx)))
        ! DX from last?
        idx = idx-1+njoined_co
        add_after = add_after .or.
     &      iocc_nonzero(iocc_xdn(2,op_co%ihpvca_occ(1:ngastp,1:2,idx)))
      end do

      if (add_before) njoined_contra = njoined_contra + 1
      if (add_after ) njoined_contra = njoined_contra + 1

      if (ntest.ge.100)
     &     write(luout,*) 'njoined_contra = ',njoined_contra

      op_contra%njoined = njoined_contra

      op_contra%dagger = op_co%dagger

      op_contra%n_occ_cls = op_co%n_occ_cls

      call init_operator(op_contra,orb_info)

      idx_contra = 0
      idx_co     = 0
      do iblk = 1, nblk

        ioffblk = (iblk-1)*njoined_contra

        ! set occupation for block
        idx_co     = idx_co    +1
        if (add_before) then
          idx_contra = idx_contra+1
          op_contra%ihpvca_occ(1:ngastp,1:2,idx_contra) =
     &         iocc_xdn(1,op_co%ihpvca_occ(1:ngastp,1:2,idx_co))
        end if
        do ijoin = 1, njoined_co-1
          idx_contra = idx_contra+1
          op_contra%ihpvca_occ(1:ngastp,1:2,idx_contra) =
     &         iocc_xdn(2,op_co%ihpvca_occ(1:ngastp,1:2,idx_co)) +
     &         iocc_xdn(1,op_co%ihpvca_occ(1:ngastp,1:2,idx_co+1))
          idx_co = idx_co + 1          
        end do
        if (add_after) then
          idx_contra = idx_contra+1
          op_contra%ihpvca_occ(1:ngastp,1:2,idx_contra) =
     &         iocc_xdn(2,op_co%ihpvca_occ(1:ngastp,1:2,idx_co))
        end if

        ! summed occupation
        op_contra%ica_occ(1:2,iblk) = 0

        do ijoin = 1, njoined_contra
          op_contra%ica_occ(1,iblk) = op_contra%ica_occ(1,iblk)+
     &         ielsum(op_contra%
     &                ihpvca_occ(1:ngastp,1,ioffblk+ijoin),ngastp)
          op_contra%ica_occ(2,iblk) = op_contra%ica_occ(2,iblk)+
     &         ielsum(op_contra%
     &                ihpvca_occ(1:ngastp,2,ioffblk+ijoin),ngastp)
        end do

      end do

      ! another loop for setting the restrictions:
      idx_contra = 0
      idx_co     = 0
      do iblk = 1, nblk

        ! set occupation for block
        idx_co     = idx_co    +1
        if (add_before) then
          idx_contra = idx_contra+1
          op_contra%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,idx_contra)=
     &          irest_xdn (  1,
     &        op_co%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,idx_co),
     &                       hpvxgas,ngas,nspin )
        end if
        do ijoin = 1, njoined_co-1
          idx_contra = idx_contra+1
          op_contra%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,idx_contra)=
     &          irest_xdn (  2,
     &        op_co%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,idx_co),
     &                       hpvxgas,ngas,nspin)
     &        + irest_xdn (  1,
     &        op_co%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,idx_co+1),
     &                       hpvxgas,ngas,nspin)
          idx_co = idx_co + 1          
        end do
        if (add_after) then
          idx_contra = idx_contra+1
          op_contra%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,idx_contra)=
     &          irest_xdn (  2,
     &        op_co%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,idx_co),
     &                       hpvxgas,ngas,nspin)
        end if

      end do

      if (ntest.ge.100) then
        write(luout,*) 'New intermediate: ',trim(op_contra%name)
        write(luout,*) 'generated ',op_contra%n_occ_cls,' blocks'
        idx = 0
        do iblk = 1, op_contra%n_occ_cls
          do ijoin = 1, op_contra%njoined
            idx = idx+1
            call wrt_occ_rstr(luout,iblk,
     &           op_contra%ihpvca_occ(1,1,idx),
     &           op_contra%igasca_restr(1,1,1,1,1,idx),ngas,nspin)
          end do
        end do
      end if

      return
      end
