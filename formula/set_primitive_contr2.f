*----------------------------------------------------------------------*
      subroutine set_primitive_contr2(contr,
     &     fac,idx_opsh,iblk_opsh,
     &     idx_op,iblk_op,use_opsh_structure,
     &     op_info)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'

      type(contraction), intent(inout) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info
      real(8), intent(in) ::
     &     fac
      integer, intent(in) ::
     &     idx_opsh, iblk_opsh, idx_op, iblk_op
      logical, intent(in) ::
     &     use_opsh_structure

      type(operator), pointer ::
     &     op, opsh
      integer ::
     &     njoined, iblk_off, iblk_off_sh, ijoin

      opsh => op_info%op_arr(idx_opsh)%op
      if (.not.use_opsh_structure)
     &     op => op_info%op_arr(idx_opsh)%op

      njoined = opsh%njoined
      if (.not.use_opsh_structure.and.njoined.ne.op%njoined) then
        write(lulog,*) 'njoined: ',opsh%njoined,op%njoined
        call quit(1,'set_primitive_contr2','inconsistency!')
      end if

      iblk_off    = (iblk_op-1)*njoined
      iblk_off_sh = (iblk_opsh-1)*njoined

      call resize_contr(contr,njoined,0,njoined,0)
      contr%fac = fac
      contr%idx_res = idx_opsh
      contr%iblk_res = iblk_off_sh+1
c      contr%dagger   = 
      contr%nvtx = njoined
      contr%svertex(1:njoined) = 1
      contr%nxarc = njoined
      contr%vertex(1:njoined)%idx_op = idx_op
      do ijoin = 1, njoined
c        if (.not.use_opsh_structure) then ! ???
c          contr%vertex(ijoin)%iblk_op = iblk_off+ijoin
c        else
        contr%vertex(ijoin)%iblk_op = iblk_op!iblk_off+1
c        end if
        contr%xarc(ijoin)%link(1) = ijoin
        contr%xarc(ijoin)%link(2) = ijoin
        if (.not.use_opsh_structure) then
          contr%xarc(ijoin)%occ_cnt =
     &           op%ihpvca_occ(1:ngastp,1:2,iblk_off+ijoin)
        else 
          contr%xarc(ijoin)%occ_cnt =
     &         opsh%ihpvca_occ(1:ngastp,1:2,iblk_off_sh+ijoin)
        end if
      end do

      call update_svtx4contr(contr)

      return
      end
