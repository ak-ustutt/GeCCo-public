*----------------------------------------------------------------------*
      subroutine clone_operator(op_clone,op_template,dagger,orb_info)
*----------------------------------------------------------------------*
*     copy all operator information from template to clone 
*     (except name and ID)
*     dagger = .true. :  build the adjoint operator
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'ifc_operators.h'
      
      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     dagger
      type(operator), intent(in) ::
     &     op_template
      type(operator), intent(inout) ::
     &     op_clone
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     n_occ_cls, ngas, nspin, njoined,
     &     id_save, ioff, iblk
      character ::
     &     name_save*(len_opname)

      ! save name and ID
      id_save = op_clone%id
      name_save = op_clone%name

      ! check to avoid old "dagger" convention
      if (op_template%dagger)
     &     call quit(1,'clone_operator',
     &     'the use of "dagger" in the operator definition is obsolete')

      ! first of all we do this:
      op_clone = op_template

      ! reset name and ID
      op_clone%id   = id_save    
      op_clone%name = name_save 
      op_clone%assoc_list(1:mxlen_melabel) = ' '

      ! however, the above statement copied the pointers as
      ! pointers, i.e. the clone points to the same arrays as
      ! the template: THIS IS HIGHLY DANGEROUS
      ! so: let's get our own space for the arrays and copy these

      call init_operator(op_clone,orb_info)

      n_occ_cls = op_template%n_occ_cls
      njoined = op_template%njoined
      ngas = orb_info%ngas
      if (associated(op_template%ihpvca_occ)) then
        if (.not.dagger) then
          op_clone%ihpvca_occ = op_template%ihpvca_occ
        else
          do iblk = 1, n_occ_cls
            ioff = (iblk-1)*njoined
            op_clone%ihpvca_occ(1:ngastp,1:2,ioff+1:ioff+njoined)
     &           = iocc_dagger_n(
     &      op_template%ihpvca_occ(1:ngastp,1:2,ioff+1:ioff+njoined),
     &           njoined)
          end do
        end if
      end if

      if (associated(op_template%ica_occ)) then
        if (.not.dagger) then
          op_clone%ica_occ = op_template%ica_occ
        else
          do iblk = 1, n_occ_cls
            op_clone%ica_occ(1,iblk) = op_template%ica_occ(2,iblk)
            op_clone%ica_occ(2,iblk) = op_template%ica_occ(1,iblk)
          end do
        end if
      end if

      if (associated(op_template%igasca_restr)) then
        if (.not.dagger) then
          op_clone%igasca_restr = op_template%igasca_restr
        else
          nspin = op_template%nspin
          do iblk = 1, n_occ_cls
            ioff = (iblk-1)*njoined
           op_clone%igasca_restr(1:2,1:ngas,1:2,1:2,
     &                                     1:nspin,ioff+1:ioff+njoined)
     &           = irest_dagger_n(
     &  op_template%igasca_restr(1:2,1:ngas,1:2,1:2,
     &                                     1:nspin,ioff+1:ioff+njoined),
     &           njoined,ngas,nspin)
          end do
        end if
      end if

      if (associated(op_template%formal_blk)) then
        op_clone%formal_blk = op_template%formal_blk
      end if

      if (ntest.ge.100) then
        write(luout,*) 'Clone operator produced the following:'
        call print_op_occ(luout,op_clone)
      end if

      return
      end
