*----------------------------------------------------------------------*
      subroutine clone_operator(op_clone,op_template,orb_info)
*----------------------------------------------------------------------*
*     copy all operator information from template to clone 
*     (except name and ID)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'

      type(operator), intent(in) ::
     &     op_template
      type(operator), intent(inout) ::
     &     op_clone
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     n_occ_cls, ngas,
     &     id_save
      character ::
     &     name_save*(len_opname)

      ! save name and ID
      id_save = op_clone%id
      name_save = op_clone%name

      ! first of all we do this:
      op_clone = op_template

      ! reset name and ID
      op_clone%id   = id_save    
      op_clone%name = name_save 

      ! however, the above statement copied the pointers as
      ! pointers, i.e. the clone points to the same arrays as
      ! the template: THIS IS HIGHLY DANGEROUS
      ! so: let's get our own space for the arrays and copy these

      n_occ_cls = op_template%n_occ_cls
      ngas = orb_info%ngas
      if (associated(op_template%ihpvca_occ)) then
        allocate(op_clone%ihpvca_occ(ngastp,2,n_occ_cls))
        op_clone%ihpvca_occ = op_template%ihpvca_occ
      else
        op_clone%ihpvca_occ => null()
      end if

      if (associated(op_template%ica_occ)) then
        allocate(op_clone%ica_occ(2,n_occ_cls))
        op_clone%ica_occ = op_template%ica_occ
      else
        op_clone%ica_occ => null()
      end if

      if (associated(op_template%igasca_restr)) then
        allocate(op_clone%igasca_restr(2,ngas,2,2,n_occ_cls))
        op_clone%igasca_restr = op_template%igasca_restr
      else
        op_clone%igasca_restr => null()
      end if

      if (associated(op_template%len_op_occ))
     &     call quit(1,'clone_operator',
     &     'not yet prepared for cloning length info')

      end
