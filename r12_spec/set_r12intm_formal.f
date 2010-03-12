*----------------------------------------------------------------------*
      subroutine set_r12intm_formal(form_out,
     &     title,label_int,label_op,nop,typ_str,op_info,orb_info)
*----------------------------------------------------------------------*
*     set the formal definition of R12 intermediates
*     label_int: intermediate
*     label_op(1:nop): labels of the operators which contribute
*     typ_str:  describes the shape of the operator
*        e.g.  V: 'gxr'
*              X: 'rxr'
*              B: 'rfxr'
*       f: 1-particle part of (Hamilton) operator
*       g: 2-particle part of (Hamilton) operator
*       x: position of open lines
*       r: R12-operator
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
c      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(formula), intent(inout), target ::
     &     form_out

      integer, intent(in) ::
     &     nop
      character*(*), intent(in) ::
     &     title,
     &     label_int,
     &     label_op(nop),
     &     typ_str
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      character(3), parameter ::
     &     opdum_scal = '_1_',
     &     opdum_x    = '_X_',
     &     opdum_f    = '_F_',
     &     opdum_g    = '_G_'

      integer ::
     &     iop, idx, nfact, n_x, n_f, n_g,
     &     idx_intm, idx_scalar, idx_x, idx_f, idx_g,
     &     idx_op(nop+2), idx_prod(nop+2)
      character ::
     &     name*(form_maxlen_label*2)
      type(formula_item), target ::
     &     flist, flist_scr
      type(formula_item), pointer ::
     &     flist_pnt
      type(operator), pointer ::
     &     op_x, op_f, op_g, op_scalar, op_int

      integer, external ::
     &     idx_oplist2, idxlist
      call quit(1,'set_r12intm_formal','call to obsolete routine')

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_r12_intm_formal')
        write(luout,*) 'setting: ',trim(label_op(1))
        write(luout,*) 'type : ',trim(typ_str)
      end if

      if (nop.lt.2)
     &     call quit(1,'set_r12intm_formal',
     &                 'nop < 2 ?? too few! phew!')
      nfact = len_trim(typ_str)
      if (nfact.ne.nop+1)
     &     call quit(1,'set_r12intm_formal',
     &                 'typ_str not compatible with nop?')

      ! get indices of operators on label_list
      idx_intm = idx_oplist2(label_int,op_info)
      if (idx_intm.lt.0)
     &     call quit(1,'set_r12intm_formal',
     &     'label not on list: '//label_int)

      iop = 0
      idx_op(1:nfact) = 0
      do idx = 1, nfact
        if (typ_str(idx:idx).eq.'x') cycle
        iop = iop+1
        idx_op(idx) = idx_oplist2(label_op(iop),op_info)
        if (idx_op(idx).lt.0)
     &     call quit(1,'set_r12intm_formal',
     &     'label not on list: '//label_op(iop))
      end do

      ! scan typ_str for needed types
      n_x = 0
      n_f = 0
      n_g = 0
      do idx = 1, nfact
        if (typ_str(idx:idx).eq.'x') n_x = n_x+1
        if (typ_str(idx:idx).eq.'f') n_f = n_f+1
        if (typ_str(idx:idx).eq.'g') n_g = n_g+1
      end do

      if (n_x.gt.1)
     &     call quit(1,'set_r12intm_formal',
     &     'more than one x in typ_str: '//trim(typ_str))
      if (n_f.gt.1)
     &     call quit(1,'set_r12intm_formal',
     &     'more than one f in typ_str: '//trim(typ_str))
      if (n_g.gt.1)
     &     call quit(1,'set_r12intm_formal',
     &     'more than one g in typ_str: '//trim(typ_str))

      ! dummy operator: a scalar as target shape for operator product
      !                 expansion
      call add_operator(opdum_scal,op_info)
      idx_scalar = idx_oplist2(opdum_scal,op_info)
      op_scalar => op_info%op_arr(idx_scalar)%op
      call set_hop(op_scalar,opdum_scal,.false.,0,0,1,.false.,orb_info)

      ! dummy operator: the co-variant counterpart of the intermediate:
      if (n_x.gt.0) then
        call add_operator(opdum_x,op_info)
        idx_x = idx_oplist2(opdum_x,op_info)
        op_x => op_info%op_arr(idx_x)%op
        op_int => op_info%op_arr(idx_intm)%op
        call set_contrav_op(op_x,opdum_x,op_int,orb_info)
      end if

      ! dummy operator: 1 particle part of H
      if (n_f.gt.0) then
        call add_operator(opdum_f,op_info)
        idx_f = idx_oplist2(opdum_f,op_info)
        op_f => op_info%op_arr(idx_f)%op
        call set_hop_p(op_f,opdum_f,.false.,
     &       1,1,0,.true.,orb_info)
      end if

      ! dummy operator: 2 particle part of H
      if (n_g.gt.0) then
        call add_operator(opdum_g,op_info)
        idx_g = idx_oplist2(opdum_g,op_info)
        op_g => op_info%op_arr(idx_g)%op
        call set_hop(op_g,opdum_g,.false.,
     &       2,2,0,.true.,orb_info)
      end if

      ! set input array for expand_op_product:
      do idx = 1, nfact
        if (typ_str(idx:idx).eq.'x') then          
          idx_prod(idx) = idx_x
        else if (typ_str(idx:idx).eq.'f') then
          idx_prod(idx) = idx_f
        else if (typ_str(idx:idx).eq.'g') then
          idx_prod(idx) = idx_g
        else
          idx_prod(idx) = idx_op(idx)
        end if
      end do

      ! set up scratch formula
      call init_formula(flist_scr)
      flist_pnt => flist_scr
      call new_formula_item(flist_pnt,
     &     command_set_target_init,idx_scalar)
      flist_pnt => flist_pnt%next

      ! expand operator product, giving the intermediate as result
      call expand_op_product(flist_pnt,idx_scalar,
     &     1d0,nfact,idx_prod,
     &     -1,-1,
     &     -1,0,.false.,
     &     op_info)

      if (ntest.ge.100) then
        write(luout,*) 'intermediate formula'
        call print_form_list(luout,flist_scr,op_info)
      end if

      ! replace f and g by their actual operator
      if (n_f.gt.0) then
        idx = idxlist(idx_f,idx_prod,nfact,1)
        call init_formula(flist)
        call set_primitive_formula(flist,idx_op(idx),
     &       1d0,idx_f,.true.,op_info)
        call expand_subexpr(flist_scr,flist,.false.,op_info)
        call dealloc_formula_list(flist)
      end if
      if (n_g.gt.0) then
        idx = idxlist(idx_g,idx_prod,nfact,1)
        call init_formula(flist)
        call set_primitive_formula(flist,idx_op(idx),
     &       1d0,idx_g,.true.,op_info)
        call expand_subexpr(flist_scr,flist,.false.,op_info)
        call dealloc_formula_list(flist)
      end if

      ! generate actual formula by taking the derivative wrt dummy
      if (n_x.gt.0) then
        call init_formula(flist)
      
        call form_deriv3(flist,flist_scr,.true.,
     &     1,idx_x,0,idx_intm,op_info)
        flist_pnt => flist
      else
        flist_pnt => flist_scr
      end if
      
      ! prepare file:
      form_out%comment = trim(title)
      write(name,'(a,".fml")') trim(form_out%label)
      call file_init(form_out%fhand,name,ftyp_sq_unf,0)
      ! and write list
      call write_form_list(form_out%fhand,flist_pnt,form_out%comment)

      if (ntest.ge.100) then
        write(luout,*) 'final formula'
        call print_form_list(luout,flist_pnt,op_info)
      end if

      call dealloc_formula_list(flist_scr)
      if (n_x.gt.0) call dealloc_formula_list(flist)

      if (n_g.gt.0) call del_operator(opdum_g,op_info)
      if (n_f.gt.0) call del_operator(opdum_f,op_info)
      if (n_x.gt.0) call del_operator(opdum_x,op_info)
      call del_operator(opdum_scal,op_info)

      return
      end
