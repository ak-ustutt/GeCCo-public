*----------------------------------------------------------------------*
      subroutine set_r12intm_formal3(form_out,
     &     title,label_int,label_op,nop,typ_str,op_info,orb_info)
*----------------------------------------------------------------------*
*     set the formal definition of R12 intermediates
*     label_int: intermediate
*     label_op(1:nop): labels of the operators which contribute
*     typ_str:  describes the shape of the operator
*        e.g.  V: 'gxr'
*              X: 'rxr'
*              B: 'rfxr'
*              C: 'rf'
*              P: 'rgxr'
*              Z: 'rxgxr'
*
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
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'def_orbinf.h'
      include 'ifc_input.h'

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
     &     opdum_f     = '_F_',
     &     opdum_g     = '_G_'

      integer ::
     &     iop, idx, nfact, n_x, n_f, n_g, n_r,
     &     idx_intm, idx_f, idx_g, nvtx, len, ivtx, ndef
      integer ::
     &     avoid(10), connect(10), navoid, nconnect
      integer, pointer ::
     &     idx_prod(:), idx_supv(:), idx_op(:), occ_def(:,:,:)
      logical ::
     &     r12fix
      character ::
     &     name*(form_maxlen_label*2)
      type(formula_item), target ::
     &     flist, flist_scr
      type(formula_item), pointer ::
     &     flist_pnt
      type(operator), pointer ::
     &     op_x, op_f, op_g, op_shape, op_int, op

      integer, external ::
     &     idx_oplist2, idxlist

      ! Check to see whether we want to fix the R12-amplitudes.
      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'set_r12_intm_formal 3')
        write(lulog,*) 'setting: ',trim(label_int)
        write(lulog,*) 'type : ',trim(typ_str)
        write(lulog,*) 'input ops: '
        do iop = 1, nop
          write(lulog,*) '"',trim(label_op(iop)),'"'
        end do
      end if

      if (nop.lt.2)
     &     call quit(1,'set_r12intm_formal3',
     &                 'nop < 2 ?? too few! phew!')
      len = len_trim(typ_str)

      ! get indices of operators on label_list
      idx_intm = idx_oplist2(label_int,op_info)
      if (idx_intm.lt.0)
     &     call quit(1,'set_r12intm_formal3',
     &     'label not on list (1): '//label_int)

      allocate(idx_op(len))

      iop = 0
      idx_op(1:len) = 0
      nvtx = 2  ! external lines at borders
      do idx = 1, len
        if (typ_str(idx:idx).eq.'x') then
          nvtx = nvtx + 2
          cycle
        end if
        nvtx = nvtx+1
        iop = iop+1
        idx_op(idx) = idx_oplist2(label_op(iop),op_info)
        if (idx_op(idx).lt.0)
     &     call quit(1,'set_r12intm_formal3',
     &     'label not on list (2): '//label_op(iop))
      end do

      allocate(idx_prod(nvtx),idx_supv(nvtx))
      
      ! scan typ_str for needed types
      n_x = 0
      n_f = 0
      n_g = 0
      n_r = 0
      do idx = 1, nvtx
        if (typ_str(idx:idx).eq.'x') n_x = n_x+1
        if (typ_str(idx:idx).eq.'f' .or.
     &      typ_str(idx:idx).eq.'F') n_f = n_f+1
        if (typ_str(idx:idx).eq.'g') n_g = n_g+1
        if (typ_str(idx:idx).eq.'r') n_r = n_r+1
      end do

      if (2*n_x+n_f+n_g+n_r+2.ne.nvtx)
     &     call quit(1,'set_r12intm_formal3',
     &     'unknown identifier in typ_str: '//trim(typ_str))
      if (n_f.gt.1)
     &     call quit(1,'set_r12intm_formal3',
     &     'more than one f in typ_str: '//trim(typ_str))
      if (n_g.gt.1)
     &     call quit(1,'set_r12intm_formal3',
     &     'more than one g in typ_str: '//trim(typ_str))

      nfact = n_f+n_g+n_r
c      if (n_x.gt.0) nfact = nfact+1
      nfact = nfact+1
c dbg
c      print *,'nfact = ',nfact
c      print *,'n_f = ',n_f
c      print *,'n_g = ',n_g
c      print *,'n_r = ',n_r
c dbg

      ! Point to the actual intermediate.
      op_int => op_info%op_arr(idx_intm)%op

      if (ntest.ge.100) then
        write(lulog,*) 'definition of intermediate:'
        call print_op_occ(lulog,op_int)
      end if

      ! dummy operator: 1 particle part of H 
      if (n_f.gt.0) then
        idx = index(typ_str,'f')
        if (idx.le.0) idx = index(typ_str,'F')
        call add_operator(opdum_f,op_info)
        idx_f = idx_oplist2(opdum_f,op_info)
        op_f => op_info%op_arr(idx_f)%op
        if(n_x.le.1)then
          allocate(occ_def(ngastp,2,4))
          if (typ_str(idx:idx).eq.'f') then ! (P/X-space only)
            ndef = 4
            occ_def = 0
            occ_def(IPART,1,1) = 1
            occ_def(IPART,2,1) = 1
            occ_def(IPART,1,2) = 1
            occ_def(IEXTR,2,2) = 1
            occ_def(IEXTR,1,3) = 1
            occ_def(IPART,2,3) = 1
            occ_def(IEXTR,1,4) = 1
            occ_def(IEXTR,2,4) = 1
          else if (typ_str(idx:idx).eq.'F') then ! (H-space only)
            ndef = 1
            occ_def = 0
            occ_def(IHOLE,1,1) = 1
            occ_def(IHOLE,2,1) = 1
          else
            call quit(1,'set_r12intm_formal3','???')
          end if
          call set_uop2(op_f,opdum_f,
     &         occ_def,ndef,1,(/0,0/),-1,orb_info)
          deallocate(occ_def)
        elseif(n_x.eq.2)then
          allocate(occ_def(ngastp,2,2))
          ndef = 2
          occ_def = 0
          occ_def(IHOLE,1,1) = 1
          occ_def(IEXTR,2,1) = 1
          occ_def(IHOLE,1,2) = 1
          occ_def(IPART,2,2) = 1
          call set_uop2(op_f,opdum_f,
     &         occ_def,ndef,1,(/0,0/),-1,orb_info)
          deallocate(occ_def)
c          call set_hop(op_f,opdum_f,.false.,
c     &         1,1,0,.true.,orb_info)
        else
          call quit(1,'set_r12intm_formal3','Too many X with F')
        endif
      end if

      ! dummy operator: 2 particle part of H
      if (n_g.gt.0) then
        call add_operator(opdum_g,op_info)
        idx_g = idx_oplist2(opdum_g,op_info)
        op_g => op_info%op_arr(idx_g)%op
        call set_hop(op_g,opdum_g,.false.,
     &       2,2,0,.true.,IEXTR,1,orb_info)
      end if

      ! set input array for expand_op_product:
      idx_prod(1) = idx_intm
      idx_supv(1) = 1
      idx_prod(nvtx) = idx_intm
      idx_supv(nvtx) = 1
      ivtx = 2
      do idx = 1, len
        if (typ_str(idx:idx).eq.'x') then
          idx_prod(ivtx:ivtx+1) = idx_intm
          idx_supv(ivtx:ivtx+1) = 1
          ivtx = ivtx+2
c          endif
        else if (typ_str(idx:idx).eq.'f'.or.
     &           typ_str(idx:idx).eq.'F') then
          idx_prod(ivtx) = idx_f
          idx_supv(ivtx) = 2
          ivtx = ivtx+1
        else if (typ_str(idx:idx).eq.'g') then
          idx_prod(ivtx) = idx_g
          idx_supv(ivtx) = 3
          ivtx = ivtx+1
        else
          idx_prod(ivtx) = idx_op(idx)
          if (idx.eq.1) idx_prod(ivtx) = -idx_prod(ivtx)
          idx_supv(ivtx) = 4+idx
          ivtx = ivtx+1
        end if
      end do

      ! set up scratch formula
      call init_formula(flist_scr)
      flist_pnt => flist_scr
      call new_formula_item(flist_pnt,
     &     command_set_target_init,idx_intm)
      flist_pnt => flist_pnt%next

      avoid = 0
      connect = 0
      navoid = 0
      nconnect = 0
      if (trim(typ_str).eq.'rgxr') then ! xrgxxrx
        avoid(1:2) = (/2,6/)
        navoid = 1
      else if (trim(typ_str).eq.'rxgxr') then ! xrxxgxxrx
        avoid(1:4) = (/2,7, 3,8/)
        navoid = 2
        connect(1:2) = (/2,8/)
        nconnect = 1
      elseif (trim(typ_str).eq.'rxfxr')then ! xrxxfxxrx
        avoid(1:6) = (/2,7, 2,5, 3,8/)
        navoid = 3
        connect(1:4) = (/2,8, 5,8/)
        nconnect = 2
      end if

c test
c      if (n_f+n_g.gt.0) then
c        nconnect = 0
c        do ivtx = 1, nvtx
c          if (idx_supv(ivtx).eq.2.or.idx_supv(ivtx).eq.3) then
c            idx = ivtx
c            exit
c          end if
c        end do
c        do ivtx = 1, nvtx
c          if (idx_supv(ivtx).gt.4) then
c            if (ivtx.lt.idx) then
c              connect(nconnect*2+1:nconnect*2+2) = (/ivtx,idx/)
c            else
c              connect(nconnect*2+1:nconnect*2+2) = (/idx,ivtx/)
c            end if
c            nconnect = nconnect + 1
c          end if
c        end do
c      end if
cc dbg
c      print *,'connect: ',connect(1:nconnect*2)
cc dbg
c test

      ! expand operator product, giving the intermediate as result
      call expand_op_product2(flist_pnt,idx_intm,
     &     1d0,nvtx,nfact,
     &     idx_prod,idx_supv,
     &     -1,-1,
     &     connect,nconnect,
     &     avoid,navoid,
     &     -1,0,
     &     .false.,op_info)

      if (ntest.ge.1000) then
        write(lulog,*) 'intermediate formula'
        call print_form_list(lulog,flist_scr,op_info)
      end if

      ! replace f and g by their actual operator
      if (n_f.gt.0) then
        idx = index(typ_str,'f')
        if (idx.le.0) idx = index(typ_str,'F')
        op   => op_info%op_arr(idx_op(idx))%op
        call form_op_replace(opdum_f,op%name,.true.,flist_scr,op_info)
      end if
      if (n_g.gt.0) then
        idx = index(typ_str,'g')
        op   => op_info%op_arr(idx_op(idx))%op
        call form_op_replace(opdum_g,op%name,.true.,flist_scr,op_info)
      end if
      
      ! prepare file:
      form_out%comment = trim(title)
      write(name,'(a,".fml")') trim(form_out%label)
      call file_init(form_out%fhand,name,ftyp_sq_unf,0)
      ! and write list
      call write_form_list(form_out%fhand,flist_scr,form_out%comment)

      if (ntest.ge.100) then
        write(lulog,*) 'final formula'
        call print_form_list(lulog,flist_pnt,op_info)
      end if

      call dealloc_formula_list(flist_scr)
      if (n_g.gt.0) call del_operator(opdum_g,op_info)
      if (n_f.gt.0) call del_operator(opdum_f,op_info)

      deallocate(idx_op,idx_prod,idx_supv)

      return
      end
