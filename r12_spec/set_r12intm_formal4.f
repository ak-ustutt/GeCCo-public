*----------------------------------------------------------------------*
      subroutine set_r12intm_formal4(form_out,
     &     title,label_int,label_op,nop,typ_str,op_info,orb_info)
*----------------------------------------------------------------------*
*     set the formal definition of R12 intermediates
*     label_int: label of intermediate
*     label_op(1:nop): labels of the operators which contribute 
*              (1 -- R12, 2 -- Hamiltonian)
*     typ_str:  describes the operator
*         X, V, B, Bp, Bh, P, Z, ... 
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
     &     ntest = 1000

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

      logical ::
     &     def_fpp, def_fhh, def_g, unknown, def_fp3f
      integer ::
     &     idx, nfact, 
     &     idx_intm, idx_r, idx_f, idx_g, idx_h, idx_rpl,
     &     nvtx, len, ivtx, ndef, njoined_int, calls
      integer ::
     &     avoid(20), connect(20), project(20), project2(20),
     &     navoid, nconnect, nproject, nproject2
      integer ::
     &     idx_prod(20), idx_supv(20)
      integer, pointer ::
     &     occ_def(:,:,:)
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

      ! preliminary: go to old code, if necessary
      if (typ_str(1:1).eq.'r'.or.typ_str(1:1).eq.'g'.or.
     &    typ_str(1:1).eq.'f') then
        call set_r12intm_formal3(form_out,
     &     title,label_int,label_op,nop,typ_str,op_info,orb_info)
        return
      end if

      ! Check to see whether we want to fix the R12-amplitudes.
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_r12_intm_formal 4')
        write(luout,*) 'setting: ',trim(label_int)
        write(luout,*) 'type : ',trim(typ_str)
        write(luout,*) 'input ops: '
        do idx = 1, nop
          write(luout,*) '"',trim(label_op(idx)),'"'
        end do
      end if

      if (nop.lt.1)
     &     call quit(1,'set_r12intm_formal',
     &                 'nop < 2 ?? too few! phew!')
      len = len_trim(typ_str)

      ! get indices of operators on label_list
      idx_intm = idx_oplist2(label_int,op_info)
      if (idx_intm.lt.0)
     &     call quit(1,'set_r12intm_formal',
     &     'label not on list (1): '//label_int)

      idx_r = idx_oplist2(label_op(1),op_info)
      if (idx_r.lt.0)
     &     call quit(1,'set_r12intm_formal',
     &     'label not on list (2): '//label_op(1))

      if (trim(typ_str).ne.'X') then
        idx_h = idx_oplist2(label_op(2),op_info)
        if (idx_h.lt.0)
     &       call quit(1,'set_r12intm_formal',
     &       'label not on list (3): '//label_op(2))
      end if

      ! Point to the actual intermediate.
      op_int => op_info%op_arr(idx_intm)%op

      if (ntest.ge.100) then
        write(luout,*) 'definition of intermediate:'
        call print_op_occ(luout,op_int)
      end if

      njoined_int = op_int%njoined

      unknown = .false.
      def_fhh = .false.
      def_fpp = .false.
      def_fp3f = .false.
      def_g   = .false.
      idx_g   = -99
      idx_f   = -99
      nconnect = 0
      navoid   = 0
      nproject = 0
      nproject2 = 0
      calls = 1
      select case(trim(typ_str))
      case('X')
        if (njoined_int.eq.1) then
          idx_prod(1:4) = (/idx_intm,-idx_r,idx_r,idx_intm/)
          idx_supv(1:4) = (/1       ,2    ,3    ,1       /)
          nvtx     = 4
          nfact    = 3
          project(1:4)  = (/2,3,2,3/)
          nproject = 1
        else if (njoined_int.eq.2) then
          idx_prod(1:6) =
     &         (/idx_intm,-idx_r,idx_intm,idx_intm,idx_r,idx_intm/)
          idx_supv(1:6) =
     &         (/1       ,2    ,1       ,1       ,3    ,1       /)
          nvtx     = 6
          nfact    = 3
          project(1:4)  = (/2,4,2,3/)
          nproject = 1
        else
          unknown = .true.
        end if
      case('X''')
        if (njoined_int.eq.1) then
          idx_prod(1:4) = (/idx_intm,-idx_r,idx_r,idx_intm/)
          idx_supv(1:4) = (/1       ,2    ,3    ,1       /)
          nvtx     = 4
          nfact    = 3
          project(1:4)  = (/2,3,1,4/)
          nproject = 1
c        else if (njoined_int.eq.2) then
        else
          unknown = .true.
        end if
      case('V')
        idx_rpl = 2
        def_g   = .true.
        if (njoined_int.eq.1) then
          idx_prod(1:4) = (/idx_intm,idx_g,idx_r,idx_intm/)
          idx_supv(1:4) = (/1       ,2    ,3    ,1       /)
          nvtx     = 4
          nfact    = 3
          project(1:4)  = (/2,3,2,3/)
          nproject = 1          
        else if (njoined_int.eq.2) then
          idx_prod(1:6) =
     &         (/idx_intm,idx_g,idx_intm,idx_intm,idx_r,idx_intm/)
          idx_supv(1:6) =
     &         (/1       ,2    ,1       ,1       ,3    ,1       /)
          nvtx     = 6
          nfact    = 3
          project(1:4)  = (/2,4,2,3/)
          nproject = 1
        else
          unknown = .true.
        end if
      case('V3')
        def_g = .true.
        idx_rpl = 2
        if(njoined_int.eq.2)then
          idx_prod(1:6) = (/idx_intm,idx_g,idx_intm,idx_intm,idx_r,
     &                      idx_intm/)
          idx_supv(1:6) = (/       1,    2,       1,       1,    3,
     &                             1/)
          nvtx = 6
          nfact = 3
          connect(1:2) = (/2,5/)
          nconnect = 1
          avoid(1:4) = (/3,5,2,4/)
          navoid = 2
        else
          unknown = .true.
        endif
      case('V4')
        def_g = .true.
        idx_rpl = 2
        if(njoined_int.eq.3)then
          idx_prod(1:9) = (/idx_intm,idx_g,idx_intm,idx_intm,idx_r,
     &                      idx_intm,idx_intm,idx_r,idx_intm/)
          idx_supv(1:9) = (/       1,    2,       1,       1,    3,
     &                             1,       1,    4,       1/)
          nvtx = 9
          nfact = 4
          connect(1:4) = (/2,5,2,8/)
          nconnect = 2
          avoid(1:10) = (/2,4,2,7,5,8,3,5,3,8/)
          navoid = 5
        else
          unknown = .true.
        endif
      case('B','Bp')
        def_fpp = .true.
        idx_rpl  = 3
        if (njoined_int.eq.1) then
          idx_prod(1:5) = (/idx_intm,-idx_r,idx_f,idx_r,idx_intm/)
          idx_supv(1:5) = (/1       ,2    ,3    ,4    ,1       /)
          nvtx     = 5
          nfact    = 4
          connect(1:4)  = (/2,3,3,4/)
          nconnect = 2
          project(1:4)  = (/2,4,1,IEXTR/) ! ... or X
          nproject = 1          
          calls = 2
          project2(1:4)  = (/2,4,1,IPART/) ! P ...
          nproject2 = 1          
        else if (njoined_int.eq.2) then
          idx_prod(1:7) = (/idx_intm,-idx_r,idx_f,idx_intm,idx_intm,
     &                                                  idx_r,idx_intm/)
          idx_supv(1:7) = (/1       ,2    ,3    ,1       ,1       ,
     &                                                  4    ,1       /)
          nvtx     = 7
          nfact    = 4
          connect(1:4)  = (/2,3,3,5/)
          nconnect = 2
          project(1:4)  = (/2,5,2,3/)
          nproject = 1
        else
          unknown = .true.
        end if
      case('Bh')
        def_fhh = .true.
        idx_rpl  = 3
        if (njoined_int.eq.1) then
          idx_prod(1:5) = (/idx_intm,-idx_r,idx_f,idx_r,idx_intm/)
          idx_supv(1:5) = (/1       ,2    ,3    ,4    ,1       /)
          nvtx     = 5
          nfact    = 4
          connect(1:4)  = (/2,3,3,4/)
          nconnect = 2
          project(1:4)  = (/2,4,2,3/)
          nproject = 1          
        else
          unknown = .true.
        end if
      case('C')
        def_fpp = .true.
        idx_rpl = 2
        if (njoined_int.eq.1) then
          idx_prod(1:4) = (/idx_intm,idx_f,idx_r,idx_intm/)
          idx_supv(1:4) = (/1       ,2    ,3    ,1       /)
          nvtx     = 4
          nfact    = 3
          connect(1:2)  = (/2,3/)
          nconnect = 1
        else
          unknown = .true.
        end if
      case('P')
        def_g = .true.
        idx_rpl = 3
        if(njoined_int.eq.2)then
          idx_prod(1:7) = (/idx_intm,-idx_r,idx_g,idx_intm,idx_intm,
     &                      idx_r,idx_intm/)
          idx_supv(1:7) = (/       1,     2,    3,       1,       1,
     &                          4,       1/)
          nvtx = 7
          nfact = 4
          connect(1:4) = (/2,3,3,6/)
          nconnect = 2
          avoid(1:4) = (/2,5,2,6/)
          navoid = 2
        else
          unknown = .true.
        endif
      case('P3G')
        def_g = .true.
        idx_rpl = 5
        if(njoined_int.eq.3)then
          idx_prod(1:9) = (/idx_intm,-idx_r,idx_intm,idx_intm,idx_g,
     &                      idx_intm,idx_intm,idx_r,idx_intm/)
          idx_supv(1:9) = (/       1,     2,       1,       1,    3,
     &                             1,       1,    4,       1/)
          nvtx = 9
          nfact = 4
          connect(1:4) = (/2,5,5,8/)
          nconnect = 2
          avoid(1:6) = (/2,7,2,8,3,8/)
          navoid = 3
        else
          unknown = .true.
        endif
      case('P3F')
        def_fp3f = .true.
        idx_rpl = 5
        if(njoined_int.eq.3)then
          idx_prod(1:9) = (/idx_intm,-idx_r,idx_intm,idx_intm,idx_f,
     &                      idx_intm,idx_intm,idx_r,idx_intm/)
          idx_supv(1:9) = (/       1,     2,       1,       1,    3,
     &                             1,       1,    4,       1/)
          nvtx = 9
          nfact = 4
          connect(1:4) = (/2,8,5,8/)
          nconnect = 2
          avoid(1:6) = (/2,7,2,5,3,8/)
          navoid = 3
        else
          unknown = .true.
        endif
      case('Z')
        def_g = .true.
        idx_rpl = 5
        if(njoined_int.eq.3)then
          idx_prod(1:9) = (/idx_intm,-idx_r,idx_intm,idx_intm,idx_g,
     &                      idx_intm,idx_intm,idx_r,idx_intm/)
          idx_supv(1:9) = (/       1,     2,       1,       1,    3,
     &                             1,       1,    4,       1/)
          nvtx = 9
          nfact = 4
          connect(1:6) = (/2,5,2,8,5,8/)
          nconnect = 3
          avoid(1:4) = (/2,7,3,8/)
          navoid = 2
        else
          unknown = .true.
        endif
      case('Z4')
        def_g = .true.
        idx_rpl = 5
        if(njoined_int.eq.4)then
          idx_prod(1:12) = (/idx_intm,-idx_r,idx_intm,idx_intm,idx_g,
     &                       idx_intm,idx_intm,idx_r,idx_intm,idx_intm,
     &                       idx_r,idx_intm/)
          idx_supv(1:12) = (/       1,     2,       1,       1,    3,
     &                              1,       1,    4,       1,       1,
     &                           5,       1/)
          nvtx = 12
          nfact = 5
          connect(1:8) = (/2,8,2,11,5,8,5,11/)
          nconnect = 4
          avoid(1:10) = (/2,7,2,10,3,8,3,11,2,5/)
          navoid = 5
        else
          unknown = .true.
        endif
      case('K4')
        def_g = .true.
        idx_rpl = 5
        if(njoined_int.eq.3)then
          idx_prod(1:9) = (/idx_intm,-idx_r,idx_intm,idx_intm,idx_g,
     &                      idx_intm,idx_intm,idx_r,idx_intm/)
          idx_supv(1:9) = (/       1,     2,       1,       1,    3,
     &                             1,       1,    4,       1/)
          nvtx = 9
          nfact = 4
          connect(1:4) = (/2,8,5,8/)
          nconnect = 2
          avoid(1:6) = (/2,7,3,8,2,5/)
          navoid = 3
        else
          unknown = .true.
        endif
      case default
        call quit(1,'set_r12intm_formal4',
     &       'unknown type: '//trim(typ_str))
      end select

      if (unknown) then
        write(luout,*) 'njoined_int = ',njoined_int
        call quit(1,'set_r12intm_formal4',
     &       'unknown shape for: '//trim(typ_str))
      end if

      ! dummy operator: 1 particle part of H 
      if (def_fhh) then
        call add_operator(opdum_f,op_info)
        idx_f = idx_oplist2(opdum_f,op_info)
        op_f => op_info%op_arr(idx_f)%op
        allocate(occ_def(ngastp,2,4))
        ndef = 1
        occ_def = 0
        occ_def(IHOLE,1,1) = 1
        occ_def(IHOLE,2,1) = 1
        call set_uop2(op_f,opdum_f,
     &       occ_def,ndef,1,orb_info)
        deallocate(occ_def)
        idx_prod(idx_rpl) = idx_f
      else if (def_fpp) then
        call add_operator(opdum_f,op_info)
        idx_f = idx_oplist2(opdum_f,op_info)
        op_f => op_info%op_arr(idx_f)%op
        allocate(occ_def(ngastp,2,4))
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
        call set_uop2(op_f,opdum_f,
     &       occ_def,ndef,1,orb_info)
        deallocate(occ_def)
        idx_prod(idx_rpl) = idx_f
      elseif(def_fp3f)then
        call add_operator(opdum_f,op_info)
        idx_f = idx_oplist2(opdum_f,op_info)
        op_f => op_info%op_arr(idx_f)%op
        allocate(occ_def(ngastp,2,2))
        ndef = 2
        occ_def = 0
        occ_def(IHOLE,1,1) = 1
        occ_def(IEXTR,2,1) = 1
        occ_def(IHOLE,1,2) = 1
        occ_def(IPART,2,2) = 1
        call set_uop2(op_f,opdum_f,
     &       occ_def,ndef,1,orb_info)
        deallocate(occ_def)
        idx_prod(idx_rpl) = idx_f
      else if (def_g) then
        call add_operator(opdum_g,op_info)
        idx_g = idx_oplist2(opdum_g,op_info)
        op_g => op_info%op_arr(idx_g)%op
        call set_hop(op_g,opdum_g,.false.,
     &       2,2,0,.true.,orb_info)
        idx_prod(idx_rpl) = idx_g
      end if

      ! set up scratch formula
      call init_formula(flist_scr)
      flist_pnt => flist_scr
      call new_formula_item(flist_pnt,
     &     command_set_target_init,idx_intm)
      flist_pnt => flist_pnt%next

      ! expand operator product, giving the intermediate as result
      call expand_op_product2(flist_pnt,idx_intm,
     &     1d0,nvtx,nfact,
     &     idx_prod,idx_supv,
     &     -1,-1,
     &     connect,nconnect,
     &     avoid,navoid,
     &     project,nproject,
     &     op_info)

      ! quick fix:
      if (calls.eq.2) then
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do        
        call expand_op_product2(flist_pnt,idx_intm,
     &     1d0,nvtx,nfact,
     &     idx_prod,idx_supv,
     &     -1,-1,
     &     connect,nconnect,
     &     avoid,navoid,
     &     project2,nproject2,
     &     op_info)
        call reorder_formula(flist_scr,op_info)
      end if

      if (ntest.ge.1000) then
        write(luout,*) 'intermediate formula'
        call print_form_list(luout,flist_scr,op_info)
      end if

      ! replace f and g by their actual operator
      if (def_fhh.or.def_fpp.or.def_fp3f) then
        op   => op_info%op_arr(idx_h)%op
        call form_op_replace(opdum_f,op%name,.true.,flist_scr,op_info)
      end if
      if (def_g) then
        op   => op_info%op_arr(idx_h)%op
        call form_op_replace(opdum_g,op%name,.true.,flist_scr,op_info)
      end if
      
      ! prepare file:
      form_out%comment = trim(title)
      write(name,'(a,".fml")') trim(form_out%label)
      call file_init(form_out%fhand,name,ftyp_sq_unf,0)
      ! and write list
      call write_form_list(form_out%fhand,flist_scr,form_out%comment)

      if (ntest.ge.100) then
        write(luout,*) 'final formula'
        call print_form_list(luout,flist_scr,op_info)
      end if

      call dealloc_formula_list(flist_scr)
      if (def_g) call del_operator(opdum_g,op_info)
      if (def_fhh.or.def_fpp.or.def_fp3f)
     &     call del_operator(opdum_f,op_info)

      return
      end
