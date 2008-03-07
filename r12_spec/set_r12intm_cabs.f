*----------------------------------------------------------------------*
      subroutine set_r12intm_cabs(form_out,
     &     title,labels,nlabels,int_type,ansatz,approx,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     future: driver for setting up all kinds of R12 intermediates
*             for a given approximation
*             (I probably need to shift some contents to subroutines)
*     
*     currently:
*     generate the CABS approximated representation of the 
*     R12 intermediates of the type
*
*     (A Q(1,2) B)^{ij}_{kl} = 
*                     (AB)^{ij}_{kl} - ( A P(1,2)  B )^{ij}_{kl}
*
*     based on set_v_intermediate() etc. by GWR
*
*     label(1) -- the intermediate to be generated
*     label(2),
*     label(3) -- the labels of the operators A, B
*     label(4) -- the (AB) product operator
*   further input for setting up more complex terms (beyond approx A)
*     label(5) -- the X intermediate
*     label(6) -- the Hamilton operator
*     label(7) -- ....
*
*     approx string: formatted string
*
*     123456789012
*     1  2   3   4
*     XX XXX XXX X
*     A
*     A' GBC HY1 S
*     B  EBC HY2 
*     C      HX
*
*     1: approximation type
*     2: assumed Brillouin condition
*     3: hybrid approximation for exchange terms
*     4: symmetrise
*
*     andreas, Jan 2008
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 1000

      type(formula), intent(inout), target ::
     &     form_out
      integer, intent(in) ::
     &     nlabels, ansatz
      character*(*), intent(in) ::
     &     title, approx, int_type,
     &     labels(nlabels)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      character(6), parameter ::
     &     op_scr_a  = '_SCR_A',
     &     op_scr_b  = '_SCR_B',
     &     op_scr_1  = '_SCR_1',
     &     op_scr_f  = '_SCR_F'
      character(7), parameter ::
     &     op_scr_xu = '_SCR_Xu',
     &     op_scr_xl = '_SCR_Xl'
      character(8), parameter ::
     &     op_scr_dum = '_SCR_C'

      logical ::
     &     symmetrise
      integer ::
     &     idx_intm, idx_a, idx_b, idx_ab, idx_opa, idx_opb,
     &     idx_dum, idx_op1, idx_h, idx_x, idx_f, idx_xl, idx_xu,
     &     ndef, idx, nblk, ioff1, ioff3, nlink
      integer, pointer ::
     &     occ_def(:,:,:), linked(:)
      character ::
     &     name*(form_maxlen_label*2)
      type(formula_item), target ::
     &     flist_a, flist_b, flist_adb, flist_a_b, flist_xdx, flist
      type(formula_item), pointer ::
     &     flist_pnt
      type(operator), pointer ::
     &     opa_pnt, opb_pnt, op_intm, opdum_pnt, op1_pnt,
     &     opxl_pnt, opxu_pnt, opf_pnt

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_r12_intm_cabs')
        write(luout,*) 'setting: ',trim(labels(1))
        write(luout,*) 'ansatz = ',ansatz
        write(luout,*) 'approx = "',approx(1:12),'"'
      end if

      symmetrise = approx(12:12).eq.'S'

      if (ntest.ge.100) write(luout,*) 'symmetrise = ',symmetrise

      if (nlabels.lt.4)
     &     call quit(1,'set_r12intm_cabs',
     &     'too few labels (<4)')

      idx_intm = idx_oplist2(labels(1),op_info)
      if (idx_intm.lt.0)
     &     call quit(1,'set_r12intm_cabs',
     &     'label not on list: '//labels(1))
      op_intm => op_info%op_arr(idx_intm)%op

      idx_a    = idx_oplist2(labels(2),op_info)
      if (idx_a.lt.0)
     &     call quit(1,'set_r12intm_cabs',
     &     'label not on list: '//labels(2))
      idx_b    = idx_oplist2(labels(3),op_info)
      if (idx_b.lt.0)
     &     call quit(1,'set_r12intm_cabs',
     &     'label not on list: '//labels(3))
      idx_ab   = idx_oplist2(labels(4),op_info)
      if (idx_ab.lt.0)
     &     call quit(1,'set_r12intm_cabs',
     &     'label not on list: '//labels(4))

      ! this should become a subroutine (up to mark "end future subroutine")

      ! make sure that for both operators only the appropriate
      ! blocks are selected 
      ! -> define two temporary operators
      call add_operator(op_scr_a,op_info)
      idx_opa = idx_oplist2(op_scr_a,op_info)
      opa_pnt => op_info%op_arr(idx_opa)%op
      ! -> and a dummy operator (the contravariant of the result)
      call add_operator(op_scr_dum,op_info)
      idx_dum = idx_oplist2(op_scr_dum,op_info)
      opdum_pnt => op_info%op_arr(idx_dum)%op
      call add_operator(op_scr_b,op_info)
      idx_opb = idx_oplist2(op_scr_b,op_info)
      opb_pnt => op_info%op_arr(idx_opb)%op
      ! -> and a scalar as preliminary result of AB expansion
      call add_operator(op_scr_1,op_info)
      idx_op1 = idx_oplist2(op_scr_1,op_info)
      op1_pnt => op_info%op_arr(idx_op1)%op
      
      ! set dummy operator
c      call set_contrav_op(opdum_pnt,op_scr_dum,op_intm,orb_info)
      call set_spacer_op(op_intm,opdum_pnt,op_scr_dum,op_info,orb_info)
      ! set scalar
c      call set_hop(op1_pnt,op_scr_1,.false.,0,0,1,.false.,orb_info)
      call get_op_shape(op_intm,op1_pnt,op_scr_1,op_info,orb_info)

      nblk = opdum_pnt%n_occ_cls

      ! set the shape operators
      if (ansatz.eq.1) ndef = 5
      if (ansatz.eq.3) ndef = 4
      if (ansatz.eq.2)
     &     call quit(1,'set_r12intm_cabs','old ansatz 2 seems obsolete')
      allocate(occ_def(ngastp,2,ndef*nblk))
      ! set the target operator shapes as appropriate for the respective
      ! ansatz
      call set_opshapes(occ_def,ndef,nblk,ansatz)

      call set_uop(opa_pnt,op_scr_a,.true.,
     &     occ_def,ndef*nblk,orb_info)
c      call set_hop(opa_pnt,op_scr_a,.false.,
c     &             2,2,2,.true.,orb_info)
      call set_uop(opb_pnt,op_scr_b,.false.,
     &     occ_def,ndef*nblk,orb_info)

      deallocate(occ_def)

      call init_formula(flist_adb)
      call new_formula_item(flist_adb,command_set_target_init,
     &     idx_op1)
      ! set up terms resulting from A P(1,2) B:
      flist_pnt => flist_adb%next
      if(int_type(1:2).eq.'V ')then
        nlink = 2
        allocate(linked(2*nlink))
        linked = (/1,3,2,3/)
      elseif(int_type(1:2).eq.'V+')then
        nlink = 2
        allocate(linked(2*nlink))
        linked = (/1,3,1,2/)
      else
        nlink = 3
        allocate(linked(2*nlink))
        linked = (/1,3,1,2,2,3/)
      endif

      do idx = 1, nblk
        ioff1 = 0
        ioff3 = 0
        if(int_type(1:2).eq.'V ')then
          ioff1 = (idx-1)*ndef                             
        elseif(int_type(1:2).eq.'V+')then
          ioff3 = (idx-1)*ndef
        endif

        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        enddo
        call expand_op_product(flist_pnt,idx_op1,
     &       -1d0,3,(/idx_opa,idx_dum,idx_opb/),
     &       (/ioff1+1,idx,ioff3+1/),(/ioff1+ndef,idx,ioff3+ndef/),
     &       linked,nlink,.true.,op_info)
        ! the ".true." in the last line:
        ! the routine also forms 'dodgy' contractions, i.e. those where
        ! the arcs form between non-matching faces of the operators. These 
        ! are needed to properly evaluate the intermediates.
      enddo
      deallocate(linked)

      if (ntest.ge.1000) then
        write(luout,*) 'result from expand_op_product'
        call print_form_list(luout,flist_adb,op_info)
      end if

      ! replace a and b by the actual operators
      ! we do this step before taking the derivatives (see (*) below)
      ! as expand_subexpr() still is buggy for results with njoined>1
      call init_formula(flist_a)
      call set_primitive_formula(flist_a,idx_a,
     &     1d0,idx_opa,.true.,op_info)
      call init_formula(flist_b)
      call set_primitive_formula(flist_b,idx_b,
     &     1d0,idx_opb,.true.,op_info)

      call expand_subexpr(flist_adb,flist_a,.true.,op_info)
      call expand_subexpr(flist_adb,flist_b,.true.,op_info)

      if (ntest.ge.1000) then
        write(luout,*) 'after replacement'
        call print_form_list(luout,flist_adb,op_info)
      end if

      call init_formula(flist_a_b)
      ! (*)
      ! obtain the actual formula by taking the derivative wrt. dummy
      call form_deriv3(flist_a_b,flist_adb,.true.,
     &     1,idx_dum,0,idx_intm,op_info)

      if (ntest.ge.1000) then
        write(luout,*) 'raw formula'
        call print_form_list(luout,flist_a_b,op_info)
      end if
      
      ! tidy up
      call dealloc_formula_list(flist_adb)

      ! add AB contribution
      ! go to end of list
      flist_pnt => flist_a_b
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! add blocks of AB
      call set_primitive_formula(flist_pnt,idx_ab,
     &     1d0,idx_intm,.false.,op_info) 

      !end future subroutine

      ! further terms
      if (int_type(1:2).eq.'B '.and.approx(1:2).ne.'A ') then
        ! preliminarily we do that here directly:

        ! only A', assuming GBC, in fact ...
        if (approx(1:2).ne.'A''')
     &       call quit(1,'set_r12intm_cabs',
     &       'only A and A'' available so far ... ')        

        if (nlabels.lt.6)
     &       call quit(1,'set_r12intm_cabs',
     &       'too few labels for A'' (must be 6)')

        idx_x    = idx_oplist2(labels(5),op_info)
        if (idx_x.lt.0)
     &       call quit(1,'set_r12intm_cabs',
     &       'label not on list: '//labels(5))
        idx_h    = idx_oplist2(labels(6),op_info)
        if (idx_h.lt.0)
     &       call quit(1,'set_r12intm_cabs',
     &       'label not on list: '//labels(6))

        ! dummy operator: 1 particle part of H
        call add_operator(op_scr_f,op_info)
        idx_f = idx_oplist2(op_scr_f,op_info)
        opf_pnt => op_info%op_arr(idx_f)%op
        allocate(occ_def(ngastp,2,1))
        ndef = 1
        occ_def = 0
        occ_def(1,1:2,1) = 1
        call set_uop(opf_pnt,op_scr_f,.false.,
     &       occ_def,ndef,orb_info)

        ! now it becomes awkward ...
        ! ... I should reprogram expand_op_product to really
        ! handle operators and targets with njoined>1
        ! dummy operator: the upper and lower vertex of x
        call add_operator(op_scr_xu,op_info)
        idx_xu = idx_oplist2(op_scr_xu,op_info)
        opxu_pnt => op_info%op_arr(idx_xu)%op
        call add_operator(op_scr_xl,op_info)
        idx_xl = idx_oplist2(op_scr_xl,op_info)
        opxl_pnt => op_info%op_arr(idx_xl)%op

        occ_def = 0
        occ_def(IPART,1,1) = 1
        occ_def(IHOLE,2,1) = 2
        call set_uop(opxu_pnt,op_scr_xu,.true.,
     &       occ_def,ndef,orb_info)
        call set_uop(opxl_pnt,op_scr_xl,.false.,
     &       occ_def,ndef,orb_info)

        deallocate(occ_def)

        call init_formula(flist_adb)
        call new_formula_item(flist_adb,command_set_target_init,
     &       idx_op1)
        flist_pnt => flist_adb%next
        ! set up terms resulting from Xu D F Xl:
        call expand_op_product(flist_pnt,idx_op1,
     &       -1d0,4,(/idx_xu,idx_dum,idx_f,idx_xl/),
     &       (/-1,-1,-1,-1/),(/-1,-1,-1,-1/),
     &       (/1,4,1,2,2,3,3,4/),4,.false.,op_info)

        if (ntest.ge.1000) then
          write(luout,*) 'Xu D F Xl formula:'
          call print_form_list(luout,flist_adb,op_info)
        end if

        ! now, we have to replace Xu/Xl -> X and F -> H

        ! replace Xu Xl by X
        ! set up X as result from X_u and X_l
        call init_formula(flist_xdx)
        call new_formula_item(flist_xdx,command_set_target_init,
     &       idx_op1)
        flist_pnt => flist_xdx%next
        call expand_op_product(flist_pnt,idx_op1,
     &       1d0,3,(/idx_xu,idx_dum,idx_xl/),
     &       (/-1,-1,-1/),(/-1,-1,-1/),
     &       (/1,3,1,2,2,3/),3,.false.,op_info)
        call init_formula(flist)
        call form_deriv3(flist,flist_xdx,.true.,
     &     1,idx_dum,0,idx_x,op_info)
        call dealloc_formula_list(flist_xdx)

        if (ntest.ge.1000) then
          write(luout,*) 'Xu .. Xl formula:'
          call print_form_list(luout,flist,op_info)
        end if

        call factor_out_subexpr(flist_adb,flist,op_info)
        call dealloc_formula_list(flist)

        if (ntest.ge.1000) then
          write(luout,*) 'XF contribution after replacing XlXu:'
          call print_form_list(luout,flist_adb,op_info)
        end if

        ! and the F -> H replacement:
        call init_formula(flist)
        call set_primitive_formula(flist,idx_h,
     &       1d0,idx_f,.true.,op_info)
        call expand_subexpr(flist_adb,flist,.false.,op_info)
        call dealloc_formula_list(flist)

        if (ntest.ge.1000) then
          write(luout,*) 'XF contribution after replacing F->H:'
          call print_form_list(luout,flist_adb,op_info)
        end if

        ! remove dummy index, append thereby to flist_a_b
        flist_pnt => flist_a_b
        do while(flist_pnt%command.ne.command_end_of_formula)
          flist_pnt => flist_pnt%next
        end do
        call form_deriv3(flist_pnt,flist_adb,.false.,
     &       1,idx_dum,0,idx_intm,op_info)

        call dealloc_formula_list(flist_adb)

        call del_operator(op_scr_f,op_info)
        call del_operator(op_scr_xu,op_info)
        call del_operator(op_scr_xl,op_info)
        
      end if

      ! add symmetrisation command if required
      flist_pnt => flist_a_b
      if (symmetrise) then
        ! go to end of list
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do
        call new_formula_item(flist_pnt,command_symmetrise,idx_intm)
      end if

      ! write result to disc
      form_out%comment = trim(title)
      write(name,'(a,".fml")') trim(form_out%label)
      call file_init(form_out%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_out%fhand,flist_a_b,form_out%comment)

      if (ntest.ge.100) then
        write(luout,*) 'final formula'
        call print_form_list(luout,flist_a_b,op_info)
      end if

      ! tidy up
      call dealloc_formula_list(flist_a)
      call dealloc_formula_list(flist_b)
      call dealloc_formula_list(flist_a_b)

      call del_operator(op_scr_a,op_info)
      call del_operator(op_scr_b,op_info)
      call del_operator(op_scr_dum,op_info)
      call del_operator(op_scr_1,op_info)

      return

      contains 

      subroutine set_opshapes(occ_def,ndef,nblk,ansatz)

      implicit none

      integer, intent(in) ::
     &     ndef, ansatz, nblk
      integer, intent(inout) ::
     &     occ_def(ngastp,2,ndef)
      integer ::
     &     idx, ioff

      occ_def(1:ngastp,2,1:ndef*nblk) = 0

      do idx = 1, nblk
        ioff = (idx-1)*ndef

        if(idx.eq.1)then
          occ_def(1,2,1:ndef) = 2
        elseif(idx.eq.2)then
          occ_def(1,2,ioff+1:ioff+ndef) = 1
          occ_def(2,2,ioff+1:ioff+ndef) = 1
        elseif(idx.eq.3)then
          occ_def(2,2,ioff+1:ioff+ndef) = 2
        else
          call quit(1,'set_opshapes','Too many blocks')
        endif

        select case (ansatz)
        case(1)
          if (ndef.ne.5)
     &         call quit(1,'set_r12intm_cabs','inconsistency (1)')
          occ_def(1:ngastp,1,ioff+1) = (/2,0,0,0/)
          occ_def(1:ngastp,1,ioff+2) = (/1,1,0,0/)
          occ_def(1:ngastp,1,ioff+3) = (/1,0,0,1/)
          occ_def(1:ngastp,1,ioff+4) = (/0,2,0,0/)
          occ_def(1:ngastp,1,ioff+5) = (/0,1,0,1/)
        case(3)
          if (ndef.ne.4)
     &         call quit(1,'set_r12intm_cabs','inconsistency (2)')
          occ_def(1:ngastp,1,ioff+1) = (/2,0,0,0/)
          occ_def(1:ngastp,1,ioff+2) = (/1,1,0,0/)
          occ_def(1:ngastp,1,ioff+3) = (/1,0,0,1/)
          occ_def(1:ngastp,1,ioff+4) = (/0,2,0,0/)
        case default
          call quit(1,'set_r12intm_cabs','inconsistency (3)')
        end select

      enddo

      end subroutine

      end
