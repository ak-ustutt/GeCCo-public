*----------------------------------------------------------------------*
      subroutine set_Xcontrib(flist,ansatz,approx,
     &     irsq,irdag,irbar,ihartree,ixmat,iham,
     &     idx_intm,idx_op,nop,op_info,orb_info)
*----------------------------------------------------------------------*
*     set X type contributions to B matrix
*     no GBC assumed:
*      B  : assemble Xbar from r^2(f+k) and r+ and rbar
*      C  :                    r^2(f+k) and r+ and rbar+
*     GBC assumed (A' and B only)
*       assemble Xbar_F from X and Fock matrix
*       if approx B: assemble Xbar_K from r+ and rbar (contains K only) 
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
     &     ntest = 00
      character(6), parameter ::
     &     op_scr_f  = '_SCR_F'

      type(formula_item), intent(inout), target ::
     &     flist
      character*(*) ::
     &     approx
      integer, intent(in) ::
     &     ansatz,nop,irsq,irdag,irbar,ixmat,ihartree,iham,
     &     idx_intm,idx_op(nop)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(formula_item), pointer ::
     &     flist_pnt, flist_pnt0

      logical ::
     &     GBC, assume_rsqbar
      integer ::
     &     idx_1, idx_2, idx_prj, idx_f, ndef,
     &     njoined_intm, njoined_op
      type(operator), pointer ::
     &     opf_pnt, op_pnt
      integer ::
     &     occ_def(ngastp,2)

      integer, external ::
     &     idx_oplist2

      if (approx(1:2).eq.'A ') return

      njoined_intm = op_info%op_arr(idx_intm)%op%njoined

      GBC = approx(4:6).eq.'GBC'.or.approx(1:2).eq.'A'''

      if (GBC.and.approx(1:2).eq.'C ')
     &     call quit(1,'set_Xcontrib','no GBC for approximation C')


      if (approx(1:1).ne.'A') then

        if (irdag.gt.nop.or.
     &      irbar.gt.nop.or.
     &      irsq.gt.nop.or.
     &      ihartree.gt.nop) then
          write(luout,*) 'idx: ',irdag,irbar,irsq,ihartree
          write(luout,*) 'nop: ',nop
          call quit(1,'set_Xcontrib',
     &         'not enough operators on input list')
        end if

        flist_pnt => flist
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do
        flist_pnt0 => flist_pnt

        ! if no hartree operator was given, we assume that
        ! the modified R.R integral was imported
        assume_rsqbar = idx_op(ihartree).le.0

        if (idx_op(irdag).le.0.or.
     &      idx_op(irbar).le.0.or.
     &      idx_op(irsq) .le.0) then
          write(luout,*) 'idx: ',idx_op(irdag),idx_op(irbar),
     &         idx_op(irsq),idx_op(ihartree)
          call quit(1,'set_Xcontrib',
     &         'operator(s) not on input list')
        end if

        !-----------------!
        ! R12^{2} * (f+k) !
        !-----------------!
        if (.not.assume_rsqbar) then

c          call quit(1,'set_Xcontrib',
c     &         'unused route. how did you get here?')

          idx_1 = idx_op(irsq)
          idx_2 = idx_op(ihartree)

          ! go to end of list
          flist_pnt => flist
          do while(associated(flist_pnt%next))
            flist_pnt => flist_pnt%next
          end do
          if (njoined_intm.eq.1) then
            call expand_op_product2(flist_pnt,idx_intm,
     &       1d0,5,3,
     &       (/idx_intm,idx_2,idx_1,idx_1,idx_intm/),
     &       (/1       ,2    ,3    ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       0,0, 
     &       (/2,3,1,0/),1,     ! def. of projector
     &       .false.,op_info)
          else if (njoined_intm.eq.2) then
            call expand_op_product2(flist_pnt,idx_intm,
     &       1d0,7,3,
     &       (/idx_intm,idx_1,idx_intm,idx_intm,idx_2,idx_1,idx_intm/),
     &       (/1       ,2    ,1       ,1       ,3    ,2    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       (/2,7, 1,5, 1,6/),3, ! avoid cross contrib. to external lines
     &       (/4,5,1,0/),1,     ! def. of projector
     &       .false.,op_info)
          else
            call quit(1,'set_Xcontrib','not prepared for this case !')
          end if
        else
          idx_1 = idx_op(irsq)
          flist_pnt => flist
          do while(associated(flist_pnt%next))
            flist_pnt => flist_pnt%next
          end do
          njoined_op   = op_info%op_arr(idx_1)%op%njoined
          if (njoined_intm.eq.njoined_op) then
            call set_primitive_formula(flist_pnt,idx_1,
     &         1d0,idx_intm,.false.,op_info) 
          else if (njoined_intm.lt.njoined_op) then
            ! generate the appropriate self-contraction
            call expand_op_product2(flist_pnt,idx_intm,
     &           1d0,4,2,
     &           (/idx_intm,idx_1,idx_1,idx_intm/),
     &           (/1       ,2   , 2, 1       /),       
     &           -1, -1,
     &           0,0,
     &           0,0,
     &           0,0,
     &           .false.,op_info)
          else
            call quit(1,'set_Xcontrib','not prepared for this case !')
          end if
        end if

        !----------------------------------!
        ! - R12+ RBAR   (C: - R12+ RBREVE) !  
        !----------------------------------!
        idx_1 = idx_op(irdag)
        idx_2 = idx_op(irbar)

        ! go to end of list
        flist_pnt => flist
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do
        idx_prj = 2
        if (ansatz.gt.1) idx_prj = 4
        if (njoined_intm.eq.1) then
c          call expand_op_product2(flist_pnt,idx_intm,
c     &       -1d0,4,3,
c     &       (/idx_intm,-idx_1,idx_2,idx_intm/),
c     &       (/1       ,2     ,3    ,1       /),       
c     &       -1, -1,
c     &       (/2,3/),1, ! force additional contraction
c     &       0,0,
c     &       (/2,3,2,idx_prj/),1, ! def. of projector
c     &       .false.,op_info)
c          ! this gives the terms with ONLY the projector as contraction
c          do while(associated(flist_pnt%next))
c            flist_pnt => flist_pnt%next
c          end do
          call expand_op_product2(flist_pnt,idx_intm,
     &       -1d0,4,3,
     &       (/idx_intm,-idx_1,idx_2,idx_intm/),
     &       (/1       ,2     ,3    ,1       /),       
     &       -1, -1,
     &       0,0,      
     &       0,0,
     &       (/2,3,2,idx_prj/),1, ! def. of projector
     &       .false.,op_info)
        else if (njoined_intm.eq.2) then
          call expand_op_product2(flist_pnt,idx_intm,
     &       -1d0,6,3,
     &       (/idx_intm,-idx_1,idx_intm,idx_intm,idx_2,idx_intm/),
     &       (/1       ,2     ,1       ,1       ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       (/2,6, 1,5/),2,    ! avoid cross contrib. to external lines
     &       (/2,5,2,idx_prj/),1, ! def. of projector
     &       .false.,op_info)
        else
          call quit(1,'set_Xcontrib','unexpected: njoined_intm>2')
        end if

        if (ntest.ge.100) then
          write(luout,*) 'Xbar contribution'
          call print_form_list(luout,flist_pnt,op_info)
        end if

      end if

      if (GBC) then

        if (iham.gt.nop.or.
     &      ixmat.gt.nop) then
          write(luout,*) 'idx: ',iham,ixmat
          write(luout,*) 'nop: ',nop
          call quit(1,'set_Xcontrib',
     &         'not enough operators on input list')
        end if

        if (idx_op(iham).le.0.or.
     &      idx_op(ixmat).le.0) then
          write(luout,*) 'idx: ',idx_op(iham),idx_op(ixmat)
          call quit(1,'set_Xcontrib',
     &         'operator(s) not on input list')
        end if

        ! dummy operator: 1 particle part of H
        call add_operator(op_scr_f,op_info)
        idx_f = idx_oplist2(op_scr_f,op_info)
        opf_pnt => op_info%op_arr(idx_f)%op
        ndef = 1
        occ_def = 0
        occ_def(1,1:2) = 1
        call set_uop(opf_pnt,op_scr_f,.false.,
     &       occ_def,ndef,orb_info)

        idx_1 = idx_op(ixmat)
        idx_2 = idx_op(iham)

        ! go to end of list
        flist_pnt => flist
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do
        idx_prj = ansatz*2
        call expand_op_product2(flist_pnt,idx_intm,
     &       -1d0,7,3,
     &       (/idx_intm,idx_1,idx_intm,idx_intm,idx_f,idx_1,idx_intm/),
     &       (/1       ,2    ,1       ,1       ,3    ,2    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       (/2,7, 1,5, 1,6/),3,  ! avoid cross contrib. to external lines
     &       0,0,
     &       .false.,op_info)

        ! F -> H replacement:
        op_pnt => op_info%op_arr(idx_2)%op
        call form_op_replace(opf_pnt%name,op_pnt%name,
     &       .true.,flist_pnt,op_info)

        if (ntest.ge.100) then
          write(luout,*) 'XF contribution:'
          call print_form_list(luout,flist_pnt0,op_info)
        end if

        call del_operator(op_scr_f,op_info)

      end if

      return
      end

      
