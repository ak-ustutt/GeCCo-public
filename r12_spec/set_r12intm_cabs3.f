*----------------------------------------------------------------------*
      subroutine set_r12intm_cabs3(form_out,
     &     title,labels,nlabels,int_type,ansatz,approx,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*    driver for setting up all kinds of R12 intermediates
*    for a given approximation
*     
*    andreas, march 2008
*
*     (based on set_v_intermediate() etc. by GWR)
*
*     expected operator labels:
*
*     label(1) -- the intermediate to be generated
*
*     label     Rbar  Rtilde Rbar+  V      V+     X     
*    ------------------------------------------------------
*       2       R12   R12    R12    G_X    R12    R12   
*       3       R12X  R12X   R12X   R12    G_X    R12   
*       4       F+K   K      F+K    H.R12  R12.H  R12^2 
*    -------------------------------------------------------
*
*     label     B (apr. A,B)  B (apr. C)  C (apr. A,B)  C (apr. C)
*    --------------------------------------------------------------
*       2       R12           R12          R12          R12  or R12C
*       3       [T12,R12]      --          H            H
*       4       R12[T12,R12]  R12[T12,R12] [T12,R12]    --
*     .............................................................
*       5       X              --          RBAR
*       6       H              --          RTILDE
*       7       RBAR          RBAR+
*       8       RTILDE        RTILDE
*       9       R12^2         R12^2 / or {R12^2}BREVE
*      10       F+K           F+K     or  --
*      11       --            RBREVE
*      12       C             C
*    --------------------------------------------------------------
*
*     approx string: formatted string
*
*     123456789012
*     1  2   3   4
*     XX XXX XXX X
*     A
*     A' GBC HY1 S
*     B  EBC HY2 
*     C  noZ HX
*
*     1: approximation type
*     2: assumed Brillouin condition (EBC is not supported)
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
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 100

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

      logical ::
     &     symmetrise
      integer ::
     &     idx_intm, idx_op(nlabels), nop, iop, iop1, iop2, nder, nopen
      integer, pointer ::
     &     occ_def(:,:,:)
      character ::
     &     name*(form_maxlen_label*2)
      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     flist_pnt
      type(formula) ::
     &     form_inp

      integer, external ::
     &     idx_oplist2, idx_formlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_r12_intm_cabs')
        write(luout,*) 'setting: ',trim(labels(1))
        write(luout,*) 'type = ',trim(int_type)
        write(luout,*) 'ansatz = ',ansatz
        write(luout,*) 'approx = "',approx(1:12),'"'
        write(luout,*) 'nlabels = ',nlabels
        do iop = 1, nlabels
          write(luout,*) '       ',trim(labels(iop))          
        end do
      end if

      symmetrise = approx(12:12).eq.'S'

      if (ntest.ge.100) write(luout,*) 'symmetrise = ',symmetrise

      nop = nlabels

      idx_op(1:nop) = -1
      do iop = 1, nop
        ! ignore empty labels
        if (ntest.ge.100)
     &       write(luout,*) 'checking label #',iop,
     &       ': ',trim(labels(iop))
        if (trim(labels(iop)).eq.'-' .or.
     &      len_trim(labels(iop)).eq.0 ) cycle
        idx_op(iop) = idx_oplist2(labels(iop),op_info)
        if (idx_op(iop).lt.0)
     &       call quit(1,'set_r12intm_cabs',
     &       'label not on list: '//labels(iop))
      end do

      if (ntest.ge.100)
     &     write(luout,*) 'idx_op: ',idx_op
      
      idx_intm = idx_op(1)

      ! initialize formula
      call init_formula(flist)
      call new_formula_item(flist,command_set_target_init,
     &     idx_intm)

      select case(trim(int_type))
      case('RB','RT','RV') ! NOT USED, NOT DEBUGGED
        call quit(1,'set_r12intm_cabs3','not debugged')
        call set_r12mod(flist,int_type,
     &       idx_intm,idx_op,nop,op_info)
      case('C')
        if (ansatz.eq.1)
     &       call quit(1,'set_r12intm_cabs3','no C for ansatz 1')
        if (approx(1:1).ne.'C') then
          call set_1contrib(flist,1d0,4,
     &       idx_intm,idx_op,nop,op_info)
          call set_C_fr(flist,approx,
     &         3,2,5,6,
     &         idx_intm,idx_op,nop,op_info,orb_info)
        else
c          if (nop.eq.2) then
c            call set_C_fr(flist,approx,
c     &         -1,2,-1,-1,
c     &         idx_intm,idx_op,nop,op_info,orb_info)
c          else
            call set_C_fr(flist,approx,
     &         3,2,-1,-1,
     &         idx_intm,idx_op,nop,op_info,orb_info)
c          end if
        end if
      case('V','V+','X')
        ! set up term arising from 1 in Q = 1 - P
        call set_1contrib(flist,1d0,4,
     &       idx_intm,idx_op,nop,op_info)
        ! set up term arising from P in Q = 1 - P
        call set_Pcontrib(flist,ansatz,
     &       2,3,
     &       idx_intm,idx_op,nop,op_info)
      case('BH')
        call set_Bhole(flist,ansatz,
     &       4,2,3,6,
     &       idx_intm,idx_op,nop,op_info,orb_info)
      case('B','BP')
        ! set up term arising from 1 in Q = 1 - P
        call set_1contrib(flist,1d0,4,
     &       idx_intm,idx_op,nop,op_info)
        if (approx(1:1).ne.'C') then
          ! set up term arising from P in Q = 1 - P
          call set_Pcontrib(flist,ansatz,
     &         2,3,
     &         idx_intm,idx_op,nop,op_info)
          ! Hartree contributions
          call set_Xcontrib(flist,ansatz,approx,
     &         9,2,7, 10, 5,6,
     &         idx_intm,idx_op,nop,op_info,orb_info)
          if (approx(8:10).ne.'HY1') then
            ! Exchange contributions
            call set_Ycontrib(flist,ansatz,approx,
     &         2,8,
     &         idx_intm,idx_op,nop,op_info)
          end if
        else
          ! Hartree contributions
          call set_Xcontrib(flist,ansatz,approx,
     &         9,2,7, 10, 5,6,
     &         idx_intm,idx_op,nop,op_info,orb_info)
          ! Exchange contributions
          if (approx(8:10).ne.'HY1') then
            call set_Ycontrib(flist,ansatz,approx,
     &         2,8,
     &         idx_intm,idx_op,nop,op_info)
          end if
        end if
        if (approx(4:6).ne.'noZ'.and.
     &      approx(4:6).ne.'GBC'.and.
     &      approx(4:6).ne.'EBC') then
          call set_Zcontrib(flist,ansatz,approx,
     &       2,11,
     &       idx_intm,idx_op,nop,op_info)
        end if
        if (ansatz.gt.1) then
          call set_RC_contrib(flist,ansatz,approx,
     &       2,12,
     &       idx_intm,idx_op,nop,op_info)
        end if

      case('PG')

        call set_p3g_contract(form_out,title,
     &       labels(1),labels(2:3),2,
     &       op_info,orb_info)

c dbg
c        stop
c dbg
      end select

      ! add symmetrisation command if required
      flist_pnt => flist
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
      call write_form_list(form_out%fhand,flist,form_out%comment)

      if (ntest.ge.100) then
        write(luout,*) 'CABS final formula'
        call print_form_list(luout,flist,op_info)
      end if

      ! tidy up
      call dealloc_formula_list(flist)

      return

      end
