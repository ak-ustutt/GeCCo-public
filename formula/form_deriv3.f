*----------------------------------------------------------------------*
      subroutine form_deriv3(fl_deriv,fl_input,init,
     &                      ncmpnd,idxder,idxmlt,idxres,
     &                      op_info)
*----------------------------------------------------------------------*
*     get the derivatives of all contractions on fl_input

*     ncmpnd: number of compounds in multi-compound operator -
*       e.g. for R12 : 'T' and 'C' are a compound op. (ncmpnd=2) 
*     idxder(1:ncmpnd) is the index of the operator with respect to which the
*     derivative has to be taken
*     idxmlt(1:ncmpnd) is the index of the operator which the derivative is 
*     multiplied with (0 if only the derivative is taken)
*     idxres is the index of the resulting operator (0, if scalar)     
*
*     andreas, end of 2006, put to gecco april 2007
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_contraction_list.h'
      include 'def_formula_item.h'

      type(formula_item), intent(inout), target ::
     &     fl_input, fl_deriv
      logical, intent(in) ::
     &     init
      integer, intent(in) ::
     &     ncmpnd, idxder(ncmpnd), idxmlt(ncmpnd), idxres
      type(operator_info) ::
     &     op_info

      type(contraction) ::
     &     contr
      type(contraction_list), target ::
     &     conder
      type(contraction_list), pointer ::
     &     cur_conder

      logical ::
     &     reo
      integer ::
     &     nvtx,
     &     nterms, nder, idx, idum, idxinp, len, icmpnd, ieqvfac

      integer, pointer ::
     &     ivtx_reo(:),occ_vtx(:,:,:)
      logical, pointer ::
     &     fix_vtx(:)
      type(formula_item), pointer ::
     &     fl_input_pnt, fl_deriv_pnt

      if (fl_input%command.ne.command_set_target_init)
     &       call quit(1,'form_deriv',
     &       'input formula definition must start with [INIT]')
      
      if (init) then
        call init_formula(fl_deriv)
        call new_formula_item(fl_deriv,command_set_target_init,idxres)
        fl_deriv_pnt => fl_deriv%next
      else
        if (fl_deriv%command.ne.command_end_of_formula)
     &      call quit(1,'form_deriv3',
     &       'init==.false., but output formula is not properly '//
     &       'initialized')
        fl_deriv_pnt => fl_deriv
      end if

      fl_input_pnt => fl_input%next

      nullify(conder%contr)
      
      ! loop over input formula list
      nterms = 0
      input_loop: do

        if (fl_input_pnt%command.ne.command_add_contribution)
     &     exit input_loop

        ! obtain derivatives for all components:
        do icmpnd = 1, ncmpnd
          call contr_deriv2(conder,nder,fl_input_pnt%contr,op_info,
     &         idxder(icmpnd),idxmlt(icmpnd),idxres)

          cur_conder => conder
          do idx = 1, nder
            nterms = nterms+1

            ! enforce law and order
            nvtx = cur_conder%contr%nvtx
            allocate(ivtx_reo(nvtx),fix_vtx(nvtx),
     &           occ_vtx(ngastp,2,nvtx))
            fix_vtx = .true.    ! "fix" all vertices -> ieqvfac will be 1
            call occvtx4contr(1,occ_vtx,cur_conder%contr,op_info)

            call topo_contr(ieqvfac,reo,ivtx_reo,
     &           cur_conder%contr,occ_vtx,fix_vtx)
            ! ieqvfac is ignored
            call canon_contr(cur_conder%contr,reo,ivtx_reo)
            deallocate(ivtx_reo,fix_vtx,occ_vtx)

            ! put result to derivative formula list
            call new_formula_item(fl_deriv_pnt,
     &           command_add_contribution,idxres)
            call copy_contr(cur_conder%contr,fl_deriv_pnt%contr)
            fl_deriv_pnt => fl_deriv_pnt%next
            if (idx.lt.nder) cur_conder => cur_conder%next
          end do

          call dealloc_contr_list(conder)
        end do

        if (.not.associated(fl_input_pnt%next))
     &       call quit(1,'form_deriv',
     &       'unexpected end of list (input)')
        fl_input_pnt => fl_input_pnt%next

      end do input_loop

      call dealloc_contr(contr)
      
      return
      end
