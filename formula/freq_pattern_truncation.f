*----------------------------------------------------------------------*
      subroutine freq_pattern_truncation(flist,order,freq_idx,ncomps,
     &                                   pop_idx,idx_tgt,op_info)
*----------------------------------------------------------------------*
*     truncate formula to terms of specific frequency pattern
*     matthias, 2008
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'def_del_list.h'
      include 'par_opnames_gen.h'
      include 'def_pert_info.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     order, idx_tgt, freq_idx(*), ncomps, pop_idx(*)

      logical ::
     &     delete, recognized, multiply
      integer ::
     &     nvtx, ivtx, op_ord, idx_op, t_max_ord, l_max_ord, op_spec,
     &     count_freq(maxcmp), ii, jj, op_ifreq,
     &     njoker(3*maxpop), pattern(maxcmp), factor

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, allocatable ::
     &     structure(:,:)

      integer, external ::
     &     ifac

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'freq_pattern_truncation')
      endif

      ! determine frequency index pattern
      if (maxval(freq_idx(1:order)).gt.maxcmp)
     &   call quit(1,'freq_pattern_truncation',
     &       'increase maxcmp or request fewer pert.op.components')
      pattern = 0
      do ii = 1,order
        pattern(freq_idx(ii)) = pattern(freq_idx(ii)) + 1
      end do

      if (ntest.ge.100) then
         write(luout,*)'pattern: ',
     &                     pattern(1:maxval(freq_idx(1:order)))
         write(luout,*)'pop_idx: ',pop_idx(1:maxval(freq_idx(1:order)))
      end if

      form_pnt => flist
      do 
        form_pnt_next => form_pnt%next
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(luout,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(luout,*) '[INIT_TARGET]'
          form_pnt%target = idx_tgt
        case(command_add_contribution)

          form_pnt%target = idx_tgt
          ! assign new result index
          ! comment: operator block should also be changed
          form_pnt%contr%idx_res = idx_tgt

          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex
          allocate(structure(maxcmp,nvtx))

          ! delete term if frequency pattern does not match freq_idx
          njoker = 0
          count_freq = 0
          structure = 0
          delete = .false.
          if (ntest.ge.100) write(luout,*)'checking ',nvtx,' vertices'
          do ivtx = 1, nvtx
c
c            print *,'ivtx = ',ivtx
c
            idx_op  = vertex(ivtx)%idx_op
            op_ord = op_info%op_arr(idx_op)%op%order
c           
c            print *,'name: ',trim(op_info%op_arr(idx_op)%op%name)
c            print *,'op_ord = ',op_ord
c
            if ((associated(op_info%op_arr(idx_op)%op%ifreq)).and.
     &          (op_ord.gt.0)) then
              op_spec = op_info%op_arr(idx_op)%op%species
c           
c            print *,'op_spec = ',op_spec
c            print *,'op_ifreq(1) = ',op_info%op_arr(idx_op)%op%ifreq(1)
c
              if (op_info%op_arr(idx_op)%op%ifreq(1).gt.0) then
                do ii = 1,op_ord
                  op_ifreq = op_info%op_arr(idx_op)%op%ifreq(ii)
                  structure(op_ifreq,ivtx) = structure(op_ifreq,ivtx)+1
                  recognized = .false.
                  do jj = 1,order
                    if (op_ifreq.eq.freq_idx(jj)) then
                      recognized = .true.
                    end if
                  end do
                  count_freq(op_ifreq) = count_freq(op_ifreq) + 1
                  if (.not.recognized) delete = .true.
                  if ((.not.recognized).and.ntest.ge.100)
     &                 write(luout,*)'non-matching vertex no.',ivtx
                end do
              else if (op_spec.ge.4) then
                if (op_spec.gt.3*maxpop+3)
     &                call quit(1,'freq_pattern_truncation',
     &                            'maxpop too small')
                ! can belong to any frequency index associated to pert op
                njoker(op_spec-3) = njoker(op_spec-3) + 1
                if (ntest.ge.100) write(luout,*)'collected joker assoc',
     &                             'iated with pert. op no. ',op_spec-3
              end if
            end if
          end do
          ! play all jokers and check if frequency pattern is ok
          ii = 0
          if (ntest.ge.100) then
            write(luout,*) 'current jokers: ', njoker
            write(luout,*) 'current counts: ', 
     &                     count_freq(1:maxval(freq_idx(1:order)))
          end if
          do while (maxval(njoker).gt.0.and.
     &              ii.lt.maxval(freq_idx(1:order)))
            ii = ii+1
            do jj = 1, nvtx
              if (sum(structure(:,jj)).eq.0) ivtx = jj
            end do
            if (count_freq(ii).lt.pattern(ii).and.
     &          njoker(pop_idx(ii)).gt.0) then
              count_freq(ii) = count_freq(ii) + 1
              njoker(pop_idx(ii)) = njoker(pop_idx(ii)) - 1
              structure(ii,ivtx) = structure(ii,ivtx)+1
            end if
          end do
          if (.not.(all(count_freq-pattern.eq.0).and.
     &        all(njoker.eq.0))) delete = .true.

          if (ntest.ge.100) write(luout,*) 'delete node?: ',delete
          if (delete) then
            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          else
            ! multiply with factor to correct for omitted ("redundant") terms
            factor = 1
            do ii = 1,maxcmp
              factor = factor * ifac(sum(structure(ii,:)))
              do jj = 1,nvtx
                factor = factor / ifac(structure(ii,jj))
              end do
            end do
            form_pnt%contr%fac = form_pnt%contr%fac * factor
          end if
        
          deallocate(structure)

        case default
          write(luout,*)'command = ',form_pnt%command
          call quit(1,'delete_non_fact','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
