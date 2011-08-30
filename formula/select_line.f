*----------------------------------------------------------------------*
      subroutine select_line(flist,idxres,idxop,ncmpnd,igastp,mode_str)
*----------------------------------------------------------------------*
*     collect those contractions on flist in which two operators
*     on the list idxop(1:ncmpnd) are contracted via a line given
*     by igastp.
*     mode_str = 'delete': delete all these terms
*     mode_str = 'keep'  : delete all other terms
*     mode_str = 'no_ext': delete terms where a listed op has external
*                          lines in the given space
*     mode_str = 'ext'   : keep only terms where some listed op has
*                          external lines in the given space
*
*     matthias, nov 2009 (adopted from form_invariant)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(inout), target ::
     &     flist
      integer, intent(in) ::
     &     ncmpnd, igastp, idxres, idxop(ncmpnd)
      character(*), intent(in) ::
     &     mode_str
      
      logical ::
     &     ok, del_mode
      integer ::
     &     nterms, iarc, idx, icmpnd
      logical, allocatable ::
     &     vtxinlist(:)

      type(formula_item), pointer ::
     &     fl_pnt, fl_pnt_next
      type(contraction), pointer ::
     &     contr

      if (trim(mode_str).eq.'delete'.or.
     &    trim(mode_str).eq.'no_ext') then
        del_mode = .true.
      else if (trim(mode_str).eq.'keep'.or.
     &         trim(mode_str).eq.'ext') then
        del_mode = .false.
      else
        call quit(1,'form_select_line',
     &            'mode must be "delete", "keep", "no_ext" or "ext"')
      end if

      fl_pnt => flist

      if (.not.associated(fl_pnt))
     &     call quit(1,'form_select_line',
     &     'empty formula list? something is buggy')

      nterms = 0
      fl_loop: do

        fl_pnt_next => fl_pnt%next
        if (fl_pnt%command.eq.command_end_of_formula)
     &       exit fl_loop

        fl_pnt%target = idxres

        if (fl_pnt%command.eq.command_add_contribution) then

          ok = .true.
          contr => fl_pnt%contr

          ! find which vertex belongs to an operator from the given list
          allocate(vtxinlist(contr%nvtx))
          vtxinlist = .false.
          cmp_loop: do idx = 1, contr%nvtx
            do icmpnd = 1, ncmpnd
              if (contr%vertex(idx)%idx_op.eq.abs(idxop(icmpnd)).and.
     &           (contr%vertex(idx)%dagger.eqv.(idxop(icmpnd).lt.0)))
     &        then
                vtxinlist(idx) = .true.
                cycle cmp_loop
              end if
            end do
          end do cmp_loop

          if (trim(mode_str).eq.'delete'
     &        .or.trim(mode_str).eq.'keep') then
            ! is there any contraction between two of these operators
            ! via a line of type igastp?
            do iarc = 1, contr%narc
              if (vtxinlist(contr%arc(iarc)%link(1)).and.
     &            vtxinlist(contr%arc(iarc)%link(2)).and.
     &            sum(contr%arc(iarc)%occ_cnt(igastp,1:2)).gt.0) then
                ok = .false.
                exit
              end if
            end do
          else
            ! is there any open line of type igastp from a listed op?
            do iarc = 1, contr%nxarc
              if (vtxinlist(contr%xarc(iarc)%link(1)).and.
     &            sum(contr%xarc(iarc)%occ_cnt(igastp,1:2)).gt.0) then
                ok = .false.
                exit
              end if
            end do
          end if

          deallocate(vtxinlist)

          if (ok.eqv.del_mode) then
            nterms = nterms+1
            contr%idx_res = idxres ! not completely OK
          else
            ! deallocate contents and re-link the list
            call delete_fl_node(fl_pnt)
            ! deallocate the node itself
            deallocate(fl_pnt)
          end if
        else 
        end if

        fl_pnt => fl_pnt_next
        if (.not.associated(fl_pnt))
     &       call quit(1,'form_select_line',
     &       'unexpected end of formula list')

      end do fl_loop

      if (ntest.ge.10) then
        write(luout,*) 'generated terms: ',nterms
      end if

      return
      end
