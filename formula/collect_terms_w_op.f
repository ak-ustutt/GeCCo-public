*----------------------------------------------------------------------*
      subroutine collect_terms_w_op(fpl_terms,form,
     &                              ncmpnd,idxop,iblkop,ntimes)
*----------------------------------------------------------------------*
*     collect all terms in formula list form which have
*       a) ntimes an entry with any of idxop(1:ncmpnd) (if iblkop.lt.0)
*       b) ntimes an entry with any of idxop(1:ncmpnd),iblkop(1:ncmpnd)
*     and store pointers to these terms on fpl_terms
*     exit, as soon as the next block occurs
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), target, intent(in) ::
     &     form
      type(formula_item_list), target, intent(inout) ::
     &     fpl_terms
      integer, intent(in) ::
     &     ncmpnd,
     &     idxop(ncmpnd), iblkop(ncmpnd), ntimes

      type(formula_item), pointer ::
     &     fl_pnt
      type(formula_item_list), pointer ::
     &     fpl_pnt
      logical ::
     &     ok
      integer ::
     &     nfound_op, nfound_blk, nvtx, ivtx, icmpnd, iblkop_cmpnd,
     &     nterms
      type(contraction), pointer ::
     &     cur_contr

      if (ntest.ge.100) then
        write(luout,*) '=============================='
        write(luout,*) ' info from collect_terms_w_op'
        write(luout,*) '=============================='
      end if

      fl_pnt => form
      fpl_pnt => fpl_terms
      if (associated(fpl_pnt%item)) then
        if (associated(fpl_pnt%next))
     &       call quit(1,'collect_terms_w_op',
     &       'formula pointer list is buggy')
      end if

      nterms = 0

      do
        if (fl_pnt%command.ne.command_add_contribution) exit

        cur_contr => fl_pnt%contr
        nvtx = cur_contr%nvtx
        nfound_op = 0
        nfound_blk = 0
        do ivtx = 1, nvtx
          do icmpnd = 1, ncmpnd
            if (cur_contr%vertex(ivtx)%idx_op.eq.
     &           abs(idxop(icmpnd)) .and.
     &         (cur_contr%vertex(ivtx)%dagger.eqv.
     &           idxop(icmpnd).lt.0)) then
              nfound_op = nfound_op+1
              iblkop_cmpnd = -1
              if (iblkop(1).gt.0) iblkop_cmpnd = iblkop(icmpnd)
              if (cur_contr%vertex(ivtx)%iblk_op.eq.iblkop_cmpnd)
     &             nfound_blk = nfound_blk+1
            end if
          end do
        end do
        if (iblkop(1).lt.0) then
          ok = nfound_op.eq.ntimes
        else
          ok = nfound_blk.eq.ntimes
        end if
c dbg
        print *,'ok = ',ok,nfound_op,nfound_blk,ntimes
c dbg
        if (ok) then
          if (.not.associated(fpl_pnt%item)) then
            fpl_pnt%item => fl_pnt
c dbg
c            print *,'(a) command = ',fpl_pnt%item%command
c dbg
          else
            allocate(fpl_pnt%next)
            fpl_pnt%next%prev => fpl_pnt
            fpl_pnt%next%next => null()
            fpl_pnt => fpl_pnt%next
            fpl_pnt%item => fl_pnt
c dbg
c            print *,'(b) command = ',fpl_pnt%item%command
c dbg
            nterms = nterms+1
          end if
        end if
        
        if (.not.associated(fl_pnt%next)) exit
        fl_pnt => fl_pnt%next
      end do
c dbg
      print *,'nterms = ',nterms
c dbg

      return
      end
