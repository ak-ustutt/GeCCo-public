*----------------------------------------------------------------------*
      subroutine set_hhat(form_hhat,
     &     nops,ops,idxhhat,idxham,idxtop)
*----------------------------------------------------------------------*
*     generate the formula for Hhat = e^{-T1}He^{T1}
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_formula.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'
      include 'ifc_baserout.h'
      include 'par_formnames_gen.h'

      integer, parameter ::
     &     ntest = 00

      type(formula), intent(inout) ::
     &     form_hhat
      integer, intent(in) ::
     &     nops, idxhhat, idxham, idxtop
      type(operator), intent(in) ::
     &     ops(nops)

      type(contraction) ::
     &     contr

      logical ::
     &     init
      integer ::
     &     luhhat,
     &     iblk_hhat, iblk_ham, idx, nt1, iblk_t1,
     &     maxvtx, maxarc
      integer ::
     &     icnt(2), icnt_part(2,4), t1_part(4), icnt_part_idx(4),
     &     iocc_hhat(ngastp,2), iocc_ham(ngastp,2),
     &     iocc_hhat_d(ngastp,2), iocc_ham_d(ngastp,2),
     &     iocc_cnt(ngastp,2), iocc_t1n(ngastp,2), iocc_scr(ngastp,2)
      character ::
     &     name*(form_maxlen_label*2)

      integer, external ::
     &     ieqfac, iblk_occ

      if (ntest.ge.100) then
        write(luout,*) '================'
        write(luout,*) 'this is set_hhat'
        write(luout,*) '================'
      end if

      if (iprlvl.ge.2)
     &     write(luout,*) 'generating e^{-T1}He^{T1}'

      maxvtx = 5
      maxarc = 4
      allocate(contr%vertex(maxvtx),contr%arc(maxarc))

      form_hhat%label = label_cchhat
      form_hhat%comment = title_cchhat

      write(name,'(a,".fml")') label_cchhat
      call file_init(form_hhat%fhand,name,ftyp_sq_unf,0)

      call file_open(form_hhat%fhand)
      luhhat = form_hhat%fhand%unit
      rewind luhhat

      write(luhhat) len_trim(title_cchhat),title_cchhat
      write(luhhat) 0,idxhhat

      ! put T1 occupation on iocc_nt1
      iocc_t1n(ipart,1) = 1
      iocc_t1n(ihole,2) = 1
      iblk_t1 = iblk_occ(iocc_t1n,.false.,ops(idxtop))
      if (iblk_t1.le.0)
     &     call quit(1,'set_hhat','T1 not on list, op = '//
     &     ops(idxtop)%name)

      contr%idx_res = idxhhat

      ! loop over occupation classes of Hhat
      do iblk_hhat = 1, ops(idxhhat)%n_occ_cls 
        contr%iblk_res = iblk_hhat
        ! get deexcitation part of Hhat
        iocc_hhat = ops(idxhhat)%ihpvca_occ(1:ngastp,1:2,iblk_hhat)
        iocc_hhat_d = iocc_xdn(2,iocc_hhat) 
        ! loop over occupation classes of H
        do iblk_ham = 1, ops(idxham)%n_occ_cls 
          iocc_ham = ops(idxham)%ihpvca_occ(1:ngastp,1:2,iblk_ham)
          
          ! iocc_hhat == iocc_hat -> just add this contrib
          if (iocc_equal(iocc_hhat,.false.,iocc_ham,.false.)) then
            contr%fac = 1d0
            contr%nvtx = 1
            contr%vertex(1)%idx_op = idxham
            contr%vertex(1)%iblk_op = iblk_ham
            contr%narc = 0
            if (ntest.ge.100) then
              call prt_contr(luout,contr,ops)
            end if
            call wrt_contr(luhhat,contr)
          else
            ! else: get deexcitation part of H
            iocc_ham_d = iocc_xdn(2,iocc_ham)
c dbg
c            print *,'deexcitation parts'
c            call wrt_occ(luout,iocc_hhat_d)
c            call wrt_occ(luout,iocc_ham_d)
c            call wrt_occ(luout,iocc_overlap(iocc_ham_d,.false.,
c     &                                    iocc_hhat_d,.false.))
c dbg
            ! deexcitation part of Hhat must be fully included in that of H:
            if (.not.iocc_equal(iocc_hhat_d,.false.,
     &                       iocc_overlap(iocc_ham_d,.false.,
     &                                    iocc_hhat_d,.false.),.false.))
     &         cycle

            ! remove dexcitation part of Hhat -> must be contraction
            iocc_cnt = iocc_add(1,iocc_ham_d,.false.,
     &                       -1,iocc_hhat_d,.false.)

            ! contraction must be non-zero
            if (ielsqsum(iocc_cnt,2*ngastp).eq.0) cycle

            ! [Hhat] = [H] + [T1^n] - [C] - [C^+] -->
            ! [T1^n] = [Hhat] - [H] + [C] + [C^+]
            iocc_t1n = iocc_add(1,iocc_hhat,.false.,-1,iocc_ham,.false.)
            iocc_scr = iocc_add(1,iocc_t1n,.false.,1,iocc_cnt,.false.)
            iocc_t1n = iocc_add(1,iocc_scr,.false.,1,iocc_cnt,.true.)
            nt1 = iocc_t1n(ipart,1)
            if (.not.(nt1.gt.0.and.
     &                iocc_t1n(ihole,2).eq.nt1.and.
     &                iocc_t1n(ipart,2).eq.0  .and.
     &                iocc_t1n(ihole,1).eq.0)) cycle

            ! the contraction must be contained in [(T1^n)^+]
            if (.not.iocc_equal(iocc_cnt,.false.,
     &                          iocc_overlap(iocc_cnt,.false.,
     &                                       iocc_t1n,.true.),.false.))
     &         cycle

            if (nt1.eq.1) then
              ! one commutator
              contr%fac = 1d0
              contr%nvtx = 2
              contr%vertex(1)%idx_op = idxham
              contr%vertex(1)%iblk_op = iblk_ham
              contr%vertex(2)%idx_op = idxtop
              contr%vertex(2)%iblk_op = iblk_t1
              contr%narc = 1
              contr%arc(1)%link(1)=1
              contr%arc(1)%link(2)=2
              contr%arc(1)%occ_cnt = iocc_cnt
              if (ntest.ge.100) then
                call prt_contr(luout,contr,ops)
              end if
              call wrt_contr(luhhat,contr)
            else
              ! two or more commutators
              
              ! holes to contract:
              icnt(1) = iocc_cnt(ihole,1)
              ! particles to contract:
              icnt(2) = iocc_cnt(ipart,2)

              init = .true.
              ! well, there is only one possibility, but this
              ! routine exists and can be used here
              do while (next_part_pair(init,.true.,icnt_part,
     &                               icnt,nt1,1,1))
                init = .false.

                t1_part(1:nt1) = 1
                icnt_part_idx(1:nt1) = icnt_part(1,1:nt1)
     &                                 + icnt_part(2,1:nt1)*(2)
                contr%fac =
     &               1d0/dble(ieqfac(t1_part,icnt_part_idx,nt1))
                contr%nvtx = nt1+1
                contr%narc = nt1
                contr%vertex(1)%idx_op = idxham
                contr%vertex(1)%iblk_op = iblk_ham
                do idx = 1, nt1
                  contr%vertex(idx+1)%idx_op = idxtop
                  contr%vertex(idx+1)%iblk_op = iblk_t1
                  contr%arc(idx)%link(1) = 1
                  contr%arc(idx)%link(2) = idx+1
                  contr%arc(idx)%occ_cnt(1:ngastp,1:2) = 0
                  contr%arc(idx)%occ_cnt(ihole,1) =
     &                 icnt_part(1,idx)
                  contr%arc(idx)%occ_cnt(ipart,2) =
     &                 icnt_part(2,idx)
                end do
                if (ntest.ge.100) then
                  call prt_contr(luout,contr,ops)
                end if
                call wrt_contr(luhhat,contr)
                
              end do
            end if

          end if

        end do
      end do

      call file_close_keep(form_hhat%fhand)
      deallocate(contr%vertex,contr%arc)

      return
      end
