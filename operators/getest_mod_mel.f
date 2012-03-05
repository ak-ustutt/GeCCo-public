*----------------------------------------------------------------------*
      subroutine getest_mod_mel(mel,iRdef,norb,icase,icaseF,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*  
*  modify list for generalize size extensivity test
*
*  andreas, feb 2012
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_operators.h'
      include 'ifc_memman.h'


      integer, intent(in) ::
     &     norb, iRdef(norb), icase, icaseF
      type(me_list), intent(in) ::
     &     mel
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      
      logical ::
     &     first, close_again, blk_buf, scalar, dagger
      integer ::
     &     idoff, idxoff, idxoff_blk, iblk, lenblk, lenprt, ifree, mmax,
     &     msamax, mscmax, idxms, ms, igam, idx_dis, ndis, nwarn, did,
     &     idum, nel, mst, lengblk,
     &     idxoff0, njoined, idx_occ, icase_cur
      integer ::
     &     msd(ngastp,2,mel%op%njoined), igamd(ngastp,2,mel%op%njoined),
     &     scr(ngastp,2,2*mel%op%njoined), occ(ngastp,2,mel%op%njoined)
      real(8) ::
     &     xnrm, xnrm_tot, xnrm_ms
      real(8), pointer ::
     &     buffer(:), curblk(:)

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      logical, external ::
     &     next_msgamdist
      real(8), external ::
     &     ddot

      ifree = mem_setmark('getest_mod_mel')

      op => mel%op
      ffop => mel%fhand
      mst = mel%mst
      dagger = mel%op%dagger

      nwarn = 0

      mmax = 0
      do iblk = 1, op%n_occ_cls
        if(op%formal_blk(iblk))cycle
        mscmax = op%ica_occ(1,iblk)
        msamax = op%ica_occ(2,iblk)
        idxms = 0
        do ms = msamax, -msamax, -2
          if (abs(mst+ms).gt.mscmax) cycle
          idxms = idxms+1
          do igam = 1, orb_info%nsym
            mmax = max(mmax,mel%len_op_gmo(iblk)%gam_ms(igam,idxms))
          end do
        end do
      end do
      ifree = mem_alloc_real(buffer,mmax,'buffer')

      close_again = .false.
      if (ffop%unit.le.0) then
        close_again = .true.
        call file_open(ffop)
      end if

      njoined = mel%op%njoined

      idx_occ = 1
      do iblk = 1, op%n_occ_cls
        if(op%formal_blk(iblk))cycle

        occ = op%ihpvca_occ(1:ngastp,1:2,idx_occ:idx_occ+njoined-1)
        if (dagger) occ = iocc_dagger_n(occ,njoined)

        blk_buf = ffop%buffered
        if (blk_buf) blk_buf = blk_buf.and.ffop%incore(iblk).gt.0

        scalar = max(op%ica_occ(1,iblk),op%ica_occ(2,iblk)).eq.0

        mscmax = op%ica_occ(1,iblk)
        msamax = op%ica_occ(2,iblk)

        ! if the Fock operator is not to be modifier, go to next blk
        if (max(mscmax,msamax).eq.0) cycle
        if (max(mscmax,msamax).eq.1) then
          if (icaseF.le.0) cycle
          icase_cur = icaseF
        else
          icase_cur = icase
        end if

        nel = msamax+op%ica_occ(1,iblk)
        idxms = 0
        idxoff = 0
        idxoff_blk = 0
        do ms = msamax, -msamax, -2
          if (abs(ms+mst).gt.mscmax) cycle
          xnrm_ms = 0d0
          idxms = idxms+1
          do igam = 1, orb_info%nsym

            ! block offest and length:
            idxoff = mel%off_op_gmo(iblk)%gam_ms(igam,idxms)
            lengblk = mel%len_op_gmo(iblk)%gam_ms(igam,idxms)
            if (lengblk.eq.0) cycle

            ! get current block
            if (.not.blk_buf) then
              idoff = ffop%length_of_record*(ffop%current_record-1)
              call get_vec(ffop,buffer,
     &             idoff+idxoff+1,idoff+idxoff+lengblk)
              curblk => buffer(1:lengblk)
            else 
c              ioff = op%off_op_gmo(iblk)%gam_ms(igam,idxms)
              ! currently: idxoff should be valid here, as well
              curblk => ffop%buffer(idxoff+1:idxoff+lengblk)
            end if

            ! reset offset within current buffer
            idxoff_blk = 0

            ndis = mel%off_op_gmox(iblk)%ndis(igam,idxms)
 
            if (icase_cur.eq.1) then ! most simple case: just ...
              curblk(1:lengblk) = 0d0  ! ... set to zero
            else if (ndis.eq.1) then
              call getest_mod_mel_blk(curblk(idxoff_blk+1),mel,
     &                   iRdef,norb,icase_cur,
     &                   iblk,igam,idxms,1,
     &                   nel,str_info,orb_info)
              idxoff_blk = idxoff_blk+lenblk
            else
              distr_loop: do idx_dis = 1, ndis
                  
                did = mel%off_op_gmox(iblk)%did(idx_dis,igam,idxms)

                call did2msgm(msd,igamd,did,
     &               op%ihpvca_occ(1,1,idx_occ),orb_info%nsym,njoined)
                lenblk =
     &               mel%len_op_gmox(iblk)%d_gam_ms(idx_dis,igam,idxms)
                if (lenblk.eq.0) cycle

                   call getest_mod_mel_blk(curblk(idxoff_blk+1),mel,
     &                     iRdef,norb,icase_cur,
     &                     iblk,igam,idxms,idx_dis,
     &                     nel,str_info,orb_info)
                idxoff_blk = idxoff_blk+lenblk
              end do distr_loop
            end if

            ! put back modified block 
            if (.not.blk_buf) then
              idoff = ffop%length_of_record*(ffop%current_record-1)
              call put_vec(ffop,buffer,
     &             idoff+idxoff+1,idoff+idxoff+lengblk)
            end if
          end do ! gam
        end do ! ms
        idx_occ = idx_occ+njoined
      end do

      if (close_again) call file_close_keep(ffop)

      ifree = mem_flushmark()

      return
      end
