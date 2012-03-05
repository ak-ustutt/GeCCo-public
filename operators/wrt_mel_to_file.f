*----------------------------------------------------------------------*
      subroutine wrt_mel_to_file(luout,mel,iblkst,iblknd,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     given: an operator definition (op) and a file handle
*         or buffer containing the matrix elements
*     print out operator blocks iblkst to iblknd as a plain list
*     simple version -- load complete symmetry block into memory
*     Modified to cater for formal operators. GWR July 2007
*     matthias, Mar 2012 (copied from wrt_mel)
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

      integer, parameter ::
     &     maxprt = 100

      integer, intent(in) ::
     &     luout, iblkst, iblknd
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
     &     idum, nel, mst,
     &     idxoff0, njoined, idx_occ
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

      ifree = mem_setmark('wrt_op2file')

      op => mel%op
      ffop => mel%fhand
      mst = mel%mst
      dagger = mel%op%dagger

      nwarn = 0

        
      mmax = 0
      do iblk = iblkst, iblknd
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

      xnrm_tot = 0d0
      idx_occ = (iblkst-1)*njoined+1
      do iblk = iblkst, iblknd
        if(op%formal_blk(iblk))cycle

        occ = op%ihpvca_occ(1:ngastp,1:2,idx_occ:idx_occ+njoined-1)
        if (dagger) occ = iocc_dagger_n(occ,njoined)

        blk_buf = ffop%buffered
        if (blk_buf) blk_buf = blk_buf.and.ffop%incore(iblk).gt.0

        scalar = max(op%ica_occ(1,iblk),op%ica_occ(2,iblk)).eq.0

        mscmax = op%ica_occ(1,iblk)
        msamax = op%ica_occ(2,iblk)
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
            lenblk = mel%len_op_gmo(iblk)%gam_ms(igam,idxms)
            if (lenblk.eq.0) cycle

            ! get current block
            if (.not.blk_buf) then
              idoff = ffop%length_of_record*(ffop%current_record-1)
              call get_vec(ffop,buffer,
     &             idoff+idxoff+1,idoff+idxoff+lenblk)
              curblk => buffer(1:lenblk)
            else
c              ioff = op%off_op_gmo(iblk)%gam_ms(igam,idxms)
              ! currently: idxoff should be valid here, as well
              curblk => ffop%buffer(idxoff+1:idxoff+lenblk)
            end if

            ! reset offset within current buffer
            idxoff_blk = 0
        
            xnrm = sqrt(ddot(lenblk,curblk,1,curblk,1))
            xnrm_tot = xnrm_tot + xnrm*xnrm
            xnrm_ms  = xnrm_ms + xnrm*xnrm

            lenprt = lenblk
            ! print out distributions
            first = .true.
            ndis = mel%off_op_gmox(iblk)%ndis(igam,idxms)
            if (ndis.eq.0) then
              write(luout,*) 'WARNING:'
              write(luout,*)
     &             ' op_info indicates no distribution at all'
              write(luout,*) ' skipping to next block'
              idxoff_blk = idxoff_blk+lenblk
              nwarn = nwarn+1
              cycle
            end if
            if (ndis.eq.1) then
              if (scalar) then
                write(luout,'(x,i10,17x,f24.14)') idxoff_blk+1,
     &             curblk(idxoff_blk+1)
              else
                call wrt_mel_blk_wi2(luout,curblk(idxoff_blk+1),
     &             mel,iblk,igam,idxms,1,
     &             nel,str_info,orb_info)
              end if
              idxoff_blk = idxoff_blk+lenblk
              cycle
            end if
            distr_loop: do idx_dis = 1, ndis
                
              did = mel%off_op_gmox(iblk)%did(idx_dis,igam,idxms)

              call did2msgm(msd,igamd,did,
     &             op%ihpvca_occ(1,1,idx_occ),orb_info%nsym,njoined)
              if (.not.dagger) then
                scr(1:ngastp,1:2,1:njoined) =
     &               msd(1:ngastp,1:2,1:njoined)
                scr(1:ngastp,1:2,njoined+1:2*njoined) =
     &               igamd(1:ngastp,1:2,1:njoined)
              else
                scr(1:ngastp,1:2,1:njoined) =
     &               iocc_dagger_n(msd,njoined)
                scr(1:ngastp,1:2,njoined+1:2*njoined) =
     &               iocc_dagger_n(igamd,njoined)
              end if

              lenblk =
     &             mel%len_op_gmox(iblk)%d_gam_ms(idx_dis,igam,idxms)
              if (lenblk.eq.0) cycle
              xnrm =
     &             sqrt(ddot(lenblk,curblk(idxoff_blk+1),1,
     &                              curblk(idxoff_blk+1),1))
              lenprt = lenblk
            
              call wrt_mel_blk_wi2(luout,curblk(idxoff_blk+1),
     &             mel,iblk,igam,idxms,idx_dis,
     &             nel,str_info,orb_info)

              idxoff_blk = idxoff_blk+lenblk
            end do distr_loop
          end do ! gam
        end do ! ms
        idx_occ = idx_occ+njoined
      end do

      if (nwarn.gt.0) then
        write(luout,*) '!!! There were ',nwarn,' warnings !!!'
        write(luout,*) 'look for "WARNING" in previous output!'
      end if

      if (close_again) call file_close_keep(ffop)

      ifree = mem_flushmark()

      return
      end
