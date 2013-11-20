*----------------------------------------------------------------------*
      subroutine orb_flip_mel(mel,str_info,orb_info)
*----------------------------------------------------------------------*
*     Retrieves relative orbital signs from the file FLIPSIGNS
*     (a list of logicals that indicates which orbitals' signs changed)
*     and flips the signs of the matrix elements accordingly.
*
*     Matthias, Nov 2013
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_operators.h'
      include 'ifc_memman.h'

      type(me_list), intent(in) ::
     &     mel
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      
      logical ::
     &     close_again, scalar, l_exist
      integer ::
     &     idoff, idxoff, idxoff_blk, iblk, lenblk, ifree, mmax,
     &     msamax, mscmax, idxms, ms, igam, idx_dis, ndis, nwarn,
     &     nel, mst, njoined, idx_occ, norb, luinp, nblk, iostatus,
     &     iorb, lenblk2
      real(8), pointer ::
     &     buffer(:)
      logical, pointer ::
     &     flipsigns_in(:), flipsigns(:)

      type(filinf) ::
     &     ffsigns
      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop

      ifree = mem_setmark('orb_flip_mel')

      op => mel%op
      ffop => mel%fhand
      mst = mel%mst
      nblk = op%n_occ_cls

      inquire(file='FLIPSIGNS',exist=l_exist)
      if (.not.l_exist) call quit(1,'orb_flip_mel',
     &                            'No FLIPSIGNS file found.')
      call file_init(ffsigns,'FLIPSIGNS',ftyp_sq_frm,0)
      call file_open(ffsigns)
      luinp = ffsigns%unit
      norb = orb_info%ntoob
c dbg
c      print *,'norb:',norb
c dbgend
      allocate(flipsigns_in(norb),flipsigns(norb))
      iostatus = 0
      iorb = 0
      do while (iostatus.ge.0.and.iorb.lt.norb)
        iorb = iorb + 1
        read(luinp,*,iostat=iostatus) flipsigns_in(iorb)
      end do
      if (iorb.ne.norb) call quit(1,'orb_flip_mel',
     &                            'Less signs than orbitals.')

      ! Reorder sign according to orbital maps
      flipsigns(1:norb) = flipsigns_in(orb_info%ireots(1:norb))
c dbg
c      print *,'flsigns_in:',flipsigns_in(1:norb)
c      print *,'flsigns   :',flipsigns(1:norb)
c dbgend

      deallocate(flipsigns_in)
      call file_close_keep(ffsigns)

      nwarn = 0

      mmax = 0
      do iblk = 1, nblk
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
      do iblk = 1, nblk
        if(op%formal_blk(iblk))cycle

        scalar = max(op%ica_occ(1,iblk),op%ica_occ(2,iblk)).eq.0

        mscmax = op%ica_occ(1,iblk)
        msamax = op%ica_occ(2,iblk)
        nel = msamax+op%ica_occ(1,iblk)
        idxms = 0
        idxoff = 0
        idxoff_blk = 0
        do ms = msamax, -msamax, -2
          if (abs(ms+mst).gt.mscmax) cycle
          idxms = idxms+1
          do igam = 1, orb_info%nsym

            ! block offest and length:
            idxoff = mel%off_op_gmo(iblk)%gam_ms(igam,idxms)
            lenblk = mel%len_op_gmo(iblk)%gam_ms(igam,idxms)
            if (lenblk.eq.0) cycle

            ! get current block
            idoff = ffop%length_of_record*(ffop%current_record-1)
            call get_vec(ffop,buffer,
     &           idoff+idxoff+1,idoff+idxoff+lenblk)

            ! reset offset within current buffer
            idxoff_blk = 0
        
            ! find distribution
            ndis = mel%off_op_gmox(iblk)%ndis(igam,idxms)
            if (ndis.eq.0) then
              write(lulog,*) 'WARNING:'
              write(lulog,*)
     &             ' op_info indicates no distribution at all'
              write(lulog,*) ' skipping to next block'
              idxoff_blk = idxoff_blk+lenblk
              nwarn = nwarn+1
              cycle
            end if
            if (ndis.eq.1) then
              if (.not.scalar)
     &           call orb_flip_blk(buffer(idxoff_blk+1),
     &                             mel,iblk,igam,idxms,1,
     &                             nel,flipsigns,str_info,orb_info)
              idxoff_blk = idxoff_blk+lenblk
              ! store current block
              call put_vec(ffop,buffer,
     &                     idoff+idxoff+1,idoff+idxoff+lenblk)
              cycle
            end if
            distr_loop: do idx_dis = 1, ndis
                
              lenblk2 =
     &             mel%len_op_gmox(iblk)%d_gam_ms(idx_dis,igam,idxms)
              if (lenblk2.eq.0) cycle
            
              call orb_flip_blk(buffer(idxoff_blk+1),
     &                          mel,iblk,igam,idxms,idx_dis,
     &                          nel,flipsigns,str_info,orb_info)
              idxoff_blk = idxoff_blk+lenblk2
            end do distr_loop

            ! store current block
            call put_vec(ffop,buffer,
     &                   idoff+idxoff+1,idoff+idxoff+lenblk)

          end do ! gam
        end do ! ms
        idx_occ = idx_occ+njoined
      end do

      if (nwarn.gt.0) then
        write(lulog,*) '!!! There were ',nwarn,' warnings !!!'
        write(lulog,*) 'look for "WARNING" in previous output!'
      end if

      if (close_again) call file_close_keep(ffop)

      deallocate(flipsigns)

      ifree = mem_flushmark()

      return
      end
