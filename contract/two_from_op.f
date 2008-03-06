      subroutine two_from_op(x2,offsets,nblk,mel,orb_info)
*-----------------------------------------------------------------------
*     Routine to extract the diagonal part of a two-electron operator.
*     GWR February 2008
*-----------------------------------------------------------------------

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'ifc_baserout.h'

      real(8), intent(out) ::
     &     x2(*)
      integer, intent(out) ::
     &     offsets(*), nblk
      type(me_list), intent(in) ::
     &     mel
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     ifree, iocc_cls, nbuff, ioff_blk, ilen_blk, ms, idxms,
     &     imo_off, isym, idxcur, idx, igas, idxbuf, x2_off, len_gam_ms,
     &     njoined
      integer ::
     &     ihpv(2), occ_temp(ngastp,2)
      integer, pointer ::
     &     ihpvgas(:), iad_gas(:), mostnd(:,:,:)

      logical ::
     &     loop(mel%op%n_occ_cls)

      real(8), pointer ::
     &     buffer(:), curblk(:)

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop

      logical, external ::
     &     occ_is_diag_blk

      ifree = mem_setmark('twodia')

      op => mel%op
      ffop => mel%fhand

      ! Initialise some things.
      nbuff = 0
      loop(1:op%n_occ_cls) = .false.
      njoined = op%njoined

      do iocc_cls = 1,op%n_occ_cls
        if(max(op%ica_occ(1,iocc_cls),
     &       op%ica_occ(2,iocc_cls)).ne.2  .or.
     &       op%formal_blk(iocc_cls)  .or. .not.
     &       occ_is_diag_blk(op%ihpvca_occ(1,1,(iocc_cls-1)*njoined+1),
     &       njoined))
     &       loop(iocc_cls) = .true.          
      enddo

      ! Allocate space if the file, ffop, is not buffered.
      if(.not.ffop%buffered) then
        nbuff = 0
        do iocc_cls = 1,op%n_occ_cls
          if (loop(iocc_cls)) cycle

          nbuff = max(nbuff, mel%len_op_occ(iocc_cls))
        enddo

        ifree = mem_alloc_real(buffer,nbuff,'buffer')
      endif

      ! Initialise the array in which the diagonal elements will be held.
      ! This is the maximum length needed if no externals (I hope).
      x2(1:4*orb_info%ntoob**2) = 0d0
      x2_off = 0

      mostnd => orb_info%mostnd
      ihpvgas => orb_info%ihpvgas
      iad_gas => orb_info%iad_gas

      nblk = 0
      do iocc_cls = 1,op%n_occ_cls
        ! Two-electron operators of the appropriate type only.
        ! Formal blocks are ignored.
        if(loop(iocc_cls))cycle
        
        ioff_blk = mel%off_op_occ(iocc_cls)
        ilen_blk = mel%len_op_occ(iocc_cls)
        if(ffop%buffered)then
          idxbuf = ffop%idxrec(iocc_cls)+1
          if(idxbuf.lt.0)
     &         call quit(1,'twodia_from_op',
     &         'idxrec inconsistent for file '//trim(ffop%name))
          buffer => ffop%buffer(idxbuf:)
        else
          call get_vec(ffop,buffer,ioff_blk+1,ioff_blk+ilen_blk)
        endif

        ! Loop over the Ms of the annihilator string (equal to that of
        ! the creator as the operators should be symmetric).
        idxms = 0
        do ms = 2, -2, -2
          idxms = idxms + 1
c          if(ms.eq.-2) idxms = idxms+1
c          imo_off = (idxms-1)*orb_info%ntoob**2

          ! Loop over irrep. of annihilator string.
          do isym = 1,orb_info%nsym

            nblk = nblk+1
            offsets(nblk) = x2_off

            ! Offset for the buffer array.
            idxcur = mel%off_op_gmo(iocc_cls)%gam_ms(isym,idxms)
     &           -ioff_blk+1
            curblk => buffer(idxcur:)

            if (mel%off_op_gmox(iocc_cls)%ndis(isym,idxms).gt.1)
     &           call quit(1,'twodia_from_op','ndis>1 occurred!')

            ! Length of the block required in x2
            len_gam_ms=mel%ld_op_gmox(iocc_cls)%d_gam_ms(1,isym,idxms)

            ! Actual loop to extract diagonal elements from each block.
            do idx = 1, len_gam_ms**2
              x2(x2_off+idx) = curblk(idx)
            enddo

            x2_off = x2_off + len_gam_ms**2

          enddo

        enddo

      enddo

      if (ntest.ge.100) then
        write(luout,*) 'blocks: ',nblk
        write(luout,*) 'offsets:',offsets(1:nblk)
        write(luout,*) 'elements:'
        write(luout,'(x,5(g18.8))') x2(1:x2_off)
      end if

      ifree = mem_flushmark()

      return
      end
