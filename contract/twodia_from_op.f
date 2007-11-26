      subroutine twodia_from_op(x2dia,ffop,op,orb_info)
*-----------------------------------------------------------------------
*     Routine to extract the diagonal part of a two-electron operator.
*     GWR November 2007
*-----------------------------------------------------------------------

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'ifc_baserout.h'
      include 'explicit.h'

      real(8), intent(inout) ::
     &     x2dia(*)
      type(operator), intent(in) ::
     &     op
      type(filinf), intent(in) ::
     &     ffop
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
     &     loop(op%n_occ_cls)

      real(8), pointer ::
     &     buffer(:), curblk(:)

      ifree = mem_setmark('twodia')

      ! Initialise some things.
      nbuff = 0
      loop(1:op%n_occ_cls) = .false.
      njoined = op%njoined

      ! Allocate space if the file, ffop, is not buffered.
      if(.not.ffop%buffered) then
        do iocc_cls = 1,op%n_occ_cls
          if(max(op%ica_occ(1,iocc_cls),
     &         op%ica_occ(2,iocc_cls)).ne.2)then
            loop(iocc_cls) = .true.
            cycle
          endif

          if(iextr.gt.0)then
            if(.not.explicit.and.max(op%ihpvca_occ(iextr,1,iocc_cls),
     &           op%ihpvca_occ(iextr,2,iocc_cls)).gt.0)then
              loop(iocc_cls) = .true.
              cycle
            endif
            if(explicit.and.op%formal_blk(iocc_cls))then
              loop(iocc_cls) = .true.
              cycle
            endif
          endif

          occ_temp(1:ngastp,1:2) = 0
          do idx = 1, njoined
            occ_temp(1:ngastp,1:2) = occ_temp(1:ngastp,1:2)+
     &           op%ihpvca_occ(1:ngastp,1:2,(iocc_cls-1)*njoined+idx)
          enddo

          if(list_cmp(occ_temp(1,1),occ_temp(1,2),ngastp)) then
            nbuff = nbuff + op%len_op_occ(iocc_cls)
          else
            loop(iocc_cls) = .true.
          endif
        enddo

        ifree = mem_alloc_real(buffer,nbuff,'buffer')
      endif

      ! Initialise the array in which the diagonal elements will be held.
      ! This is the maximum length needed if no externals (I hope).
      x2dia(1:4*orb_info%ntoob**2) = 0d0
      x2_off = 0

      mostnd => orb_info%mostnd
      ihpvgas => orb_info%ihpvgas
      iad_gas => orb_info%iad_gas

      do iocc_cls = 1,op%n_occ_cls
        ! Two-electron operators of the appropriate type only.
        if(loop(iocc_cls))cycle

        ioff_blk = op%off_op_occ(iocc_cls)
        ilen_blk = op%len_op_occ(iocc_cls)
        if(ffop%buffered)then
          idxbuf = ffop%idxrec(iocc_cls)+1
          if(idxbuf.lt.0)
     &         call quit(1,'twodia_from_op',
     &         'idxrec inconsistent for file '//trim(ffop%name))
          buffer => ffop%buffer(idxbuf:)
        else
          call get_vec(ffop,buffer,ioff_blk+1,ioff_blk+ilen_blk)
        endif

        ! Ascertain the nature of the operator.
        ihpv = 0
        ihpv(1) = idxlist(1,op%ihpvca_occ(1,1,iocc_cls),ngastp,1)
        ihpv(2) = idxlist(1,op%ihpvca_occ(1,ihpv(1)+1,iocc_cls),
     &       ngastp-ihpv(1)+1,1)
        if(ihpv(1).le.0)then
          ihpv(1) = idxlist(2,op%ihpvca_occ(1,1,iocc_cls),ngastp,1)
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

            ! Offset for the buffer array.
            idxcur = op%off_op_gmo(iocc_cls)%gam_ms(isym,idxms)
     &           -ioff_blk+1
            curblk => buffer(idxcur:)

            ! Length of the block required in x2dia
            len_gam_ms=int(sqrt(dble(op%
     &           len_op_gmo(iocc_cls)%gam_ms(isym,idxms))))

            ! Actual loop to extract diagonal elements from each block.
            do idx = 1, len_gam_ms
              x2dia(x2_off+idx) = curblk((idx-1)*len_gam_ms+idx)
            enddo

            x2_off = x2_off + len_gam_ms

          enddo

        enddo

      enddo

      ifree = mem_flushmark()

      return
      end
