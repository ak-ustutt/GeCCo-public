*----------------------------------------------------------------------*
      subroutine onedia_from_op(x1dia,ffop,op,orb_info)
*----------------------------------------------------------------------*
*     extract the diagonal of the rank-one part of operator op
*----------------------------------------------------------------------*

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

      real(8), intent(out) ::
     &     x1dia(*)
      type(operator), intent(in) ::
     &     op
      type(filinf), intent(in) ::
     &     ffop
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     ifree, iocc_cls, nbuff, ioff_blk, ilen_blk, idx, idxcur, 
     &     most, mond, ms, idxms, isym, imo, imo_off, len, igas, idxbuf,
     &     ihpv

      integer, pointer ::
     &     ihpvgas(:), iad_gas(:), mostnd(:,:,:)

      real(8), pointer ::
     &     buffer(:), curblk(:)

      ifree = mem_setmark('onedia')

      if (.not.ffop%buffered) then
        do iocc_cls = 1, op%n_occ_cls
          ! see below for explanations
          if (max(op%ica_occ(1,iocc_cls),op%ica_occ(2,iocc_cls)).ne.1)
     &         cycle
c          if (iextr.gt.0.and.max(op%ihpvca_occ(iextr,1,iocc_cls),
c     &                         op%ihpvca_occ(iextr,2,iocc_cls)).gt.0)
c     &         cycle
          if(iextr.gt.0)then
            if(.not.explicit.and.max(op%ihpvca_occ(iextr,1,iocc_cls),
     &           op%ihpvca_occ(iextr,2,iocc_cls)).gt.0)
     &           cycle
            if(explicit.and.op%formal_blk(iocc_cls)) cycle
          endif
          if (list_cmp(op%ihpvca_occ(1,1,iocc_cls),
     &                 op%ihpvca_occ(1,2,iocc_cls),ngastp)) then
            nbuff = nbuff + op%len_op_occ(iocc_cls)
          end if
        end do
        ifree = mem_alloc_real(buffer,nbuff,'buffer')
      end if

      x1dia(1:2*(orb_info%ntoob+orb_info%caborb)) = 0d0

      mostnd => orb_info%mostnd
      ihpvgas => orb_info%ihpvgas
      iad_gas => orb_info%iad_gas

      do iocc_cls = 1, op%n_occ_cls
        ! 1-electron operators only ....
        if (max(op%ica_occ(1,iocc_cls),op%ica_occ(2,iocc_cls)).ne.1)
     &       cycle
        ! ... but the only normal ones (i.e. no 
        ! external/auxiliary orbitals)
c        if (iextr.gt.0.and.max(op%ihpvca_occ(iextr,1,iocc_cls),
c     &                         op%ihpvca_occ(iextr,2,iocc_cls)).gt.0)
c     &       cycle
          if(iextr.gt.0)then
            if(.not.explicit.and.max(op%ihpvca_occ(iextr,1,iocc_cls),
     &           op%ihpvca_occ(iextr,2,iocc_cls)).gt.0)
     &           cycle
            if(explicit.and.op%formal_blk(iocc_cls)) cycle
          endif  

        ! diagonal: so C and A must have same occ
        if (.not.list_cmp(op%ihpvca_occ(1,1,iocc_cls),
     &                    op%ihpvca_occ(1,2,iocc_cls),ngastp)) cycle

        ioff_blk = op%off_op_occ(iocc_cls)
        ilen_blk = op%len_op_occ(iocc_cls)
        if (ffop%buffered) then
          idxbuf = ffop%idxrec(iocc_cls)+1
          if (idxbuf.le.0)
     &         call quit(1,'onedia_from_op',
     &         'idxrec inconsistent for file '//trim(ffop%name))
          buffer => ffop%buffer(idxbuf:)
        else
          call get_vec(ffop,buffer,ioff_blk+1,ioff_blk+ilen_blk)
        end if

        ihpv = idxlist(1,op%ihpvca_occ(1,1,iocc_cls),ngastp,1)

        do ms = 1, -1, -2
          idxms =1
          if (ms.eq.-1) idxms = 2
          imo_off = (idxms-1)*(orb_info%ntoob+orb_info%caborb)

          do isym = 1, orb_info%nsym

            idxcur = op%off_op_gmo(iocc_cls)%gam_ms(isym,idxms)
     &               -ioff_blk+1
            curblk => buffer(idxcur:)
            idx = 0
            do igas = 1, orb_info%ngas
              if (ihpvgas(igas).ne.ihpv) cycle
              if (iad_gas(igas).ne.2) cycle
              most = mostnd(1,isym,igas)
              mond = mostnd(2,isym,igas)
              len = mond - most + 1              
              if (len.le.0) cycle

              do imo = most, mond
                idx = idx+1
                x1dia(imo+imo_off) = curblk((idx-1)*len+idx)
              end do

            end do

          end do

        end do

      end do
      
      ifree = mem_flushmark()

      if (ntest.ge.100) then
        write(luout,*) 'extracted diagonal: '
        idx = 0
        do ms = 1, -1, -2
          do imo = 1, orb_info%ntoob+orb_info%caborb
            idx = idx+1
            write(luout,'(x,i2,"/2",i5,2x,g12.6)') ms, imo, x1dia(idx)
          end do
        end do
      end if

      return
      end
