*----------------------------------------------------------------------*
      subroutine btran_one(ffao,ffcmo,add_dref,me_dens,orb_info)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      type(filinf), intent(inout) ::
     &     ffao, ffcmo
      type(me_list), intent(inout) ::
     &     me_dens
      type(orbinf), intent(in), target ::
     &     orb_info
      logical, intent(in) ::
     &     add_dref


      logical ::
     &     closeit, close_ffmo
      integer ::
     &     ifree, nsym, ngas, nblk, ncmo, nmo, nao, nhlf,
     &     isym, igas, iblk, idxst, idxnd, idxms, norb, njoined,
     &     hpvx_a, hpvx_c, iblkoff, ijoin

      integer, pointer ::
     &     nbas(:), ntoobs(:), mostnd(:,:,:),
     &     ica_occ(:,:), hpvx_occ(:,:,:),
     &     iad_gas(:), hpvx_gas(:), idxcmo(:,:)
      real(8), pointer ::
     &     cmo(:), xao(:), xmo(:), xop(:), xhlf(:)
      type(filinf), pointer ::
     &     ffmo
      type(operator), pointer ::
     &     opdef

      integer, external ::
     &     max_dis_blk, idxlist

      ifree = mem_setmark('btran_one')

      opdef => me_dens%op
      ffmo  => me_dens%fhand

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'btran_one')
        write(luout,*) 'backtransform one-particle part of: ',
     &       trim(me_dens%label)
      end if

      if (me_dens%gamt.ne.1)
     &     call quit(1,'btran_one',
     &                 'adapt for non-totally symmetric density')

      ! set up some dimensions
      nsym = orb_info%nsym
      ngas = orb_info%ngas
      nbas => orb_info%nbas
      ntoobs => orb_info%ntoobs
      mostnd => orb_info%mostnd
      iad_gas => orb_info%iad_gas
      hpvx_gas => orb_info%ihpvgas

      allocate(idxcmo(nsym,ngas))

      nao  = 0
      do isym = 1, nsym
        nao  = nao  + nbas(isym)*nbas(isym)
      end do
      ncmo = 0
      nmo  = 0
      nhlf = 0
      do igas = 1, ngas
        do isym = 1, nsym
          idxcmo(isym,igas) = ncmo+1
          norb = mostnd(2,isym,igas)-mostnd(1,isym,igas)+1
          ncmo = ncmo + nbas(isym)*norb
          nhlf = max(nhlf,ncmo)
        end do
      end do

      nblk = opdef%n_occ_cls
      njoined = opdef%njoined
      ica_occ => opdef%ica_occ
      hpvx_occ => opdef%ihpvca_occ

      do iblk = 1, nblk
        ! only one-particle part is of interest:
        if (ica_occ(1,iblk).ne.1 .or.
     &      ica_occ(2,iblk).ne.1) cycle
        nmo = max(nmo,max_dis_blk(-1,me_dens,iblk,orb_info))
      end do

      ifree = mem_alloc_real(cmo,ncmo,'CMO')
      ifree = mem_alloc_real(xao,nao,'AObuf')
      ifree = mem_alloc_real(xmo,nmo,'MObuf')
      ifree = mem_alloc_real(xhlf,nhlf,'HLFbuf')

      closeit = .false.
      if (ffcmo%unit.lt.0) then
        call file_open(ffcmo)
        closeit = .true.
      end if

      call get_vec(ffcmo,cmo,1,ncmo)

      if (closeit)
     &     call file_close_keep(ffcmo)

      close_ffmo = .false.
      if (ffmo%unit.lt.0) then
        call file_open(ffmo)
        close_ffmo = .false.
      end if

      ! init to zero
      xao(1:nao) = 0d0 

      ! add reference contribution, if requested
      if (add_dref) then
        call make_refmat(xao,cmo,orb_info)
      end if

      ! loop over all blocks of density
      do iblk = 1, nblk
        ! process only one-particle density part
        if (ica_occ(1,iblk).ne.1 .or.
     &      ica_occ(2,iblk).ne.1) cycle

        iblkoff = (iblk-1)*njoined

        if (ntest.ge.100) then
          call wrt_occ_n(luout,hpvx_occ(1,1,iblkoff+1),njoined)
        end if

        ! get HPVX of C and A
        do ijoin = 1, njoined
          hpvx_c = idxlist(1,hpvx_occ(1,1,iblkoff+ijoin),ngastp,1)
          if (hpvx_c.gt.0) exit
        end do
        do ijoin = 1, njoined
          hpvx_a = idxlist(1,hpvx_occ(1,2,iblkoff+ijoin),ngastp,1)
          if (hpvx_a.gt.0) exit
        end do

        do idxms = 1, 2
        
          ! load density, or reassign buffer pointer
          idxst = me_dens%off_op_gmo(iblk)%gam_ms(1,idxms) + 1
          idxnd = idxst-1 
          do isym = 1, nsym
            idxnd = idxnd + me_dens%len_op_gmo(iblk)%gam_ms(isym,idxms)
          end do
          if (ffmo%buffered.and.ffmo%incore(iblk).gt.0) then
            xop => ffmo%buffer(idxst:)
          else
            call get_vec(ffmo,xmo,idxst,idxnd)
            xop => xmo
          end if

          call btran_one_blk(xao,xop,cmo,xhlf,
     &         me_dens%gamt,idxcmo,hpvx_c,hpvx_a,
     &         nbas,mostnd,iad_gas,hpvx_gas,ngas,nsym)

          if (ntest.ge.100) then
            write(luout,*) 'XAO after block,idxms: ',iblk, idxms
            call wr_blkmat(xao,nbas,nbas,nsym,0)
          end if

        end do
      end do

      if (close_ffmo)
     &     call file_close_keep(ffmo)

      closeit = .false.
      if (ffao%unit.lt.0) then
        call file_open(ffao)
        closeit = .true.
      end if

      call put_vec(ffao,xao,1,nao)

      if (closeit)
     &     call file_close_keep(ffao)

      deallocate(idxcmo)

      ifree = mem_flushmark('btran_one')

      return
      end
