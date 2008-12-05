*----------------------------------------------------------------------*
      subroutine tran_one(me_mo,ffao,ffcmo,orb_info)
*----------------------------------------------------------------------*
*     transform one-particle quantity on file ffao using CMOs from
*     file ffcmo and store result on ME-list me_mo
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'multd2h.h'
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
     &     me_mo
      type(orbinf), intent(in), target ::
     &     orb_info


      logical ::
     &     closeit, close_ffmo, is_xmo
      integer ::
     &     ifree, nsym, ngas, nspin, nblk, ncmo, nmo, nao, nhlf,
     &     isym, jsym, igas, ispin, iblk, idxst, idxnd,
     &     idxms, norb, njoined,
     &     hpvx_a, hpvx_c, iblkoff, ijoin, len_blk,
     &     me_sym, ipass

      integer, pointer ::
     &     nbas(:), nxbas(:), ntoobs(:), mostnd(:,:,:),
     &     ica_occ(:,:), hpvx_occ(:,:,:),
     &     iad_gas(:), hpvx_gas(:,:), idxcmo(:,:,:)
      real(8), pointer ::
     &     cmo(:), xao(:), xmo(:), xop(:), xhlf(:)
      type(filinf), pointer ::
     &     ffmo
      type(operator), pointer ::
     &     opdef
      real(8) ::
     &     xref

      integer, external ::
     &     max_dis_blk, idxlist

      ifree = mem_setmark('tran_one')

      opdef => me_mo%op
      ffmo  => me_mo%fhand

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'tran_one')
        write(luout,*) 'transform one-particle part of: ',
     &       trim(me_mo%label)
      end if

      me_sym = me_mo%gamt

      ! set up some dimensions
      nsym = orb_info%nsym
      ngas = orb_info%ngas
      nspin = orb_info%nspin
      nbas => orb_info%nbas
      nxbas => orb_info%nxbas
      ntoobs => orb_info%ntoobs
      mostnd => orb_info%mostnd
      iad_gas => orb_info%iad_gas
      hpvx_gas => orb_info%ihpvgas

      allocate(idxcmo(nsym,ngas,2))

      nao  = 0
      do isym = 1, nsym
        jsym = multd2h(isym,me_sym)
        nao  = nao  + (nbas(isym)+nxbas(isym))*
     &                (nbas(jsym)+nxbas(jsym))
      end do
      ncmo = 0
      nmo  = 0
      nhlf = 0
      do ipass = 1, 2
        do igas = 1, ngas
          is_xmo = hpvx_gas(igas,1).eq.IEXTR
          do isym = 1, nsym
            jsym = multd2h(isym,me_sym)
            idxcmo(isym,igas,ipass) = ncmo+1
            if (     is_xmo.and.ipass.eq.1) cycle
            if (.not.is_xmo.and.ipass.eq.2) cycle
            norb = mostnd(2,isym,igas)-mostnd(1,isym,igas)+1
            if (.not.is_xmo) ncmo = ncmo + nbas(isym)*norb
            if (     is_xmo) ncmo = ncmo + (nbas(isym)+nxbas(isym))*norb
            nhlf = max(nhlf,(nxbas(jsym)+nbas(jsym))*norb)
          end do
        end do
      end do

      nblk = opdef%n_occ_cls
      njoined = opdef%njoined
      ica_occ => opdef%ica_occ
      hpvx_occ => opdef%ihpvca_occ

      do iblk = 1, nblk
        ! only zero- and one-particle part is of interest:
        if (ica_occ(1,iblk).gt.1 .or.
     &      ica_occ(2,iblk).gt.1) cycle
        nmo = max(nmo,max_dis_blk(-1,me_mo,iblk,orb_info))
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

      ! read from file
      closeit = .false.
      if (ffao%unit.lt.0) then
        call file_open(ffao)
        closeit = .true.
      end if

      call get_vec(ffao,xao,1,nao)

      if (closeit)
     &     call file_close_keep(ffao)

      close_ffmo = .false.
      if (ffmo%unit.lt.0) then
        call file_open(ffmo)
        close_ffmo = .false.
      end if

      ! loop over all blocks of operator
      do iblk = 1, nblk
        ! process zero-particle term: <0|V|0>
        if (max(ica_occ(1,iblk),ica_occ(2,iblk)).eq.0 .and.
     &      me_mo%len_op_gmo(iblk)%gam_ms(1,1).gt.0) then
          call calc_xref(xref,xao,cmo,xhlf,nsym,nbas,nxbas,
     &                   ngas,hpvx_gas(:,1),idxcmo(:,:,1),
     &                   mostnd)
          idxst = me_mo%off_op_gmo(iblk)%gam_ms(1,1) + 1
          call put_vec(ffmo,xref,idxst,idxst)
          cycle
        end if

        ! process one-particle part
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

          ispin = 1
          if (nspin.eq.2.and.idxms.eq.2) ispin = 2
   
          ! buffer pointer
          idxst = me_mo%off_op_gmo(iblk)%gam_ms(1,idxms) + 1
          idxnd = idxst-1 
          do isym = 1, nsym
            idxnd = idxnd + me_mo%len_op_gmo(iblk)%gam_ms(isym,idxms)
          end do
          if (ffmo%buffered.and.ffmo%incore(iblk).gt.0) then
            xop => ffmo%buffer(idxst:)
          else
            xop => xmo
          end if
          len_blk = idxnd-idxst+1
          xop(1:len_blk) = 0d0

          call tran_one_blk(xop,xao,cmo,xhlf,
     &         me_sym,idxcmo,hpvx_c,hpvx_a,
     &         nbas,nxbas,mostnd,iad_gas,hpvx_gas(1,ispin),ngas,nsym)

          if (.not.ffmo%buffered.or.ffmo%incore(iblk).eq.0) then
            call put_vec(ffmo,xop,idxst,idxnd)
          end if

        end do
      end do

      if (close_ffmo)
     &     call file_close_keep(ffmo)

      deallocate(idxcmo)

      ifree = mem_flushmark('tran_one')

      return
      end
