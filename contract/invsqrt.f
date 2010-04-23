      subroutine invsqrt(mel_inp,mel_inv,nocc_cls,half,
     &     op_info,orb_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     Routine to calculate S^(-0.5) of density matrices.
*     The ME list is split into matrices that can either contain
*     elements from a single distribution or from several distributions.
*     half = true: only makes a half transform U*s^(-0.5)
*     half = false: full transform U*s^(-0.5)*U^+
*      
*     works also for sums of density matrices of the structure:
*       /0 0 0 0\
*       \0 0 x 0/  (not necessarily igastp=3)
*       /0 0 y 0\ 
*       \0 0 y 0/ 
*       /0 0 x 0\
*       \0 0 0 0/
*
*     matthias, dec 2009 (adopted from invert.f)
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'
      include 'hpvxseq.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 10

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(me_list), intent(in) ::
     &     mel_inp, mel_inv
      integer, intent(in) ::
     &     nocc_cls
      logical, intent(in) ::
     &     half

      type(graph), pointer ::
     &     graphs(:)

      logical ::
     &     bufin, bufout, first, ms_fix, fix_success, onedis, transp,
     &     logdum
c      logical ::
c     &     loop(nocc_cls)
      integer ::
     &     ifree, nbuff, idxmsa, iocc_cls,
     &     msmax, msa, igama, idx, jdx, ngam,
     &     ioff, njoined,
     &     idxdis, lenca, iblkoff, ncblk, nablk,
     &     msc, igamc, idxmsc, mscmax,
     &     ndis, ndim, off_col, off_line,
     &     idxa1, idxa2, idxc1, idxc2, icol, iline,
     &     msmax_sub, ms1, ms2, msc1, msc2, msa1, msa2,
     &     gamc1, gamc2, gama1, gama2, igam, na1, na2, nc1, nc2, ij,
     &     maxbuf, ngraph, ioff2, idxmsa2, ndis2, idxdis2, idx2, nel,
     &     nsing, ising, itrip, ntrip, icnt_sv, icnt_sv0
      real(8) ::
     &     fac
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:), scratch(:,:),
     &     sing(:,:), trip(:,:)

      integer, pointer ::
     &     hpvx_csub(:),hpvx_asub(:),
     &     occ_csub(:), occ_asub(:),
     &     graph_csub(:), graph_asub(:),
     &     msdis_c(:),  msdis_a(:),
     &     idxmsdis_c(:),  idxmsdis_a(:),
     &     msdis_c2(:),  msdis_a2(:),
     &     idxmsdis_c2(:),  idxmsdis_a2(:),
     &     gamdis_c(:), gamdis_a(:),
     &     len_str(:),
     &     hpvx_occ(:,:,:), ca_occ(:,:), idx_graph(:,:,:),
     &     ldim_opin_c(:), ldim_opin_a(:),
     &     istr_csub(:), istr_asub(:),
     &     istr_csub_flip(:), istr_asub_flip(:),
     &     flipmap_c(:), flipmap_a(:),
     &     flmap(:,:),
     &     iocc(:,:,:), idx_g(:,:,:), msdst(:,:,:), igamdst(:,:,:),
     &     idorb(:), idspn(:), idspc(:), lexlscr(:,:)

      type(filinf), pointer ::
     &     ffinp, ffinv
      type(operator), pointer ::
     &     op_inv, op_inp

      integer, external ::
     &     ielprd, idx_msgmdst2, idx_str_blk3, msa2idxms4op, idxcount,
     &     idxlist
      logical, external ::
     &     next_tupel_ca

      if (ntest.ge.100) write(luout,*) 'entered invsqrt'

      ffinp => mel_inp%fhand
      ffinv => mel_inv%fhand
      op_inp => mel_inp%op
      op_inv => mel_inv%op
      njoined = op_inp%njoined

      ngraph = str_info%ngraph
      hpvx_occ => op_inp%ihpvca_occ
      idx_graph => mel_inp%idx_graph
      ca_occ => op_inp%ica_occ
      graphs => str_info%g

      ms_fix = .false.
      if(mel_inp%fix_vertex_ms.or.mel_inv%fix_vertex_ms)then
        ms_fix = mel_inp%fix_vertex_ms.and.mel_inv%fix_vertex_ms
        if(.not.ms_fix) call quit(1,'invsqrt',
     &                            'fix ms or not?')
      endif

      ! Check whether files are buffered.
      bufin = .false.
      if(ffinp%buffered) bufin = .true.
      bufout = .false.
      if(ffinv%buffered) bufout = .true.

      ifree = mem_setmark('invsqrt')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Allocations made to maximum block length to save time.
      if(.not.bufin)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls
          if(op_inp%formal_blk(iocc_cls))
     &         cycle

          nbuff = nbuff+mel_inp%len_op_occ(iocc_cls)
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
        call get_vec(ffinp,buffer_in,1,nbuff)
      else
        if(ntest.ge.100)
     &       write(luout,*)'Invert: input not incore'
        buffer_in => ffinp%buffer(1:)
      endif

      if(.not.bufout)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls
          nbuff = nbuff + mel_inv%len_op_occ(iocc_cls)
        enddo
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
        buffer_out(1:nbuff) = 0d0
      else
        if(ntest.ge.100)
     &       write(luout,*)'Invert: output not incore'
        buffer_out => ffinv%buffer(1:)
      endif

      icnt_sv  = 0 ! we will count the
      icnt_sv0 = 0 ! number of singular values below threshold

      ! Loop over occupation class.
      iocc_loop: do iocc_cls = 1, nocc_cls
        if(op_inp%formal_blk(iocc_cls)) cycle iocc_loop
        iblkoff = (iocc_cls-1)*njoined

        if (ntest.ge.10) write(luout,*) 'current occ_cls: ',iocc_cls
        if (mel_inp%len_op_occ(iocc_cls).eq.1) then
          ioff = mel_inp%off_op_gmo(iocc_cls)%gam_ms(1,1)
          buffer_out(ioff+1) = buffer_in(ioff+1)
          cycle
        end if 

        ifree = mem_setmark('invsqrt_blk')

        call get_num_subblk(ncblk,nablk,
     &       hpvx_occ(1,1,iblkoff+1),njoined)

        ! for simplicity we distinguish the following two types:
        !  /0 0 0 0\       /0 0 0 0\  i.e. occupation classes  with
        !  \0 0 x 0/       \0 0 x 0/  several distributions per
        !  /0 0 y 0\  and  /0 0 0 0\  per Ms(A)/GAMMA(A) block
        !  \0 0 y 0/       \0 0 0 0/  and those with only single
        !  /0 0 x 0\       /0 0 x 0\  distributions
        !  \0 0 0 0/       \0 0 0 0/
        if (ncblk.eq.1.and.nablk.eq.1) then
          onedis = .true.
        else if (ncblk.eq.2.and.nablk.eq.2) then
          onedis = .false.
        else
          call quit(1,'invsqrt','not adapted for this case')
        end if

        allocate(hpvx_csub(ncblk),hpvx_asub(nablk),
     &           occ_csub(ncblk), occ_asub(nablk),
     &           graph_csub(ncblk), graph_asub(nablk),
     &           msdis_c(ncblk),  msdis_a(nablk),
     &           idxmsdis_c(ncblk),  idxmsdis_a(nablk),
     &           msdis_c2(ncblk),  msdis_a2(nablk),
     &           idxmsdis_c2(ncblk),  idxmsdis_a2(nablk),
     &           gamdis_c(ncblk), gamdis_a(nablk),
     &           len_str(ncblk+nablk),
     &           istr_csub(ncblk),istr_asub(nablk),
     &           istr_csub_flip(ncblk),istr_asub_flip(nablk),
     &           ldim_opin_c(ncblk),ldim_opin_a(nablk))

        ! set HPVX and OCC info
        call condense_occ(occ_csub, occ_asub,
     &                    hpvx_csub,hpvx_asub,
     &                    hpvx_occ(1,1,iblkoff+1),njoined,hpvxblkseq)
        ! do the same for the graph info
        call condense_occ(graph_csub, graph_asub,
     &                    hpvx_csub,hpvx_asub,
     &                    idx_graph(1,1,iblkoff+1),njoined,hpvxblkseq)

        ! set flip maps
        call strmap_man_flip(
     &       maxbuf,
     &       graph_csub,ncblk,
     &       str_info,strmap_info,orb_info)
        ifree = mem_alloc_int(flipmap_c,maxbuf,'flipmap_c')
        call strmap_man_flip(
     &       maxbuf,
     &       graph_asub,nablk,
     &       str_info,strmap_info,orb_info)
        ifree = mem_alloc_int(flipmap_a,maxbuf,'flipmap_a')

        ! simple case: only single distributions:
        if (onedis) then

          ! we also need to distinguish:
          ! /0 0 0 0\     /0 0 0 0\
          ! \0 0 x 0/     \0 0 0 0/ i.e. if creations come first,
          ! /0 0 0 0\ and /0 0 x 0\ we need to transpose the
          ! \0 0 0 0/     \0 0 x 0/ scratch matrix
          ! /0 0 x 0\     /0 0 0 0\ (important if result is non-symmetric!)
          ! \0 0 0 0/     \0 0 0 0/
          transp = .false.
          do ij = 1, njoined
            if (sum(hpvx_occ(1:ngastp,1,iblkoff+ij)).gt.0) then
              transp = .true.
              exit
            else if (sum(hpvx_occ(1:ngastp,2,iblkoff+ij)).gt.0) then
              exit
            end if
          end do
          if (transp.and.ntest.ge.100) then
            write(luout,*) 'Using transposed scratch matrix!'
          end if

          ! Loop over Ms of annihilator string.
          idxmsa = 0
          msmax = op_inp%ica_occ(1,iocc_cls)
          mscmax = op_inp%ica_occ(2,iocc_cls)
          if (msmax.ne.mscmax) call quit(1,'invsqrt',
     &            'need particle conserving operator')
          msa_loop : do msa = msmax, -msmax, -2

            idxmsa = idxmsa+1
            msc = msa + mel_inp%mst
            idxmsc = (msmax-msc)/2 + 1
      
            ! Loop over Irrep of annihilator string.
            igama_loop: do igama =1, ngam

              igamc = multd2h(igama,mel_inp%gamt)
              ndis = mel_inp%off_op_gmox(iocc_cls)%ndis(igamc,idxmsc)

              ndim = int(sqrt(dble(mel_inv%
     &           len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa))))
              if (ndim.gt.0.and.ndis.ne.1)
     &                call quit(1,'invsqrt','cannot handle this')
              if (ndim.eq.0) cycle igama_loop

              if (ntest.ge.10)
     &           write(luout,'(a,3i8)') 'msa, gama, ndim:',msa,igama,
     &           int(sqrt(dble(mel_inp%len_op_gmo(iocc_cls)%
     &                         gam_ms(igama,idxmsa))))
              if (ntest.ge.100)
     &           write(luout,*) ' len = ',
     &             mel_inp%len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa),
     &             ' ndis = ',ndis

              ioff = mel_inv%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)

              call get_flipmap_blk(flipmap_c,
     &            ncblk,occ_csub,ndim,graph_csub,idxmsc,igamc,
     &            strmap_info,ngam,ngraph)

              ! single distribution can simply be read in as simple matrix
              allocate(scratch(ndim,ndim))
              if (transp) then
                do idx = 1,ndim
                  do jdx = 1,ndim
                    scratch(jdx,idx) = buffer_in(ioff+(idx-1)*ndim+jdx)
                  enddo
                enddo
              else
                do idx = 1,ndim
                  do jdx = 1,ndim
                    scratch(idx,jdx) = buffer_in(ioff+(idx-1)*ndim+jdx)
                  enddo
                enddo
              end if

              if (msc.eq.0) then
                ! here a splitting into "singlet" and "triplet" blocks is needed:

c dbg
c                write(luout,*) 'flmap:'
c                do icol = 1, ndim
c                  write(luout,'(i4,3i6)') icol,flipmap_c(icol)
c                end do
c dbgend
                nsing = ndim
                do icol = 1, ndim
                  if (abs(flipmap_c(icol)).eq.icol) nsing = nsing
     &                                  + sign(1,flipmap_c(icol))
                end do
                nsing = nsing/2
                ntrip = ndim - nsing
c dbg
c                print *,'nsing: ',nsing
c dbgend

                ! do the pre-diagonalization
                allocate(sing(nsing,nsing),trip(ntrip,ntrip))
                call spinsym_traf(1,ndim,scratch,flipmap_c,nsing,
     &                            sing,trip,.false.)

                ! calculate T^(-0.5) for both blocks
                call invsqrt_mat(nsing,sing,half,icnt_sv,icnt_sv0)
                call invsqrt_mat(ntrip,trip,half,icnt_sv,icnt_sv0)

                ! partial undo of pre-diagonalization: Upre*T^(-0.5)
                call spinsym_traf(2,ndim,scratch,flipmap_c,nsing,
     &                            sing,trip,half)
                deallocate(sing,trip)
              else

                ! calculate S^(-0.5)
                call invsqrt_mat(ndim,scratch,half,icnt_sv,icnt_sv0)

              end if

              ! write to output buffer
              if (transp) then
                do idx = 1,ndim
                  do jdx = 1,ndim
                    buffer_out((idx-1)*ndim+jdx+ioff) = scratch(jdx,idx)
                  enddo
                enddo
              else
                do idx = 1,ndim
                  do jdx = 1,ndim
                    buffer_out((idx-1)*ndim+jdx+ioff) = scratch(idx,jdx)
                  enddo
                enddo
              end if

              deallocate(scratch)

            enddo igama_loop
          enddo msa_loop

          deallocate(hpvx_csub,hpvx_asub,
     &             occ_csub, occ_asub,
     &             graph_csub, graph_asub,
     &             msdis_c,  msdis_a,
     &             idxmsdis_c,  idxmsdis_a,
     &             msdis_c2,  msdis_a2,
     &             idxmsdis_c2,  idxmsdis_a2,
     &             gamdis_c, gamdis_a,
     &             len_str,
     &             istr_csub,istr_asub,
     &             istr_csub_flip,istr_asub_flip,
     &             ldim_opin_c,ldim_opin_a)
          ifree = mem_flushmark('invsqrt_blk')
          cycle iocc_loop
        end if

        ! Here comes the complicated part for densities with 3 vertices
        ! and multiple distributions per Ms(A)/GAMMA(A) block
        ! Here also elements from different Ms(A)/GAMMA(A) blocks couple
        ! since we need to sort the whole occupation class into
        ! vectors of certain A1/C1 tuples
        ! i.e. the matrix must be resorted from A | C to A1/C1 | A2/C2

        ! determine number of c/a-operators of A1/C1 tuple
        na1 = sum(hpvx_occ(1:ngastp,2,iblkoff+1))
        nc1 = sum(hpvx_occ(1:ngastp,1,iblkoff+2))
        na2 = sum(hpvx_occ(1:ngastp,2,iblkoff+2))
        nc2 = sum(hpvx_occ(1:ngastp,1,iblkoff+3))
        msmax_sub = na1 + nc1
        ! must be symmetric and so on and so forth
        if (njoined.ne.3.or.msmax_sub.ne.na2+nc2.or.na1.ne.nc2.or.
     &      na2.ne.nc1)
     &     call quit(1,'invsqrt','too difficult operator!')
        msmax = op_inp%ica_occ(1,iocc_cls)
        mscmax = op_inp%ica_occ(2,iocc_cls)
        if (msmax.ne.mscmax.or.msmax.ne.na1+na2.or.mscmax.ne.nc1+nc2)
     &            call quit(1,'invsqrt','need part.cons.op')

        ! occupation and so on of A2/C2 vector
        nel = na2 + nc2
        allocate(idorb(nel),idspn(nel),idspc(nel),lexlscr(nel,3),
     &           iocc(ngastp,2,2),idx_g(ngastp,2,2),msdst(ngastp,2,2),
     &           igamdst(ngastp,2,2))
        iocc = 0
        idx_g = 0
        iocc(1:ngastp,2,1) = hpvx_occ(1:ngastp,2,iblkoff+2)
        iocc(1:ngastp,1,2) = hpvx_occ(1:ngastp,1,iblkoff+3)
        idx_g(1:ngastp,2,1) = idx_graph(1:ngastp,2,iblkoff+2)
        idx_g(1:ngastp,1,2) = idx_graph(1:ngastp,1,iblkoff+3)


        ! loop over Ms/Gamma combinations of A1/C1 tuple
        ! which are the decoupled blocks of the A1/C1 | A2/C2 matrix
        ! Ms/Gamma of C2/A2 is then already defined
        ! For now we consider densities: Ms(total)=0, Gamma(total) = 1
        do ms1 = msmax_sub, -msmax_sub, -2
         ms2 = -ms1 ! particle conserving operator
         do igam = 1, ngam ! irrep of the second dimension must be the same

          ! First we need the dimension of the A1/C1 | A2/C2 block
          ndim = 0
          do msa1 = na1, -na1, -2
           msc1 = msa1 + ms1
           if (abs(msc1).gt.nc1) cycle
           msdis_c(1) = msc1
           msdis_a(1) = msa1
           do gama1 = 1, ngam
            gamc1 = multd2h(gama1,igam)
            gamdis_c(1) = gamc1
            gamdis_a(1) = gama1
            msdis_c(2) = msa1
            msdis_a(2) = msc1
            gamdis_c(2) = gama1
            gamdis_a(2) = gamc1
            call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncblk)
            call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nablk)

            call set_len_str(len_str,ncblk,nablk,
     &                       graphs,
     &                       graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                       graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                       hpvxseq,.false.)
            ndim = ndim + len_str(1)*len_str(3)
           end do
          end do

          if (ndim.eq.0) cycle

          if (ntest.ge.10) 
     &       write(luout,'(a,3i8)'),'ms1, igam, ndim:',ms1,igam,ndim

          allocate(scratch(ndim,ndim),flmap(ndim,3))
          scratch = 0d0
          flmap(1:ndim,3) = 1

          ! read in current Ms1/GAMMA block
          ! loop over Ms(A1)/GAMMA(A1) --> Ms(C1)/GAMMA(C1) already defined
          off_line = 0
          do msa1 = na1, -na1, -2
           msc1 = msa1 + ms1
           if (abs(msc1).gt.nc1) cycle
           msdis_c(1) = msc1
           msdis_a(1) = msa1
           do gama1 = 1, ngam
            gamc1 = multd2h(gama1,igam)
            gamdis_c(1) = gamc1
            gamdis_a(1) = gama1
            ! loop over Ms(C2)/GAMMA(C2) --> Ms(A2)/GAMMA(A2) already defined
            off_col = 0
            do msc2 = nc2, -nc2, -2
             msa2 = msc2 - ms2
             if (abs(msa2).gt.na2) cycle
             msdis_c(2) = msc2
             msdis_a(2) = msa2
             msa = msa1 + msa2
             idxmsa = (msmax-msa)/2 + 1
             idxmsa2 = msa2idxms4op(-msa,ms1+ms2,msmax,msmax)
             do gamc2 = 1, ngam
              gama2 = multd2h(gamc2,igam)
              gamdis_c(2) = gamc2
              gamdis_a(2) = gama2
              igama = multd2h(gama1,gama2)

              ! determine the distribution in question
              call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncblk)
              call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nablk)

              call set_len_str(len_str,ncblk,nablk,
     &                         graphs,
     &                         graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                         graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                         hpvxseq,.false.)
              lenca = ielprd(len_str,ncblk+nablk)
              if (lenca.eq.0) cycle

              ndis = mel_inp%off_op_gmox(iocc_cls)%ndis(igama,idxmsa)
              idxdis = 1
              if (ndis.gt.1)
     &           idxdis =
     &               idx_msgmdst2(
     &                iocc_cls,idxmsa,igama,
     &                occ_csub,idxmsdis_c,gamdis_c,ncblk,
     &                occ_asub,idxmsdis_a,gamdis_a,nablk,
     &                .false.,-1,-1,mel_inp,ngam)
              if (lenca.ne.mel_inp%len_op_gmox(iocc_cls)%
     &             d_gam_ms(idxdis,igama,idxmsa))
     &           call quit(1,'invsqrt','inconsistency!')

              call set_op_ldim_c(ldim_opin_c,ldim_opin_a,
     &             hpvx_csub,hpvx_asub,
     &             len_str,ncblk,nablk,.false.)

              ioff = mel_inp%off_op_gmox(iocc_cls)%
     &               d_gam_ms(idxdis,igama,idxmsa)

              ! now for spin-flipped counterpart:
              msdis_c2(1:ncblk) = -msdis_c(1:ncblk)
              msdis_a2(1:nablk) = -msdis_a(1:nablk)
              call ms2idxms(idxmsdis_c2,msdis_c2,occ_csub,ncblk)
              call ms2idxms(idxmsdis_a2,msdis_a2,occ_asub,nablk)

              ndis2 = mel_inp%off_op_gmox(iocc_cls)%ndis(igama,idxmsa2)
              idxdis2 = 1
              if (ndis2.gt.1)
     &             idxdis2 =
     &                 idx_msgmdst2(
     &                  iocc_cls,idxmsa2,igama,
     &                  occ_csub,idxmsdis_c2,gamdis_c,ncblk,
     &                  occ_asub,idxmsdis_a2,gamdis_a,nablk,
     &                  .false.,-1,-1,mel_inp,ngam)

              ioff2 = mel_inp%off_op_gmox(iocc_cls)%
     &               d_gam_ms(idxdis2,igama,idxmsa2)

              call get_flipmap_blk(flipmap_c,
     &            ncblk,occ_csub,len_str,graph_csub,idxmsdis_c,gamdis_c,
     &            strmap_info,ngam,ngraph)
              call get_flipmap_blk(flipmap_a,
     &            nablk,occ_asub,len_str(ncblk+1),
     &                                   graph_asub,idxmsdis_a,gamdis_a,
     &            strmap_info,ngam,ngraph)

              ! assemble distribution of A2/C2 tuple
              ! inelegant: we are currently using next_tupel_ca for setting up
              ! final flipmap. But there should be a way to do this using
              ! only the flipmaps from above. 
              msdst = 0
              igamdst = 1
              do idx = 1, ngastp
                if (iocc(idx,2,1).ne.0) then
                  msdst(idx,2,1) = msdis_a(2)
                  igamdst(idx,2,1) = gamdis_a(2)
                end if
                if (iocc(idx,1,2).ne.0) then
                  msdst(idx,1,2) = msdis_c(2)
                  igamdst(idx,1,2) = gamdis_c(2)
                end if
              end do

c dbg
c              write(luout,'(a,4i4)') 'ms  : ',msa1,msc1,msa2,msc2
c              write(luout,'(a,4i4)') 'gam : ',gama1,gamc1,gama2,gamc2
c              write(luout,'(a,2i4)') 'msa, igama: ',msa,igama
c              write(luout,'(a,2i4)') 'dist, len: ',idxdis, lenca
c              write(luout,'(a,1i4)') 'dist2: ',idxdis2
c              write(luout,'(a,2i4)') 'off_line/col: ',off_line,off_col
c              print *,'len1: ',len_str(1)*len_str(3)
c              print *,'flipmap_c: len=',len_str(1:ncblk)
c              print '(10i6)',flipmap_c(1:sum(len_str(1:ncblk)))
c              print *,'flipmap_a:'
c              print '(10i6)',
c     &             flipmap_a(1:sum(len_str(ncblk+1:ncblk+nablk)))
c dbgend

              ! copy all required elements of this distribution
              ! to their block in the A1/C1 | A2/C2 matrix
              first = .true.
              iline = off_line
              do idxc1 = 1, len_str(1)
                istr_csub(1) = idxc1-1
                istr_csub_flip(1) = abs(flipmap_c(idxc1))-1
                do idxa1 = 1, len_str(3)
                  istr_asub(1) = idxa1-1
                  istr_asub_flip(1) = abs(flipmap_a(idxa1))-1
                  iline = iline + 1
                  icol = off_col
                  do idxa2 = 1, len_str(4)
                    istr_asub(2) = idxa2-1
                    istr_asub_flip(2) = 
     &                       abs(flipmap_a(len_str(3)+idxa2))-1
                    do idxc2 = 1, len_str(2)
                      istr_csub(2) = idxc2-1
                      istr_csub_flip(2) =
     &                       abs(flipmap_c(len_str(1)+idxc2))-1
                      icol = icol + 1
                      idx = ioff + idx_str_blk3(istr_csub,istr_asub,
     &                       ldim_opin_c,ldim_opin_a,ncblk,nablk)
                      scratch(iline,icol) = buffer_in(idx)
                      if (ms1.eq.0.and.iline.eq.icol) then
                        idx2 = ioff2 + idx_str_blk3(istr_csub_flip,
     &                         istr_asub_flip,
     &                         ldim_opin_c,ldim_opin_a,ncblk,nablk)
                        flmap(icol,1) = idx
                        flmap(icol,2) = idx2
                        logdum = next_tupel_ca(idorb,idspn,idspc,
     &                   nel,2,iocc,idx_g,
     &                   msdst,igamdst,first,
     &                   str_info%igas_restr,
     &                   orb_info%mostnd,orb_info%igamorb,
     &                   orb_info%nsym,orb_info%ngas,
     &                   orb_info%ngas_hpv,orb_info%idx_gas,
     &                   hpvxseq,lexlscr)
                        first = .false.
                        if (.not.logdum) call quit(1,'invsqrt',
     &                       'no next tuple found!')
                        if (mod(idxcount(2,idspn,nel,1),4).ne.0)
     &                        flmap(icol,3) = -1
c dbg
c                        write(luout,'(i4,x,2i4,x,2i4)')idx,idorb,idspn
c dbgend
                      end if
                    end do
                  end do
                end do
              end do

              off_col = off_col + len_str(2)*len_str(4)
             end do
            end do
            off_line = off_line + len_str(1)*len_str(3)
           end do
          end do

          if (ms1.eq.0) then
            ! here a splitting into "singlet" and "triplet" blocks is needed:

            do icol = 1, ndim
              idx = idxlist(flmap(icol,1),flmap(1:ndim,2),ndim,1)
              if (idx.eq.-1) call quit(1,'invsqrt','idx not found!')
              flmap(icol,3) = flmap(icol,3)*idx
            end do
c dbg
c            write(luout,*) 'flmap:'
c            do icol = 1, ndim
c              write(luout,'(i4,3i6)') icol,flmap(icol,1:3)
c            end do
c dbgend
            nsing = ndim
            do icol = 1, ndim
              if (abs(flmap(icol,3)).eq.icol) nsing = nsing
     &                              + sign(1,flmap(icol,3))
            end do
            nsing = nsing/2
            ntrip = ndim - nsing
c dbg
c            print *,'nsing: ',nsing
c dbgend

            ! do the pre-diagonalization
            allocate(sing(nsing,nsing),trip(ntrip,ntrip))
            call spinsym_traf(1,ndim,scratch,flmap(1:ndim,3),nsing,
     &                        sing,trip,.false.)

            ! calculate T^(-0.5) for both blocks
            call invsqrt_mat(nsing,sing,half,icnt_sv,icnt_sv0)
            call invsqrt_mat(ntrip,trip,half,icnt_sv,icnt_sv0)

            ! partial undo of pre-diagonalization: Upre*T^(-0.5)
            call spinsym_traf(2,ndim,scratch,flmap(1:ndim,3),nsing,
     &                        sing,trip,half)
            deallocate(sing,trip)
          else

            ! calculate S^(-0.5)
            call invsqrt_mat(ndim,scratch,half,icnt_sv,icnt_sv0)

          end if

          ! write to output buffer
          ! loop over Ms(A1)/GAMMA(A1) --> Ms(C1)/GAMMA(C1) already defined
          off_line = 0
          do msa1 = na1, -na1, -2
           msc1 = msa1 + ms1
           if (abs(msc1).gt.nc1) cycle
           msdis_c(1) = msc1
           msdis_a(1) = msa1
           do gama1 = 1, ngam
            gamc1 = multd2h(gama1,igam)
            gamdis_c(1) = gamc1
            gamdis_a(1) = gama1
            ! loop over Ms(C2)/GAMMA(C2) --> Ms(A2)/GAMMA(A2) already defined
            off_col = 0
            do msc2 = nc2, -nc2, -2
             msa2 = msc2 - ms2
             if (abs(msa2).gt.na2) cycle
             msdis_c(2) = msc2
             msdis_a(2) = msa2
             msa = msa1 + msa2
             idxmsa = (msmax-msa)/2 + 1
             do gamc2 = 1, ngam
              gama2 = multd2h(gamc2,igam)
              gamdis_c(2) = gamc2
              gamdis_a(2) = gama2
              igama = multd2h(gama1,gama2)

              ! determine the distribution in question
              call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncblk)
              call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nablk)

              call set_len_str(len_str,ncblk,nablk,
     &                         graphs,
     &                         graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                         graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                         hpvxseq,.false.)
              if (ielprd(len_str,ncblk+nablk).eq.0) cycle

              ndis = mel_inp%off_op_gmox(iocc_cls)%ndis(igama,idxmsa)
              idxdis = 1
              if (ndis.gt.1)
     &           idxdis =
     &               idx_msgmdst2(
     &                iocc_cls,idxmsa,igama,
     &                occ_csub,idxmsdis_c,gamdis_c,ncblk,
     &                occ_asub,idxmsdis_a,gamdis_a,nablk,
     &                .false.,-1,-1,mel_inp,ngam)

              call set_op_ldim_c(ldim_opin_c,ldim_opin_a,
     &             hpvx_csub,hpvx_asub,
     &             len_str,ncblk,nablk,.false.)

              ioff = mel_inp%off_op_gmox(iocc_cls)%
     &               d_gam_ms(idxdis,igama,idxmsa)

              ! copy all required elements of this distribution
              ! to their block in the A1/C1 | A2/C2 matrix
              iline = off_line
              do idxc1 = 1, len_str(1)
                istr_csub(1) = idxc1-1
                do idxa1 = 1, len_str(3)
                  istr_asub(1) = idxa1-1
                  iline = iline + 1
                  icol = off_col
                  do idxa2 = 1, len_str(4)
                    istr_asub(2) = idxa2-1
                    do idxc2 = 1, len_str(2)
                      istr_csub(2) = idxc2-1
                      icol = icol + 1
                      idx = ioff + idx_str_blk3(istr_csub,istr_asub,
     &                       ldim_opin_c,ldim_opin_a,ncblk,nablk)
                      buffer_out(idx) = scratch(iline,icol)
                    end do
                  end do
                end do
              end do

              off_col = off_col + len_str(2)*len_str(4)
             end do
            end do
            off_line = off_line + len_str(1)*len_str(3)
           end do
          end do

          deallocate(scratch,flmap)

         end do
        end do

        deallocate(hpvx_csub,hpvx_asub,
     &           occ_csub, occ_asub,
     &           graph_csub, graph_asub,
     &           msdis_c,  msdis_a,
     &           idxmsdis_c,  idxmsdis_a,
     &           msdis_c2,  msdis_a2,
     &           idxmsdis_c2,  idxmsdis_a2,
     &           gamdis_c, gamdis_a,
     &           len_str,
     &           istr_csub,istr_asub,
     &           istr_csub_flip,istr_asub_flip,
     &           ldim_opin_c,ldim_opin_a,
     &           idorb,idspn,idspc,lexlscr,
     &           iocc,idx_g,msdst,
     &           igamdst)

        ifree = mem_flushmark('invsqrt_blk')

      enddo iocc_loop

      if (ntest.ge.10) write(luout,'(i8,a,i8,a)') icnt_sv0,' out of ',
     &        icnt_sv,' singular values were below threshold'
 

      if(.not.bufout)then
        call put_vec(ffinv,buffer_out,1,nbuff)
      endif  

      ifree = mem_flushmark('invsqrt')

      return
      end
