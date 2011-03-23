      subroutine invsqrt(mel_inp,mel_inv,nocc_cls,half,sgrm,
     &     op_info,orb_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     Routine to calculate S^(-0.5) of density matrices.
*     The ME list is split into matrices that can either contain
*     elements from a single distribution or from several distributions.
*     half = true: only makes a half transform U*s^(-0.5)
*     half = false: also does half transform and in addition
*                   returns projector matrix U*1s*U^+ on input list
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
     &     half,sgrm

      type(graph), pointer ::
     &     graphs(:)

      logical ::
     &     bufin, bufout, first, ms_fix, fix_success, onedis, transp,
     &     logdum, sing_remove, normalize
      logical, pointer ::
     &     blk_used(:)
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
     &     nsing, ising, itrip, ntrip, icnt_sv, icnt_sv0,
     &     bins(17), nalph, nbeta, nc1mx, na1mx, jocc_cls, jblkoff,
     &     ncblk2, nablk2, len2(4), off_line2, off_col2, off_colmax,
c dbg
c     &     ipass,
c dbgend
     &     off_linmax, maxbuf_tmp
      real(8) ::
     &     fac, xmax, xmin
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:), scratch(:,:), scratch2(:,:),
     &     sing(:,:), trip(:,:), sing2(:,:), trip2(:,:),
     &     palph(:), pbeta(:), proj_sing(:,:), norm(:), scratch3(:,:)
c dbg
c      integer, pointer ::
c     &     matrix(:,:)
c dbgend

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
     &     idorb(:), idspn(:), idspc(:), lexlscr(:,:),
     &     hpvx_csub2(:),hpvx_asub2(:),
     &     occ_csub2(:), occ_asub2(:),
     &     graph_csub2(:), graph_asub2(:),
     &     iocc2(:,:,:),idx_g2(:,:,:)

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
c dbg
c      ipass = 0
c dbgend

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

      normalize = .false.
c dbg
      if (normalize) write(luout,*)  'Normalizing Metric!'
c dbgend
      icnt_sv  = 0 ! we will count the
      icnt_sv0 = 0 ! number of singular values below threshold
      xmax = 0d0   ! largest excluded singular value
      xmin = 1234567890d0   ! smallest included singular value
      bins = 0 ! binning for singular values:
               ! >10E0,>10E-1,...,>10E-15,0
      allocate(blk_used(nocc_cls))
      blk_used(1:nocc_cls) = .false.

      if (.not.half.and.max(iprlvl,ntest).ge.3) write(luout,*)
     &         'Input list will be overwritten by projector.'

      ! Loop over occupation class.
      iocc_loop: do iocc_cls = 1, nocc_cls
        if(op_inp%formal_blk(iocc_cls)) cycle iocc_loop
        if (blk_used(iocc_cls)) cycle iocc_loop
        blk_used(iocc_cls) = .true.
        iblkoff = (iocc_cls-1)*njoined

        if (ntest.ge.10) write(luout,*) 'current occ_cls: ',iocc_cls
        if (mel_inp%len_op_occ(iocc_cls).eq.1) then
          ioff = mel_inp%off_op_gmo(iocc_cls)%gam_ms(1,1)
          buffer_out(ioff+1) = buffer_in(ioff+1)
          call invsqrt_mat(1,buffer_out(ioff+1),buffer_in(ioff+1),
     &                     half,icnt_sv,icnt_sv0,xmax,xmin,bins)
c          if (.not.half) buffer_in(ioff+1) = 1d0
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
c dbg
c        print *,'maxbuf flipmap_c:',maxbuf
        maxbuf = maxbuf + 5
c dbgend
        ifree = mem_alloc_int(flipmap_c,maxbuf,'flipmap_c')
        call strmap_man_flip(
     &       maxbuf,
     &       graph_asub,nablk,
     &       str_info,strmap_info,orb_info)
c dbg
c        print *,'maxbuf flipmap_a:',maxbuf
        maxbuf = maxbuf + 5
c dbgend
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
              if (.not.half) allocate(scratch2(ndim,ndim))
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

              ! normalization
              if (normalize) then
                allocate(norm(ndim))
                do idx = 1, ndim
                  norm(idx) = sqrt(abs(scratch(idx,idx)))
                end do
                call dmdiagm(ndim,scratch,norm,.true.,.true.)
                call dmdiagm(ndim,scratch,norm,.false.,.true.)
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
                if (.not.half) allocate(sing2(nsing,nsing),
     &                                  trip2(ntrip,ntrip))
                call spinsym_traf(1,ndim,scratch,flipmap_c,nsing,
     &                            sing,trip,.false.)

                ! calculate T^(-0.5) for both blocks
                call invsqrt_mat(nsing,sing,sing2,half,icnt_sv,icnt_sv0,
     &                           xmax,xmin,bins)
                call invsqrt_mat(ntrip,trip,trip2,half,icnt_sv,icnt_sv0,
     &                           xmax,xmin,bins)

                ! partial undo of pre-diagonalization: Upre*T^(-0.5)
                call spinsym_traf(2,ndim,scratch,flipmap_c,nsing,
     &                            sing,trip,.true.)
                if (.not.half) then
                  ! full undo of pre-diagonalization for projector
                  call spinsym_traf(2,ndim,scratch2,flipmap_c,nsing,
     &                              sing2,trip2,.false.)
                  deallocate(sing2,trip2)
                end if
                deallocate(sing,trip)
              else

                ! calculate S^(-0.5)
                call invsqrt_mat(ndim,scratch,scratch2,
     &                           half,icnt_sv,icnt_sv0,xmax,xmin,bins)

              end if

              ! "undo" normalization
              if (normalize) then ! X = N*X'
                call dmdiagm(ndim,scratch,norm,.false.,.true.)
                if (.not.half) then ! P = N*P'*N^-1
                  call dmdiagm(ndim,scratch2,norm,.false.,.true.)
                  call dmdiagm(ndim,scratch2,norm,.true.,.false.)
                end if
                deallocate(norm)
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

              if (.not.half) then
                ! write projector to input buffer
                if (transp) then
                  do idx = 1,ndim
                    do jdx = 1,ndim
                      buffer_in((idx-1)*ndim+jdx+ioff) 
     &                       = scratch2(jdx,idx)
                    enddo
                  enddo
                else
                  do idx = 1,ndim
                    do jdx = 1,ndim
                      buffer_in((idx-1)*ndim+jdx+ioff) 
     &                       = scratch2(idx,jdx)
                    enddo
                  enddo
                end if

                deallocate(scratch2)
              end if

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
        na1mx = sum(hpvx_occ(1:ngastp,2,iblkoff+1))
        nc1mx = sum(hpvx_occ(1:ngastp,1,iblkoff+2))
        na2 = sum(hpvx_occ(1:ngastp,2,iblkoff+2))
        nc2 = sum(hpvx_occ(1:ngastp,1,iblkoff+3))
        msmax_sub = na1mx + nc1mx
        ! must be symmetric and so on and so forth
        if (njoined.ne.3.or.msmax_sub.ne.na2+nc2.or.na1mx.ne.nc2.or.
     &      na2.ne.nc1mx)
     &     call quit(1,'invsqrt','too difficult operator!')
        msmax = op_inp%ica_occ(1,iocc_cls)
        mscmax = op_inp%ica_occ(2,iocc_cls)
        if (msmax.ne.mscmax.or.msmax.ne.na1mx+na2
     &      .or.mscmax.ne.nc1mx+nc2)
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
          sing_remove = .false.
          ! loops over coupling blocks. Must be in correct order!
          jocc_cls = iocc_cls
          jblkoff = (jocc_cls-1)*njoined
          na1 = na1mx
          nc1 = nc1mx
          blk_loop: do while(min(na1,nc1).ge.0.and.na1+nc1.ge.abs(ms1))
           na2 = nc1mx
           nc2 = na1mx
           do while(min(na2,nc2).ge.0.and.na2+nc2.ge.abs(ms1))
            ! exit if there is no block like this
            if (na1.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+1)).or.
     &         nc1.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+2)).or.
     &         na2.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+2)).or.
     &         nc2.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+3))) exit blk_loop
            ! skip off-diagonal blocks
            if (na1.ne.nc2) then
              na2 = na2 - 1
              nc2 = nc2 - 1
              jocc_cls = jocc_cls + 1
              if (na2+nc2.lt.abs(ms1)) !jump to next line
     &                  jocc_cls = jocc_cls + min(na1mx,nc1mx)
     &                  -(na1mx+nc1mx-abs(ms1))/2
              jblkoff = (jocc_cls-1)*njoined
              cycle
            end if
            ! exception for pure inactive block:
            if (na1.eq.0.and.nc1.eq.0.and.igam.eq.1) then
              ndim = ndim + 1
              exit blk_loop
            end if
            call get_num_subblk(ncblk2,nablk2,
     &           hpvx_occ(1,1,jblkoff+1),njoined)
            allocate(hpvx_csub2(ncblk2),hpvx_asub2(nablk2),
     &               occ_csub2(ncblk2), occ_asub2(nablk2),
     &               graph_csub2(ncblk2), graph_asub2(nablk2))
            ! set HPVX and OCC info
            call condense_occ(occ_csub2, occ_asub2,
     &                      hpvx_csub2,hpvx_asub2,
     &                      hpvx_occ(1,1,jblkoff+1),njoined,hpvxblkseq)
            ! do the same for the graph info
            call condense_occ(graph_csub2, graph_asub2,
     &                      hpvx_csub2,hpvx_asub2,
     &                      idx_graph(1,1,jblkoff+1),njoined,hpvxblkseq)

c           ! First we need the dimension of the A1/C1 | A2/C2 block
c           ndim = 0
            do msa1 = na1, -na1, -2
             msc1 = msa1 + ms1
             if (abs(msc1).gt.nc1) cycle
             msdis_c(1) = msc1
             msdis_a(1) = msa1
             do gama1 = 1, ngam
              gamc1 = multd2h(gama1,igam)
              if (nablk2.eq.1.and.na2.eq.0) then ! diag.: nablk2=ncblk2
                if (msc1.ne.0.or.gamc1.ne.1) cycle
                gamdis_c(1) = gama1
                msdis_c(1) = msa1
                gamdis_a(1) = gama1
              else if (nablk2.eq.1.and.na1.eq.0) then
                if (msa1.ne.0.or.gama1.ne.1) cycle
                gamdis_a(1) = gamc1
                msdis_a(1) = msc1
                gamdis_c(1) = gamc1
              else
                gamdis_c(1) = gamc1
                gamdis_a(1) = gama1
                msdis_c(2) = msa1
                msdis_a(2) = msc1
                gamdis_c(2) = gama1
                gamdis_a(2) = gamc1
              end if
              call ms2idxms(idxmsdis_c,msdis_c,occ_csub2,ncblk2)
              call ms2idxms(idxmsdis_a,msdis_a,occ_asub2,nablk2)

              call set_len_str(len_str,ncblk2,nablk2,
     &                       graphs,
     &                       graph_csub2,idxmsdis_c,gamdis_c,hpvx_csub2,
     &                       graph_asub2,idxmsdis_a,gamdis_a,hpvx_asub2,
     &                       hpvxseq,.false.)
              if (nablk2.eq.1) then
                ndim = ndim + len_str(1)
              else
                ndim = ndim + len_str(1)*len_str(3)
              end if
             end do
            end do

            ! remove singles only for one specific block for now
            if (sing_remove) call quit(1,'invsqrt','unexpected block')
            if (sgrm.and.ms1.eq.0.and.igam.eq.1
     &                    .and.na1.eq.nc1.and.na1.eq.1) then
              sing_remove = .true.
              allocate(palph(ndim),pbeta(ndim),proj_sing(ndim,ndim))
              palph = 0d0
              pbeta = 0d0
              nalph = 0
              nbeta = 0
              ! get pseudo-singles vectors that should be removed
              iline = 0
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
                call ms2idxms(idxmsdis_c,msdis_c,occ_csub2,ncblk2)
                call ms2idxms(idxmsdis_a,msdis_a,occ_asub2,nablk2)

                call set_len_str(len_str,ncblk2,nablk2,
     &                       graphs,
     &                       graph_csub2,idxmsdis_c,gamdis_c,hpvx_csub2,
     &                       graph_asub2,idxmsdis_a,gamdis_a,hpvx_asub2,
     &                       hpvxseq,.false.)
                if (len_str(1).ne.len_str(3)) call quit(1,'invsqrt',
     &              'this should not happen.')
                ! simple square matrix, the diagonal elements
                ! contribute to the pseudo-singles component
                do idxc1 = 1, len_str(1)
                 do idxa1 = 1, len_str(3)
                  iline = iline + 1
                  if (idxc1.ne.idxa1) cycle
ctest
                  palph(iline) = 1d0/2d0 !Nact=2 hard-coded
ctest                  if (msa1.eq.1) then
ctest                    palph(iline) = 1d0
ctest                    nalph = nalph + 1
ctest                  else if (msa1.eq.-1) then
ctest                    pbeta(iline) = 1d0
ctest                    nbeta = nbeta + 1
ctest                  else
ctest                    call quit(1,'invsqrt','this should not happen(2).')
ctest                  end if
                 end do
                end do
               end do
              end do
c              ! normalize vectors
c              if (nalph.gt.0) palph = palph/sqrt(dble(nalph))
c              if (nbeta.gt.0) pbeta = pbeta/sqrt(dble(nbeta))
c              ! get projector
c              proj_sing = 0d0
c              do idx = 1, ndim
c                proj_sing(idx,idx) = 1d0
c              end do
c              call dgemm('n','t',ndim,ndim,1,
c     &                   -1d0,palph,ndim,
c     &                   palph,ndim,
c     &                   1d0,proj_sing,ndim)
c              call dgemm('n','t',ndim,ndim,1,
c     &                   -1d0,pbeta,ndim,
c     &                   pbeta,ndim,
c     &                   1d0,proj_sing,ndim)
c              if (ntest.ge.100) then
c                write(luout,*) 'projector for removing singles:'
c                call wrtmat2(proj_sing,ndim,ndim,ndim,ndim)
c              end if
            end if

            deallocate(hpvx_csub2,hpvx_asub2,occ_csub2,
     &               occ_asub2,graph_csub2,graph_asub2)
            na2 = na2 - 1
            nc2 = nc2 - 1
            jocc_cls = jocc_cls + 1
            if (na2+nc2.lt.abs(ms1)) !jump to next line
     &                jocc_cls = jocc_cls + min(na1mx,nc1mx)
     &                -(na1mx+nc1mx-abs(ms1))/2
            jblkoff = (jocc_cls-1)*njoined
           end do
           na1 = na1 - 1
           nc1 = nc1 - 1
          end do blk_loop

          if (ndim.eq.0) cycle

          if (ntest.ge.10)
     &       write(luout,'(a,3i8)'),'ms1, igam, ndim:',ms1,igam,ndim

          allocate(scratch(ndim,ndim),flmap(ndim,3))
          scratch = 0d0
          if (.not.half.or.sing_remove) allocate(scratch2(ndim,ndim))
c dbg
c          allocate(matrix(ndim,ndim))
c dbgend

          flmap(1:ndim,1:3) = 1
          ! loops over coupling blocks. Must be in correct order!
          jocc_cls = iocc_cls
          jblkoff = (jocc_cls-1)*njoined
          na1 = na1mx
          nc1 = nc1mx
          off_line2 = 0
          do while(min(na1,nc1).ge.0.and.na1+nc1.ge.abs(ms1))
           na2 = nc1mx
           nc2 = na1mx
           off_col2 = 0
           off_linmax = 0
           do while(min(na2,nc2).ge.0.and.na2+nc2.ge.abs(ms1))
            ! exit if there is no block like this
            if (na1.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+1)).or.
     &          nc1.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+2)).or.
     &          na2.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+2)).or.
     &          nc2.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+3))) exit
            blk_used(jocc_cls) = .true.
            ! exception for pure inactive block:
            if (na1+nc1+na2+nc2.eq.0.and.igam.eq.1) then
              ioff = mel_inp%off_op_gmox(jocc_cls)%
     &                 d_gam_ms(1,1,1)
              scratch(off_line2+1,off_col2+1) = buffer_in(ioff+1)
              flmap(off_col2+1,1:2) = ioff+1
              flmap(off_col2+1,3) = 1
              exit
            end if
c dbg
c            print *,'na1,nc1,na2,nc2,jocc_cls:',na1,nc1,na2,nc2,
c     &          jocc_cls
c            print *,'off_line2, off_col2: ',off_line2,off_col2
c dbgend
            msmax = op_inp%ica_occ(1,jocc_cls)
            call get_num_subblk(ncblk2,nablk2,
     &           hpvx_occ(1,1,jblkoff+1),njoined)
            allocate(hpvx_csub2(ncblk2),hpvx_asub2(nablk2),
     &               occ_csub2(ncblk2), occ_asub2(nablk2),
     &               graph_csub2(ncblk2), graph_asub2(nablk2),
     &               iocc2(ngastp,2,2),idx_g2(ngastp,2,2))
            iocc2 = 0
            idx_g2 = 0
            iocc2(1:ngastp,2,1) = hpvx_occ(1:ngastp,2,jblkoff+2)
            iocc2(1:ngastp,1,2) = hpvx_occ(1:ngastp,1,jblkoff+3)
            idx_g2(1:ngastp,2,1) = idx_graph(1:ngastp,2,jblkoff+2)
            idx_g2(1:ngastp,1,2) = idx_graph(1:ngastp,1,jblkoff+3)
            ! set HPVX and OCC info
            call condense_occ(occ_csub2, occ_asub2,
     &                      hpvx_csub2,hpvx_asub2,
     &                      hpvx_occ(1,1,jblkoff+1),njoined,hpvxblkseq)
            ! do the same for the graph info
            call condense_occ(graph_csub2, graph_asub2,
     &                      hpvx_csub2,hpvx_asub2,
     &                      idx_graph(1,1,jblkoff+1),njoined,hpvxblkseq)

            ! read in current Ms1/GAMMA block
            ! loop over Ms(A1)/GAMMA(A1) --> Ms(C1)/GAMMA(C1) already defined
            off_line = off_line2
            off_colmax = 0
            do msa1 = na1, -na1, -2
             msc1 = msa1 + ms1
             if (abs(msc1).gt.nc1) cycle
             msdis_c(1) = msc1
             msdis_a(1) = msa1
c dbg
c             print *,'msa1,msc1: ',msa1,msc1
c dbgend
             do gama1 = 1, ngam
              gamc1 = multd2h(gama1,igam)
              if (na1.eq.0.and.gama1.ne.1.or.
     &            nc1.eq.0.and.gamc1.ne.1) cycle
              gamdis_c(1) = gamc1
              gamdis_a(1) = gama1
              ! loop over Ms(C2)/GAMMA(C2) --> Ms(A2)/GAMMA(A2) already defined
              off_col = off_col2
c dbg
c              print *,'gama1,gamc1: ',gama1,gamc1
c dbgend
              do msc2 = nc2, -nc2, -2
               msa2 = msc2 - ms2
               if (abs(msa2).gt.na2) cycle
               msa = msa1 + msa2
               idxmsa = (msmax-msa)/2 + 1
               idxmsa2 = msa2idxms4op(-msa,ms1+ms2,msmax,msmax)
c dbg
c               print *,'msc2,msa2: ',msc2,msa2
c dbgend
               do gamc2 = 1, ngam
                gama2 = multd2h(gamc2,igam)
                if (na2.eq.0.and.gama2.ne.1.or.
     &              nc2.eq.0.and.gamc2.ne.1) cycle
                igama = multd2h(gama1,gama2)
c dbg
c                print *,'gamc2,gama2: ',gamc2,gama2
c                print *,'off_line, off_col: ',off_line,off_col
c dbgend
                if (nablk2.eq.1.and.na1.eq.0) then
                  msdis_a(1) = msa2
                  gamdis_a(1) = gama2
                else if (nablk2.eq.2) then
                  msdis_a(2) = msa2
                  gamdis_a(2) = gama2
                end if
                if (ncblk2.eq.1.and.nc1.eq.0) then
                  msdis_c(1) = msc2
                  gamdis_c(1) = gamc2
                else if (ncblk2.eq.2) then
                  msdis_c(2) = msc2
                  gamdis_c(2) = gamc2
                end if

                ! determine the distribution in question
                call ms2idxms(idxmsdis_c,msdis_c,occ_csub2,ncblk2)
                call ms2idxms(idxmsdis_a,msdis_a,occ_asub2,nablk2)

                call set_len_str(len_str,ncblk2,nablk2,
     &                       graphs,
     &                       graph_csub2,idxmsdis_c,gamdis_c,hpvx_csub2,
     &                       graph_asub2,idxmsdis_a,gamdis_a,hpvx_asub2,
     &                       hpvxseq,.false.)
                lenca = ielprd(len_str,ncblk2+nablk2)
                len2(1:4) = 1
                if (ncblk2.eq.1.and.nc2.eq.0) len2(1) = len_str(1)
                if (ncblk2.eq.1.and.nc1.eq.0) len2(2) = len_str(1)
                if (ncblk2.eq.2) then
                  len2(1) = len_str(1)
                  len2(2) = len_str(2)
                end if
                if (nablk2.eq.1.and.na2.eq.0) len2(3) =len_str(ncblk2+1)
                if (nablk2.eq.1.and.na1.eq.0) len2(4) =len_str(ncblk2+1)
                if (nablk2.eq.2) then
                  len2(3) = len_str(ncblk2+1)
                  len2(4) = len_str(ncblk2+2)
                end if
                if (lenca.eq.0) cycle

                ndis = mel_inp%off_op_gmox(jocc_cls)%ndis(igama,idxmsa)
                idxdis = 1
                if (ndis.gt.1)
     &             idxdis =
     &                 idx_msgmdst2(
     &                  jocc_cls,idxmsa,igama,
     &                  occ_csub2,idxmsdis_c,gamdis_c,ncblk2,
     &                  occ_asub2,idxmsdis_a,gamdis_a,nablk2,
     &                  .false.,-1,-1,mel_inp,ngam)
c dbg
c                print *,'igama,idxmsa: ',igama,idxmsa
c                print *,'lenca,ndis,idxdis,len:',lenca,ndis,idxdis,
c     &           mel_inp%len_op_gmox(jocc_cls)%
c     &               d_gam_ms(idxdis,igama,idxmsa)
c dbgend
                if (lenca.ne.mel_inp%len_op_gmox(jocc_cls)%
     &               d_gam_ms(idxdis,igama,idxmsa))
     &             call quit(1,'invsqrt','inconsistency!')

                call set_op_ldim_c(ldim_opin_c,ldim_opin_a,
     &               hpvx_csub2,hpvx_asub2,
     &               len_str,ncblk2,nablk2,.false.)

                ioff = mel_inp%off_op_gmox(jocc_cls)%
     &                 d_gam_ms(idxdis,igama,idxmsa)

                ! now for spin-flipped counterpart:
                msdis_c2(1:ncblk2) = -msdis_c(1:ncblk2)
                msdis_a2(1:nablk2) = -msdis_a(1:nablk2)
                call ms2idxms(idxmsdis_c2,msdis_c2,occ_csub2,ncblk2)
                call ms2idxms(idxmsdis_a2,msdis_a2,occ_asub2,nablk2)

                ndis2 =mel_inp%off_op_gmox(jocc_cls)%ndis(igama,idxmsa2)
                idxdis2 = 1
                if (ndis2.gt.1)
     &               idxdis2 =
     &                   idx_msgmdst2(
     &                    jocc_cls,idxmsa2,igama,
     &                    occ_csub2,idxmsdis_c2,gamdis_c,ncblk2,
     &                    occ_asub2,idxmsdis_a2,gamdis_a,nablk2,
     &                    .false.,-1,-1,mel_inp,ngam)

                ioff2 = mel_inp%off_op_gmox(jocc_cls)%
     &                 d_gam_ms(idxdis2,igama,idxmsa2)

                ! set flip maps (if not existing yet)
                call strmap_man_flip(
     &               maxbuf_tmp,graph_csub2,ncblk2,
     &               str_info,strmap_info,orb_info)
                if (maxbuf_tmp.gt.maxbuf) call quit(1,'invsqrt',
     &             'add more to maxbuf')
                call strmap_man_flip(
     &               maxbuf_tmp,graph_asub2,nablk2,
     &               str_info,strmap_info,orb_info)
                if (maxbuf_tmp.gt.maxbuf) call quit(1,'invsqrt',
     &             'add more to maxbuf')

                call get_flipmap_blk(flipmap_c,
     &              ncblk2,occ_csub2,len_str,
     &              graph_csub2,idxmsdis_c,gamdis_c,
     &              strmap_info,ngam,ngraph)
                call get_flipmap_blk(flipmap_a,
     &              nablk2,occ_asub2,len_str(ncblk2+1),
     &              graph_asub2,idxmsdis_a,gamdis_a,
     &              strmap_info,ngam,ngraph)

                ! assemble distribution of A2/C2 tuple
                ! inelegant: we are currently using next_tupel_ca for setting up
                ! final flipmap. But there should be a way to do this using
                ! only the flipmaps from above. 
                msdst = 0
                igamdst = 1
                do idx = 1, ngastp
                  if (iocc2(idx,2,1).ne.0) then
                    msdst(idx,2,1) = msdis_a(2)
                    igamdst(idx,2,1) = gamdis_a(2)
                  end if
                  if (iocc2(idx,1,2).ne.0) then
                    msdst(idx,1,2) = msdis_c(2)
                    igamdst(idx,1,2) = gamdis_c(2)
                  end if
                end do

c dbg
c                write(luout,'(a,4i4)') 'ms  : ',msa1,msc1,msa2,msc2
c                write(luout,'(a,4i4)') 'gam : ',gama1,gamc1,gama2,gamc2
c                write(luout,'(a,2i4)') 'msa, igama: ',msa,igama
c                write(luout,'(a,2i4)') 'dist, len: ',idxdis, lenca
c                write(luout,'(a,1i4)') 'dist2: ',idxdis2
c                write(luout,'(a,2i4)') 'off_line/col: ',off_line,off_col
c                print *,'len1: ',len_str(1)*len_str(3)
c                print *,'flipmap_c: len=',len_str(1:ncblk2)
c                print '(10i6)',flipmap_c(1:sum(len_str(1:ncblk2)))
c                print *,'flipmap_a:'
c                print '(10i6)',
c     &               flipmap_a(1:sum(len_str(ncblk2+1:ncblk2+nablk2)))
c dbgend

                ! copy all required elements of this distribution
                ! to their block in the A1/C1 | A2/C2 matrix
                first = .true.
                iline = off_line
                do idxc1 = 1, len2(1)
                  if (nc1.ne.0) then
                   istr_csub(1) = idxc1-1
                   istr_csub_flip(1) = abs(flipmap_c(idxc1))-1
                  end if
                  do idxa1 = 1, len2(3)
                    if (na1.ne.0) then
                     istr_asub(1) = idxa1-1
                     istr_asub_flip(1) = abs(flipmap_a(idxa1))-1
                    end if
                    iline = iline + 1
                    icol = off_col
                    do idxa2 = 1, len2(4)
                      if (na1.eq.0.and.na2.ne.0) then
                       istr_asub(1) = idxa2-1
                       istr_asub_flip(1) = 
     &                          abs(flipmap_a(idxa2))-1
c     &                          abs(flipmap_a(len2(3)+idxa2))-1
                      else if (na2.ne.0) then
                       istr_asub(2) = idxa2-1
                       istr_asub_flip(2) =
     &                          abs(flipmap_a(len2(3)+idxa2))-1
                      end if
                      do idxc2 = 1, len2(2)
                        if (nc1.eq.0.and.nc2.ne.0) then
                         istr_csub(1) = idxc2-1
                         istr_csub_flip(1) =
     &                          abs(flipmap_c(idxc2))-1
c     &                          abs(flipmap_c(len2(1)+idxc2))-1
                        else if (nc2.ne.0) then
                         istr_csub(2) = idxc2-1
                         istr_csub_flip(2) =
     &                          abs(flipmap_c(len2(1)+idxc2))-1
                        end if
                        icol = icol + 1
                        idx = ioff + idx_str_blk3(istr_csub,istr_asub,
     &                         ldim_opin_c,ldim_opin_a,ncblk2,nablk2)
                        scratch(iline,icol) = buffer_in(idx)
c dbg
c                        matrix(iline,icol) = idx
c dbgend
                        if (ms1.eq.0.and.iline.eq.icol) then
                          idx2 = ioff2 + idx_str_blk3(istr_csub_flip,
     &                           istr_asub_flip,
     &                           ldim_opin_c,ldim_opin_a,ncblk2,nablk2)
                          flmap(icol,1) = idx
                          flmap(icol,2) = idx2
                          logdum = next_tupel_ca(idorb,idspn,idspc,
     &                     nel,2,iocc2,idx_g2,
     &                     msdst,igamdst,first,
     &                     str_info%igas_restr,
     &                     orb_info%mostnd,orb_info%igamorb,
     &                     orb_info%nsym,orb_info%ngas,
     &                     orb_info%ngas_hpv,orb_info%idx_gas,
     &                     hpvxseq,lexlscr)
                          first = .false.
                          if (.not.logdum) call quit(1,'invsqrt',
     &                         'no next tuple found!')
                          if (mod(idxcount(2,idspn,nel,1),4).ne.0)
     &                          flmap(icol,3) = -1
c dbg
c                          write(luout,'(i4,x,4i4,x,4i4)')idx,idorb,idspn
c dbgend
                        end if
                      end do
                    end do
                  end do
                end do

                off_col = off_col + len2(2)*len2(4)
                off_colmax = max(off_colmax,off_col)
               end do
              end do
              off_line = off_line + len2(1)*len2(3)
              off_linmax = max(off_linmax,off_line)
             end do
            end do

            deallocate(hpvx_csub2,hpvx_asub2,occ_csub2,
     &              occ_asub2,graph_csub2,graph_asub2,
     &              iocc2,idx_g2)
            na2 = na2 - 1
            nc2 = nc2 - 1
            jocc_cls = jocc_cls + 1
            if (na2+nc2.lt.abs(ms1)) !jump to next line
     &                jocc_cls = jocc_cls + min(na1mx,nc1mx)
     &                -(na1mx+nc1mx-abs(ms1))/2
            jblkoff = (jocc_cls-1)*njoined
            off_col2 = off_colmax
           end do
           na1 = na1 - 1
           nc1 = nc1 - 1
           off_line2 = off_linmax
          end do

c dbg
c          print *,'index matrix:'
c          write(*,'(5x,18i4)') (icol, icol=1,ndim)
c          do iline = 1, ndim
c            write(*,'(i4,x,18i4)') iline, matrix(iline,1:ndim)
c          end do
c dbgend

          if (sing_remove) then
            ! project out pseudo-singles components
c dbg
c           if (ipass.eq.0) then
c dbgend
            if (max(iprlvl,ntest).ge.3) write(luout,*)
     &           'Projecting out conventional singles.'
c dbg
c           else
c            if (max(iprlvl,ntest).ge.3) write(luout,*)
c     &           'Projecting out triples (keep pseudo-doubles).'
c           end if
c           ipass = ipass+1
c dbgend
ctest
            print *,'...with crazy new projector (CAS(2,2) only).'
ctestend
            ! remove vanishing excitations from proj. vectors
            do iline = 1, ndim
ctest
              pbeta(iline) = scratch(iline,iline) !=overlap with singles
c              ! ad hoc: only spectator excitations (4-det. case)
c              if (iline.eq.2.or.iline.eq.3
c     &            .or.iline.eq.6.or.iline.eq.7) then
c                palph(iline) = 0d0
c                pbeta(iline) = scratch(1,2) !0d0
c              end if
ctest              if (scratch(iline,iline).lt.1d-14) then
ctest                if (palph(iline).gt.1d-14) then
ctest                  palph(iline) = 0d0
ctest                  nalph = nalph - 1
ctest                end if
ctest                if (pbeta(iline).gt.1d-14) then
ctest                  pbeta(iline) = 0d0
ctest                  nbeta = nbeta - 1
ctest                end if
ctest              end if
            end do
ctest            ! normalize vectors
ctest            if (nalph.gt.0) palph = palph/sqrt(dble(nalph))
ctest            if (nbeta.gt.0) pbeta = pbeta/sqrt(dble(nbeta))
c dbg
            print *,'palph:',palph
            print *,'pbeta:',pbeta
c dbgend
            ! get projector
            proj_sing = 0d0
            do idx = 1, ndim
              proj_sing(idx,idx) = 1d0
            end do
ctest
c dbg
c           if (ipass.eq.1) then
c dbgend
            call dgemm('n','t',ndim,ndim,1,
     &                 -1d0,palph,ndim,
     &                 pbeta,ndim,
     &                 1d0,proj_sing,ndim)
c dbg
c           else
c            call dgemm('n','t',ndim,ndim,1,
c     &                 1d0,palph,ndim,
c     &                 pbeta,ndim,
c     &                 0d0,proj_sing,ndim)
c           end if
c dbgend
ctest            call dgemm('n','t',ndim,ndim,1,
ctest     &                 -1d0,palph,ndim,
ctest     &                 palph,ndim,
ctest     &                 1d0,proj_sing,ndim)
ctest            call dgemm('n','t',ndim,ndim,1,
ctest     &                 -1d0,pbeta,ndim,
ctest     &                 pbeta,ndim,
ctest     &                 1d0,proj_sing,ndim)
            if (ntest.ge.100) then
              write(luout,*) 'projector for removing singles:'
              call wrtmat2(proj_sing,ndim,ndim,ndim,ndim)
            end if
            if (ntest.ge.100) then
              write(luout,*) 'matrix before removing singles:'
              call wrtmat2(scratch,ndim,ndim,ndim,ndim)
            end if
ctest
            scratch2(1:ndim,1:ndim) = scratch(1:ndim,1:ndim)
ctest            call dgemm('n','n',ndim,ndim,ndim,
ctest     &             1d0,proj_sing,ndim,
ctest     &                 scratch,ndim,
ctest     &             0d0,scratch2,ndim)
            call dgemm('n','n',ndim,ndim,ndim,
     &             1d0,scratch2,ndim,
     &                 proj_sing,ndim,
     &             0d0,scratch,ndim)
            if (ntest.ge.100) then
              write(luout,*) 'matrix after removing singles:'
              call wrtmat2(scratch,ndim,ndim,ndim,ndim)
            end if
c            deallocate(palph,pbeta,proj_sing)
            deallocate(palph,pbeta)
          end if
c dbg
c          if (.not.sing_remove.and.iocc_cls.eq.8)
c     &      scratch(1:ndim,1:ndim) = 0d0
c dbgend

          ! normalization
          if (normalize) then
            allocate(norm(ndim))
            do idx = 1, ndim
              norm(idx) = sqrt(abs(scratch(idx,idx)))
            end do
            call dmdiagm(ndim,scratch,norm,.true.,.true.)
            call dmdiagm(ndim,scratch,norm,.false.,.true.)
          end if

          if (ms1.eq.0) then
            ! here a splitting into "singlet" and "triplet" blocks is needed:

c dbg
c            write(luout,*) 'flmap:'
c            do icol = 1, ndim
c              write(luout,'(i4,2i6)') icol,flmap(icol,1:2)
c            end do
c dbgend
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
            if (.not.half) allocate(sing2(nsing,nsing),
     &                              trip2(ntrip,ntrip))
            call spinsym_traf(1,ndim,scratch,flmap(1:ndim,3),nsing,
     &                        sing,trip,.false.)

            ! calculate T^(-0.5) for both blocks
            call invsqrt_mat(nsing,sing,sing2,half,icnt_sv,icnt_sv0,
     &                       xmax,xmin,bins)
            call invsqrt_mat(ntrip,trip,trip2,half,icnt_sv,icnt_sv0,
     &                       xmax,xmin,bins)

            ! partial undo of pre-diagonalization: Upre*T^(-0.5)
            call spinsym_traf(2,ndim,scratch,flmap(1:ndim,3),nsing,
     &                        sing,trip,.true.)
            if (.not.half) then
              ! full undo of pre-diagonalization for projector
              call spinsym_traf(2,ndim,scratch2,flmap(1:ndim,3),nsing,
     &                          sing2,trip2,.false.)
              deallocate(sing2,trip2)
            end if
            deallocate(sing,trip)
          else

            ! calculate S^(-0.5)
            call invsqrt_mat(ndim,scratch,scratch2,
     &                       half,icnt_sv,icnt_sv0,xmax,xmin,bins)

          end if

          ! "undo" normalization
          if (normalize) then ! X = N*X'
            call dmdiagm(ndim,scratch,norm,.false.,.true.)
            if (.not.half) then ! P = N*P'*N^-1
              call dmdiagm(ndim,scratch2,norm,.false.,.true.)
              call dmdiagm(ndim,scratch2,norm,.true.,.false.)
            end if
            deallocate(norm)
          end if

c not needed for symmetric projector
          ! multiply trafo (and proj.) matrix with proj. for singles
          if (sing_remove) then
            if (ntest.ge.10)
     &       write(luout,*) 'Multiplying X (&P) with singles projector'
            allocate(scratch3(ndim,ndim))
            call dgemm('n','n',ndim,ndim,ndim,
     &                 1d0,proj_sing,ndim,
     &                 scratch,ndim,
     &                 0d0,scratch3,ndim)
            scratch(1:ndim,1:ndim) = scratch3(1:ndim,1:ndim)
            if (ntest.ge.100) then
              write(luout,*) 'Trafo matrix:'
              call wrtmat2(scratch,ndim,ndim,ndim,ndim)
            end if
            if (.not.half) then
              call dgemm('n','n',ndim,ndim,ndim,
     &                   1d0,proj_sing,ndim,
     &                   scratch2,ndim,
     &                   0d0,scratch3,ndim)
              scratch2(1:ndim,1:ndim) = scratch3(1:ndim,1:ndim)
              if (ntest.ge.100) then
                write(luout,*) 'Projector matrix:'
                call wrtmat2(scratch2,ndim,ndim,ndim,ndim)
              end if
            end if
            deallocate(scratch3,proj_sing)
          end if

          ! write to output buffer
          ! loops over coupling blocks. Must be in correct order!
          jocc_cls = iocc_cls
          jblkoff = (jocc_cls-1)*njoined
          na1 = na1mx
          nc1 = nc1mx
          off_line2 = 0
          do while(min(na1,nc1).ge.0.and.na1+nc1.ge.abs(ms1))
           na2 = nc1mx
           nc2 = na1mx
           off_col2 = 0
           off_linmax = 0
           do while(min(na2,nc2).ge.0.and.na2+nc2.ge.abs(ms1))
            ! exit if there is no block like this
            if (na1.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+1)).or.
     &          nc1.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+2)).or.
     &          na2.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+2)).or.
     &          nc2.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+3))) exit
            ! exception for pure inactive block:
            if (na1+nc1+na2+nc2.eq.0.and.igam.eq.1) then
              ioff = mel_inp%off_op_gmox(jocc_cls)%
     &                 d_gam_ms(1,1,1)
              buffer_out(ioff+1) = scratch(off_line2+1,off_col2+1)
              if (.not.half) ! copy projector to input buffer
     &           buffer_in(ioff+1) = scratch2(off_line2+1,off_col2+1)
              exit
            end if
            msmax = op_inp%ica_occ(1,jocc_cls)
            call get_num_subblk(ncblk2,nablk2,
     &           hpvx_occ(1,1,jblkoff+1),njoined)
            allocate(hpvx_csub2(ncblk2),hpvx_asub2(nablk2),
     &               occ_csub2(ncblk2), occ_asub2(nablk2),
     &               graph_csub2(ncblk2), graph_asub2(nablk2))
            ! set HPVX and OCC info
            call condense_occ(occ_csub2, occ_asub2,
     &                      hpvx_csub2,hpvx_asub2,
     &                      hpvx_occ(1,1,jblkoff+1),njoined,hpvxblkseq)
            ! do the same for the graph info
            call condense_occ(graph_csub2, graph_asub2,
     &                      hpvx_csub2,hpvx_asub2,
     &                      idx_graph(1,1,jblkoff+1),njoined,hpvxblkseq)

            ! loop over Ms(A1)/GAMMA(A1) --> Ms(C1)/GAMMA(C1) already defined
            off_line = off_line2
            off_colmax = 0
            do msa1 = na1, -na1, -2
             msc1 = msa1 + ms1
             if (abs(msc1).gt.nc1) cycle
             msdis_c(1) = msc1
             msdis_a(1) = msa1
             do gama1 = 1, ngam
              gamc1 = multd2h(gama1,igam)
              if (na1.eq.0.and.gama1.ne.1.or.
     &            nc1.eq.0.and.gamc1.ne.1) cycle
              gamdis_c(1) = gamc1
              gamdis_a(1) = gama1
              ! loop over Ms(C2)/GAMMA(C2) --> Ms(A2)/GAMMA(A2) already defined
              off_col = off_col2
              do msc2 = nc2, -nc2, -2
               msa2 = msc2 - ms2
               if (abs(msa2).gt.na2) cycle
               msa = msa1 + msa2
               idxmsa = (msmax-msa)/2 + 1
               do gamc2 = 1, ngam
                gama2 = multd2h(gamc2,igam)
                if (na2.eq.0.and.gama2.ne.1.or.
     &              nc2.eq.0.and.gamc2.ne.1) cycle
                igama = multd2h(gama1,gama2)
                if (nablk2.eq.1.and.na1.eq.0) then
                  msdis_a(1) = msa2
                  gamdis_a(1) = gama2
                else if (nablk2.eq.2) then
                  msdis_a(2) = msa2
                  gamdis_a(2) = gama2
                end if
                if (ncblk2.eq.1.and.nc1.eq.0) then
                  msdis_c(1) = msc2
                  gamdis_c(1) = gamc2
                else if (ncblk2.eq.2) then
                  msdis_c(2) = msc2
                  gamdis_c(2) = gamc2
                end if

                ! determine the distribution in question
                call ms2idxms(idxmsdis_c,msdis_c,occ_csub2,ncblk2)
                call ms2idxms(idxmsdis_a,msdis_a,occ_asub2,nablk2)

                call set_len_str(len_str,ncblk2,nablk2,
     &                       graphs,
     &                       graph_csub2,idxmsdis_c,gamdis_c,hpvx_csub2,
     &                       graph_asub2,idxmsdis_a,gamdis_a,hpvx_asub2,
     &                       hpvxseq,.false.)
                len2(1:4) = 1
                if (ncblk2.eq.1.and.nc2.eq.0) len2(1) = len_str(1)
                if (ncblk2.eq.1.and.nc1.eq.0) len2(2) = len_str(1)
                if (ncblk2.eq.2) then
                  len2(1) = len_str(1)
                  len2(2) = len_str(2)
                end if
                if (nablk2.eq.1.and.na2.eq.0) len2(3) =len_str(ncblk2+1)
                if (nablk2.eq.1.and.na1.eq.0) len2(4) =len_str(ncblk2+1)
                if (nablk2.eq.2) then
                  len2(3) = len_str(ncblk2+1)
                  len2(4) = len_str(ncblk2+2)
                end if
                if (ielprd(len_str,ncblk2+nablk2).eq.0) cycle

                ndis = mel_inp%off_op_gmox(jocc_cls)%ndis(igama,idxmsa)
                idxdis = 1
                if (ndis.gt.1)
     &             idxdis =
     &                 idx_msgmdst2(
     &                  jocc_cls,idxmsa,igama,
     &                  occ_csub2,idxmsdis_c,gamdis_c,ncblk2,
     &                  occ_asub2,idxmsdis_a,gamdis_a,nablk2,
     &                  .false.,-1,-1,mel_inp,ngam)

                call set_op_ldim_c(ldim_opin_c,ldim_opin_a,
     &               hpvx_csub2,hpvx_asub2,
     &               len_str,ncblk2,nablk2,.false.)

                ioff = mel_inp%off_op_gmox(jocc_cls)%
     &                 d_gam_ms(idxdis,igama,idxmsa)

                ! copy all required elements of this distribution
                ! to their block in the A1/C1 | A2/C2 matrix
                iline = off_line
                do idxc1 = 1, len2(1)
                  if (nc1.ne.0) istr_csub(1) = idxc1-1
                  do idxa1 = 1, len2(3)
                    if (na1.ne.0) istr_asub(1) = idxa1-1
                    iline = iline + 1
                    icol = off_col
                    do idxa2 = 1, len2(4)
                      if (na1.eq.0.and.na2.ne.0) then
                       istr_asub(1) = idxa2-1
                      else if (na2.ne.0) then
                       istr_asub(2) = idxa2-1
                      end if
                      do idxc2 = 1, len2(2)
                        if (nc1.eq.0.and.nc2.ne.0) then
                         istr_csub(1) = idxc2-1
                        else if (nc2.ne.0) then
                         istr_csub(2) = idxc2-1
                        end if
                        icol = icol + 1
                        idx = ioff + idx_str_blk3(istr_csub,istr_asub,
     &                         ldim_opin_c,ldim_opin_a,ncblk2,nablk2)
                        buffer_out(idx) = scratch(iline,icol)
                        if (.not.half) ! copy projector to input buffer
     &                     buffer_in(idx) = scratch2(iline,icol)
                      end do
                    end do
                  end do
                end do

                off_col = off_col + len2(2)*len2(4)
                off_colmax = max(off_colmax,off_col)
               end do
              end do
              off_line = off_line + len2(1)*len2(3)
              off_linmax = max(off_linmax,off_line)
             end do
            end do

            deallocate(hpvx_csub2,hpvx_asub2,occ_csub2,
     &              occ_asub2,graph_csub2,graph_asub2)
            na2 = na2 - 1
            nc2 = nc2 - 1
            jocc_cls = jocc_cls + 1
            if (na2+nc2.lt.abs(ms1)) !jump to next line
     &                jocc_cls = jocc_cls + min(na1mx,nc1mx)
     &                -(na1mx+nc1mx-abs(ms1))/2
            jblkoff = (jocc_cls-1)*njoined
            off_col2 = off_colmax
           end do
           na1 = na1 - 1
           nc1 = nc1 - 1
           off_line2 = off_linmax
          end do

          deallocate(scratch,flmap)
          if (.not.half.or.sing_remove) deallocate(scratch2)
c dbg
c          deallocate(matrix)
c dbgend

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
      deallocate(blk_used)

      if (ntest.ge.5) then
        write(luout,'(x,i8,a,i8,a)') icnt_sv-icnt_sv0,
     &        ' out of ',icnt_sv,' singular values were included'
        write(luout,'(x,a,E19.10)') 
     &        'The  largest excluded singular value is ',xmax
        write(luout,'(x,a,E19.10)')
     &        'The smallest included singular value is ',xmin
        write(luout,'(x,a)') 'Singular value histogram'
        write(luout,'(x,a)') '------------------------'
        write(luout,'(x,a,x,i16)') ' 1E+00',bins(1)
        do idx = 2, 16
          write(luout,'(x,a,i2.2,x,i16)') ' 1E-',idx-1,bins(idx)
        end do
        write(luout,'(x,a,x,i16)') '     0',bins(17)
        write(luout,'(x,a)') '------------------------'
      end if
 

      if(.not.bufout)then
        call put_vec(ffinv,buffer_out,1,nbuff)
      endif  
      if(.not.bufin.and..not.half)then
        ! return projector matrix on input list
        call put_vec(ffinp,buffer_in,1,nbuff)
      endif

      ifree = mem_flushmark('invsqrt')

      return
      end
