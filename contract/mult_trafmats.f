      subroutine mult_trafmats(mel_inp,mel_inv,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*     multiplies two (transformation) matrices of the types:
*      
*       /0 0 0 0\
*       \0 0 x 0/  (not necessarily igastp=3)
*       /0 0 y 0\ 
*       \0 0 y 0/ 
*       /0 0 x 0\
*       \0 0 0 0/
*
*     While the first matrix "mel_inp" may have off-diagonal blocks,
*     the second ("mel_inv") is assumed to have diagonal blocks only.
*     The result is written to mel_inp.
*
*     matthias, nov 2013
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'hpvxseq.h'
      include 'multd2h.h'
      include 'ifc_input.h'
      include 'routes.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(me_list), intent(in) ::
     &     mel_inp, mel_inv

      integer, parameter ::
     &     maxrank = 5
      type(graph), pointer ::
     &     graphs(:)

      logical ::
     &     bufin, buf2, onedis, transp, sgrm
      logical, pointer ::
     &     blk_used(:)
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
     &     ngraph, ioff2, idx2, nel,
     &     nc1mx, na1mx, jocc_cls, jblkoff,
     &     ncblk2, nablk2, len2(4), off_line2, off_col2, off_colmax,
     &     off_linmax,
     &     rankdim(maxrank), rankoff(maxrank), nrank,
     &     rdim, gno, project,
     &     nocc_cls2, blk_corresp(maxrank), iocc_cls2, iblkoff2,
     &     nocc_cls, nbuff_in
      real(8), pointer ::
     &     buffer_in(:), buffer_2(:), scratch(:,:), scratch2(:,:),
     &     scratch3(:,:)

      integer, pointer ::
     &     msdis_c(:),  msdis_a(:),
     &     idxmsdis_c(:),  idxmsdis_a(:),
     &     gamdis_c(:), gamdis_a(:),
     &     len_str(:),
     &     hpvx_occ(:,:,:), idx_graph(:,:,:),
     &     ldim_opin_c(:), ldim_opin_a(:),
     &     istr_csub(:), istr_asub(:),
     &     hpvx_csub2(:),hpvx_asub2(:),
     &     occ_csub2(:), occ_asub2(:),
     &     graph_csub2(:), graph_asub2(:),
     &     iocc_inp(:,:,:), iocc_2(:,:,:)
c dbg
c      integer, pointer ::
c     &     matrix(:,:)
c dbgend

      type(filinf), pointer ::
     &     ffinp, ffinv
      type(operator), pointer ::
     &     op_inv, op_inp

      integer, external ::
     &     ielprd, idx_msgmdst2, idx_str_blk3
      logical, external ::
     &     iocc_equal_n

      if (ntest.ge.100) write(lulog,*) 'entered mult_trafmats'
      call get_argument_value('method.MR','project',ival=project)
      call get_argument_value('method.MR','GNO',ival=gno)
      sgrm = project.eq.1.and.gno.eq.0 ! only diagonal blocks in mel_inp

      ffinp => mel_inp%fhand
      ffinv => mel_inv%fhand
      op_inp => mel_inp%op
      op_inv => mel_inv%op
      njoined = op_inp%njoined
      nocc_cls = op_inp%n_occ_cls
      nocc_cls2 = op_inv%n_occ_cls

      ngraph = str_info%ngraph
      hpvx_occ => op_inp%ihpvca_occ
      idx_graph => mel_inp%idx_graph
      graphs => str_info%g

      ! Check whether files are buffered.
      bufin = .false.
      if(ffinp%buffered) bufin = .true.
      buf2 = .false.
      if(ffinv%buffered) buf2 = .true.

      ifree = mem_setmark('mult_trafmats')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Allocations made to maximum block length to save time.
      if(.not.bufin)then
        nbuff_in = 0
        do iocc_cls = 1, nocc_cls
          if(op_inp%formal_blk(iocc_cls))
     &         cycle
          nbuff_in = nbuff_in+mel_inp%len_op_occ(iocc_cls)
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff_in,'buffer_in')
        call get_vec(ffinp,buffer_in,1,nbuff_in)
      else
        if(ntest.ge.100)
     &       write(lulog,*)'mult_trafmats: input not incore'
        buffer_in => ffinp%buffer(1:)
      endif

      if(.not.buf2)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls2
          if(op_inv%formal_blk(iocc_cls))
     &         cycle
          nbuff = nbuff + mel_inv%len_op_occ(iocc_cls)
        enddo
        ifree= mem_alloc_real(buffer_2,nbuff,'buffer_2')
        call get_vec(ffinv,buffer_2,1,nbuff)
      else
        if(ntest.ge.100)
     &       write(lulog,*)'mult_trafmats: output not incore'
        buffer_2 => ffinv%buffer(1:)
      endif

      allocate(blk_used(nocc_cls))
      blk_used(1:nocc_cls) = .false.

      ! Loop over occupation class.
      iocc_loop: do iocc_cls = 1, nocc_cls
        iblkoff = (iocc_cls-1)*njoined
        if(op_inp%formal_blk(iocc_cls)) cycle iocc_loop
        if (blk_used(iocc_cls)) cycle iocc_loop
        blk_used(iocc_cls) = .true.
        iocc_inp => op_inp%ihpvca_occ(1:ngastp,1:2,
     &                       iblkoff+1:iblkoff+njoined)

        if (ntest.ge.10) write(lulog,*) 'current occ_cls: ',iocc_cls
        ! only one element? easy!
        ! (also regularization is never needed in this case)
        if (mel_inp%len_op_occ(iocc_cls).eq.1) then
          ! exception for purely inactive block: no change
          if (op_inp%ica_occ(1,iocc_cls).eq.0.and.
     &        op_inp%ica_occ(2,iocc_cls).eq.0) cycle iocc_loop
          ioff = mel_inp%off_op_gmo(iocc_cls)%gam_ms(1,1)
          ! find corresponding block
          do iocc_cls2 = 1, nocc_cls2
            iblkoff2 = (iocc_cls2-1)*njoined
            iocc_2 => op_inv%ihpvca_occ(1:ngastp,1:2,
     &                       iblkoff2+1:iblkoff2+njoined)
            if (iocc_equal_n(iocc_inp,.true.,iocc_2,.true.,njoined))
     &         exit
            if (iocc_cls2.eq.nocc_cls2)
     &         call quit(1,'mult_trafmats','1) No corresponding block!')
          end do
          if (ntest.ge.10) write(lulog,*)
     &       'corresponding block:',iocc_cls2
          ioff2 = mel_inv%off_op_gmo(iocc_cls2)%gam_ms(1,1)
          buffer_in(ioff+1) = buffer_in(ioff+1) * buffer_2(ioff2+1)
          cycle iocc_loop
        end if 

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
          call quit(1,'mult_trafmats','not adapted for this case')
        end if

        allocate(msdis_c(ncblk),msdis_a(nablk),
     &           idxmsdis_c(ncblk),idxmsdis_a(nablk),
     &           gamdis_c(ncblk), gamdis_a(nablk),
     &           len_str(ncblk+nablk),
     &           istr_csub(ncblk),istr_asub(nablk),
     &           ldim_opin_c(ncblk),ldim_opin_a(nablk))

        ! simple case: only single distributions:
        if (onedis) then

          ! we also need to distinguish:
          ! /0 0 0 0\     /0 0 0 0\
          ! \0 0 x 0/     \0 0 0 0/ i.e. if creators come first,
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
            write(lulog,*) 'Using transposed scratch matrix!'
          end if

          ! find corresponding block
          do iocc_cls2 = 1, nocc_cls2
            iblkoff2 = (iocc_cls2-1)*njoined
            iocc_2 => op_inv%ihpvca_occ(1:ngastp,1:2,
     &                       iblkoff2+1:iblkoff2+njoined)
            if (iocc_equal_n(iocc_inp,.true.,iocc_2,.true.,njoined))
     &         exit
            if (iocc_cls2.eq.nocc_cls2)
     &         call quit(1,'mult_trafmats','2) No corresponding block!')
          end do
          if (ntest.ge.10) write(lulog,*)
     &       'corresponding block:',iocc_cls2

          ! Loop over Ms of annihilator string.
          idxmsa = 0
          msmax = op_inp%ica_occ(1,iocc_cls)
          mscmax = op_inp%ica_occ(2,iocc_cls)
          if (msmax.ne.mscmax) call quit(1,'mult_trafmats',
     &            'need particle conserving operator')
          msa_loop : do msa = msmax, -msmax, -2

            idxmsa = idxmsa+1
            msc = msa + mel_inp%mst
            idxmsc = (msmax-msc)/2 + 1
      
            ! Loop over Irrep of annihilator string.
            igama_loop: do igama =1, ngam

              igamc = multd2h(igama,mel_inp%gamt)
              ndis = mel_inp%off_op_gmox(iocc_cls)%ndis(igamc,idxmsc)

              ndim = int(sqrt(dble(mel_inp%
     &           len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa))))
              if (ndim.gt.0.and.ndis.ne.1)
     &                call quit(1,'mult_trafmats','cannot handle this')
              if (ndim.eq.0) cycle igama_loop

              if (ntest.ge.10)
     &           write(lulog,'(a,3i8)') 'msa, gama, ndim:',msa,igama,
     &           ndim
              if (ntest.ge.100)
     &           write(lulog,*) ' len = ',
     &             mel_inp%len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa),
     &             ' ndis = ',ndis

              ioff = mel_inp%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)
              ioff2 = mel_inv%off_op_gmo(iocc_cls2)%gam_ms(igama,idxmsa)

              ! single distribution can simply be read in as simple matrix
              allocate(scratch3(ndim,ndim),scratch2(ndim,ndim))
              if (transp) then
                do idx = 1,ndim
                  do jdx = 1,ndim
                    scratch3(jdx,idx) = buffer_in(ioff+(idx-1)*ndim+jdx)
                    scratch2(jdx,idx) = buffer_2(ioff2+(idx-1)*ndim+jdx)
                  enddo
                enddo
              else
                do idx = 1,ndim
                  do jdx = 1,ndim
                    scratch3(idx,jdx) = buffer_in(ioff+(idx-1)*ndim+jdx)
                    scratch2(idx,jdx) = buffer_2(ioff2+(idx-1)*ndim+jdx)
                  enddo
                enddo
              end if

              if (ntest.ge.100) then
                write(lulog,*) 'trafo matrix 1:'
                call wrtmat2(scratch3,ndim,ndim,ndim,ndim)
                write(lulog,*) 'trafo matrix 2:'
                call wrtmat2(scratch2,ndim,ndim,ndim,ndim)
              end if

              allocate(scratch(ndim,ndim))
              ! multiplication
              call dgemm('n','n',ndim,ndim,ndim,
     &                   1d0,scratch3,ndim,
     &                   scratch2,ndim,
     &                   0d0,scratch,ndim)
              deallocate(scratch2,scratch3)

              if (ntest.ge.100) then
                write(lulog,*) 'product of trafo matrices:'
                call wrtmat2(scratch,ndim,ndim,ndim,ndim)
              end if

              ! write to output buffer
              if (transp) then
                do idx = 1,ndim
                  do jdx = 1,ndim
                    buffer_in((idx-1)*ndim+jdx+ioff) = scratch(jdx,idx)
                  enddo
                enddo
              else
                do idx = 1,ndim
                  do jdx = 1,ndim
                    buffer_in((idx-1)*ndim+jdx+ioff) = scratch(idx,jdx)
                  enddo
                enddo
              end if

              deallocate(scratch)

            enddo igama_loop
          enddo msa_loop

          deallocate(msdis_c,msdis_a,
     &               idxmsdis_c,idxmsdis_a,
     &               gamdis_c,gamdis_a,
     &               len_str,
     &               istr_csub,istr_asub,
     &               ldim_opin_c,ldim_opin_a)

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
     &     call quit(1,'mult_trafmats','too difficult operator!')
        msmax = op_inp%ica_occ(1,iocc_cls)
        mscmax = op_inp%ica_occ(2,iocc_cls)
        if (msmax.ne.mscmax.or.msmax.ne.na1mx+na2
     &      .or.mscmax.ne.nc1mx+nc2)
     &            call quit(1,'mult_trafmats','need part.cons.op')

        ! loop over Ms/Gamma combinations of A1/C1 tuple
        ! which are the decoupled blocks of the A1/C1 | A2/C2 matrix
        ! Ms/Gamma of C2/A2 is then already defined
        ! For now we consider densities: Ms(total)=0, Gamma(total) = 1
        nel = na2 + nc2
        do ms1 = msmax_sub, -msmax_sub, -2
         ms2 = -ms1 ! particle conserving operator
         do igam = 1, ngam ! irrep of the second dimension must be the same

          ! First we need the dimension of the A1/C1 | A2/C2 block
          ndim = 0
          nrank = 0
          ! loops over coupling blocks. Must be in correct order!
          jocc_cls = iocc_cls
          jblkoff = (jocc_cls-1)*njoined
          na1 = na1mx
          nc1 = nc1mx
          blk_loop: do while(min(na1,nc1).ge.0.and.na1+nc1.ge.abs(ms1))
           na2 = nc1mx
           nc2 = na1mx
           rdim = 0
           do while(min(na2,nc2).ge.0.and.na2+nc2.ge.abs(ms1))
            ! no off-diagonal blocks in SepO (separate orthog.)
            if (sgrm.and.na1.ne.nc2) then
              na2 = na2 - 1
              nc2 = nc2 - 1
              cycle
            end if
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
            if (na1.eq.0.and.nc1.eq.0) then
              if (igam.eq.1) ndim = ndim + 1
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
                rdim = rdim + len_str(1)
              else
                rdim = rdim + len_str(1)*len_str(3)
              end if
             end do
            end do
            ndim = ndim + rdim

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
           if (project.gt.0.and.rdim.gt.0) then
             nrank = nrank + 1
             if (nrank.gt.maxrank) call quit(1,'mult_trafmats',
     &                                       'increase maxrank')
             do idx = nrank, 2, -1
               rankdim(idx) = rankdim(idx-1)
               rankoff(idx) = rankoff(idx-1)
             end do
             rankdim(1) = rdim
             rankoff(1) = 0
             if (nrank.ge.2) rankoff(1) = rankoff(2) + rankdim(2)
           end if
          end do blk_loop

          if (ndim.eq.0) cycle
          if (project.eq.0) then
            nrank = 1
            rankdim(1) = ndim
            rankoff(1) = 0
          else if (ndim.eq.rankoff(1)+rankdim(1)+1) then !inactive blk.
             nrank = nrank + 1
             if (nrank.gt.maxrank) call quit(1,'mult_trafmats',
     &                                       'increase maxrank')
             do idx = nrank, 2, -1
               rankdim(idx) = rankdim(idx-1)
               rankoff(idx) = rankoff(idx-1)
             end do
             rankdim(1) = 1
             rankoff(1) = 0
             if (nrank.ge.2) rankoff(1) = rankoff(2) + rankdim(2)
          else if (ndim.ne.rankoff(1)+rankdim(1)) then
            call quit(1,'mult_trafmats','dimensions don''t add up!')
          end if

          if (ntest.ge.10)
     &       write(lulog,'(a,3i8)') 'ms1, igam, ndim:',ms1,igam,ndim
          if (ntest.ge.100)
     &       write(lulog,'(a,5i8)') 'dim. per rank:',rankdim(1:nrank)

          allocate(scratch3(ndim,ndim),scratch2(ndim,ndim))
          scratch3 = 0d0
          scratch2 = 0d0

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
            ! no off-diagonal blocks in SepO (separate orthog.)
            if (sgrm.and.na1.ne.nc2) then
              na2 = na2 - 1
              nc2 = nc2 - 1
              off_col2 = off_colmax
              cycle
            end if
            ! exit if there is no block like this
            if (na1.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+1)).or.
     &          nc1.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+2)).or.
     &          na2.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+2)).or.
     &          nc2.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+3))) exit
            blk_used(jocc_cls) = .true.
            ! exception for pure inactive block: assume 2nd matrix is 1
            if (na1+nc1+na2+nc2.eq.0) then
              if (igam.ne.1) exit
              ioff = mel_inp%off_op_gmox(jocc_cls)%
     &                 d_gam_ms(1,1,1)
              scratch3(off_line2+1,off_col2+1) = buffer_in(ioff+1)
              scratch2(off_line2+1,off_col2+1) = 1d0
              exit
            end if
c dbg
c            print *,'na1,nc1,na2,nc2,jocc_cls:',na1,nc1,na2,nc2,
c     &          jocc_cls
c            print *,'off_line2, off_col2: ',off_line2,off_col2
c dbgend

            ! find corresponding block (if diagonal)
            if (na1.eq.nc2) then
              iocc_inp => op_inp%ihpvca_occ(1:ngastp,1:2,
     &                         jblkoff+1:jblkoff+njoined)
              do iocc_cls2 = 1, nocc_cls2
                iblkoff2 = (iocc_cls2-1)*njoined
                iocc_2 => op_inv%ihpvca_occ(1:ngastp,1:2,
     &                           iblkoff2+1:iblkoff2+njoined)
                if (iocc_equal_n(iocc_inp,.true.,iocc_2,.true.,njoined))
     &             exit
                if (iocc_cls2.eq.nocc_cls2)
     &             call quit(1,'mult_trafmats',
     &                       '3) No corresponding block!')
              end do
              if (ntest.ge.10) write(lulog,*)
     &           'corresponding block:',iocc_cls2
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
     &                 idx_msgmdst2(.true.,
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
     &             call quit(1,'mult_trafmats','inconsistency!')

                call set_op_ldim_c(ldim_opin_c,ldim_opin_a,
     &               hpvx_csub2,hpvx_asub2,
     &               len_str,ncblk2,nablk2,.false.)

                ioff = mel_inp%off_op_gmox(jocc_cls)%
     &                 d_gam_ms(idxdis,igama,idxmsa)

                if (na1.eq.nc2)
     &             ioff2 = mel_inv%off_op_gmox(iocc_cls2)%
     &                 d_gam_ms(idxdis,igama,idxmsa)

c dbg
c                write(lulog,'(a,4i4)') 'ms  : ',msa1,msc1,msa2,msc2
c                write(lulog,'(a,4i4)') 'gam : ',gama1,gamc1,gama2,gamc2
c                write(lulog,'(a,2i4)') 'msa, igama: ',msa,igama
c                write(lulog,'(a,2i4)') 'dist, len: ',idxdis, lenca
c                write(lulog,'(a,2i4)') 'off_line/col: ',off_line,off_col
c                print *,'len1: ',len_str(1)*len_str(3)
c dbgend

                ! copy all required elements of this distribution
                ! to their block in the A1/C1 | A2/C2 matrix
                iline = off_line
                do idxc1 = 1, len2(1)
                  if (nc1.ne.0) then
                   istr_csub(1) = idxc1-1
                  end if
                  do idxa1 = 1, len2(3)
                    if (na1.ne.0) then
                     istr_asub(1) = idxa1-1
                    end if
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
                        scratch3(iline,icol) = buffer_in(idx)

                        if (na1.eq.nc2) then
                          idx2 = ioff2 + idx - ioff
                          scratch2(iline,icol) = buffer_2(idx2)
                        end if
c dbg
c                        matrix(iline,icol) = idx
c dbgend
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
c dbg
c          print *,'index matrix:'
c          write(*,'(5x,18i4)') (icol, icol=1,ndim)
c          do iline = 1, ndim
c            write(*,'(i4,x,18i4)') iline, matrix(iline,1:ndim)
c          end do
c dbgend

          ! multiplication:
          ! / X2  X21\ / Z2  0  \   / X2*Z2  X21*Z1\
          ! |        |*|        | = |              |
          ! \ X12 X1 / \ 0   Z1 /   \ X12*Z2 X1*Z1 /
          if (ntest.ge.100) then
            write(lulog,*) 'trafo matrix 1:'
            call wrtmat2(scratch3,ndim,ndim,ndim,ndim)
            write(lulog,*) 'trafo matrix 2:'
            call wrtmat2(scratch2,ndim,ndim,ndim,ndim)
          end if

          allocate(scratch(ndim,ndim))
          call dgemm('n','n',ndim,ndim,ndim,
     &               1d0,scratch3,ndim,
     &               scratch2,ndim,
     &               0d0,scratch,ndim)
          deallocate(scratch2,scratch3)

          if (ntest.ge.100) then
            write(lulog,*) 'product of trafo matrices:'
            call wrtmat2(scratch,ndim,ndim,ndim,ndim)
          end if

          ! write to input buffer
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
            ! no off-diagonal blocks in SepO (separate orthog.)
            if (sgrm.and.na1.ne.nc2) then
              na2 = na2 - 1
              nc2 = nc2 - 1
              off_col2 = off_colmax
              cycle
            end if
            ! exit if there is no block like this
            if (na1.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+1)).or.
     &          nc1.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+2)).or.
     &          na2.ne.sum(hpvx_occ(1:ngastp,2,jblkoff+2)).or.
     &          nc2.ne.sum(hpvx_occ(1:ngastp,1,jblkoff+3))) exit
            ! exception for pure inactive block:
            if (na1+nc1+na2+nc2.eq.0) then
              if (igam.ne.1) exit
              ioff = mel_inp%off_op_gmox(jocc_cls)%
     &                 d_gam_ms(1,1,1)
              buffer_in(ioff+1) = scratch(off_line2+1,off_col2+1)
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
     &                 idx_msgmdst2(.true.,
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
                        buffer_in(idx) = scratch(iline,icol)
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

          deallocate(scratch)
c dbg
c          deallocate(matrix)
c dbgend

         end do
        end do

        deallocate(msdis_c,  msdis_a,
     &           idxmsdis_c,  idxmsdis_a,
     &           gamdis_c, gamdis_a,
     &           len_str,
     &           istr_csub,istr_asub,
     &           ldim_opin_c,ldim_opin_a)

      enddo iocc_loop
      deallocate(blk_used)

      if(.not.bufin)then
        call put_vec(ffinp,buffer_in,1,nbuff_in)
      endif

      ifree = mem_flushmark('mult_trafmats')

      return
      end
