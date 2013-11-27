      subroutine invsqrt(mel_inp,mel_inv,nocc_cls,half,
     &     get_u,mel_u,pass_spc,mel_spc,lmodspc,
     &     op_info,orb_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     Routine to calculate S^(-0.5) of density matrices.
*     The ME list is split into matrices that can either contain
*     elements from a single distribution or from several distributions.
*     half = true: only makes a half transform U*s^(-0.5)
*     half = false: also does half transform and in addition
*                   returns projector matrix U*1s*U^+ on input list
*     get_u= true: return unitary matrix U on mel_u
*     lmodspc=true: special mode: no sequential projections; symmetrize;
*                   just diagonal blocks; no spin splitting for ms=0;
*                   return diagonalized list on inp and unitary on inv
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
      include 'ifc_input.h'
      include 'routes.h'

      integer, parameter ::
     &     ntest = 5

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(me_list), intent(in) ::
     &     mel_inp, mel_inv, mel_u, mel_spc
      integer, intent(in) ::
     &     nocc_cls
      logical, intent(in) ::
     &     half, get_u, pass_spc, lmodspc

      integer, parameter ::
     &     maxrank = 5
      type(graph), pointer ::
     &     graphs(:)

      logical ::
     &     bufin, bufout, first, ms_fix, fix_success, onedis, transp,
     &     logdum, sgrm, bufu, svdonly, reg_tik, first2
      logical, pointer ::
     &     blk_used(:), blk_redundant(:)
c      logical ::
c     &     loop(nocc_cls)
      integer ::
     &     ifree, nbuff, idxmsa, iocc_cls, iexc_cls,
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
     &     nalph, nbeta, nc1mx, na1mx, jocc_cls, jblkoff,
     &     ncblk2, nablk2, len2(4), off_line2, off_col2, off_colmax,
c dbg
c     &     ipass,
c dbgend
     &     off_linmax, maxbuf_tmp,
     &     rankdim(maxrank), rankoff(maxrank), nrank,
     &     irank, jrank, idxst, idxnd, rdim, idxst2, idxnd2, rdim2,
     &     icnt_cur, ih, ip, iexc, tocc_cls, gno, project,
     &     krank, idxst3, idxnd3, rdim3
      real(8) ::
     &     fac, xmax, xmin, xdum, omega2
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:), scratch(:,:), scratch2(:,:),
     &     sing(:,:), trip(:,:), sing2(:,:), trip2(:,:),
     &     palph(:), pbeta(:), proj(:,:), norm(:), scratch3(:,:),
     &     scratch4(:,:), sing3(:,:), trip3(:,:), buffer_u(:),
     &     svs(:)
      real(8), target ::
     &     xdummy(1,1)
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
     &     iocc2(:,:,:),idx_g2(:,:,:),
     &     bins(:,:), ex2occ_cls(:)

      type(filinf), pointer ::
     &     ffinp, ffinv, ffu
      type(operator), pointer ::
     &     op_inv, op_inp, op_u, op_t

      integer, external ::
     &     ielprd, idx_msgmdst2, idx_str_blk3, msa2idxms4op, idxcount,
     &     idxlist, idx_oplist2
      logical, external ::
     &     next_tupel_ca, occ_is_diag_blk

      if (ntest.ge.100) write(lulog,*) 'entered invsqrt'
c dbg
c      ipass = 0
c dbgend
      call get_argument_value('method.MR','project',ival=project)
      call get_argument_value('method.MR','GNO',ival=gno)
      sgrm = project.eq.1.and.gno.eq.0.or.lmodspc ! only diagonal blocks
      if (lmodspc.and.(get_u.or..not.half))
     &   call quit(1,'invsqrt','lmodspc only with get_u=F, half=T')

      reg_tik = tikhonov.ne.0d0
      omega2 = tikhonov**2

      ffinp => mel_inp%fhand
      ffinv => mel_inv%fhand
      ffu => mel_u%fhand
      op_inp => mel_inp%op
      op_inv => mel_inv%op
      op_u => mel_u%op
      njoined = op_inp%njoined

      ngraph = str_info%ngraph
      hpvx_occ => op_inp%ihpvca_occ
      idx_graph => mel_inp%idx_graph
      ca_occ => op_inp%ica_occ
      graphs => str_info%g

      ms_fix = .false.
      if(mel_inp%fix_vertex_ms.or.mel_inv%fix_vertex_ms
     &   .or.mel_u%fix_vertex_ms)then
        ms_fix = mel_inp%fix_vertex_ms.and.mel_inv%fix_vertex_ms
     &           .and.mel_u%fix_vertex_ms
        if(.not.ms_fix) call quit(1,'invsqrt',
     &                            'fix ms or not?')
      endif

      ! Check whether files are buffered.
      bufin = .false.
      if(ffinp%buffered) bufin = .true.
      bufout = .false.
      if(ffinv%buffered) bufout = .true.
      bufu = .false.
      if (get_u.and.ffu%buffered) bufu = .true.

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
     &       write(lulog,*)'Invert: input not incore'
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
     &       write(lulog,*)'Invert: output not incore'
        buffer_out => ffinv%buffer(1:)
      endif

      if(get_u.and..not.bufu)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls
          nbuff = nbuff + mel_u%len_op_occ(iocc_cls)
        enddo
        ifree= mem_alloc_real(buffer_u,nbuff,'buffer_u')
        buffer_u(1:nbuff) = 0d0
      else if (get_u) then
        if(ntest.ge.100)
     &       write(lulog,*)'Invert: output (2) not incore'
        buffer_u => ffu%buffer(1:)
      endif

      icnt_sv  = 0 ! we will count the
      icnt_sv0 = 0 ! number of singular values below threshold
      xmax = 0d0   ! largest excluded singular value
      xmin = 1234567890d0   ! smallest included singular value
      allocate(blk_used(nocc_cls),blk_redundant(nocc_cls),
     &         bins(17,nocc_cls),ex2occ_cls(nocc_cls))
      bins = 0 ! binning for singular values:
               ! >10E0,>10E-1,...,>10E-15,0
      blk_used(1:nocc_cls) = .false.
      blk_redundant(1:nocc_cls) = .true.
      iexc_cls = 0
      tocc_cls = 0

      ! (assuming name 'T' for cluster op.
      op_t => op_info%op_arr(idx_oplist2('T',op_info))%op
      if (is_keyword_set('method.MRCI').gt.0)
     &   op_t => op_info%op_arr(idx_oplist2('C',op_info))%op

      if (.not.half.and.max(iprlvl,ntest).ge.3) write(lulog,*)
     &         'Input list will be overwritten by projector.'

      ! Loop over occupation class.
      iocc_loop: do iocc_cls = 1, nocc_cls
        iblkoff = (iocc_cls-1)*njoined
        if (occ_is_diag_blk(hpvx_occ(1,1,iblkoff+1),njoined))
     &     tocc_cls = tocc_cls + 1 ! increment occ.cls of corresp. Op.
        if(op_inp%formal_blk(iocc_cls)) cycle iocc_loop
        if (blk_used(iocc_cls)) cycle iocc_loop
        blk_used(iocc_cls) = .true.
        iexc_cls = iexc_cls + 1
        ex2occ_cls(iexc_cls) = tocc_cls

        if (ntest.ge.10) write(lulog,*) 'current occ_cls: ',iocc_cls
        ! only one element? easy!
        ! (also regularization is never needed in this case)
        if (mel_inp%len_op_occ(iocc_cls).eq.1) then
          icnt_cur = icnt_sv - icnt_sv0
          ioff = mel_inp%off_op_gmo(iocc_cls)%gam_ms(1,1)
          buffer_out(ioff+1) = buffer_in(ioff+1)
          if (lmodspc) then
            call mat_svd_traf(1,buffer_out(ioff+1),buffer_in(ioff+1),
     &                       icnt_sv,icnt_sv0,xmax,xmin,
     &                       bins(1,iexc_cls))
          else if (get_u) then
            call invsqrt_mat(1,buffer_out(ioff+1),buffer_in(ioff+1),
     &                       half,buffer_u(ioff+1),get_u,
     &                       xdum,icnt_sv,icnt_sv0,xmax,xmin,
     &                       bins(1,iexc_cls))
          else
            call invsqrt_mat(1,buffer_out(ioff+1),buffer_in(ioff+1),
     &                       half,xdummy,get_u, !buffer_u: dummy
     &                       xdum,icnt_sv,icnt_sv0,xmax,xmin,
     &                       bins(1,iexc_cls))
          end if
          if (sgrm.and.icnt_cur.lt.icnt_sv-icnt_sv0)
     &       blk_redundant(iocc_cls) = .false.
          cycle iocc_loop
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
          icnt_cur = icnt_sv - icnt_sv0

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
     &           write(lulog,'(a,3i8)') 'msa, gama, ndim:',msa,igama,
     &           int(sqrt(dble(mel_inp%len_op_gmo(iocc_cls)%
     &                         gam_ms(igama,idxmsa))))
              if (ntest.ge.100)
     &           write(lulog,*) ' len = ',
     &             mel_inp%len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa),
     &             ' ndis = ',ndis

              ioff = mel_inv%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)

              call get_flipmap_blk(flipmap_c,
     &            ncblk,occ_csub,ndim,graph_csub,idxmsc,igamc,
     &            strmap_info,ngam,ngraph)

              ! single distribution can simply be read in as simple matrix
              allocate(scratch(ndim,ndim),svs(ndim))
              if (.not.half.or.lmodspc) then
                allocate(scratch2(ndim,ndim))
              else
                scratch2 => xdummy
              end if
              if (get_u) then
                allocate(scratch3(ndim,ndim))
              else
                scratch3 => xdummy
              end if
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

              if (msc.eq.0.and..not.lmodspc) then
                ! here a splitting into "singlet" and "triplet" blocks is needed:

c dbg
c                write(lulog,*) 'flmap:'
c                do icol = 1, ndim
c                  write(lulog,'(i4,3i6)') icol,flipmap_c(icol)
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
                if (.not.half) then
                  allocate(sing2(nsing,nsing),trip2(ntrip,ntrip))
                else
                  sing2 => xdummy
                  trip2 => xdummy
                end if
                if (get_u) then
                  allocate(sing3(nsing,nsing),trip3(ntrip,ntrip))
                else
                  sing3 => xdummy
                  trip3 => xdummy
                end if
                call spinsym_traf(1,ndim,scratch,flipmap_c,nsing,
     &                            sing,trip,.false.)

                ! calculate T^(-0.5) for both blocks
                call invsqrt_mat(nsing,sing,sing2,half,sing3,get_u,
     &                           svs,icnt_sv,icnt_sv0,
     &                           xmax,xmin,bins(1,iexc_cls))
                call invsqrt_mat(ntrip,trip,trip2,half,trip3,get_u,
     &                           svs(nsing+min(1,ntrip)),!avoid segfault
     &                           icnt_sv,icnt_sv0,
     &                           xmax,xmin,bins(1,iexc_cls))

                ! partial undo of pre-diagonalization: Upre*T^(-0.5)
                call spinsym_traf(2,ndim,scratch,flipmap_c,nsing,
     &                            sing,trip,.true.)
                if (.not.half) then
                  ! full undo of pre-diagonalization for projector
                  call spinsym_traf(2,ndim,scratch2,flipmap_c,nsing,
     &                              sing2,trip2,.false.)
                  deallocate(sing2,trip2)
                end if
                if (get_u) then
                  ! partial undo of pre-diagonalization: Upre*U
                  call spinsym_traf(2,ndim,scratch3,flipmap_c,nsing,
     &                              sing3,trip3,.true.)
                  deallocate(sing3,trip3)
                end if
                deallocate(sing,trip)
              else

               if (lmodspc) then
                call mat_svd_traf(ndim,scratch,scratch2,
     &                           icnt_sv,icnt_sv0,xmax,xmin,
     &                           bins(1,iexc_cls))
               else
                ! calculate S^(-0.5)
                call invsqrt_mat(ndim,scratch,scratch2,
     &                           half,scratch3,get_u,svs,
     &                           icnt_sv,icnt_sv0,xmax,xmin,
     &                           bins(1,iexc_cls))
               end if

              end if

              ! Tikhonov regularization?
              if (reg_tik.and..not.lmodspc)
     &           call regular_tikhonov(ndim,ndim,scratch,svs,omega2)

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

              deallocate(scratch,svs)

              if (.not.half.or.lmodspc) then
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

              if (get_u) then
                ! write unitary matrix to buffer
                if (transp) then
                  do idx = 1,ndim
                    do jdx = 1,ndim
                      buffer_u((idx-1)*ndim+jdx+ioff)
     &                       = scratch3(jdx,idx)
                    enddo
                  enddo
                else
                  do idx = 1,ndim
                    do jdx = 1,ndim
                      buffer_u((idx-1)*ndim+jdx+ioff)
     &                       = scratch3(idx,jdx)
                    enddo
                  enddo
                end if

                deallocate(scratch3)
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
          if (sgrm.and.icnt_cur.lt.icnt_sv-icnt_sv0)
     &       blk_redundant(iocc_cls) = .false.
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
             if (nrank.gt.maxrank) call quit(1,'invsqrt',
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
             if (nrank.gt.maxrank) call quit(1,'invsqrt',
     &                                       'increase maxrank')
             do idx = nrank, 2, -1
               rankdim(idx) = rankdim(idx-1)
               rankoff(idx) = rankoff(idx-1)
             end do
             rankdim(1) = 1
             rankoff(1) = 0
             if (nrank.ge.2) rankoff(1) = rankoff(2) + rankdim(2)
          else if (ndim.ne.rankoff(1)+rankdim(1)) then
            call quit(1,'invsqrt','dimensions don''t add up!')
          end if

          if (ntest.ge.10)
     &       write(lulog,'(a,3i8)') 'ms1, igam, ndim:',ms1,igam,ndim
          if (ntest.ge.100)
     &       write(lulog,'(a,5i8)') 'dim. per rank:',rankdim(1:nrank)

          allocate(scratch(ndim,ndim),flmap(ndim,3),svs(ndim))
          scratch = 0d0
          svs = 0d0
          if (.not.half.or.lmodspc) then
            allocate(scratch2(ndim,ndim))
            ! initialize: important for off-dia. blks for GNO/seq.orth.
            scratch2(1:ndim,1:ndim) = 0d0
          end if
          if (get_u) allocate(scratch3(ndim,ndim))
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
            ! exception for pure inactive block:
            if (na1+nc1+na2+nc2.eq.0) then
              if (igam.ne.1) exit
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
     &                   idx_msgmdst2(.true.,
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
                    msdst(idx,2,1) = msdis_a(nablk2)
                    igamdst(idx,2,1) = gamdis_a(nablk2)
                  end if
                  if (iocc2(idx,1,2).ne.0) then
                    msdst(idx,1,2) = msdis_c(ncblk2)
                    igamdst(idx,1,2) = gamdis_c(ncblk2)
                  end if
                end do

c dbg
c                write(lulog,'(a,4i4)') 'ms  : ',msa1,msc1,msa2,msc2
c                write(lulog,'(a,4i4)') 'gam : ',gama1,gamc1,gama2,gamc2
c                write(lulog,'(a,2i4)') 'msa, igama: ',msa,igama
c                write(lulog,'(a,2i4)') 'dist, len: ',idxdis, lenca
c                write(lulog,'(a,1i4)') 'dist2: ',idxdis2
c                write(lulog,'(a,2i4)') 'off_line/col: ',off_line,off_col
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
     &                     na2+nc2,2,iocc2,idx_g2,
     &                     msdst,igamdst,first,
     &                     str_info%igas_restr,
     &                     orb_info%mostnd,orb_info%igamorb,
     &                     orb_info%nsym,orb_info%ngas,
     &                     orb_info%ngas_hpv,orb_info%idx_gas,
     &                     hpvxseq,lexlscr)
                          first = .false.
                          if (.not.logdum) call quit(1,'invsqrt',
     &                         'no next tuple found!')
c                          if (mod(idxcount(2,idspn,na2+nc2,1),4).ne.0)
                          if (mod(na2+nc2
     &                            -idxcount(2,idspn,na2+nc2,1),4).ne.0)
     &                          flmap(icol,3) = -1
c dbg
c                          write(lulog,'(i8,x,4i4,x,4i4)')idx,idorb,idspn
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
c            if (ntest.ge.100) then
c              write(lulog,*) 'initial overlap matrix:'
c              call wrtmat2(scratch,ndim,ndim,ndim,ndim)
c            end if
c dbgend
c dbg
c          print *,'index matrix:'
c          write(*,'(5x,18i4)') (icol, icol=1,ndim)
c          do iline = 1, ndim
c            write(*,'(i4,x,18i4)') iline, matrix(iline,1:ndim)
c          end do
c dbgend

          ! loop over blocks that should be orthogonalized separately
          do irank = 1, nrank
           icnt_cur = icnt_sv - icnt_sv0
           rdim = rankdim(irank)
           idxst = rankoff(irank) + 1
           idxnd = rankoff(irank) + rdim

           ! build projector and apply to current block
           if (irank.ge.2) then
            if (ntest.ge.10)
     &        write(lulog,'(x,a,i8)') 'next rank:',irank
            select case(project)
            case (1,2)
             if (project.eq.1) then
              ! write tensors "p" on off-diagonal blocks of scratch
              do jrank = 1, irank-1
                rdim2 = rankdim(jrank)
                idxst2 = rankoff(jrank) + 1
                idxnd2 = rankoff(jrank) + rdim2
                scratch(idxst:idxnd,idxst2:idxnd2) = 0d0
                call get_vec2remove(scratch(idxst:idxnd,idxst2:idxnd2),
     &                             scratch(idxst2:idxnd2,idxst2:idxnd2),
     &                             rdim,rdim2,ms1,igam,irank-jrank,
     &                             na1mx-nrank+irank,nc1mx-nrank+irank,
     &                             gno,pass_spc,mel_spc,na1mx,nc1mx,
     &                             orb_info,str_info,strmap_info)
              end do
             else ! if project.eq.2
              ! multiply off-diagonal blocks with lower trafo matrices
              do jrank = 1, irank-1
               rdim2 = rankdim(jrank)
               idxst2 = rankoff(jrank) + 1
               idxnd2 = rankoff(jrank) + rdim2
               call dgemm('t','n',rdim,rdim2,rdim2,
     &                   1d0,scratch(idxst2:idxnd2,idxst:idxnd),rdim2,
     &                   scratch(idxst2:idxnd2,idxst2:idxnd2),rdim2,
     &                   0d0,scratch(idxst:idxnd,idxst2:idxnd2),rdim)
c dbg
               ! transpose (only needed for double-check below)
               do idx = 0, rdim2-1
                 scratch(idxst2+idx,idxst:idxnd) 
     &           = scratch(idxst:idxnd,idxst2+idx)
               end do
c dbgend
              end do
              ! singular value decompose off-diagonal block
              ! and overwrite by nonredundant part of W^+ (left-hand)
              call svd_get_left(rdim,ndim-idxnd,
     &                          scratch(idxst:idxnd,idxnd+1:ndim))
             end if
             allocate(proj(rdim,rdim))
             if (project.eq.1) then
              allocate(scratch4(rdim,rdim))
              ! form outer product of these vectors
              call dgemm('n','t',rdim,rdim,ndim-idxnd,
     &                   1d0,scratch(idxst:idxnd,idxnd+1:ndim),rdim,
     &                   scratch(idxst:idxnd,idxnd+1:ndim),rdim,
     &                   0d0,scratch4,rdim)
              ! assemble projector: 1 - sum(p*p^T)*S
              call dgemm('n','n',rdim,rdim,rdim,
     &                   -1d0,scratch4,rdim,
     &                   scratch(idxst:idxnd,idxst:idxnd),rdim,
     &                   0d0,proj,rdim)
             else ! if project.eq.2
              ! form outer product -W^+ * W
              call dgemm('n','t',rdim,rdim,ndim-idxnd,
     &                   -1d0,scratch(idxst:idxnd,idxnd+1:ndim),rdim,
     &                   scratch(idxst:idxnd,idxnd+1:ndim),rdim,
     &                   0d0,proj,rdim)
             end if
             do idx = 1, rdim
               proj(idx,idx) = proj(idx,idx) + 1d0
             end do
             if (ntest.ge.100) then
               write(lulog,*) 'Projector for removing lower-rank exc.:'
               call wrtmat2(proj,rdim,rdim,rdim,rdim)
             end if
             ! Apply projector to current metric block
             ! (a) Q^+*S
             !     This step is not strictly necessary, but numerically
             !     it helps to keep the matrix S*Q symmetric
             if (project.eq.2) allocate(scratch4(rdim,rdim))
             call dgemm('t','n',rdim,rdim,rdim,
     &                  1d0,proj,rdim,
     &                  scratch(idxst:idxnd,idxst:idxnd),rdim,
     &                  0d0,scratch4,rdim)
             ! (b) (Q^+*S)*Q
             call dgemm('n','n',rdim,rdim,rdim,
     &                  1d0,scratch4,rdim,proj,rdim,
     &                  0d0,scratch(idxst:idxnd,idxst:idxnd),rdim)
             deallocate(scratch4)
c dbg
             if (project.eq.2) then
               ! check that off-diagonal blocks are projected out
               call dgemm('t','t',rdim,ndim-idxnd,rdim,
     &                    1d0,proj,rdim,
     &                    scratch(idxnd+1:ndim,idxst:idxnd),ndim-idxnd,
     &                    0d0,scratch(idxst:idxnd,idxnd+1:ndim),rdim)
               if (ntest.ge.100) then
                 write(lulog,*) 'Projected off-diagonal block:'
                 call wrtmat2(scratch(idxst:idxnd,idxnd+1:ndim),
     &                        rdim,ndim-idxnd,rdim,ndim-idxnd)
               end if
               if (any(abs(scratch(idxst:idxnd,idxnd+1:ndim)).gt.1d-12))
     &           call warn('invsqrt','off-dia block not projected out?')
             end if
c dbgend
             if (gno.eq.1.or.project.eq.2) then
               ! zero off-diagonal blocks (assuming the projector works)
               scratch(idxst:idxnd,idxnd+1:ndim) = 0d0
               scratch(idxnd+1:ndim,idxst:idxnd) = 0d0
             end if

            case (3)
             ! Gram-Schmidt step: substract lower-rank components from
             ! the diagonal block and the blocks above
             do jrank = irank, nrank ! target block: row j, column i
              rdim2 = rankdim(jrank)
              idxst2 = rankoff(jrank) + 1
              idxnd2 = rankoff(jrank) + rdim2
              do krank = 1, irank-1 ! input blocks: row k, columns j and i
               rdim3 = rankdim(krank)
               idxst3 = rankoff(krank) + 1
               idxnd3 = rankoff(krank) + rdim3
               call dgemm('t','n',rdim2,rdim,rdim3,
     &                  -1d0,scratch(idxst3:idxnd3,idxst2:idxnd2),rdim3,
     &                  scratch(idxst3:idxnd3,idxst:idxnd),rdim3,
     &                  1d0,scratch(idxst2:idxnd2,idxst:idxnd),rdim2)
              end do
             end do
             ! compute the blocks below (intermediately stored on the right)
             ! for the formation of the trafo matrix
             do jrank = 1, irank-1 ! target block: row i, column j
              rdim2 = rankdim(jrank)
              idxst2 = rankoff(jrank) + 1
              idxnd2 = rankoff(jrank) + rdim2
              scratch(idxst:idxnd,idxst2:idxnd2) = 0d0 !initialize
              do krank = jrank, irank-1 ! input blocks: (k,i) and (j,k)
               rdim3 = rankdim(krank)
               idxst3 = rankoff(krank) + 1
               idxnd3 = rankoff(krank) + rdim3
               call dgemm('t','t',rdim,rdim2,rdim3,
     &                  -1d0,scratch(idxst3:idxnd3,idxst:idxnd),rdim3,
     &                  scratch(idxst2:idxnd2,idxst3:idxnd3),rdim2,
     &                  1d0,scratch(idxst:idxnd,idxst2:idxnd2),rdim)
              end do
             end do

            case default
              call quit(1,'invsqrt','unknown case for keyword project')
            end select
           end if

          ! core step: singular value decomposition
          if (ms1.eq.0.and.rdim.gt.1.and..not.lmodspc) then
            ! here a splitting into "singlet" and "triplet" blocks is needed:

c dbg
c            write(lulog,*) 'flmap:'
c            do icol = idxst, idxnd
c              write(lulog,'(i4,2i6)') icol,flmap(icol,1:2)
c            end do
c dbgend
            do icol = idxst, idxnd
              idx = idxlist(flmap(icol,1),flmap(idxst:idxnd,2),rdim,1)
              if (idx.eq.-1) call quit(1,'invsqrt','idx not found!')
              flmap(icol,3) = flmap(icol,3)*idx
            end do
c dbg
c            write(lulog,*) 'flmap:'
c            do icol = idxst, idxnd
c              write(lulog,'(i4,3i6)') icol,flmap(icol,1:3)
c            end do
c dbgend
            nsing = rdim
            do icol = idxst, idxnd
              if (abs(flmap(icol,3)).eq.icol-idxst+1) nsing = nsing
     &                              + sign(1,flmap(icol,3))
            end do
            nsing = nsing/2
            ntrip = rdim - nsing

            ! do the pre-diagonalization
            allocate(sing(nsing,nsing),trip(ntrip,ntrip))
            if (.not.half) then
              allocate(sing2(nsing,nsing),trip2(ntrip,ntrip))
            else
              sing2 => xdummy
              trip2 => xdummy
            end if
            if (get_u) then
              allocate(sing3(nsing,nsing),trip3(ntrip,ntrip))
            else
              sing3 => xdummy
              trip3 => xdummy
            end if
            call spinsym_traf(1,rdim,
     &                        scratch(idxst:idxnd,idxst:idxnd),
     &                        flmap(idxst:idxnd,3),nsing,
     &                        sing,trip,.false.)

            ! calculate T^(-0.5) for both blocks
            call invsqrt_mat(nsing,sing,sing2,half,sing3,get_u,
     &                       svs(idxst),icnt_sv,icnt_sv0,
     &                       xmax,xmin,bins(1,iexc_cls))
            call invsqrt_mat(ntrip,trip,trip2,half,trip3,get_u,
     &                       svs(idxst-1+nsing+min(1,ntrip)),!avoid segfault
     &                       icnt_sv,icnt_sv0,
     &                       xmax,xmin,bins(1,iexc_cls))

            ! partial undo of pre-diagonalization: Upre*T^(-0.5)
            call spinsym_traf(2,rdim,
     &                        scratch(idxst:idxnd,idxst:idxnd),
     &                        flmap(idxst:idxnd,3),nsing,
     &                        sing,trip,.true.)
            if (.not.half) then
              ! full undo of pre-diagonalization for projector
              call spinsym_traf(2,rdim,
     &                          scratch2(idxst:idxnd,idxst:idxnd),
     &                          flmap(idxst:idxnd,3),nsing,
     &                          sing2,trip2,.false.)
              deallocate(sing2,trip2)
            end if
            if (get_u) then
              ! partial undo of pre-diagonalization: Upre*U
              call spinsym_traf(2,rdim,
     &                          scratch3(idxst:idxnd,idxst:idxnd),
     &                          flmap(idxst:idxnd,3),nsing,
     &                          sing3,trip3,.true.)
              deallocate(sing3,trip3)
            end if
            deallocate(sing,trip)

          else if (lmodspc) then
            call mat_svd_traf(rdim,scratch(idxst:idxnd,idxst:idxnd),
     &                       scratch2(idxst:idxnd,idxst:idxnd),
     &                       icnt_sv,icnt_sv0,xmax,xmin,
     &                       bins(1,iexc_cls))
          else if (.not.half.and.get_u) then
            ! calculate S^(-0.5)
            call invsqrt_mat(rdim,scratch(idxst:idxnd,idxst:idxnd),
     &                       scratch2(idxst:idxnd,idxst:idxnd),
     &                       half,
     &                       scratch3(idxst:idxnd,idxst:idxnd),get_u,
     &                       svs(idxst),icnt_sv,icnt_sv0,xmax,xmin,
     &                       bins(1,iexc_cls))
          else if (.not.half) then
            ! calculate S^(-0.5)
            call invsqrt_mat(rdim,scratch(idxst:idxnd,idxst:idxnd),
     &                       scratch2(idxst:idxnd,idxst:idxnd),
     &                       half,xdummy,get_u, !scratch3: dummy
     &                       svs(idxst),icnt_sv,icnt_sv0,xmax,xmin,
     &                       bins(1,iexc_cls))
          else if (get_u) then
            ! calculate S^(-0.5)
            call invsqrt_mat(rdim,scratch(idxst:idxnd,idxst:idxnd),
     &                       xdummy,half, !scratch2: dummy
     &                       scratch3(idxst:idxnd,idxst:idxnd),get_u,
     &                       svs(idxst),icnt_sv,icnt_sv0,xmax,xmin,
     &                       bins(1,iexc_cls))
          else
            ! calculate S^(-0.5)
            call invsqrt_mat(rdim,scratch(idxst:idxnd,idxst:idxnd),
     &                       xdummy,half, !scratch2: dummy
     &                       xdummy,get_u, !scratch3: dummy
     &                       svs(idxst),icnt_sv,icnt_sv0,xmax,xmin,
     &                       bins(1,iexc_cls))
          end if

           select case(project)
           case (1,2)
            ! apply projector again: X = Q*U*s^(-0.5)
            if (irank.ge.2) then
              allocate(scratch4(rdim,rdim))
              call dgemm('n','n',rdim,rdim,rdim,
     &                   1d0,proj,rdim,
     &                   scratch(idxst:idxnd,idxst:idxnd),rdim,
     &                   0d0,scratch4,rdim)
              scratch(idxst:idxnd,idxst:idxnd) = scratch4
              if (ntest.ge.100) then
                write(lulog,*) 'Trafo matrix:'
                call wrtmat2(scratch(idxst:idxnd,idxst:idxnd),
     &                       rdim,rdim,rdim,rdim)
              end if
              if (.not.half) then
                call dgemm('n','n',rdim,rdim,rdim,
     &                     1d0,proj,rdim,
     &                     scratch2(idxst:idxnd,idxst:idxnd),rdim,
     &                     0d0,scratch4,rdim)
                scratch2(idxst:idxnd,idxst:idxnd) = scratch4
                if (ntest.ge.100) then
                  write(lulog,*) 'Projector matrix:'
                  call wrtmat2(scratch2(idxst:idxnd,idxst:idxnd),
     &                         rdim,rdim,rdim,rdim)
                end if
              end if
              deallocate(scratch4,proj)
            end if

           case (3)
            ! multiply off-diagonal blocks of trafo matrix with dia.blks
            ! and set upper off-diagonal blocks to zero
            do jrank = 1, irank-1 ! target: (j,i) input: (i,j) and (i,i)
             rdim2 = rankdim(jrank)
             idxst2 = rankoff(jrank) + 1
             idxnd2 = rankoff(jrank) + rdim2
             call dgemm('t','n',rdim2,rdim,rdim,
     &                  1d0,scratch(idxst:idxnd,idxst2:idxnd2),rdim,
     &                  scratch(idxst:idxnd,idxst:idxnd),rdim,
     &                  0d0,scratch(idxst2:idxnd2,idxst:idxnd),rdim2)
             scratch(idxst:idxnd,idxst2:idxnd2) = 0d0
            end do
            do jrank = irank+1, nrank ! target: (i,j) input: (i,i) and (j,i)
             rdim2 = rankdim(jrank)
             idxst2 = rankoff(jrank) + 1
             idxnd2 = rankoff(jrank) + rdim2
             call dgemm('t','t',rdim,rdim2,rdim,
     &                  1d0,scratch(idxst:idxnd,idxst:idxnd),rdim,
     &                  scratch(idxst2:idxnd2,idxst:idxnd),rdim2,
     &                  0d0,scratch(idxst:idxnd,idxst2:idxnd2),rdim)
            end do

           end select

           if (sgrm.and.icnt_cur.lt.icnt_sv-icnt_sv0)
c     &        blk_redundant(iocc_cls+min(na1mx,nc1mx)+irank-nrank)
c     &        blk_redundant(iocc_cls+irank-1)
     &        blk_redundant(iocc_cls+nrank-irank)
     &        = .false.
          end do

          ! Tikhonov regularization?
          if (reg_tik.and..not.lmodspc)
     &       call regular_tikhonov(ndim,ndim,scratch,svs,omega2)

c dbg
c            if (ntest.ge.100) then
c              write(lulog,*) 'final transformation matrix:'
c              call wrtmat2(scratch,ndim,ndim,ndim,ndim)
c            end if
c dbgend

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
              buffer_out(ioff+1) = scratch(off_line2+1,off_col2+1)
              if (.not.half.or.lmodspc) ! copy projector to input buffer
     &           buffer_in(ioff+1) = scratch2(off_line2+1,off_col2+1)
              if (get_u)
     &           buffer_u(ioff+1) = scratch3(off_line2+1,off_col2+1)
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
                        buffer_out(idx) = scratch(iline,icol)
                        if (.not.half.or.lmodspc) ! copy projector to input buffer
     &                     buffer_in(idx) = scratch2(iline,icol)
                        if (get_u)
     &                     buffer_u(idx) = scratch3(iline,icol)
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

          deallocate(scratch,flmap,svs)
          if (.not.half.or.lmodspc) deallocate(scratch2)
          if (get_u) deallocate(scratch3)
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
        write(lulog,'(x,77("="))')
        write(lulog,'(x,a)')
     &       'Singular value histogram (by excitation classes)'
        if (lmodspc) then
          write(lulog,'(x,a,x,14i10)') 'class:', (idx,idx=1,iexc_cls)
        else
          write(lulog,'(x,a,x,14i10)') 'n_h = ',
     &         op_t%ihpvca_occ(IHOLE,2,ex2occ_cls(1:iexc_cls))
          write(lulog,'(x,a,x,14i10)') 'n_p = ',
     &         op_t%ihpvca_occ(IPART,1,ex2occ_cls(1:iexc_cls))
        end if
        write(lulog,'(x,77("-"))')
        write(lulog,'(x,a,x,14i10)') ' 1E+00',bins(1,1:iexc_cls)
c        if (lmodspc) then
c          do idx = 2, 8
c            write(lulog,'(x,a,i2.2,x,14i10)') ' 1E-',idx-1,
c     &                                        bins(idx,1:iexc_cls)
c          end do
c          write(lulog,'(x,a,x,14i10)') '     0',bins(9,1:iexc_cls)
c          do idx = 10, 16
c            write(lulog,'(x,a,i2.2,x,14i10)') '-1E-',17-idx,
c     &                                        bins(idx,1:iexc_cls)
c          end do
c          write(lulog,'(x,a,x,14i10)') '-1E+00',bins(17,1:iexc_cls)
c          write(lulog,'(x,77("="))')
c          write(lulog,'(x,i8,a,i8,a)') icnt_sv-icnt_sv0,
c     &          ' out of ',icnt_sv,' eigenvalues were included'
c          write(lulog,'(x,a,E19.10)')
c     &          'The  largest excluded eigenvalue is ',xmax
c          write(lulog,'(x,a,E19.10)')
c     &          'The smallest included eigenvalue is ',xmin
c        else
        do idx = 2, 16
          write(lulog,'(x,a,i2.2,x,14i10)') ' 1E-',idx-1,
     &                                      bins(idx,1:iexc_cls)
        end do
        write(lulog,'(x,a,x,14i10)') '     0',bins(17,1:iexc_cls)
        write(lulog,'(x,77("="))')
        write(lulog,'(x,i8,a,i8,a)') icnt_sv-icnt_sv0,
     &        ' out of ',icnt_sv,' singular values were included'
        write(lulog,'(x,a,E19.10)')
     &        'The  largest excluded singular value is ',xmax
        write(lulog,'(x,a,E19.10)')
     &        'The smallest included singular value is ',xmin
c        end if
      end if
      deallocate(bins,ex2occ_cls)
 
      if (sgrm.and.any(blk_redundant(1:nocc_cls))) then
        ! Print out which blocks are redundant
        write(lulog,'(x,a)') 'There are redundant blocks:'
c        write(lulog,'(x,a)') 'There are redundant blocks in T:'
        write(lulog,'(x,a,26i3)') 'Block #  :',(idx,idx=1,nocc_cls)
        write(lulog,'(x,a,26(2x,L1))') 'Redundant?',
     &                               (blk_redundant(idx),idx=1,nocc_cls)
        call get_argument_value('method.MR','svdonly',lval=svdonly)
        if (svdonly.and..not.lmodspc) then
          ! Excplicitly print restrictions for input file
          write(lulog,*)
          write(lulog,'(x,a)') 'Copy the following into the input file:'
          ih = -1
          ip = -1
          iexc = 0
          first = .false.
          do iocc_cls = 1, nocc_cls
            first2 = .false.
            if (ih.ne.op_t%ihpvca_occ(IHOLE,2,iocc_cls)) then
              ih = op_t%ihpvca_occ(IHOLE,2,iocc_cls)
              first = .true.
              first2 = .true.
            end if
            if (ip.ne.op_t%ihpvca_occ(IPART,1,iocc_cls)) then
              ip = op_t%ihpvca_occ(IPART,1,iocc_cls)
              first = .true.
              first2 = .true.
            end if
            if (first2) then
              ! how many blocks belong to this excitation class?
              jocc_cls = nocc_cls
              do while(op_t%ihpvca_occ(IHOLE,2,jocc_cls).ne.ih.or.
     &                 op_t%ihpvca_occ(IPART,1,jocc_cls).ne.ip)
                jocc_cls = jocc_cls - 1
              end do
              nrank = jocc_cls - iocc_cls
              irank = 0
            else
              irank = irank + 1
            end if
            iexc = op_t%ica_occ(1,iocc_cls)
            if (blk_redundant(iocc_cls+nrank-2*irank).and.first) then
              first = .false.
              write(lulog,'(x,a,i1,a,i1,a,i1,a,i1,a,i1,a)')
     &          'MR excrestr=(',ih,',',ih,',',ip,',',ip,',1,',iexc-1,')'
            else if (.not.blk_redundant(iocc_cls+nrank-2*irank).and.
     &               .not.first) then
              write(lulog,'(x,a,i4)') 'Watch out: Non-redundant block:',
     &                                iocc_cls
              call warn('invsqrt',
     &                  'non-redundant block beyond redundant one?')
            end if
          end do
          write(lulog,*)
        end if
      end if
      deallocate(blk_redundant)

      if(.not.bufout)then
        call put_vec(ffinv,buffer_out,1,nbuff)
      endif  
      if(.not.bufin.and.(.not.half.or.lmodspc))then
        ! return projector matrix on input list
        call put_vec(ffinp,buffer_in,1,nbuff)
      endif
      if (.not.bufu.and.get_u) call put_vec(ffu,buffer_u,1,nbuff)

      ifree = mem_flushmark('invsqrt')

      return
      end
