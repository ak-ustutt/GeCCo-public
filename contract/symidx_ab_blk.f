*----------------------------------------------------------------------*
      subroutine symidx_ab_blk(idxlist_out,
     &     idxlist_in,nlist,
     &     mel,iblk,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     core routine for spin symmetrizer (of idxlist)
*
*     andreas, nov 2008
*     
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
     &     ntest = 000

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(me_list), intent(inout) ::
     &     mel
      integer, intent(in) ::
     &     iblk, nlist, idxlist_in(nlist)
      integer, intent(inout) ::
     &     idxlist_out(nlist)

      logical ::
     &     first, ms_fix, ms_fix_ok
      integer ::
     &     nocc_cls, njoined, iel,
     &     ifree, nblk, nbuff, idxmsa, idxmsc, idxdis_1, idxdis_1_x,
     &     idxdis_2, ioff_1, ioff_2, ioff0_1, ioff0_2, ioff,
     &     msa, msc, igama, igamc, idxa, idxc, ngam, lena, lenc,
     &     msa2, msc2,  idxmsa2, idxmsc2,
     &     iblkoff, ncblk, nablk, msc_max, msa_max,
     &     istr, idx1, idx2, icmp, ngraph, maxbuf,
     &     asign, csign, gsign, imap
      real(8) ::
     &     fac_off, fac_dia, value

      type(graph), pointer ::
     &     graphs(:)
      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffin, ffout

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
     &     ldim_op_c(:), ldim_op_a(:),
     &     ldim_optr_c(:), ldim_optr_a(:),
     &     istr_csub(:), istr_asub(:),
     &     istr_csub_flip(:), istr_asub_flip(:)

      integer, pointer ::
     &     flipmap_c(:), flipmap_a(:)

      real(8), external ::
     &     ddot
      logical, external ::
     &     iocc_equal_n, next_msgamdist2
      integer, external ::
     &     ielprd, idx_msgmdst2, idx_str_blk3, std_spsign_msdis,
     &     list_in_bounds, idxlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'symidx_ab_blk')
      end if

      ifree = mem_setmark('symidx_ab_blk')

      if (mel%mst.ne.0) then
        call quit(1,'symidx_ab_blk',
     &       'I should only be called if the total '//
     &       'Ms-projection equals 0')
      end if

      ms_fix = mel%fix_vertex_ms

      op  => mel%op

      ioff0_1 = mel%off_op_occ(iblk)

      njoined  = op%njoined
      iblkoff = (iblk-1)*njoined

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ngraph = str_info%ngraph

      hpvx_occ => op%ihpvca_occ
      idx_graph => mel%idx_graph
      ca_occ => op%ica_occ
      graphs => str_info%g

      call get_num_subblk(ncblk,nablk,
     &     op%ihpvca_occ(1,1,iblkoff+1),njoined)
c dbg
c      call wrt_occ_n(6,op%ihpvca_occ(1,1,iblkoff+1),njoined)
c      print *,'ncblk,nablk: ',ncblk,nablk
c dbg

      allocate(hpvx_csub(ncblk),hpvx_asub(nablk),
     &         occ_csub(ncblk), occ_asub(nablk),
     &         graph_csub(ncblk), graph_asub(nablk),
     &         msdis_c(ncblk),  msdis_a(nablk),
     &         idxmsdis_c(ncblk),  idxmsdis_a(nablk),
     &         msdis_c2(ncblk),  msdis_a2(nablk),
     &         idxmsdis_c2(ncblk),  idxmsdis_a2(nablk),
     &         gamdis_c(ncblk), gamdis_a(nablk),
     &         len_str(ncblk+nablk),
     &         istr_csub(ncblk),istr_asub(nablk),
     &         istr_csub_flip(ncblk),istr_asub_flip(nablk),
     &         ldim_op_c(ncblk),ldim_op_a(nablk),
     &         ldim_optr_c(ncblk),ldim_optr_a(nablk))

      ! set HPVX and OCC info
      call condense_occ(occ_csub, occ_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  hpvx_occ(1,1,iblkoff+1),njoined,hpvxblkseq)
      ! do the same for the graph info
      call condense_occ(graph_csub, graph_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  idx_graph(1,1,iblkoff+1),njoined,hpvxblkseq)

      ! set flip maps
      call strmap_man_flip(
     &     maxbuf,
     &     graph_csub,ncblk,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(flipmap_c,maxbuf,'flipmap_c')
      call strmap_man_flip(
     &     maxbuf,
     &     graph_asub,nablk,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(flipmap_a,maxbuf,'flipmap_a')

      ! Loop over Ms of annihilator string of output array
      idxmsa = 0
      msc_max = ca_occ(1,iblk)
      msa_max = ca_occ(2,iblk)

      idxmsa = 0
      msa_loop : do msa = msa_max, -msa_max, -2

        msc = msa + mel%mst ! <- should be zero (checked above)

        if (abs(msc).gt.msc_max) cycle msa_loop

        msa2 = -msa
        msc2 = msa2 + mel%mst 

        idxmsa = idxmsa+1
        idxmsc = (msc_max-msc)/2 + 1
        idxmsa2 = (msa_max-msa2)/2 + 1
        idxmsc2 = (msc_max-msc2)/2 + 1

        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam          

          igamc = multd2h(igama,mel%gamt)

          if (ntest.ge.100)
     &         write(luout,*) 'MS(A), GAMMA(A): ',msa,igama,' len = ',
     &           mel%len_op_gmo(iblk)%gam_ms(igama,idxmsa)

          if (mel%len_op_gmo(iblk)%
     &         gam_ms(igama,idxmsa).le.0) cycle

          first = .true.
          idxdis_1_x = 0
          distr_loop: do

            if (.not.next_msgamdist2(first,
     &            msdis_c,msdis_a,gamdis_c,gamdis_a,
     &            ncblk, nablk,
     &            occ_csub,occ_asub,
     &            msc,msa,igamc,igama,ngam,
     &            ms_fix,ms_fix_ok)) exit
            first = .false.

            if (ms_fix.and..not.ms_fix_ok) cycle

            msdis_c2(1:ncblk) = -msdis_c(1:ncblk)
            msdis_a2(1:nablk) = -msdis_a(1:nablk)

            gsign =       std_spsign_msdis(msdis_c,occ_csub,ncblk)
            gsign = gsign*std_spsign_msdis(msdis_a,occ_asub,nablk)
            gsign = gsign*std_spsign_msdis(msdis_c2,occ_csub,ncblk)
            gsign = gsign*std_spsign_msdis(msdis_a2,occ_asub,nablk)

            call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncblk)
            call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nablk)
            call ms2idxms(idxmsdis_c2,msdis_c2,occ_csub,ncblk)
            call ms2idxms(idxmsdis_a2,msdis_a2,occ_asub,nablk)

            call set_len_str(len_str,ncblk,nablk,
     &                         graphs,
     &                         graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                         graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                         hpvxseq,.false.)

            lenc = ielprd(len_str,ncblk)
            lena = ielprd(len_str(ncblk+1),nablk)
c dbg
c            print *,'len_str: ',len_str
c dbg

            if (lenc.eq.0.or.lena.eq.0) cycle

            idxdis_1_x = idxdis_1_x+1
c test
            idxdis_1 = 1
            if (mel%off_op_gmox(iblk)%ndis(igama,idxmsa).gt.1)
     &           idxdis_1 =
     &               idx_msgmdst2(
     &                iblk,idxmsa,igama,
     &                occ_csub,idxmsdis_c,gamdis_c,ncblk,
     &                occ_asub,idxmsdis_a,gamdis_a,nablk,
     &                .false.,-1,-1,mel,ngam)
c            print *,'idxdis_1: ',idxdis_1,idxdis_1_x
            if (idxdis_1.ne.idxdis_1_x) print *,'!!!OHA!!!'
c test

            if (ntest.ge.1000)
     &         write(luout,*) 'dist: ',idxdis_1,' len = ',
     &             mel%len_op_gmox(iblk)%
     &             d_gam_ms(idxdis_1,igama,idxmsa)

            ioff_1 = mel%off_op_gmox(iblk)%
     &             d_gam_ms(idxdis_1,igama,idxmsa) !- ioff0_1

            if (.not.list_in_bounds(idxlist_in,nlist,
     &                              ioff_1+1,ioff_1+lenc*lena)) cycle

            idxdis_2 = 1
            if (mel%off_op_gmox(iblk)%ndis(igama,idxmsa2).gt.1)
     &           idxdis_2 =
     &               idx_msgmdst2(
     &                iblk,idxmsa2,igama,
     &                occ_csub,idxmsdis_c2,gamdis_c,ncblk,
     &                occ_asub,idxmsdis_a2,gamdis_a,nablk,
     &                .false.,-1,-1,mel,ngam)

            ioff_2 = mel%off_op_gmox(iblk)%
     &             d_gam_ms(idxdis_2,igama,idxmsa2) !- ioff0_1

            if (msa.eq.0.and.idxdis_1.gt.idxdis_2) cycle distr_loop
            if (ntest.ge.1000) then
              write(luout,*) 'idxdis_1,idxdis_2: ',idxdis_1,idxdis_2
              write(luout,*) 'ioff_1,ioff_2:     ',ioff_1,ioff_2
            end if

            call get_flipmap_blk(flipmap_c,
     &           ncblk,occ_csub,len_str,graph_csub,idxmsdis_c,gamdis_c,
     &           strmap_info,ngam,ngraph)
            call get_flipmap_blk(flipmap_a,
     &           nablk,occ_asub,len_str(ncblk+1),
     &                                  graph_asub,idxmsdis_a,gamdis_a,
     &           strmap_info,ngam,ngraph)

            call set_op_ldim_c(ldim_op_c,ldim_op_a,
     &           hpvx_csub,hpvx_asub,
     &           len_str,ncblk,nablk,.false.)
            
c dbg
c            print *,'ldim_op_c: ',ldim_op_c
c            print *,'ldim_op_a: ',ldim_op_a
c            print *,'flipmap_c: len=',len_str(1:ncblk)
c            print '(10i6)',flipmap_c(1:sum(len_str(1:ncblk)))
c            print *,'flipmap_a:'
c            print '(10i6)',
c     &           flipmap_a(1:sum(len_str(ncblk+1:ncblk+nablk)))
c dbg

            idxa_loop: do idxa = 1, lena
              istr = idxa-1
              ioff = 1
              asign = gsign
              do icmp = 1, nablk
                istr_asub(icmp) = mod(istr,len_str(ncblk+icmp))
                imap = flipmap_a(ioff+istr_asub(icmp))
                istr_asub_flip(icmp) = abs(imap)-1
                asign = asign*sign(1,imap)
                istr = istr/len_str(ncblk+icmp)
                ioff = ioff + len_str(ncblk+icmp)
              end do
              idxc_loop: do idxc = 1, lenc
                istr = idxc-1
                ioff = 1
                csign = 1
                do icmp = 1, ncblk
                  istr_csub(icmp) = mod(istr,len_str(icmp))
                  imap = flipmap_c(ioff+istr_csub(icmp))
                  istr_csub_flip(icmp) = abs(imap)-1
                  csign = csign*sign(1,imap)
                  istr = istr/len_str(icmp)
                  ioff = ioff + len_str(icmp)
                end do

                idx1 = ioff_1 + idx_str_blk3(
     &               istr_csub,istr_asub,
     &               ldim_op_c,ldim_op_a,
     &               ncblk,nablk)

                iel = idxlist(idx1,idxlist_in,nlist,1)
                if (iel.le.0) cycle

                idx2 = ioff_2 + idx_str_blk3(
     &               istr_csub_flip,istr_asub_flip,
     &               ldim_op_c,ldim_op_a,
     &               ncblk,nablk)

                idxlist_out(iel) = csign*asign*idx2

              end do idxc_loop
            end do idxa_loop

          end do distr_loop
          
        end do igama_loop
          
      end do msa_loop

      deallocate(hpvx_csub,hpvx_asub,
     &         occ_csub, occ_asub,
     &         graph_csub, graph_asub,
     &         msdis_c,  msdis_a,
     &         idxmsdis_c,  idxmsdis_a,
     &         msdis_c2,  msdis_a2,
     &         idxmsdis_c2,  idxmsdis_a2,
     &         gamdis_c, gamdis_a,
     &         len_str,
     &         istr_csub,istr_asub,
     &         istr_csub_flip,istr_asub_flip,
     &         ldim_op_c,ldim_op_a,
     &         ldim_optr_c,ldim_optr_a)

      ifree = mem_flushmark('symidx_ab_blk')

      return
      end
