*----------------------------------------------------------------------*
      subroutine dia_from_blk(buffer_out,buffer_inp,
     &     meinp,meout,iblkinp,iblkout,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     writes diagonal elements of buffer_inp (from meinp) to
*     buffer_out (to meout) for a given pair of operator blocks.
*     currently, the output block must be a single vertex block.
*     
*     matthias, fall 2009 (adapted from symmetrise_blk1blk2)
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

      integer, parameter ::
     &     ntest = 000

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(me_list), intent(in), target ::
     &     meinp, meout
      integer, intent(in) ::
     &     iblkinp, iblkout
      real(8), intent(inout) ::
     &     buffer_inp(*), buffer_out(*)

      logical ::
     &     first, ms_fix, fix_success
      integer ::
     &     njinp, njout, idxmsa, idxmsc, idxdis_1, idxdis_2,
     &     ioff_1, ioff_2, ioffinp, ioffout, msa, msc, 
     &     igama, igamc, idxa, idxc, ngam, lena, lenc, iblkoff,
     &     ncblk, nablk, msc_max, msa_max, istr, idx1, idx2,
     &     icmp, ncablk2, igamca, msca, msca_max, idxmsca
      real(8) ::
     &     fac

      type(graph), pointer ::
     &     graphs(:)
      type(operator), pointer ::
     &     opinp, opout
      type(filinf), pointer ::
     &     ffin, ffout

      integer, pointer ::
     &     hpvx_csub(:),hpvx_asub(:),
     &     occ_csub(:), occ_asub(:),
     &     graph_csub(:), graph_asub(:),
     &     msdis_c(:),  msdis_a(:),
     &     idxmsdis_c(:),  idxmsdis_a(:),
     &     gamdis_c(:), gamdis_a(:),
     &     len_str(:),
     &     hpvx_occ(:,:,:), ca_occ(:,:), idx_graph(:,:,:),
     &     ldim_op_c(:), ldim_op_a(:),
     &     ldim_op_c2(:), ldim_op_a2(:),
     &     istr_csub(:), istr_asub(:),
     &     istr_csub2(:), istr_asub2(:),
     &     hpvx_csub2(:),hpvx_asub2(:),
     &     occ_csub2(:),occ_asub2(:),
     &     graph_csub2(:), graph_asub2(:),
     &     msdis_c2(:),  msdis_a2(:),
     &     idxmsdis_c2(:),  idxmsdis_a2(:),
     &     gamdis_c2(:), gamdis_a2(:),
     &     len_str2(:)

      logical, external ::
     &     next_msgamdist2
      integer, external ::
     &     ielprd, idx_msgmdst2, idx_str_blk3

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'dia_from_blk')
      end if

      opinp => meinp%op
      opout => meout%op

      ioffinp = meinp%off_op_occ(iblkinp)
      ioffout = meout%off_op_occ(iblkout)

      njinp = opinp%njoined
      njout = opout%njoined
      if (njout.ne.1) call quit(1,'dia_from_blk',
     &         'nj>1 not yet available for output operator')
      iblkoff = (iblkout-1)*njout
      ms_fix = meout%fix_vertex_ms
      if (ms_fix.ne.meinp%fix_vertex_ms) call quit(1,'dia_from_blk',
     &         'inconsistent ms_fix')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      hpvx_occ => opout%ihpvca_occ
      idx_graph => meout%idx_graph
      ca_occ => opout%ica_occ
      graphs => str_info%g

      ! sign for splitting middle vertex (due to normal ordering):
      ! \  ...  /     /0 b c d\           /a b c d\
      ! /e b c d\ --> \e 0 0 0/       --> \e f g h/
      ! \e b c d/     /e 0 0 0\            ^-------- e odd : -1 
      ! /  ...  \     \0 b c d/                      e even: +1
      fac = dble(1-2*mod(hpvx_occ(1,2,iblkoff+1),2))

      call get_num_subblk(ncblk,nablk,
     &     opout%ihpvca_occ(1,1,iblkoff+1),njout)
      ncablk2 = ncblk+nablk

      allocate(hpvx_csub(ncblk),hpvx_asub(nablk),
     &         occ_csub(ncblk), occ_asub(nablk),
     &         graph_csub(ncblk), graph_asub(nablk),
     &         msdis_c(ncblk),  msdis_a(nablk),
     &         idxmsdis_c(ncblk),  idxmsdis_a(nablk),
     &         gamdis_c(ncblk), gamdis_a(nablk),
     &         len_str(ncblk+nablk),
     &         istr_csub(ncblk),istr_asub(nablk),
     &         ldim_op_c(ncblk),ldim_op_a(nablk),
     &         ldim_op_c2(ncablk2),ldim_op_a2(ncablk2),
     &         istr_csub2(ncablk2),istr_asub2(ncablk2),
     &         hpvx_csub2(ncablk2),hpvx_asub2(ncablk2),
     &         occ_csub2(ncablk2),occ_asub2(ncablk2),
     &         graph_csub2(ncablk2), graph_asub2(ncablk2),
     &         msdis_c2(ncablk2),  msdis_a2(ncablk2),
     &         idxmsdis_c2(ncablk2),  idxmsdis_a2(ncablk2),
     &         gamdis_c2(ncablk2), gamdis_a2(ncablk2),
     &         len_str2(2*ncablk2))

      ! set HPVX and OCC info
      call condense_occ(occ_csub, occ_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  hpvx_occ(1,1,iblkoff+1),njout,hpvxblkseq)
      ! do the same for the graph info
      call condense_occ(graph_csub, graph_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  idx_graph(1,1,iblkoff+1),njout,hpvxblkseq)

      ! Loop over Ms of annihilator string of output array
      idxmsa = 0
      msc_max = ca_occ(1,iblkout)
      msa_max = ca_occ(2,iblkout)
      msca_max = msc_max + msa_max

      idxmsa = 0
      msa_loop : do msa = msa_max, -msa_max, -2

        msc = msa + meout%mst
        msca = msa + msc

        if (abs(msc).gt.msc_max) cycle msa_loop
        idxmsa = idxmsa+1
        idxmsc = (msc_max-msc)/2 + 1
        idxmsca = (msca_max-msca)/2 + 1

        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam          

          igamc = multd2h(igama,meout%gamt)
          igamca = multd2h(igama,igamc)

          if (ntest.ge.100)
     &         write(luout,*) 'MS(A), GAMMA(A): ',msa,igama,' len = ',
     &           meout%len_op_gmo(iblkout)%gam_ms(igama,idxmsa)
          if (ntest.ge.1000)
     &         write(luout,*) 'ms(A)2, gamma(A)2: ',msca,igamca,
     &           ' len2 = ',
     &           meinp%len_op_gmo(iblkinp)%gam_ms(igamca,idxmsca)

          if (meout%len_op_gmo(iblkout)%
     &         gam_ms(igama,idxmsa).le.0) cycle

          ! loop over distributions of output operator
          first = .true.
          idxdis_1 = 0
          distr_loop: do

            if (.not.next_msgamdist2(first,
     &            msdis_c,msdis_a,gamdis_c,gamdis_a,
     &            ncblk, nablk,
     &            occ_csub,occ_asub,
     &            msc,msa,igamc,igama,ngam,
     &            ms_fix,fix_success)) exit
            first = .false.
            if(ms_fix.and..not.fix_success)cycle distr_loop

            call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncblk)
            call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nablk)

            ! get ms, graph, hpvx strings of input operator,
            ! assuming that output operator represents the diagonal elements
            call diag_condensed_occ(occ_csub2,occ_asub2,
     &                              occ_csub,occ_asub,
     &                              hpvx_occ(1,1,iblkoff+1),
     &                              njout,hpvxseq)
            call diag_condensed_occ(hpvx_csub2,hpvx_asub2,
     &                              hpvx_csub,hpvx_asub,
     &                              hpvx_occ(1,1,iblkoff+1),
     &                              njout,hpvxseq)
            call diag_condensed_occ(graph_csub2,graph_asub2,
     &                              graph_csub,graph_asub,
     &                              hpvx_occ(1,1,iblkoff+1),
     &                              njout,hpvxseq)
            call diag_condensed_occ(msdis_c2,msdis_a2,
     &                              msdis_c,msdis_a,
     &                              hpvx_occ(1,1,iblkoff+1),
     &                              njout,hpvxseq)
            call diag_condensed_occ(gamdis_c2,gamdis_a2,
     &                              gamdis_c,gamdis_a,
     &                              hpvx_occ(1,1,iblkoff+1),
     &                              njout,hpvxseq)

            call ms2idxms(idxmsdis_c2,msdis_c2,occ_csub2,ncablk2)
            call ms2idxms(idxmsdis_a2,msdis_a2,occ_asub2,ncablk2)

            call set_len_str(len_str,ncblk,nablk,
     &                         graphs,
     &                         graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                         graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                         hpvxseq,.false.)

            call set_len_str(len_str2,ncablk2,ncablk2,
     &                         graphs,
     &                         graph_csub2,idxmsdis_c2,gamdis_c2,
     &                         hpvx_csub2,
     &                         graph_asub2,idxmsdis_a2,gamdis_a2,
     &                         hpvx_asub2,
     &                         hpvxseq,.false.)

            lenc = ielprd(len_str,ncblk)
            lena = ielprd(len_str(ncblk+1),nablk)

            if (lenc.eq.0.or.lena.eq.0) cycle

            idxdis_1 = idxdis_1+1


            if (ntest.ge.1000)
     &         write(luout,*) 'dist: ',idxdis_1,' len = ',
     &             meout%len_op_gmox(iblkout)%
     &             d_gam_ms(idxdis_1,igama,idxmsa)

            ioff_1 = meout%off_op_gmox(iblkout)%
     &             d_gam_ms(idxdis_1,igama,idxmsa) - ioffout

            ! find corresponding distribution of input operator
            idxdis_2 = 1
            if (meinp%off_op_gmox(iblkinp)%ndis(igamca,idxmsca).gt.1)
     &           idxdis_2 =
     &               idx_msgmdst2(
     &                iblkinp,idxmsca,igamca,
     &                occ_csub2,idxmsdis_c2,gamdis_c2,ncablk2,
     &                occ_asub2,idxmsdis_a2,gamdis_a2,ncablk2,
     &                .false.,meinp,ngam)

            ioff_2 = meinp%off_op_gmox(iblkinp)%
     &             d_gam_ms(idxdis_2,igamca,idxmsca) - ioffinp

            if (ntest.ge.1000) then
              write(luout,*) 'idxdis_1,idxdis_2: ',idxdis_1,idxdis_2
              write(luout,*) 'ioff_1,ioff_2:     ',ioff_1,ioff_2
            end if

            call set_op_ldim_c(ldim_op_c,ldim_op_a,
     &           hpvx_csub,hpvx_asub,
     &           len_str,ncblk,nablk,.false.)
            call set_op_ldim_c(ldim_op_c2,ldim_op_a2,
     &           hpvx_csub2,hpvx_asub2,
     &           len_str2,ncablk2,ncablk2,.false.)

            if (ntest.ge.1000) then
              write(luout,*) 'len_str(C):',len_str(1:ncblk)
              write(luout,*) 'len_str(A):',len_str(ncblk+1:ncblk+nablk)
              write(luout,*) 'ldim_op_c: ',ldim_op_c(1:ncblk)
              write(luout,*) 'ldim_op_a: ',ldim_op_a(1:nablk)
              write(luout,*) 'ldim_op_c2: ',ldim_op_c2(1:ncablk2)
              write(luout,*) 'ldim_op_a2: ',ldim_op_a2(1:ncablk2)
            end if

            ! loop over output me-list elements
            idxc_loop: do idxc = 1, lenc
              istr = idxc-1
              do icmp = 1, ncblk
                istr_csub(icmp) = mod(istr,len_str(icmp)) !+1
                istr = istr/len_str(icmp)
              end do
              idxa_loop: do idxa = 1, lena
                istr = idxa-1
                do icmp = 1, nablk
                  istr_asub(icmp) = mod(istr,len_str(ncblk+icmp))!+1
                  istr = istr/len_str(ncblk+icmp)
                end do

                ! find indices for corresponding diagonal element of input mel
                call diag_condensed_occ(istr_csub2,istr_asub2,
     &                                  istr_csub,istr_asub,
     &                                  hpvx_occ(1,1,iblkoff+1),
     &                                  njout,hpvxseq)


                idx1 = ioff_1 + idx_str_blk3(istr_csub,istr_asub,
     &               ldim_op_c,ldim_op_a,
     &               ncblk,nablk)
                idx2 = ioff_2 + idx_str_blk3(istr_csub2,istr_asub2,
     &               ldim_op_c2,ldim_op_a2,
     &               ncablk2,ncablk2)

                ! write diagonal element of input operator to output list
                buffer_out(idx1) = fac*buffer_inp(idx2)

              end do idxa_loop
            end do idxc_loop

          end do distr_loop
          
        end do igama_loop
          
      end do msa_loop

      deallocate(hpvx_csub,hpvx_asub,
     &         occ_csub, occ_asub,
     &         graph_csub, graph_asub,
     &         msdis_c,  msdis_a,
     &         idxmsdis_c,  idxmsdis_a,
     &         gamdis_c, gamdis_a,
     &         len_str,
     &         istr_csub,istr_asub,
     &         ldim_op_c,ldim_op_a,
     &         ldim_op_c2,ldim_op_a2,
     &         hpvx_csub2,hpvx_asub2,
     &         occ_csub2, occ_asub2,
     &         graph_csub2, graph_asub2,
     &         msdis_c2,  msdis_a2,
     &         idxmsdis_c2,  idxmsdis_a2,
     &         gamdis_c2, gamdis_a2,
     &         len_str2)

      return
      end
