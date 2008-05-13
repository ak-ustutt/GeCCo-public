*----------------------------------------------------------------------*
      subroutine symmetrise_blk1blk2(buffer_out1,buffer_out2,
     &     buffer_in1,buffer_in2,fac,facd,
     &     mel,iblk_1,iblk_2,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     core routine for symmetrizer
*
*     andreas, may 2008
*     
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
      type(me_list), intent(inout) ::
     &     mel
      integer, intent(in) ::
     &     iblk_1, iblk_2
      real(8), intent(inout) ::
     &     buffer_in1(*), buffer_out1(*),
     &     buffer_in2(*), buffer_out2(*)
      real(8), intent(in) ::
     &     fac, facd

      logical ::
     &     first
      integer ::
     &     nocc_cls, njoined,
     &     ifree, nblk, nbuff, idxmsa, idxmsc, idxdis_1,
     &     idxdis_2, ioff_1, ioff_2, ioff0_1, ioff0_2,
     &     msa, msc, igama, igamc, idxa, idxc, ngam, lena, lenc,
     &     iblkoff, ncblk, nablk, msc_max, msa_max
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
     &     gamdis_c(:), gamdis_a(:),
     &     len_str(:),
     &     hpvx_occ(:,:,:), ca_occ(:,:), idx_graph(:,:,:)

      real(8), external ::
     &     ddot
      logical, external ::
     &     iocc_equal_n, next_msgamdist2
      integer, external ::
     &     ielprd, idx_msgmdst2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'symmetrise_blk1blk2')
      end if

c dbg
c      print *,'buffer_in1: ',buffer_in1(1:mel%len_op_occ(iblk_1))
c      print *,'buffer_in2: ',buffer_in2(1:mel%len_op_occ(iblk_2))
c dbg
      op  => mel%op

      ioff0_1 = mel%off_op_occ(iblk_1)
      ioff0_2 = mel%off_op_occ(iblk_2)

      njoined  = op%njoined
      iblkoff = (iblk_1-1)*njoined

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      hpvx_occ => op%ihpvca_occ
      idx_graph => mel%idx_graph
      ca_occ => op%ica_occ
      graphs => str_info%g

      call get_num_subblk(ncblk,nablk,
     &     op%ihpvca_occ(1,1,iblkoff+1),njoined)

      allocate(hpvx_csub(ncblk),hpvx_asub(nablk),
     &         occ_csub(ncblk), occ_asub(nablk),
     &         graph_csub(ncblk), graph_asub(nablk),
     &         msdis_c(ncblk),  msdis_a(nablk),
     &         idxmsdis_c(ncblk),  idxmsdis_a(nablk),
     &         gamdis_c(ncblk), gamdis_a(nablk),
     &         len_str(ncblk+nablk))

      ! set HPVX and OCC info
      call condense_occ(occ_csub, occ_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  hpvx_occ(1,1,iblkoff+1),njoined,hpvxblkseq)
      ! do the same for the graph info
      call condense_occ(graph_csub, graph_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  idx_graph(1,1,iblkoff+1),njoined,hpvxblkseq)

      ! Loop over Ms of annihilator string of output array
      idxmsa = 0
      msc_max = ca_occ(1,iblk_1)
      msa_max = ca_occ(2,iblk_1)

      idxmsa = 0
      msa_loop : do msa = msa_max, -msa_max, -2

        msc = msa + mel%mst

        if (abs(msc).gt.msc_max) cycle msa_loop
        idxmsa = idxmsa+1
        idxmsc = (msc_max-msc)/2 + 1

        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam          

          igamc = multd2h(igama,mel%gamt)

          if (ntest.ge.100)
     &         write(luout,*) 'MS(A), GAMMA(A): ',msa,igama,' len = ',
     &           mel%len_op_gmo(iblk_1)%gam_ms(igama,idxmsa)

          if (mel%len_op_gmo(iblk_1)%
     &         gam_ms(igama,idxmsa).le.0) cycle

          first = .true.
          idxdis_1 = 0
          distr_loop: do

            if (.not.next_msgamdist2(first,
     &            msdis_c,msdis_a,gamdis_c,gamdis_a,
     &            ncblk, nablk,
     &            occ_csub,occ_asub,
     &            msc,msa,igamc,igama,ngam)) exit

            call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncblk)
            call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nablk)

            call set_len_str(len_str,ncblk,nablk,
     &                         graphs,
     &                         graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                         graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                         hpvxseq,.false.)

            lenc = ielprd(len_str,ncblk)
            lena = ielprd(len_str(ncblk+1),nablk)

            if (lenc.eq.0.or.lena.eq.0) cycle

            idxdis_1 = idxdis_1+1

            first = .false.

            if (ntest.ge.1000)
     &         write(luout,*) 'dist: ',idxdis_1,' len = ',
     &             mel%len_op_gmox(iblk_1)%
     &             d_gam_ms(idxdis_1,igama,idxmsa)

            ioff_1 = mel%off_op_gmox(iblk_1)%
     &             d_gam_ms(idxdis_1,igama,idxmsa) - ioff0_1

            idxdis_2 = 1
            if (mel%off_op_gmox(iblk_2)%ndis(igamc,idxmsc).gt.1)
     &           idxdis_2 =
     &               idx_msgmdst2(
     &                iblk_2,idxmsa,igama,
     &                occ_csub,idxmsdis_c,gamdis_c,ncblk,
     &                occ_asub,idxmsdis_a,gamdis_a,nablk,
     &                .true.,mel,ngam)

            ioff_2 = mel%off_op_gmox(iblk_2)%
     &             d_gam_ms(idxdis_2,igamc,idxmsc) - ioff0_2

            if (ntest.ge.1000) then
              write(luout,*) 'idxdis_1,idxdis_2: ',idxdis_1,idxdis_2
              write(luout,*) 'ioff_1,ioff_2:     ',ioff_1,ioff_2
            end if

            idxc_loop: do idxc = 1, lenc
              idxa_loop: do idxa = 1, lena

                value =fac*buffer_in1(ioff_1+(idxc-1)*lena+idxa) +
     &                 fac*buffer_in2(ioff_2+(idxa-1)*lenc+idxc)

                buffer_out1(ioff_1+(idxc-1)*lena+idxa) = value
                buffer_out2(ioff_2+(idxa-1)*lenc+idxc) = value

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
     &         len_str)

      return
      end
