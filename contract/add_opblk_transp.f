      subroutine add_opblk_transp(xnorm2,type,fac,
     &     me_in,me_out,tra_in,tra_out,
     &     iblk_in,iblk_out,
     &     op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to add a transposed list:
*
*     OUT_{I_C,I_A} = OUT_{I_C,I_A} + IN_{I_A,I_C} 
*
*     andreas, march 2008
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
     &     ntest = 00

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      logical, intent(in) ::
     &     tra_in, tra_out
      type(me_list), intent(inout) ::
     &     me_in, me_out
      integer, intent(in) ::
     &     iblk_in, iblk_out, type
      real(8), intent(in) ::
     &     fac
      real(8), intent(out) ::
     &     xnorm2

      logical ::
     &     bufin, bufout, open_close_in, open_close_out, same, first,
     &     ms_fix, fix_success
      integer ::
     &     nocc_cls, njoined,
     &     ifree, nblk, nbuff, idxmsa, idxmsc, idxdis,
     &     idxdis_in, ioff_out, ioff_in, ioff0_out, ioff0_in,
     &     msa, msc, igama, igamc, idxa, idxc, ngam, lena, lenc,
     &     iblkoff, iblkoff_in, ncblk, nablk, msc_max, msa_max,
     &     istr, idx_in, idx_out, icmp
      real(8) ::
     &     fac_off, fac_dia, value

      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      type(graph), pointer ::
     &     graphs(:)
      type(operator), pointer ::
     &     op_in, op_out
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
     &     ldim_opin_c(:), ldim_opin_a(:),
     &     ldim_opout_c(:), ldim_opout_a(:),
     &     istr_csub(:), istr_asub(:),
     &     dis_map_c(:), dis_map_a(:)

      real(8), external ::
     &     ddot
      logical, external ::
     &     iocc_equal_n, next_msgamdist2
      integer, external ::
     &     ielprd, idx_msgmdst2, idx_str_blk3

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'add_opblk_transp')
        write(luout,*) 'input list:  ',trim(me_in%label),
     &       ' op: ',trim(me_in%op%name)
        write(luout,*) 'output list: ',trim(me_out%label),
     &       ' op: ',trim(me_out%op%name)
      end if

      ffin  => me_in%fhand
      ffout => me_out%fhand

      ms_fix = .false.
      if(me_in%fix_vertex_ms.or.me_out%fix_vertex_ms)then
        ms_fix = me_in%fix_vertex_ms.and.me_out%fix_vertex_ms
        if(.not.ms_fix) call quit(1,'add_opblk_transp',
     &                            'fix ms or not?')
      endif
        
      ! in and out on same file?
      same = associated(ffin,ffout)
      
      op_in  => me_in%op
      op_out => me_out%op

      njoined  = op_in%njoined
      iblkoff = (iblk_out-1)*njoined
      iblkoff_in = (iblk_in-1)*njoined

      if (njoined.ne.op_out%njoined.or.
     &    .not.iocc_equal_n(op_in%ihpvca_occ(1,1,iblkoff_in+1),tra_in,
     &                      op_out%ihpvca_occ(1,1,iblkoff+1),  tra_out,
     &                      njoined)) then
        if (njoined.ne.op_out%njoined) then
          write(luout,*) njoined, op_out%njoined
        else
          write(luout,*) 'IN, OUT: ',tra_in, tra_out
          call wrt_occ_n(luout,op_in%ihpvca_occ(1,1,iblkoff_in+1),
     &                                                          njoined)
          call wrt_occ_n(luout,op_out%ihpvca_occ(1,1,iblkoff_in+1),
     &                                                          njoined)
        end if
        call quit(1,'add_opblk_transp',
     &     'Input and output list do not have compatible operators')
      end if

      open_close_in  = ffin%unit.le.0
      open_close_out = ffout%unit.le.0

      if (open_close_in ) call file_open(ffin )
      if (open_close_out.and..not.same) call file_open(ffout)

      ! Check whether files are buffered.
      bufin = ffin%buffered
      bufout = ffout%buffered

      if (bufin.or.bufout)
     &   call quit(1,'add_opblk_tranposed',
     &               'buffered files are not yet debugged')
      ifree = mem_setmark('add_transp')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ioff0_in = me_in%off_op_occ(iblk_in)
      ioff0_out = me_out%off_op_occ(iblk_out)

      ! Allocations made to full block length
      if(.not.bufin)then
        nbuff = me_in%len_op_occ(iblk_in)

        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
        call get_vec(ffin,buffer_in,ioff0_in+1,ioff0_in+nbuff)

      else
        buffer_in => ffin%buffer(ioff0_in+1:)
      endif

      if(.not.bufout.or.same)then
        nbuff = me_out%len_op_occ(iblk_out)
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
        call get_vec(ffout,buffer_out,ioff0_out+1,ioff0_out+nbuff)
      else
        buffer_out => ffout%buffer(ioff0_out+1:)
      endif

      hpvx_occ => op_out%ihpvca_occ
      idx_graph => me_out%idx_graph
      ca_occ => op_out%ica_occ
      graphs => str_info%g

      call get_num_subblk(ncblk,nablk,
     &     op_out%ihpvca_occ(1,1,iblkoff+1),njoined)

      allocate(hpvx_csub(ncblk),hpvx_asub(nablk),
     &         occ_csub(ncblk), occ_asub(nablk),
     &         graph_csub(ncblk), graph_asub(nablk),
     &         msdis_c(ncblk),  msdis_a(nablk),
     &         idxmsdis_c(ncblk),  idxmsdis_a(nablk),
     &         gamdis_c(ncblk), gamdis_a(nablk),
     &         len_str(ncblk+nablk),
     &         istr_csub(ncblk),istr_asub(nablk),
     &         ldim_opin_c(ncblk),ldim_opin_a(nablk),
     &         ldim_opout_c(ncblk),ldim_opout_a(nablk),
     &         dis_map_c(ncblk),dis_map_a(nablk))

      ! set HPVX and OCC info
      call condense_occ(occ_csub, occ_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  hpvx_occ(1,1,iblkoff+1),njoined,hpvxblkseq)
      ! do the same for the graph info
      call condense_occ(graph_csub, graph_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  idx_graph(1,1,iblkoff+1),njoined,hpvxblkseq)

      ! set map as some complicated changes in the sequence may happen
      call set_dis_tra_map(dis_map_c,dis_map_a,
     &     hpvx_csub,hpvx_asub,ncblk,nablk)

      ! Loop over Ms of annihilator string of output array
      idxmsa = 0
      msc_max = ca_occ(1,iblk_out)
      msa_max = ca_occ(2,iblk_out)

      idxmsa = 0
      msa_loop : do msa = msa_max, -msa_max, -2

        msc = msa + me_out%mst

        if (abs(msc).gt.msc_max) cycle msa_loop
        idxmsa = idxmsa+1
        idxmsc = (msc_max-msc)/2 + 1

        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam

          igamc = multd2h(igama,me_out%gamt)

          if (me_out%len_op_gmo(iblk_out)%
     &         gam_ms(igama,idxmsa).le.0) cycle

          first = .true.
          idxdis = 0
          distr_loop: do

            if (.not.next_msgamdist2(first,
     &            msdis_c,msdis_a,gamdis_c,gamdis_a,
     &            ncblk, nablk,
     &            occ_csub,occ_asub,
     &            msc,msa,igamc,igama,ngam,ms_fix,fix_success)) exit
            first = .false.
            if(ms_fix.and..not.fix_success)cycle distr_loop

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

            idxdis = idxdis+1

            ioff_out = me_out%off_op_gmox(iblk_out)%
     &             d_gam_ms(idxdis,igama,idxmsa) - ioff0_out

            idxdis_in = 1
            if (me_in%off_op_gmox(iblk_in)%ndis(igamc,idxmsc).gt.1)
     &           idxdis_in =
     &               idx_msgmdst2(
     &                iblk_in,idxmsa,igama,
     &                occ_csub,idxmsdis_c,gamdis_c,ncblk,
     &                occ_asub,idxmsdis_a,gamdis_a,nablk,
     &                .true.,dis_map_c,dis_map_a,me_in,ngam)

            ioff_in = me_in%off_op_gmox(iblk_in)%
     &             d_gam_ms(idxdis_in,igamc,idxmsc) - ioff0_in

            if (me_out%len_op_gmox(iblk_out)%
     &             d_gam_ms(idxdis,igama,idxmsa).le.0) cycle

            call set_op_ldim_c(ldim_opin_c,ldim_opin_a,
     &           hpvx_csub,hpvx_asub,
     &           len_str,ncblk,nablk,tra_in)
            call set_op_ldim_c(ldim_opout_c,ldim_opout_a,
     &           hpvx_csub,hpvx_asub,
     &           len_str,ncblk,nablk,tra_out)

            idxc_loop: do idxc = 1, lenc
              istr = idxc-1
              do icmp = 1, ncblk
                istr_csub(icmp) = mod(istr,len_str(icmp))
                istr = istr/len_str(icmp)
              end do
              idxa_loop: do idxa = 1, lena
                istr = idxa-1
                do icmp = 1, nablk
                  istr_asub(icmp) = mod(istr,len_str(ncblk+icmp))
                  istr = istr/len_str(icmp)
                end do

                idx_in  = ioff_in + idx_str_blk3(istr_csub,istr_asub,
     &               ldim_opin_c,ldim_opin_a,
     &               ncblk,nablk)
                idx_out = ioff_out + idx_str_blk3(istr_csub,istr_asub,
     &               ldim_opout_c,ldim_opout_a,
     &               ncblk,nablk)                

                buffer_out(idx_out) =
     &                   buffer_out(idx_out) +
     &               fac*buffer_in (idx_in)

              end do idxa_loop
            end do idxc_loop

          end do distr_loop
          
        end do igama_loop
          
      end do msa_loop

      ! update norm^2
      if (type.eq.1) then
        xnorm2 = ddot(me_out%len_op_occ(iblk_out),
     &       buffer_out,1,buffer_out,1)
      else
        xnorm2 = buffer_out(1)
      end if

      if(.not.bufout)then
        call put_vec(ffout,buffer_out,ioff0_out+1,ioff0_out+nbuff)
      endif  

      if (open_close_in ) call file_close_keep(ffin)
      if (open_close_out.and..not.same) call file_close_keep(ffout)

      deallocate(hpvx_csub,hpvx_asub,
     &         occ_csub, occ_asub,
     &         graph_csub, graph_asub,
     &         msdis_c,  msdis_a,
     &         idxmsdis_c,  idxmsdis_a,
     &         gamdis_c, gamdis_a,
     &         len_str,
     &         istr_csub,istr_asub,
     &         ldim_opin_c,ldim_opin_a,
     &         ldim_opout_c,ldim_opout_a,
     &         dis_map_c, dis_map_a)

      ifree = mem_flushmark('add_transp')

      return
      end
