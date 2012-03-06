      subroutine add_unity(fac,ssign,mel_out,iblkout,orb_info,str_info)
*----------------------------------------------------------------------*
*
*     Add (antisymmetrised) unit operator
*     ssign = 1: include sign flips for hole lines
*           = 0: total Ms value as additional factor (for Sz operator)
*           =-1: only spacial part is diagonal (for S+, S- operators)
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'
      include 'hpvxseq.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      real(8), intent(in)::
     &     fac
      type(me_list), intent(inout) ::
     &     mel_out
      type(strinf), intent(in), target ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     iblkout, ssign

      logical ::
     &     bufout,closeit, first, ms_fix, fix_success
      integer ::
     &     ifree, nbuff, njoined, join_off, iblk,
     &     idxmsa, msmax, msa, msc, igama, igamc,
     &     idxc, ngam, len_blk, ioff_blk, lenc, lena, icmp,
     &     ihpv, ica, ncblk, nablk, ioff_1, idx1, istr, idxdis_1,
     &     iocc(ngastp,2), idx_graph(ngastp,2), ms_sym_sign
      integer, pointer ::
     &     hpvx_csub(:),hpvx_asub(:),
     &     occ_csub(:), occ_asub(:),
     &     graph_csub(:), graph_asub(:),
     &     msdis_c(:),  msdis_a(:),
     &     idxmsdis_c(:),  idxmsdis_a(:),
     &     gamdis_c(:), gamdis_a(:),
     &     len_str(:), istr_csub(:),
     &     ldim_op_c(:), ldim_op_a(:),
     &     mapca(:), diag_idx(:), diag_ca(:)
      real(8) ::
     &     xsign, totsign

      type(graph), pointer ::
     &     graphs(:)
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      type(filinf), pointer ::
     &     ffout
      type(operator), pointer ::
     &     opout

      logical, external ::
     &     occ_is_diag_blk, next_msgamdist2, 
     &     nondia_blk, nondia_distr
      integer, external ::
     &     msa2idxms4op, ielprd, idx_str_blk3

      ffout => mel_out%fhand
      opout => mel_out%op
      ms_sym_sign = ssign
      if (ssign.eq.0) ms_sym_sign = 1

      if (ntest.ge.100) then
        write(luout,*) '=========================='
        write(luout,*) ' add_unity messing around'
        write(luout,*) '=========================='
        write(luout,*) ' mode = ',ssign
        write(luout,*) ' fac = ',fac
        write(luout,*) ' list = ',trim(mel_out%label)
        write(luout,*) ' ffout: ',trim(ffout%name),
     &                   ' rec: ',ffout%current_record
        write(luout,*) ' opout: ',opout%name(1:len_trim(opout%name)),
     &       '  block: ',iblkout
      end if

      closeit = .false.
      njoined = opout%njoined

      join_off=(iblkout-1)*njoined
      ms_fix = mel_out%fix_vertex_ms

      if (ngastp.ne.4) call quit(1,'add_unity','only for ngastp=4')
      ! handle njoined>1 operator only if occ can be merged to njoined=1
      iocc = 0
      idx_graph = 0
      do iblk = 1, njoined
        do ihpv = 1, ngastp
          do ica = 1, 2
            if (opout%ihpvca_occ(ihpv,ica,join_off+iblk).ne.0) then
              if (iocc(ihpv,ica).ne.0)
     &            call quit(1,'add_unity','cannot handle this operator')
              iocc(ihpv,ica) = opout%ihpvca_occ(ihpv,ica,join_off+iblk)
              idx_graph(ihpv,ica) = mel_out%idx_graph(ihpv,ica,
     &                               join_off+iblk)
            end if
          end do
        end do
      end do

      ! check whether the out operator is a diagonal block:
      graphs => str_info%g

      if (.not.
     &     occ_is_diag_blk(iocc(1:ngastp,1:2),1))
     &     return ! else we just tacitly return

      bufout = .false.
      if(ffout%buffered) bufout = ffout%incore(iblkout).gt.0

      ifree = mem_setmark('add_unity')

      ngam = orb_info%nsym

      ! Find total block length.
      ioff_blk = mel_out%off_op_occ(iblkout)
      len_blk  = mel_out%len_op_occ(iblkout)

      if(len_blk.gt.ifree)
     &     call quit(1,'add_unity','insufficient space')
      
      ! Allocations made to maximum block length to save time.
      if(.not.bufout)then
        nbuff=len_blk
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
      else
        if(ntest.ge.100)
     &       write(luout,*)'Add_unity: Not incore'
        buffer_out => ffout%buffer(ioff_blk+1:)
      endif

      if(.not.bufout)then
        if(ffout%unit.le.0)then
          call file_open(ffout)
          closeit = .true.
        endif
        call get_vec(ffout,buffer_out,ioff_blk+1,ioff_blk+len_blk)
      endif  

      ! due to normal ordering, hole lines introduce sign flips:
      ! (for ssign.ne.1: maybe we will have to think about this again)
      if (njoined.eq.1) then
        xsign = dble(1-2*mod(iocc(1,1),2))
        if (ssign.ne.1) xsign = 1d0 ! maybe we will have to rethink this
      else if (.not.all(iocc(1:ngastp,1:2).eq.0)) then
        ! here we have to think again ... have a look at dia_from_blk!
        call quit(1,'add_unity','for njoined>1, think about the sign')
      else
        xsign = 1d0
      end if
      totsign = xsign

      call get_num_subblk(ncblk,nablk,iocc(1,1),1)

      allocate(hpvx_csub(ncblk),hpvx_asub(nablk),
     &         occ_csub(ncblk), occ_asub(nablk),
     &         graph_csub(ncblk), graph_asub(nablk),
     &         msdis_c(ncblk),  msdis_a(nablk),
     &         idxmsdis_c(ncblk),  idxmsdis_a(nablk),
     &         gamdis_c(ncblk), gamdis_a(nablk),
     &         len_str(ncblk+nablk),
     &         istr_csub(ncblk),
     &         ldim_op_c(ncblk),ldim_op_a(nablk))

      ! set HPVX and OCC info
      call condense_occ(occ_csub, occ_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  iocc(1,1),1,hpvxblkseq)
      ! do the same for the graph info
      call condense_occ(graph_csub, graph_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  idx_graph(1,1),1,hpvxblkseq)

      if (mel_out%diag_type.ne.0) then
        allocate(mapca(ncblk),diag_idx(ncblk),diag_ca(ncblk))
        if (nondia_blk(mapca,diag_idx,diag_ca,
     &                 iocc(1,1),1,ncblk,
     &                 nablk,mel_out%diag_type))
     &       call quit(1,'add_unity','non-diagonal block!')
      end if

      ! Loop over Ms of annihilator string.
      idxmsa = 0
c      msmax = 2
      msmax = opout%ica_occ(2,iblkout)
      if (msmax.ne.opout%ica_occ(1,iblkout)) call quit(1,'add_unity',
     &        'number of creators and annihilators must be the same')
      if (ssign.eq.-1.and.msmax.ne.1) call quit(1,'add_unity',
     &        'mode=-1 for one-particle operators only')
      msa_loop : do msa = msmax, -msmax, -2

cmh        idxmsa = idxmsa+1
        msc = msa + mel_out%mst
        idxmsa = msa2idxms4op(msa,mel_out%mst,msmax,msmax)
        ! Usually have mst=0 operators => Ms(c)=Ms(a)
      
        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam
          igamc = multd2h(igama,mel_out%gamt)

          if (mel_out%len_op_gmo(iblkout)%
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
            if (mel_out%diag_type.ne.0) then
              ! skip non-diagonal distributions ...
              if (nondia_distr(mapca,diag_idx,diag_ca,
     &              msdis_c,msdis_a,gamdis_c,gamdis_a,
     &              ncblk,mel_out%msdiag,
     &              mel_out%gamdiag)) cycle distr_loop
            end if

            call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncblk)
            call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nablk)
c dbg
c      print *,'ncblk,nablk:',ncblk,nablk
c      write(luout,'(a,10i4)') 'len_str:',len_str(1:ncblk+nablk)
c      write(luout,'(a,10i4)') 'graph_csub:',graph_csub(1:ncblk)
c      write(luout,'(a,10i4)') 'idxmsdis_c:',idxmsdis_c(1:ncblk)
c      write(luout,'(a,10i4)') 'gamdis_c:',gamdis_c(1:ncblk)
c      write(luout,'(a,10i4)') 'hpvx_csub:',hpvx_csub(1:ncblk)
c      write(luout,'(a,10i4)') 'graph_asub:',graph_asub(1:nablk)
c      write(luout,'(a,10i4)') 'idxmsdis_a:',idxmsdis_a(1:nablk)
c      write(luout,'(a,10i4)') 'gamdis_a:',gamdis_a(1:nablk)
c      write(luout,'(a,10i4)') 'hpvx_asub:',hpvx_asub(1:nablk)
c dbgend

            call set_len_str(len_str,ncblk,nablk,
     &                       graphs,
     &                       graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                       graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                       hpvxseq,.false.)

            lenc = ielprd(len_str,ncblk)
            lena = ielprd(len_str(ncblk+1),nablk)

            if (lenc.eq.0.or.lena.eq.0) cycle

            idxdis_1 = idxdis_1+1

            ! ms-dst and gamma-dst should also be diagonal:
            if (nablk.ne.ncblk) cycle distr_loop
            if (.not.all(msdis_c(1:ncblk)
     &                   -ms_sym_sign*msdis_a(1:ncblk).eq.0).or.
     &          .not.all(gamdis_c(1:ncblk)-gamdis_a(1:ncblk).eq.0))
     &           cycle distr_loop
            if (lenc.ne.lena) call quit(1,'add_unity','inconsistency')
            if (ssign.eq.0)
     &               totsign = xsign*0.5d0*dble(sum(msdis_c(1:ncblk)))

            if (ntest.ge.100) then
              write(luout,*) 'msc,msa,igamc,igama: ',
     &             msc,msa,igamc,igama
              write(luout,*) 'idxdis: ',idxdis_1
            end if

            ioff_1 = mel_out%off_op_gmox(iblkout)%
     &             d_gam_ms(idxdis_1,igama,idxmsa) - ioff_blk

            call set_op_ldim_c(ldim_op_c,ldim_op_a,
     &           hpvx_csub,hpvx_asub,
     &           len_str,ncblk,nablk,.false.)

            ! loop over diagonal me-list elements
            idxc_loop: do idxc = 1, lenc
              istr = idxc-1
              do icmp = 1, ncblk
                istr_csub(icmp) = mod(istr,len_str(icmp)) !+1
                istr = istr/len_str(icmp)
              end do

              idx1 = ioff_1 + idx_str_blk3(istr_csub,istr_csub,
     &             ldim_op_c,ldim_op_a,
     &             ncblk,nablk)

              buffer_out(idx1)=totsign*fac + buffer_out(idx1)

            end do idxc_loop

          end do distr_loop

        enddo igama_loop
          
      enddo msa_loop

      if (mel_out%diag_type.ne.0) deallocate(mapca,diag_idx,diag_ca)

      deallocate(hpvx_csub,hpvx_asub,
     &         occ_csub, occ_asub,
     &         graph_csub, graph_asub,
     &         msdis_c,  msdis_a,
     &         idxmsdis_c,  idxmsdis_a,
     &         gamdis_c, gamdis_a,
     &         len_str, istr_csub,
     &         ldim_op_c,ldim_op_a)

      if(.not.bufout)then
        call put_vec(ffout,buffer_out,ioff_blk+1,ioff_blk+len_blk)
      endif  

      call touch_file_rec(ffout)

      if(closeit)
     &     call file_close_keep(ffout)

      ifree = mem_flushmark('add_unity')

      return
      end
