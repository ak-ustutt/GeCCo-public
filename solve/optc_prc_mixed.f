*----------------------------------------------------------------------*
      subroutine optc_prc_mixed(me_grd,me_special,nspecial,
     &                            name_opt,w_shift,
     &                            nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                            orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     experimental routine for testing mixed preconditioner
*     if nincore>1: xbuf1 contains the gradient vector on entry
*                   and should contain the preconditioned gradient
*                   on exit
*                   xbuf2 contains DIA for op. blocks with version=1
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'hpvxseq.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_reorder_info.h'
      include 'multd2h.h'
c      include 'ifc_input.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nincore, lenbuf, nspecial
      type(me_list) ::
     &     me_grd
      type(me_list_array) ::
     &     me_special(nspecial)
      character(*), intent(in) ::
     &     name_opt
      real(8), intent(in) ::
     &     w_shift
      real(8), intent(inout), target ::
     &     xbuf1(lenbuf), xbuf2(lenbuf), xbuf3(lenbuf)
      type(orbinf), intent(in), target ::
     &     orb_info
      type(strinf),intent(in), target ::
     &     str_info
      type(strinf),intent(in), target ::
     &     strmap_info
      type(operator_info), intent(inout) ::
     &     op_info

      real(8) ::
     &     xdum
      logical ::
     &     ldum1, ldum2, ms_fix, fix_success
      integer ::
     &     ngrd, nblk_grd, nbmat, nxmat, nfmat,
     &     idxmsa, msa, msc, gama, gamc,  ngam,
     &     msc_max, msa_max, idxdis,  idx, idoff_grd,
     &     ncsub, nasub, nidx_p,
     &     nidx_cstr, nidx_astr, ms_cstr, ms_astr, gam_cstr, gam_astr,
     &     graph_cstr, graph_astr, len_cstr, len_astr, len_grd, ifree,
     &     idx_grd, idx_b, idx_x, mode,
     &     ngas, nspin, njoined, njoined_tmp, mstotal, gamtotal,
     &     iblk, iblk_off, idxms_bx, gam_bx, ld_bx, iblk_b,
     &     occ_b(ngastp,2),
     &     occ_blk_tmp(ngastp,2,2),
     &     rst_blk_tmp(2,orb_info%ngas,2,2,orb_info%nspin,2),
     &     rst_blk_reo(2,orb_info%ngas,2,2,orb_info%nspin,2)
      integer, target ::
     &     occ_blk_reo(ngastp,2,2),
     &     occ_blk_dag(ngastp,2,2),
     &     graph_blk_dag(ngastp,2,2)
      logical ::
     &     first, beyond_A, ca_reverse

      type(operator), pointer ::
     &     op_grd_reo
      type(me_list), pointer ::
     &     me_bmat, me_xmat, me_fmat, me_grd_reo
      type(graph), pointer ::
     &     graphs(:)
      type(reorder_info) ::
     &     reo_info

      integer ::
     &     f_hpvx(ngastp,2), msdst(ngastp,2), igamdst(ngastp,2),
     &     ioff_xsum(2*ngastp), nstr(2*ngastp), nincr(2*ngastp)

      integer, pointer ::
     &     hpvx_occ(:,:,:), idx_graph(:,:,:),
     &     ca_occ(:,:), occ_blk(:,:,:), rst_blk(:,:,:,:,:,:),
     &     igasca_restr(:,:,:,:,:,:),
     &     igas_restr(:,:,:,:,:),
     &     occ_blk_pnt(:,:,:), graph_blk_pnt(:,:,:),
     &     off_grd_d_gam_ms(:,:,:), len_grd_d_gam_ms(:,:,:),
     &     off_bx_gam_ms(:,:),
     &     mostnd(:,:,:), igamorb(:), idx_gas(:), ngas_hpv(:),
     &     hpvx_csub(:), hpvx_asub(:), occ_csub(:), occ_asub(:),
     &     graph_csub(:), graph_asub(:), msdis_c(:), msdis_a(:),
     &     idxmsdis_c(:), idxmsdis_a(:), gamdis_c(:), gamdis_a(:),
     &     len_str(:)
      real(8), pointer ::
     &     f_dia(:), xmat(:), bmat(:), xgrd_pnt(:), xgrd_buf(:),
     &     scrbuf(:)

      integer, external ::
     &     idxlist, iblk_occ, idx_oplist2, idx_mel_list, imltlist
      logical,external ::
     &     next_msgamdist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'mixed preconditioner')
      end if

      if (nincore.ne.3)
     &     call quit(1,'optc_prc_mixed',
     &     'currently: special preconditioning for nincore==3, only')

      if (nspecial.lt.1)
     &     call quit(1,'optc_prc_mixed',
     &     'nspecial must be >= 1')

      me_bmat => me_special(1)%mel
      nbmat = me_bmat%len_op
      allocate(bmat(nbmat))
      ! should be open
      call vec_from_da(me_bmat%fhand,1,bmat,nbmat)

      beyond_A = me_bmat%op%name.ne.op_b_inv

      if (beyond_A) then
        if (nspecial.lt.3)
     &       call quit(1,'optc_prc_mixed','nspecial must be >= 3')

        me_xmat => me_special(2)%mel
        nxmat = me_xmat%len_op

        allocate(xmat(nxmat))
        ! should be open as well
        call vec_from_da(me_xmat%fhand,1,xmat,nxmat)

        me_fmat => me_special(3)%mel
        nfmat = me_fmat%len_op

        allocate(f_dia(2*orb_info%ntoob))

        ! xtract diagonal
        call onedia_from_op(f_dia,me_fmat,.false.,orb_info)
        
      else
        allocate(xmat(1),f_dia(1))
      end if
      
      ngam = orb_info%nsym
      ngas = orb_info%ngas
      nspin = orb_info%nspin
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      idx_gas => orb_info%idx_gas
      ngas_hpv => orb_info%ngas_hpv

      nblk_grd = me_grd%op%n_occ_cls
      njoined  = me_grd%op%njoined
      hpvx_occ => me_grd%op%ihpvca_occ
      igasca_restr => me_grd%op%igasca_restr
      idx_graph => me_grd%idx_graph
      ca_occ => me_grd%op%ica_occ
      mstotal = me_grd%mst
      gamtotal = me_grd%gamt
      graphs => str_info%g
      igas_restr => str_info%igas_restr

      ms_fix = me_grd%fix_vertex_ms
      do idx = 1, nspecial
        if(ms_fix.or.me_special(idx)%mel%fix_vertex_ms)then
          ms_fix = ms_fix.and.me_special(idx)%mel%fix_vertex_ms
          if(.not.ms_fix) call quit(1,'optc_prc_mixed',
     &                              'fix ms or not?')
        endif
      enddo

c      if (njoined.ne.2)
c     &     call quit(1,'optc_prc_special',
c     &     'strange -- expected njoined==2')

      ! loops over occupations of GRD

      ! a) first use DIA for blocks with version = 1
      do iblk = 1, nblk_grd
        if (me_grd%op%blk_version(iblk).ne.1) cycle
        idx_grd = me_grd%off_op_occ(iblk)+1
        len_grd = me_grd%len_op_occ(iblk)
        if (ntest.ge.100) then
          write(luout,*) 'iblk = ',iblk
          write(luout,*) 'gradient vector before:'
          write(luout,*) xbuf1(idx_grd:idx_grd+len_grd)
        end if
        call diavc(xbuf1(idx_grd),xbuf1(idx_grd),1d0,
     &             xbuf2(idx_grd),w_shift,len_grd)
        if (ntest.ge.100) then
          write(luout,*) 'gradient vector afterwards:'
          write(luout,*) xbuf1(idx_grd:idx_grd+len_grd)
        end if
      end do

      ! b) then use BLK for blocks with version <> 1
      do iblk = 1, nblk_grd
        if (me_grd%op%blk_version(iblk).eq.1) cycle

        if (ntest.ge.100) then
          idx_grd = me_grd%off_op_occ(iblk)+1
          len_grd = me_grd%len_op_occ(iblk)
          write(luout,*) 'iblk = ',iblk
          write(luout,*) 'gradient vector before:'
          write(luout,*) xbuf1(idx_grd+1:idx_grd+len_grd)
        end if

        iblk_off = (iblk-1)*njoined

        occ_blk =>
     &       hpvx_occ(1:ngastp,1:2,iblk_off+1:iblk_off+njoined)
        rst_blk =>
     &       igasca_restr(1:ngastp,1:ngas,1:2,1:2,1:nspin,
     &                             iblk_off+1:iblk_off+njoined)
c        len_grd_d_gam_ms =>
c     &       me_grd%len_op_gmox(iblk)%d_gam_ms
        off_grd_d_gam_ms =>
     &       me_grd%off_op_gmox(iblk)%d_gam_ms

        if (ntest.ge.100) then
          write(luout,*) 'Now caring for GRD block: '
          call wrt_occ_n(luout,occ_blk,njoined)
        end if

        call get_num_subblk(ncsub,nasub,
     &       hpvx_occ(1,1,iblk_off+1),njoined)

        if (ntest.ge.100)
     &       write(luout,*) 'ncsub, nasub: ',ncsub, nasub

        mode = 1
        nidx_p = 1
        ca_reverse = .false.
        ! find the relevant block of B
        if (njoined.eq.2) then
          iblk_b = 1
        else if (njoined.eq.1) then
          if (name_opt.eq.op_cex.or.
     &        name_opt.eq.op_rp) then
            mode = 1
            nidx_p = occ_blk(IPART,1,1)
            occ_b = 0
            occ_b(IPART,1) = 1
            occ_b(IPART,2) = 1
          else if (
     &        name_opt.eq.op_cexbar.or.
     &        name_opt.eq.op_lp) then
            mode = 1
            nidx_p = occ_blk(IPART,2,1)
            occ_b = 0
            occ_b(IPART,1) = 1
            occ_b(IPART,2) = 1
            ca_reverse = .true.
          else if (name_opt.eq.op_cexx) then
            mode = 2
            nidx_p = occ_blk(IPART,1,1)
            occ_b = 0
            occ_b(IPART,1) = 2
            occ_b(IPART,2) = 2
          else
            occ_b(1:ngastp,1) = occ_blk(1:ngastp,1,1)
            occ_b(1:ngastp,2) = occ_blk(1:ngastp,1,1)
          end if
          iblk_b = iblk_occ(occ_b,.false.,me_bmat%op,1)
          if (iblk_b.le.0)
     &         call quit(1,'optc_prc_mixed',
     &         'did not find an appropriate block of B (precond)')
        else
          call quit(1,'optc_prc_mixed',
     &         'gradient -- njoined > 2 ??')
        end if

        off_bx_gam_ms => me_bmat%off_op_gmo(iblk_b)%gam_ms

        ! special case: scalar:
        if (ncsub.eq.0.and.nasub.eq.0) then
          idx_grd = off_grd_d_gam_ms(1,1,1)+1
          idx_b   = off_bx_gam_ms(1,1)+1
          xbuf1(idx_grd) = xbuf1(idx_grd)/bmat(idx_b)
          cycle
        end if

        if (ncsub.ne.1.and.ncsub.ne.2.and.nasub.ne.1) then
          write(luout,*) 'ncsub, nasub: ',ncsub,nasub
          call quit(1,'optc_prc_mixed','this is not what I expected')
        end if

        if (njoined.eq.1.and.mode.ne.nidx_p) then

          call init_reo_info(reo_info)

          call set_reo_special(1,reo_info,mode)

          ! cheat get_reo_info: append zero block, so both
          ! reordered and original occupation have njoined==2
          occ_blk_tmp = 0
          rst_blk_tmp = 0
          occ_blk_tmp(1:ngastp,1:2,1) = occ_blk(1:ngastp,1:2,1)
          rst_blk_tmp(1:2,1:ngas,1:2,1:2,1:nspin,1)
     &         = rst_blk(1:2,1:ngas,1:2,1:2,1:nspin,1)

          ! post-processing
          occ_blk_reo = occ_blk_tmp
          call get_reo_info(ldum1,ldum2,
     &         occ_blk_reo,occ_blk_tmp,
     &         rst_blk_reo,rst_blk_tmp,
     &         2,
     &         reo_info,str_info,orb_info)

          ! set up operator and list for re-ordered block
          call add_operator('GRD_REO',op_info)
          idx = idx_oplist2('GRD_REO',op_info)
          op_grd_reo => op_info%op_arr(idx)%op
          call set_uop2(op_grd_reo,'GRD_REO',
     &         occ_blk_reo,1,2,(/.true.,.true./),orb_info)
c          call add_me_list('L_GRD_REO',op_info)
          call define_me_list('L_GRD_REO','GRD_REO',
     &         me_grd%absym,me_grd%casym,me_grd%gamt,
     &         me_grd%s2,me_grd%mst,.false.,
     &         1,1,
     &         op_info,orb_info,str_info,strmap_info)
          idx = idx_mel_list('L_GRD_REO',op_info)
          me_grd_reo => op_info%mel_arr(idx)%mel

          len_grd = me_grd_reo%len_op
          ifree = mem_setmark('grd_reo')
          ifree = mem_alloc_real(xgrd_buf,len_grd,'reo_buf')
          xgrd_pnt => xgrd_buf

          ! make the buffer pointer point to the buffer
          if (me_grd_reo%fhand%buffered)
     &         call quit(1,'optc_prc_mixed',
     &         'OK, adapt the hand-made buffering in this routine ...')
          me_grd_reo%fhand%buffered = .true.
          allocate(me_grd_reo%fhand%incore(1))
          me_grd_reo%fhand%incore(1) = 1
          me_grd_reo%fhand%buffer => xgrd_pnt

          if (me_grd%fhand%buffered)
     &         call quit(1,'optc_prc_mixed',
     &         'OK, adapt the hand-made buffering in this routine ...')
          me_grd%fhand%buffered = .true.
          allocate(me_grd%fhand%incore(nblk_grd))
          me_grd%fhand%incore(1:nblk_grd) = 0
          me_grd%fhand%incore(iblk) = 1
          me_grd%fhand%buffer => xbuf1
          ! we must resort
          idoff_grd = 0
          call reo_op_wmaps_c(
     &         .false.,xdum,0,
     &         me_grd,me_grd_reo,
     &         .false.,.false.,
     &         iblk,1,
     &         idoff_grd,0,
     &         reo_info,
     &         str_info,strmap_info,orb_info)

          call get_num_subblk(ncsub,nasub,
     &       occ_blk_reo,2)
          occ_blk_pnt => occ_blk_reo
          idx_graph => me_grd_reo%idx_graph
          graph_blk_pnt => idx_graph(1:ngastp,1:2,1:2)
          njoined_tmp = 2
          off_grd_d_gam_ms =>
     &         me_grd_reo%off_op_gmox(1)%d_gam_ms
        else
c          if (.not.ca_reverse) then
            occ_blk_pnt =>
     &         hpvx_occ(1:ngastp,1:2,iblk_off+1:iblk_off+njoined)
            graph_blk_pnt => idx_graph(1:ngastp,1:2,
     &         iblk_off+1:iblk_off+njoined)
c          else
c            occ_blk_dag(1:ngastp,1,1:njoined) =
c     &           hpvx_occ(1:ngastp,2,iblk_off+1:iblk_off+njoined)
c            occ_blk_dag(1:ngastp,2,1:njoined) =
c     &           hpvx_occ(1:ngastp,1,iblk_off+1:iblk_off+njoined)
c            graph_blk_dag(1:ngastp,1,1:njoined) =
c     &           idx_graph(1:ngastp,2,iblk_off+1:iblk_off+njoined)
c            graph_blk_dag(1:ngastp,2,1:njoined) =
c     &           idx_graph(1:ngastp,1,iblk_off+1:iblk_off+njoined)
c            occ_blk_pnt   => occ_blk_dag
c            graph_blk_pnt => graph_blk_dag
c          end if
          njoined_tmp = njoined
          xgrd_pnt => xbuf1
        end if

        ! uses new ncsub/nasub if reordered
        allocate(hpvx_csub(ncsub),hpvx_asub(nasub),
     &           occ_csub(ncsub), occ_asub(nasub),
     &           graph_csub(ncsub), graph_asub(nasub),
     &           msdis_c(ncsub),  msdis_a(nasub),
     &           idxmsdis_c(ncsub),  idxmsdis_a(nasub),
     &           gamdis_c(ncsub), gamdis_a(nasub),
     &           len_str(ncsub+nasub))

        ! the same for reordered and un-reordered case...
        msc_max = ca_occ(1,iblk)
        msa_max = ca_occ(2,iblk)

        ! set HPVX and OCC info
        call condense_occ(occ_csub, occ_asub,
     &                    hpvx_csub,hpvx_asub,
     &                    occ_blk_pnt,njoined_tmp,hpvxblkseq)
        ! do the same for the graph info
        call condense_occ(graph_csub, graph_asub,
     &                    hpvx_csub,hpvx_asub,
     &                    graph_blk_pnt,
     &                                njoined_tmp,hpvxblkseq)

        ! loop over MS(A)
        idxmsa = 0
        msa_loop: do msa = msa_max, -msa_max, -2

          msc = msa + mstotal
          
          if (abs(msc).gt.msc_max) cycle msa_loop
          idxmsa = idxmsa+1

          gama_loop: do gama = 1, ngam
            
            gamc = multd2h(gama,gamtotal)

            if (ntest.ge.100)
     &           write(luout,*) 'current MS(A), MS(C), GAM(A), GAM(C):',
     &                           msa, msc, gama, gamc

            idxdis = 0
            first = .true.
            distr_loop: do

              if (.not.next_msgamdist2(first,
     &             msdis_c,msdis_a,gamdis_c,gamdis_a,
     &             ncsub,nasub,
     &             occ_csub,occ_asub,
     &             msc,msa,gamc,gama,ngam,
     &             ms_fix,fix_success)) exit distr_loop

              first = .false.
c              if(ms_fix.and..not.fix_success)cycle distr_loop
              
              call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncsub)
              call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nasub)

              call set_len_str(len_str,ncsub,nasub,
     &                         graphs,
     &                         graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                         graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                         hpvxseq,.false.)

              if (idxlist(0,len_str,ncsub+nasub,1).gt.0)
     &             cycle distr_loop

              idxdis = idxdis+1

c test -- special insert 
              if (.false..and.njoined.eq.1) then

                if (mode.eq.nidx_p) then
                  idx = idxlist(IPART,hpvx_csub,ncsub,1)
                  if (idx.ne.1) stop '???'
                  idxms_bx  = idxmsdis_c(idx)
                  gam_bx = gamdis_c(idx)
                  ld_bx  = len_str(idx)
                  idx_b = off_bx_gam_ms(gam_bx,idxms_bx)+1
                  len_astr = len_str(2)
                  idx_grd = off_grd_d_gam_ms(idxdis,gama,idxmsa)+1
                  call optc_prc_test(xbuf1(idx_grd),
     &               xbuf2,xbuf3,bmat(idx_b),xmat(idx_b),ld_bx,len_astr)
                  cycle
                else
                  stop 'I cannot do this ...'
                end if
                
              end if
  
              ! Gamma and Ms of B and X
              idx = idxlist(IHOLE,hpvx_csub,ncsub,1)
              if (idx.le.0.and.njoined.eq.2)
     &             call quit(1,'optc_prc_special','no HOLE??')
              if (njoined.eq.1.and..not.ca_reverse)
     &             idx = imltlist(IPART,hpvx_csub,ncsub,1)
              if (njoined.eq.1.and.ca_reverse)
     &             idx = imltlist(IPART,hpvx_asub,nasub,1)
              if (idx.le.0.or.idx.gt.2)
     &             call quit(1,'optc_prc_special','strange')
              if (.not.ca_reverse) then
                idxms_bx  = idxmsdis_c(idx)
                gam_bx = gamdis_c(idx)
                ld_bx  = len_str(idx)
              else
                idxms_bx  = idxmsdis_a(idx)
                gam_bx = gamdis_a(idx)
                ld_bx  = len_str(ncsub+idx)
              end if

              if (ld_bx.le.0) cycle

              ! occupation, Gamma and Ms of C-string
              if (njoined.eq.2) then
                idx = idxlist(IPART,hpvx_csub,ncsub,1)
              else
                if (.not.ca_reverse)
     &               idx = imltlist(IPART,hpvx_csub,ncsub,1)
                if (     ca_reverse)
     &               idx = imltlist(IPART,hpvx_asub,nasub,1)
                if (idx.eq.2) then
                  idx = 1
                else
                  idx = 0
                end if
              end if
              if (idx.le.0) then
                nidx_cstr = 0
                ms_cstr = 0
                gam_cstr = 1
                len_cstr = 1
                graph_cstr = 1 ! dummy
              else
                if (.not.ca_reverse) then
                  nidx_cstr = occ_csub(idx)
                  ms_cstr   = msdis_c(idx)
                  gam_cstr  = gamdis_c(idx)
                  len_cstr = len_str(idx)
                  graph_cstr = graph_csub(idx)
                else
                  nidx_cstr = occ_asub(idx)
                  ms_cstr   = msdis_a(idx)
                  gam_cstr  = gamdis_a(idx)
                  len_cstr = len_str(ncsub+idx)
                  graph_cstr = graph_asub(idx)
                end if
              end if

              if (len_cstr.le.0) cycle

              ! occupation, Gamma and Ms of A-string
              ! idx = 1,  we know that for sure
              if (.not.ca_reverse) then
                nidx_astr = occ_asub(1)
                ms_astr   = msdis_a(1)
                gam_astr  = gamdis_a(1)
                len_astr = len_str(ncsub+1)
                graph_astr = graph_asub(1)
              else
                nidx_astr = occ_csub(1)
                ms_astr   = msdis_c(1)
                gam_astr  = gamdis_c(1)
                len_astr = len_str(1)
                graph_astr = graph_csub(1)
              end if

              if (len_astr.le.0) cycle

              ! offset of GRD
              idx_grd = off_grd_d_gam_ms(idxdis,gama,idxmsa)+1

              ! offset of appropriate block of B (and X)
              idx_b = off_bx_gam_ms(gam_bx,idxms_bx)+1
              idx_x = idx_b

              idx_b = off_bx_gam_ms(gam_bx,idxms_bx)+1
              idx_x = idx_b
              if (.not.beyond_A) idx_x = 1

              if (ld_bx.lt.len_cstr*len_astr) then
                scrbuf => xbuf3
              else
                allocate(scrbuf(ld_bx*ld_bx))
              end if

              call optc_prc_special2_inner
     &             (xgrd_pnt(idx_grd), beyond_A,njoined.eq.1
     &                 .and.nidx_cstr.gt.0, w_shift,
     &              xbuf2,scrbuf,bmat(idx_b),xmat(idx_x),f_dia,
     &              ld_bx,len_cstr, len_astr,
     &              nidx_cstr,ms_cstr,gam_cstr,
     &               igas_restr(1,1,1,1,graph_cstr),
     &               mostnd(1,1,idx_gas(IPART)),ngas_hpv(IPART),
     &              nidx_astr,ms_astr,gam_astr,
     &               igas_restr(1,1,1,1,graph_astr),
     &               mostnd(1,1,idx_gas(IHOLE)),ngas_hpv(IHOLE),
     &              igamorb,ngam,ngas)
              
              if (ld_bx.ge.len_cstr*len_astr) then
                deallocate(scrbuf)
              end if

            end do distr_loop

          end do gama_loop

        end do msa_loop

        if (njoined.eq.1.and.mode.ne.nidx_p) then

          call set_reo_special(2,reo_info,mode)

          ! post-processing
          occ_blk_tmp = occ_blk_reo
          call get_reo_info(ldum1,ldum2,
     &         occ_blk_tmp,occ_blk_reo,
     &         rst_blk_tmp,rst_blk_reo,
     &         2,
     &         reo_info,str_info,orb_info)

          ! we must sort back, as well
          idoff_grd = 0
          call reo_op_wmaps_c(
     &         .false.,xdum,0,
     &         me_grd_reo,me_grd,
     &         .false.,.false.,
     &         1,iblk,
     &         0,idoff_grd,
     &         reo_info,
     &         str_info,strmap_info,orb_info)

          ! clean up
          deallocate(me_grd%fhand%incore,
     &               me_grd_reo%fhand%incore)

          me_grd%fhand%buffered = .false.
          me_grd%fhand%incore => null()
          me_grd%fhand%buffer => null()

          me_grd_reo%fhand%buffered = .false.
          me_grd_reo%fhand%incore => null()
          me_grd_reo%fhand%buffer => null()

          call del_me_list('L_GRD_REO',op_info)
          call del_operator('GRD_REO',op_info)
          call dealloc_reo_info(reo_info)
        end if

        if (ntest.ge.100) then
          write(luout,*) 'gradient vector afterwards:'
          write(luout,*) xbuf1(idx_grd+1:idx_grd+len_grd)
        end if

      end do

      deallocate(f_dia,xmat,bmat)

      return
      end
