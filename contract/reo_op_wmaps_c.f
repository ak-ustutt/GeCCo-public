*----------------------------------------------------------------------*
      subroutine reo_op_wmaps_c(fac,
     &     update,xret,type_xret,
     &     me_opori,me_opreo,
     &     tra_opori, tra_opreo,
     &     iblkopori,iblkopreo,
     &     idoffopori,idoffopreo,
     &     reo_info,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     initial driver routine for reordering list ME(Opori) 
*     to list ME(Opreo)
*
*     andreas, june 2008
*
*----------------------------------------------------------------------*
      implicit none

      include 'routes.h'
      include 'contr_times.h'

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'multd2h.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_filinf.h'
      include 'def_strmapinf.h'
      include 'def_reorder_info.h'
      include 'def_contraction_info.h'
      include 'ifc_memman.h'
      include 'ifc_operators.h'
      include 'hpvxseq.h'

      integer, parameter ::
     &     ntest = 000

      logical, intent(in) ::
     &     update
      real(8), intent(inout), target ::
     &     xret(1)
      real(8), intent(in) ::
     &     fac
      integer, intent(in) ::
     &     type_xret,
     &     iblkopori, iblkopreo, 
     &     idoffopori,idoffopreo
      logical, intent(in) ::
     &     tra_opori, tra_opreo
      type(me_list), intent(in) ::
     &     me_opori, me_opreo
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(reorder_info), intent(in), target ::
     &     reo_info

      logical ::
     &     bufopori, bufopreo, first,
     &     nonzero, ms_fix, fix_success, ltmp
      integer ::
     &     njoined_opori, njoined_opreo, ngam,
     &     mstopori,mstopreo,
     &     igamtopori,igamtopreo,
     &     nc_opori, na_opori, nc_opreo, na_opreo,
     &     ifree, lenscr, lenblock,
     &     idxst_opori, idxst_opreo,
     &     idxopori, idxopreo,
     &     lenopori, lenopreo,
     &     msa_max, msc_max, msa, msc, igama, igamc,
     &     idxms, idxms_tra, idxdis, idxdis_tra, lenmap
      integer ::
     &     ncblk_opori, nablk_opori, 
     &     ncblk_opreo, nablk_opreo, 
     &     ncblk_oporiopreo_0, nablk_oporiopreo_0,
     &     ncblk_reo12,    nablk_reo12
      type(filinf), pointer ::
     &     ffopori,ffopreo
      type(operator), pointer ::
     &     opori, opreo
      integer, pointer ::
     &     cinfo_oporic(:,:),cinfo_oporia(:,:),
     &     cinfo_opreoc(:,:),cinfo_opreoa(:,:),
     &     occ_opori_blk1(:,:,:), graph_opori_blk1(:,:,:),
     &     occ_opreo_blk2(:,:,:), graph_opreo_blk2(:,:,:),
     &     dis_map_ca(:), dis_map_ac(:)

      real(8) ::
     &     xnrm
      real(8) ::
     &     cpu, sys, cpu0, sys0, cpu00, sys00

      real(8), pointer ::
     &     xopori(:), xopreo(:)
      real(8), pointer ::
     &     xbf1(:), xbf2(:)


      integer, pointer ::
     &     gm_dis_c (:), gm_dis_a (:),
     &     ms_dis_c (:), ms_dis_a (:),
     &     idxms_dis_c (:), idxms_dis_a (:),
     &     lstropori(:)

      integer, pointer ::
     &     ndis_opori(:,:), d_gam_ms_opori(:,:,:), gam_ms_opori(:,:),
     &     ndis_opreo(:,:), d_gam_ms_opreo(:,:,:), gam_ms_opreo(:,:)

      integer, pointer ::
     &     cinfo_reo12c(:,:), cinfo_reo12a(:,:),
     &     cinfo_op0c(:,:), cinfo_op0a(:,:),
     &     map_info_reo1c(:), map_info_reo1a(:),
     &     map_info_reo2c(:), map_info_reo2a(:)

      type(graph), pointer ::
     &     graphs(:)

      integer, external ::
     &     idxlist, max_dis_blk, idx_msgmdst2, msa2idxms4op
      logical, external ::
     &     next_msgamdist2
      real(8), external ::
     &     ddot

      if (ntest.gt.0) then
        call write_title(lulog,wst_dbg_subr,
     &       'reo_op_wmaps_c at work')
      end if

      opori => me_opori%op
      opreo => me_opreo%op

      ffopori => me_opori%fhand
      ffopreo => me_opreo%fhand

      ms_fix = .false.
      if(me_opori%fix_vertex_ms.or.me_opreo%fix_vertex_ms)then
        ms_fix = me_opori%fix_vertex_ms.and.me_opreo%fix_vertex_ms
        if(.not.ms_fix) call quit(1,'reo_op_wmaps_c',
     &                            'fix ms or not?')
      endif

      ngam = orb_info%nsym

      if (ntest.ge.10) then
        write(lulog,*) 'list1:   ',trim(me_opori%label),
     &       ' transp:',tra_opori
        write(lulog,*) 'list2:   ',trim(me_opreo%label),' transp:',
     &       tra_opreo
        write(lulog,*) 'ffopori:   ',
     &       ffopori%name(1:len_trim(ffopori%name))
        write(lulog,*) 'ffopreo:   ',
     &       ffopreo%name(1:len_trim(ffopreo%name))
        if (type_xret.ne.0)
     &       write(lulog,*) 'xret on entry = ',xret(1)
        write(lulog,*) 'opori: ',trim(opori%name),
     &       ' block ',iblkopori
        write(lulog,*) 'opreo: ',trim(opreo%name),
     &       ' block ',iblkopreo
      end if

      if (tra_opreo) call quit(1,'reo_op_wmaps_c',
     &     'savety trap for tra_opreo==.true. (does it work?)')
      ! in fact, a few things in connection with complicated operators
      !  (including a dis_map) need to be fixed

      ! flag whether any block block was non-zero
      nonzero = .true. ! not used currently

      njoined_opori  = opori%njoined
      njoined_opreo  = opreo%njoined

      idxst_opori = me_opori%off_op_occ(iblkopori) + 1
      lenopori    = me_opori%len_op_occ(iblkopori)
      idxst_opreo = me_opreo%off_op_occ(iblkopreo) + 1
      lenopreo    = me_opreo%len_op_occ(iblkopreo)

      mstopori = me_opori%mst
      mstopreo = me_opreo%mst
      igamtopori = me_opori%gamt
      igamtopreo = me_opreo%gamt

      if (mstopori.ne.mstopreo)
     &     call quit(1,'reo_op_wmaps_c','inconsistent MS')
      if (igamtopori.ne.igamtopreo)
     &     call quit(1,'reo_op_wmaps_c','inconsistent symmetries')

      if (me_opori%op%formal_blk(iblkopori).or.
     &    me_opreo%op%formal_blk(iblkopreo)) then
        write(lulog,*) me_opori%op%formal_blk(iblkopori),
     &                 me_opreo%op%formal_blk(iblkopreo)
        write(lulog,*) 'opori: ',trim(opori%name),
     &       ' block ',iblkopori
        write(lulog,*) 'opreo: ',trim(opreo%name),
     &       ' block ',iblkopreo
        call quit(1,'reo_op_wmaps_c','called for formal block')
      end if

! just accept that in some strange cases this happens not being an error
!      if (lenopori.le.0.or.lenopreo.le.0) then
!        write(lulog,*)
!     &       trim(opori%name),' ',
!     &       trim(opreo%name),' '
!        write(lulog,*) 'lenopori, lenopreo: ',
!     &                  lenopori, lenopreo
!        call quit(1,'reo_oporiopreo_wmaps_c',
!     &     'zero length for operator?')
!      end if

      ifree = mem_setmark('reo1')

      ltmp = ffopori%buffered
      if (ltmp) ltmp = ffopori%incore(iblkopori).ge.0
      if (ltmp) then
        bufopori = .true.
        xopori => ffopori%buffer(idxst_opori:)
      else
        bufopori = .false.
        ! LOWER incore requirements:
        ! check for length of operator
        ! if length < 2 * blocks * da_reclen, try to alloc
        !   full block
        ! else get largest symmetry block
        ifree = mem_alloc_real(xbf1,lenopori,'xbf1')
        xopori => xbf1
        call get_vec(ffopori,xopori,idoffopori+idxst_opori,
     &                          idoffopori+idxst_opori-1+lenopori)
      end if
      ltmp = ffopreo%buffered
      if (ltmp) ltmp = ffopreo%incore(iblkopreo).ge.0
      if (ltmp) then
        bufopreo = .true.
        xopreo => ffopreo%buffer(idxst_opreo:)
      else
        bufopreo = .false.
        ! LOWER incore requirements:
        ! see above
        ifree = mem_alloc_real(xbf2,lenopreo,'xbf2')
        xopreo => xbf2
        if (update)
     &       call get_vec(ffopreo,xopreo,idoffopreo+idxst_opreo,
     &                          idoffopreo+idxst_opreo-1+lenopreo)
      end if

      if (.not.update)
     &     xopreo(1:lenopreo) = 0d0

      if (ntest.ge.100) write(lulog,*) ' bufopori/2: ',bufopori,bufopreo

      if (ntest.ge.1000) then
        ! this will work if all blocks incore, only:
        write(lulog,*) 'operator 1 (',trim(opori%name),
     &                    ',list=',trim(me_opori%label),')'
        call wrt_mel_buf(lulog,5,xopori,me_opori,iblkopori,iblkopori,
     &                  str_info,orb_info)
        if (update) then
          write(lulog,*) 'operator 2 (',trim(opreo%name),
     &                    ',list=',trim(me_opreo%label),')'
          call wrt_mel_buf(lulog,5,xopreo,me_opreo,iblkopreo,iblkopreo,
     &                  str_info,orb_info)
        end if
      end if

      if (.not.tra_opori) then
        occ_opori_blk1 => opori%ihpvca_occ(1:ngastp,1:2,
     &                       njoined_opori*(iblkopori-1)+1:
     &                       njoined_opori*(iblkopori-1)+njoined_opori)
        graph_opori_blk1 => me_opori%idx_graph(1:ngastp,1:2,
     &                       njoined_opori*(iblkopori-1)+1:
     &                       njoined_opori*(iblkopori-1)+njoined_opori)
      else
        allocate(occ_opori_blk1(ngastp,2,njoined_opori),
     &         graph_opori_blk1(ngastp,2,njoined_opori) )
        occ_opori_blk1 = iocc_dagger_n(
     &                       opori%ihpvca_occ(1:ngastp,1:2,
     &                       njoined_opori*(iblkopori-1)+1:
     &                       njoined_opori*(iblkopori-1)+njoined_opori),
     &                          njoined_opori)
        graph_opori_blk1 = iocc_dagger_n(
     &                     me_opori%idx_graph(1:ngastp,1:2,
     &                       njoined_opori*(iblkopori-1)+1:
     &                       njoined_opori*(iblkopori-1)+njoined_opori),
     &                          njoined_opori)
      end if

      if (.not.tra_opreo) then
        occ_opreo_blk2 => opreo%ihpvca_occ(1:ngastp,1:2,
     &                       njoined_opreo*(iblkopreo-1)+1:
     &                       njoined_opreo*(iblkopreo-1)+njoined_opreo)
        graph_opreo_blk2 => me_opreo%idx_graph(1:ngastp,1:2,
     &                       njoined_opreo*(iblkopreo-1)+1:
     &                       njoined_opreo*(iblkopreo-1)+njoined_opreo)
      else
        allocate(occ_opreo_blk2(ngastp,2,njoined_opreo),
     &         graph_opreo_blk2(ngastp,2,njoined_opreo) )
        occ_opreo_blk2 = iocc_dagger_n(
     &                     opreo%ihpvca_occ(1:ngastp,1:2,
     &                       njoined_opreo*(iblkopreo-1)+1:
     &                       njoined_opreo*(iblkopreo-1)+njoined_opreo),
     &                          njoined_opreo)
        graph_opreo_blk2 = iocc_dagger_n(
     &                     me_opreo%idx_graph(1:ngastp,1:2,
     &                       njoined_opreo*(iblkopreo-1)+1:
     &                       njoined_opreo*(iblkopreo-1)+njoined_opreo),
     &                          njoined_opreo)
      end if

      call get_num_subblk(ncblk_opori,nablk_opori,
     &     occ_opori_blk1,njoined_opori)
      call get_num_subblk(ncblk_opreo,nablk_opreo,
     &     occ_opreo_blk2,njoined_opreo)

      allocate(
     &       gm_dis_c(ncblk_opori), gm_dis_a(nablk_opori),
     &       ms_dis_c(ncblk_opori), ms_dis_a(nablk_opori),
     &       idxms_dis_c(ncblk_opori), idxms_dis_a(nablk_opori),
     &       lstropori(ncblk_opori+nablk_opori),
     &       dis_map_ca(ncblk_opori), dis_map_ac(nablk_opori))
      
      graphs => str_info%g

      ndis_opori => me_opori%off_op_gmox(iblkopori)%ndis
      gam_ms_opori => me_opori%off_op_gmo(iblkopori)%gam_ms
      d_gam_ms_opori => me_opori%off_op_gmox(iblkopori)%d_gam_ms
      ndis_opreo => me_opreo%off_op_gmox(iblkopreo)%ndis
      gam_ms_opreo => me_opreo%off_op_gmo(iblkopreo)%gam_ms
      d_gam_ms_opreo => me_opreo%off_op_gmox(iblkopreo)%d_gam_ms

      allocate(cinfo_oporic(ncblk_opori,3),cinfo_oporia(nablk_opori,3),
     &         cinfo_opreoc(ncblk_opreo,3),cinfo_opreoa(nablk_opreo,3))

      ! set HPVX and OCC info
      call condense_occ(cinfo_oporic(1,1), cinfo_oporia(1,1),
     &                  cinfo_oporic(1,3), cinfo_oporia(1,3),
     &                  occ_opori_blk1,njoined_opori,hpvxblkseq)
      ! do the same for the graph info
      call condense_occ(cinfo_oporic(1,2), cinfo_oporia(1,2),
     &                  cinfo_oporic(1,3), cinfo_oporia(1,3),
     &                  graph_opori_blk1,njoined_opori,hpvxblkseq)

      call set_dis_tra_map(dis_map_ca,dis_map_ac,
     &     cinfo_oporic(1,3),cinfo_oporia(1,3),ncblk_opori,nablk_opori)

      ! and the same for OPREO
      call condense_occ(cinfo_opreoc(1,1), cinfo_opreoa(1,1),
     &                  cinfo_opreoc(1,3), cinfo_opreoa(1,3),
     &                  occ_opreo_blk2,njoined_opreo,hpvxblkseq)
      ! ...
      call condense_occ(cinfo_opreoc(1,2), cinfo_opreoa(1,2),
     &                  cinfo_opreoc(1,3), cinfo_opreoa(1,3),
     &                  graph_opreo_blk2,njoined_opreo,hpvxblkseq)

      call sum_occ(nc_opori,cinfo_oporic,ncblk_opori)
      call sum_occ(na_opori,cinfo_oporia,nablk_opori)
      call sum_occ(nc_opreo,cinfo_opreoc,ncblk_opreo)
      call sum_occ(na_opreo,cinfo_opreoa,nablk_opreo)

      if (nc_opori.ne.nc_opreo .or. na_opori.ne.na_opreo) then
        write(lulog,*) 'NC (ORI/REO): ',nc_opori,nc_opreo
        write(lulog,*) 'NA (ORI/REO): ',na_opori,na_opreo
        call quit(1,'reo_op_wmaps_c','inconsistent particle numbers')
      end if

      msa_max = na_opori
      msc_max = nc_opori
      
      cinfo_reo12c => reo_info%cinfo_reo_c
      cinfo_reo12a => reo_info%cinfo_reo_a
      cinfo_op0c => reo_info%cinfo_opreo0c
      cinfo_op0a => reo_info%cinfo_opreo0a
      map_info_reo1c  => reo_info%map_reo1c
      map_info_reo1a  => reo_info%map_reo1a
      map_info_reo2c  => reo_info%map_reo2c
      map_info_reo2a  => reo_info%map_reo2a
      ncblk_oporiopreo_0  = reo_info%ncblk_reo0
      nablk_oporiopreo_0  = reo_info%nablk_reo0
      ncblk_reo12  = reo_info%ncblk_reo
      nablk_reo12  = reo_info%nablk_reo
      call strmap_man_c(2,lenmap,
     &     cinfo_op0c(1,2),ncblk_oporiopreo_0,
     &     cinfo_reo12c(1,2),ncblk_reo12,
     &     cinfo_oporic(1,2),ncblk_opori,map_info_reo1c,
     &     str_info,strmap_info,orb_info)
      call strmap_man_c(2,lenmap,
     &     cinfo_reo12a(1,2),nablk_reo12,
     &     cinfo_op0a(1,2),nablk_oporiopreo_0,
     &     cinfo_oporia(1,2),nablk_opori,map_info_reo1a,
     &     str_info,strmap_info,orb_info)
      call strmap_man_c(2,lenmap,
     &     cinfo_op0c(1,2),ncblk_oporiopreo_0,
     &     cinfo_reo12c(1,2),ncblk_reo12,
     &     cinfo_opreoc(1,2),ncblk_opreo,map_info_reo2c,
     &     str_info,strmap_info,orb_info)
      call strmap_man_c(2,lenmap,
     &     cinfo_reo12a(1,2),nablk_reo12,
     &     cinfo_op0a(1,2),nablk_oporiopreo_0,
     &     cinfo_opreoa(1,2),nablk_opreo,map_info_reo2a,
     &     str_info,strmap_info,orb_info)

      idxms = 0
      ms_loop: do msa = msa_max, -msa_max, -2

        msc = msa + mstopori
        if (abs(msc).gt.msc_max) cycle ms_loop

        idxms = idxms+1

        gam_loop: do igama = 1, ngam
          
          igamc = multd2h(igama,igamtopori)

          if (.not.tra_opori) then
            idxopori = gam_ms_opori(igama,idxms) + 1
     &           - idxst_opori+1
          else
            idxms_tra = msa2idxms4op(msc,mstopori,msa_max,msc_max)
            idxopori = gam_ms_opori(igamc,idxms_tra) + 1
     &           - idxst_opori+1
          end if
          if (.not.tra_opori) then
            idxopreo = gam_ms_opreo(igama,idxms) + 1
     &             - idxst_opreo+1
          else
            idxms_tra = msa2idxms4op(msc,mstopreo,msa_max,msc_max)
            idxopreo = gam_ms_opreo(igamc,idxms_tra) + 1
     &             - idxst_opreo+1            
          end if

          ! LOWER incore requirements:
          !  load here Opori, Opreo, OporiOpreo block

          idxdis = 0            

          first = .true.
          distr_loop: do
                  
            if (.not.next_msgamdist2(first,
     &            ms_dis_c,ms_dis_a,gm_dis_c,gm_dis_a,
     &            ncblk_opori, nablk_opori,
     &            cinfo_oporic,cinfo_oporia,
     &            msc,msa,igamc,igama,ngam,
     &            ms_fix,fix_success))
     &           exit distr_loop
            first = .false.                            
            if(ms_fix.and..not.fix_success)cycle distr_loop

c dbg
c            print *,'igama, igamc: ',igama, igamc
c            print *,'distr: '
c            print *,'GM(C): ',gm_dis_c(1:ncblk_opori)
c            print *,'GM(A): ',gm_dis_a(1:nablk_opori)
c dbg

            call ms2idxms(idxms_dis_c,ms_dis_c,cinfo_oporic,ncblk_opori)
            call ms2idxms(idxms_dis_a,ms_dis_a,cinfo_oporia,nablk_opori)

            call set_len_str(lstropori,ncblk_opori,nablk_opori,
     &                  graphs,
     &                  cinfo_oporic(1,2),idxms_dis_c,
     &                                 gm_dis_c,cinfo_oporic(1,3),
     &                  cinfo_oporia(1,2),idxms_dis_a,
     &                                 gm_dis_a,cinfo_oporia(1,3),
     &                  hpvxseq,.false.)

            if ( ncblk_opori+nablk_opori.gt.0 .and.
     &             idxlist(0,lstropori,
     &                          ncblk_opori+nablk_opori,1).gt.0)
     &             cycle
            idxdis = idxdis + 1
              
            if (.not.tra_opori) then
              idxopori = d_gam_ms_opori(idxdis,igama,idxms)+1
     &                   - idxst_opori+1
            else
              ! if we transpose, the sequence is different
              ! should be changed in order to pass sequentially
              ! through ORI list
              idxdis_tra =
     &             idx_msgmdst2(.true.,
     &                   iblkopori,idxms,igama,
     &                   cinfo_oporic,idxms_dis_c,
     &                              gm_dis_c,ncblk_opori,
     &                   cinfo_oporia,idxms_dis_a,
     &                              gm_dis_a,nablk_opori,
     &                   tra_opori,dis_map_ca,dis_map_ac,me_opori,ngam)
              idxopori = d_gam_ms_opori(idxdis_tra,igama,idxms)+1
     &                   - idxst_opori+1
            end if

c dbg
c          write(lulog,*) 'input block '
c          write(lulog,'(x,5g15.8)')    xbf12tmp(1:lblk_oporiopreotmp)
c          call wrt_mel_buf(lulog,5,xoporiopreo,me_oporiopreo,
c     &         iblkoporiopreo,iblkoporiopreo,str_info,orb_info)
c dbg
            call reo_blk_wmaps_c(fac,xopreo,xopori(idxopori),
     &                   lenopreo,lenopori-idxopori+1, ! just for checks
     &                   reo_info%sign_reo,
     &                   tra_opreo,tra_opori,
     &                   msc,msa,igamc,igama,
     &                   ms_dis_c,ms_dis_a,gm_dis_c,gm_dis_a,
     &                   ncblk_opori,nablk_opori,
     &                   cinfo_oporic,cinfo_oporia,
     &                   lstropori,
     &                   me_opreo,iblkopreo,
     &                   ncblk_opreo,nablk_opreo,
     &                   cinfo_opreoc,cinfo_opreoa,
     &                   reo_info%ncblk_reo,reo_info%nablk_reo,
     &                   reo_info%cinfo_reo_c,reo_info%cinfo_reo_a,
     &                   reo_info%ncblk_reo0,reo_info%nablk_reo0,
     &                   reo_info%cinfo_opreo0c,reo_info%cinfo_opreo0a,
     &                   reo_info%map_reo1c,reo_info%map_reo1a,
     &                   reo_info%map_reo2c,reo_info%map_reo2a,
     &                   ngam,str_info,strmap_info)
c dbg
c          write(lulog,*) 'reordered operator (',trim(oporiopreo%name),')'
c          call wrt_mel_buf(lulog,5,xoporiopreo,me_oporiopreo,
c     &         iblkoporiopreo,iblkoporiopreo,str_info,orb_info)
c dbg
            
          end do distr_loop

        end do gam_loop
      end do ms_loop
          
      if (ntest.ge.1000) then
        write(lulog,*) 'operator 2 on exit (',trim(opreo%name),')'
          call wrt_mel_buf(lulog,5,xopreo,me_opreo,
     &         iblkopreo,iblkopreo,str_info,orb_info)
      end if

      if (type_xret.eq.2) then
        xret(1) = xopreo(1)
      else if (type_xret.eq.1) then
        xret(1) = ddot(lenopreo,xopreo,1,xopreo,1)
      end if

      ! put result to disc
      if (.not.bufopreo) then
        call put_vec(ffopreo,xopreo,idoffopreo+idxst_opreo,
     &                    idoffopreo+idxst_opreo-1+lenopreo)
      end if

      deallocate(
     &     gm_dis_c , gm_dis_a ,
     &     ms_dis_c , ms_dis_a ,
     &     idxms_dis_c, idxms_dis_a,
     &     lstropori,
     &     cinfo_oporic, cinfo_oporia,
     &     cinfo_opreoa, cinfo_opreoc,
     &     dis_map_ca, dis_map_ac
     &     )

      if (tra_opori) deallocate(occ_opori_blk1,graph_opori_blk1)
      if (tra_opreo) deallocate(occ_opreo_blk2,graph_opreo_blk2)

      ifree = mem_flushmark()

      if (ntest.ge.100) then
        if (type_xret.ne.0)
     &       write(lulog,*) 'xret on exit = ',xret(1)
      end if

      return
      end


