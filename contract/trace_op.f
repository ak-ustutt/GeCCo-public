*----------------------------------------------------------------------*
      subroutine trace_op(xfac,casign,
     &     update,xret,type_xret,
     &     me_op,me_trop,me_troptmp,
     &     tra_op, tra_trop,
     &     iblkop,iblktrop,iblktroptmp,
     &     idoffop,idofftrop,
     &     cnt_info,reo_info,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     calculate (partial) trace of operator:
*
*        OPtr(IC,IA) = sum_J OP(IC*JC,IA*JA)
*
*     andreas, may 2008
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
      real(8), intent(in) ::
     &     xfac, casign
      real(8), intent(inout), target ::
     &     xret(1)
      type(contraction_info), target ::
     &     cnt_info
      integer, intent(in) ::
     &     type_xret,
     &     iblkop, iblktrop, iblktroptmp,
     &     idoffop,idofftrop
      logical, intent(in) ::
     &     tra_op, tra_trop
      type(me_list), intent(in) ::
     &     me_op, me_trop, me_troptmp
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(reorder_info), intent(in), target ::
     &     reo_info

      logical ::
     &     bufop, buftrop, 
     &     first1, first2, first3, first4, first5,
     &     reo_trop, nonzero, ms_fix, fix_success
      integer ::
     &     mstop,msttrop,
     &     igamtop,igamttrop,
     &     nc_op, na_op,
     &     nc_trop, na_trop,
     &     nc_troptmp, na_troptmp,
     &     nc_cnt, na_cnt,
     &     nsym, ifree, lenscr, lenblock,
     &     idxst_op, idxst_trop,
     &     idxop, idxtrop,
     &     lenop, lentrop,
     &     mscmx_a, mscmx_c, msc_ac, msc_a, msc_c,
     &     igamc_ac, igamc_a, igamc_c,
     &     idxms, idxdis, lenmap, lblk_troptmp,
     &     idxdis_trop
      integer ::
     &     ncblk_op, nablk_op, 
     &     ncblk_trop, nablk_trop, ncblk_troptmp, nablk_troptmp, 
     &     ncblk_cnt, nablk_cnt,
     &     ncblk_trop_0, nablk_trop_0,
     &     ncblk_reo12,    nablk_reo12
      type(filinf), pointer ::
     &     ffop,ffop2,fftrop
      type(operator), pointer ::
     &     op, op2, trop, troptmp
      integer, pointer ::
     &     cinfo_opc(:,:),cinfo_opa(:,:),
     &     cinfo_tropc(:,:),
     &     cinfo_tropa(:,:),
     &     cinfo_troptmpc(:,:),
     &     cinfo_troptmpa(:,:),
     &     cinfo_cntc(:,:),cinfo_cnta(:,:),
     &     map_info_c(:),
     &     map_info_a(:)

      real(8) ::
     &     xnrm
      real(8) ::
     &     cpu, sys, cpu0, sys0, cpu00, sys00

      real(8), pointer ::
     &     xop(:), xtrop(:), xscr(:)
      real(8), pointer ::
     &     xbf(:), xbftr(:), xbftrtmp(:), xtropblk(:)

      integer ::
     &     msbnd(2,2), igambnd(2,2),
     &     ms_op_tr_a(2), ms_op_tr_c(2),
c dbg mh igam_op_tr_a(2) --> igam_op_tr_a(3)
     &     igam_op_tr_a(3), igam_op_tr_c(2)

      integer, pointer ::
     &     gmopdis_c(:), gmopdis_a(:),
     &     gmc_dis_c (:), gmc_dis_a (:),
     &     gmtropdis_c (:), gmtropdis_a (:),
     &     msopdis_c(:), msopdis_a(:),
     &     msc_dis_c (:), msc_dis_a (:),
     &     mstropdis_c(:), mstropdis_a(:),
     &     idxmsopdis_c(:), idxmsopdis_a(:),
     &     idxmsc_dis_c (:), idxmsc_dis_a (:),
     &     idxmstropdis_c(:), idxmstropdis_a(:),
     &     lstrcnt(:),
     &     lstrop(:),lstrtroptmp(:)

      integer, pointer ::
     &     map_tropcnta(:), map_tropcntc(:)

      integer, pointer ::
     &     ndis_op(:,:), d_gam_ms_op(:,:,:), gam_ms_op(:,:),
     &     ndis_troptmp(:,:), d_gam_ms_trop(:,:,:),
     &     gam_ms_trop(:,:),
     &     len_gam_ms_troptmp(:,:), len_d_gam_ms_troptmp(:,:,:)

      integer, pointer ::
     &     cinfo_reo12c(:,:), cinfo_reo12a(:,:),
     &     cinfo_trop_0c(:,:), cinfo_trop_0a(:,:),
     &     map_info_reo1c(:), map_info_reo1a(:),
     &     map_info_reo2c(:), map_info_reo2a(:)

      type(graph), pointer ::
     &     graphs(:)

      integer, external ::
     &     ielsum, ielprd, idx_msgmdst2, get_lenmap, idxlist,
     &     max_dis_blk
      logical, external ::
     &     next_dist, next_msgamdist2, msa2idxms4op, 
     &     next_msgamdist_diag
      real(8), external ::
     &     ddot

      if (ntest.gt.0) then
        call write_title(luout,wst_dbg_subr,
     &       'trace_op at work')
      end if

      op => me_op%op
      trop => me_trop%op
      troptmp => me_troptmp%op

      ffop => me_op%fhand
      fftrop => me_trop%fhand

      ms_fix=.false.
      if(me_op%fix_vertex_ms.or.me_trop%fix_vertex_ms)then
        ms_fix = me_op%fix_vertex_ms.and.me_trop%fix_vertex_ms
        if(.not.ms_fix) call quit(1,'trace_op',
     &                            'fix ms or not?')
      endif
 
      if (ntest.ge.10) then
        write(luout,*) 'list(in):   ',trim(me_op%label),
     &       ' transp:',tra_op
        write(luout,*) 'list(out):  ',trim(me_trop%label),
     &       ' transp:',tra_trop
        write(luout,*) 'ffop:   ',ffop%name(1:len_trim(ffop%name))
        write(luout,*) 'fftrop:',
     &       fftrop%name(1:len_trim(fftrop%name))
        write(luout,*) 'xfac = ',xfac
        write(luout,*) 'casign = ',casign
        if (type_xret.ne.0)
     &       write(luout,*) 'xret on entry = ',xret(1)
        write(luout,*) 'op: ',trim(op%name),
     &       ' block ',iblkop
        if (iblktrop.gt.0) then
          write(luout,*) 'trop: ',trim(trop%name),
     &       ' block ',iblktrop
        else
          write(luout,*) 'trop: scalar'
        end if
      end if

      ! flag whether non-zero contribution to trop occurred
      nonzero = .false.

      reo_trop = reo_info%n_op_reo.gt.0
      if (reo_trop.and..not.associated(reo_info%map_reo1c))
     &     call quit(1,'contr_trop_wmaps_c',
     &     'reo_info is not consistent')
      if (ntest.ge.10) write(luout,*) 'reo_trop: ',reo_trop
      if (ntest.ge.10 .and. reo_trop) then
        write(luout,*) 'troptmp: ',trim(troptmp%name),
     &       ' block ',iblktroptmp
      end if

      ! set pointers to cnt_info entries
      ncblk_op = cnt_info%ncblk_op1 
      nablk_op = cnt_info%nablk_op1 
      ncblk_cnt = cnt_info%ncblk_cnt  
      nablk_cnt = cnt_info%nablk_cnt  
      ncblk_trop = cnt_info%ncblk_op1op2
      nablk_trop = cnt_info%nablk_op1op2
      ncblk_troptmp = cnt_info%ncblk_op1op2tmp
      nablk_troptmp = cnt_info%nablk_op1op2tmp

      cinfo_opc => cnt_info%cinfo_op1c
      cinfo_opa => cnt_info%cinfo_op1a
      cinfo_tropc => cnt_info%cinfo_op1op2c
      cinfo_tropa => cnt_info%cinfo_op1op2a
      cinfo_troptmpc => cnt_info%cinfo_op1op2tmpc
      cinfo_troptmpa => cnt_info%cinfo_op1op2tmpa
      cinfo_cntc => cnt_info%cinfo_cntc
      cinfo_cnta => cnt_info%cinfo_cnta

      map_info_c => cnt_info%map_info_1c
      map_info_a => cnt_info%map_info_1a

      allocate(
     &       gmopdis_c(ncblk_op), gmopdis_a(nablk_op),
     &       gmtropdis_c(ncblk_trop), gmtropdis_a(nablk_trop),
     &       gmc_dis_c(ncblk_cnt), gmc_dis_a(nablk_cnt),
     &       msopdis_c(ncblk_op), msopdis_a(nablk_op),
     &       mstropdis_c(ncblk_trop), mstropdis_a(nablk_trop),
     &       msc_dis_c(ncblk_cnt), msc_dis_a(nablk_cnt),
     &       idxmsopdis_c(ncblk_op), idxmsopdis_a(nablk_op),
     &       idxmstropdis_c(ncblk_trop), idxmstropdis_a(nablk_trop),
     &       idxmsc_dis_c(ncblk_cnt), idxmsc_dis_a(nablk_cnt),
     &       lstrcnt(ncblk_cnt+nablk_cnt),
     &       lstrop(ncblk_op+nablk_op),
     &       lstrtroptmp(ncblk_troptmp+nablk_troptmp)
     &       )
      
      idxst_op = me_op%off_op_occ(iblkop) + 1
      lenop    = me_op%len_op_occ(iblkop)
      if (iblktrop.gt.0) then
        ! refers to reordered trop (if that makes a difference)
        idxst_trop = me_trop%off_op_occ(iblktrop) + 1
        lentrop    = me_trop%len_op_occ(iblktrop)
      else
        idxst_trop = 1
        lentrop = 1
      end if

      mstop = me_op%mst
      msttrop = me_trop%mst
      igamtop = me_op%gamt
      igamttrop = me_trop%gamt

      if (igamtop.ne.igamttrop)
     &     call quit(1,'trace_op','inconsistent symmetries')

      if (lenop.le.0.or.lentrop.le.0) then
        write(luout,*)
     &       trim(op%name),' ',
     &       trim(trop%name)
        write(luout,*) 'lenop, lentrop: ',
     &                  lenop, lentrop
        call quit(1,'trace_op',
     &     'zero length for operator?')
      end if

      ifree = mem_setmark('trace_op')

      if (ffop%buffered.and.ffop%incore(iblkop).gt.0) then
        bufop = .true.
        xop => ffop%buffer(idxst_op:)
      else
        bufop = .false.
        ! LOWER incore requirements:
        ! check for length of operator
        ! if length < 2 * blocks * da_reclen, try to alloc
        !   full block
        ! else get largest symmetry block
        ifree = mem_alloc_real(xbf,lenop,'xbf')
        xop => xbf
        call get_vec(ffop,xop,idoffop+idxst_op,
     &                          idoffop+idxst_op-1+lenop)
      end if

      if (ntest.ge.100) write(luout,*) ' bufop: ',bufop

      ! get result vector as well (as we update)
      ! refers to reordered trop
      if (iblktrop.gt.0) then
        if (fftrop%buffered.and.fftrop%incore(iblktrop).gt.0) then
          buftrop = .true.
          xtrop => fftrop%buffer(idxst_trop:)
        else
          buftrop = .false.
          ! LOWER incore requirements:
          ! see above
          ifree = mem_alloc_real(xbftr,lentrop,'xbftr')
          xtrop => xbftr
          if (update) then
            ! read from disc
            call get_vec(fftrop,xtrop,idofftrop+idxst_trop,
     &                             idofftrop+idxst_trop-1+lentrop)
          else
            ! init with zero
            xtrop(1:lentrop) = 0d0
          end if
        end if
        if (ntest.ge.100) write(luout,*) ' buftrop: ',buftrop
      else
        buftrop = .true.
        xtrop => xret
        if (ntest.ge.100) write(luout,*) ' result is scalar '
      end if

      if (ntest.ge.1000) then
        ! this will work if all blocks incore, only:
        write(luout,*) 'operator(in) (',trim(op%name),
     &                    ',list=',trim(me_op%label),')'
        call wrt_mel_buf(luout,5,xop,me_op,iblkop,iblkop,
     &                  str_info,orb_info)
        if (iblktrop.gt.0) then
          write(luout,*) 'operator(out) on entry (',trim(trop%name),
     &                                ',list=',trim(me_trop%label),')'

          call wrt_mel_buf(luout,5,xtrop,me_trop,
     &                    iblktrop,iblktrop,
     &                    str_info,orb_info)
        end if
      end if

      if (reo_trop) then
        call quit(1,'trace_op','reo-route is not yet debugged')
        lblk_troptmp =
     &       max_dis_blk(0,me_troptmp,iblktroptmp,orb_info)
        ifree = mem_alloc_real(xbftrtmp,lblk_troptmp,'xbftrtmp')
      end if

      graphs => str_info%g

      ndis_op => me_op%off_op_gmox(iblkop)%ndis
      gam_ms_op => me_op%off_op_gmo(iblkop)%gam_ms
      d_gam_ms_op => me_op%off_op_gmox(iblkop)%d_gam_ms
      ndis_troptmp => me_troptmp%off_op_gmox(iblktroptmp)%ndis
      gam_ms_trop => me_trop%off_op_gmo(iblktrop)%gam_ms
      len_gam_ms_troptmp =>
     &                   me_troptmp%len_op_gmo(iblktroptmp)%gam_ms
      d_gam_ms_trop => me_trop%off_op_gmox(iblktrop)%d_gam_ms
      len_d_gam_ms_troptmp =>
     &                  me_troptmp%len_op_gmox(iblktroptmp)%d_gam_ms

      call sum_occ(nc_op,cinfo_opc,ncblk_op)
      call sum_occ(na_op,cinfo_opa,nablk_op)
      call sum_occ(nc_trop,cinfo_tropc,ncblk_trop)
      call sum_occ(na_trop,cinfo_tropa,nablk_trop)
      call sum_occ(nc_troptmp,cinfo_troptmpc,ncblk_troptmp)
      call sum_occ(na_troptmp,cinfo_troptmpa,nablk_troptmp)
c dbg
      if (na_trop.ne.na_troptmp)
     &     call quit(1,'contr_trop_wmaps_c','unexpected 1a')
      if (nc_trop.ne.nc_troptmp) then
        write(luout,*) 'TROP (C)   : ',nc_trop,
     &       ' <- ',cinfo_tropc(1:ncblk_trop,1)
        write(luout,*) 'TROPTMP (C): ',nc_troptmp,
     &       ' <- ',cinfo_troptmpc(1:ncblk_troptmp,1)
        call quit(1,'contr_trop_wmaps_c','unexpected 1b')
      end if
c dbg
      call sum_occ(nc_cnt,cinfo_cntc,ncblk_cnt)
      call sum_occ(na_cnt,cinfo_cnta,nablk_cnt)

      ! set up maps (if necessary)
      call strmap_man_c(lenmap,
     &     cinfo_cntc(1,2),ncblk_cnt,
     &     cinfo_tropc(1,2),ncblk_trop,
     &     cinfo_opc(1,2),ncblk_op,map_info_c,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_tropcntc,lenmap,'opmap_c')
      call strmap_man_c(lenmap,
     &     cinfo_cnta(1,2),nablk_cnt,
     &     cinfo_tropa(1,2),nablk_trop,
     &     cinfo_opa(1,2),nablk_op,map_info_a,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_tropcnta,lenmap,'opmap_a')

      if (reo_trop) then
        cinfo_reo12c => reo_info%cinfo_reo_c
        cinfo_reo12a => reo_info%cinfo_reo_a
        cinfo_trop_0c => reo_info%cinfo_opreo0c
        cinfo_trop_0a => reo_info%cinfo_opreo0a
        map_info_reo1c  => reo_info%map_reo1c
        map_info_reo1a  => reo_info%map_reo1a
        map_info_reo2c  => reo_info%map_reo2c
        map_info_reo2a  => reo_info%map_reo2a
        ncblk_trop_0  = reo_info%ncblk_reo0
        nablk_trop_0  = reo_info%nablk_reo0
        ncblk_reo12  = reo_info%ncblk_reo
        nablk_reo12  = reo_info%nablk_reo
        call strmap_man_c(lenmap,
     &     cinfo_trop_0c(1,2),ncblk_trop_0,
     &     cinfo_reo12c(1,2),ncblk_reo12,
     &     cinfo_troptmpc(1,2),ncblk_troptmp,map_info_reo1c,
     &     str_info,strmap_info,orb_info)
        call strmap_man_c(lenmap,
     &     cinfo_reo12a(1,2),nablk_reo12,
     &     cinfo_trop_0a(1,2),nablk_trop_0,
     &     cinfo_troptmpa(1,2),nablk_troptmp,map_info_reo1a,
     &     str_info,strmap_info,orb_info)
        call strmap_man_c(lenmap,
     &     cinfo_trop_0c(1,2),ncblk_trop_0,
     &     cinfo_reo12c(1,2),ncblk_reo12,
     &     cinfo_tropc(1,2),ncblk_trop,map_info_reo2c,
     &     str_info,strmap_info,orb_info)
        call strmap_man_c(lenmap,
     &     cinfo_reo12a(1,2),nablk_reo12,
     &     cinfo_trop_0a(1,2),nablk_trop_0,
     &     cinfo_tropa(1,2),nablk_trop,map_info_reo2a,
     &     str_info,strmap_info,orb_info)
      end if

      ! minimum Ms(A) for ...
      msbnd(1,1) = -na_op ! operator 
      msbnd(1,2) = -na_trop ! trace
      ! maximum Ms(A) for ...
      msbnd(2,1) = -msbnd(1,1)
      msbnd(2,2) = -msbnd(1,2)
      ! max |Ms| for ...
      mscmx_a = na_cnt!+nc_cnt ! C(A)
      mscmx_c = nc_cnt!+na_cnt ! C(C)
      ! minimum IRREP for operators
      igambnd(1,1) = 1
      igambnd(1,2) = 1
      ! maximum IRREP
      nsym = orb_info%nsym
      igambnd(2,1) = nsym
      igambnd(2,2) = nsym
      ! loop Ms-cases of (Op(A),TrOp(A))
      first1 = .true.
      ms_loop: do
        if (first1) then
          first1 = .false.
          ! initial Ms distribution
          ms_op_tr_a(1:2) = msbnd(2,1:2)
        else
          ! next Ms distribution
          if (.not.next_dist(ms_op_tr_a,2,msbnd,-2)) exit
        end if

        ms_op_tr_c(1) = ms_op_tr_a(1) + mstop
        ms_op_tr_c(2) = ms_op_tr_a(2) + msttrop
        if (abs(ms_op_tr_c(1)).gt.nc_op) cycle ms_loop
        if (abs(ms_op_tr_c(2)).gt.nc_trop) cycle ms_loop

        msc_c = ms_op_tr_c(1) - ms_op_tr_c(2)
        msc_a = ms_op_tr_a(1) - ms_op_tr_a(2)
        if (msc_c.ne.msc_a)
     &       call quit(1,'trace_op','trap 1')

        if (mscmx_a.lt.abs(msc_a)) cycle ms_loop

        ! loop IRREP cases of (Op(A),TrOp(A))
        first2 = .true.
        gam_loop: do
          if (first2) then
            first2 = .false.
            ! initial IRREP distribution
            igam_op_tr_a(1:2) = igambnd(1,1:2)
          else
            ! next IRREP distribution
            if (.not.next_dist(igam_op_tr_a,2,igambnd,+1)) exit
          end if
          ! set up start addresses
          ! need to be modified, if more than one distribution
          ! exists, see below
          idxms = msa2idxms4op(ms_op_tr_a(1),mstop,na_op,nc_op)
          idxop = gam_ms_op(igam_op_tr_a(1),idxms) + 1
     &             - idxst_op+1
          idxms = msa2idxms4op(ms_op_tr_a(2),msttrop,na_trop,nc_trop)
          ! relevant for case where no reordering necessary
          ! then we have: troptmp == trop
          if (iblktrop.gt.0)
     &           idxtrop = gam_ms_trop(igam_op_tr_a(2),idxms) + 1
     &                - idxst_trop+1
          if (reo_trop)
     &           lblk_troptmp=len_gam_ms_troptmp(igam_op_tr_a(2),idxms)
          if (iblktrop.eq.0) idxtrop = 1

          igam_op_tr_c(1) = multd2h(igam_op_tr_a(1),igamtop)
          igam_op_tr_c(2) = multd2h(igam_op_tr_a(2),igamttrop)

          igamc_a = multd2h(igam_op_tr_a(1),igam_op_tr_a(2))
          igamc_c = multd2h(igam_op_tr_c(1),igam_op_tr_c(2))

          ! loop over distributions of current Ms and IRREP 
          ! of Aex and Cex (=trOP(A),trOP(C)) over ngastypes
          first3 = .true.
          caex_loop: do
            if (.not.next_msgamdist2(first3,
     &          mstropdis_c,mstropdis_a,gmtropdis_c,gmtropdis_a,
     &          ncblk_trop, nablk_trop,
     &          cinfo_tropc,cinfo_tropa,
     &          ms_op_tr_c(2),ms_op_tr_a(2),
     &          igam_op_tr_c(2),igam_op_tr_a(2),nsym,
     &          ms_fix,fix_success))
     &          exit caex_loop
            first3 = .false.
            if(ms_fix.and..not.fix_success)cycle caex_loop

            call ms2idxms(idxmstropdis_c,mstropdis_c,
     &           cinfo_tropc,ncblk_trop)
            call ms2idxms(idxmstropdis_a,mstropdis_a,
     &           cinfo_tropa,nablk_trop)

            call set_len_str(lstrtroptmp,ncblk_trop,nablk_trop,
     &           graphs,
     &           cinfo_tropc(1,2),idxmstropdis_c,
     &             gmtropdis_c,cinfo_tropc(1,3),
     &           cinfo_tropa(1,2),idxmstropdis_a,
     &             gmtropdis_a,cinfo_tropa(1,3),
     &           hpvxseq,.false.)
                
            ! test C and A separately to avoid overflow
            if ( ncblk_trop+nablk_trop.gt.0 .and.
     &               idxlist(0,lstrtroptmp,
     &                          ncblk_trop+nablk_trop,1).gt.0)
     &               cycle

            idxms =  msa2idxms4op(ms_op_tr_a(2),msttrop,na_trop,nc_trop)
            if (iblktroptmp.gt.0.and.
     &           ndis_troptmp(igam_op_tr_a(2),idxms).gt.1) then
              idxdis =
     &             idx_msgmdst2(
     &                   iblktroptmp,idxms,igam_op_tr_a(2),
     &                   cinfo_troptmpc,idxmstropdis_c,
     &                              gmtropdis_c,ncblk_troptmp,
     &                   cinfo_troptmpa,idxmstropdis_a,
     &                              gmtropdis_a,nablk_troptmp,
     &                   tra_trop,me_troptmp,nsym)
              idxdis_trop = idxdis

              ! relevant for case w/o reordering
              ! then we have troptmp == trop
              idxtrop = 
     &             d_gam_ms_trop(idxdis,igam_op_tr_a(3),idxms)+1
     &                   - idxst_trop+1

              if (reo_trop)
     &             lblk_troptmp =
     &             len_d_gam_ms_troptmp(idxdis,igam_op_tr_a(3),idxms)

            end if

            if (.not.reo_trop) then
              ! direct update of result block
              xtropblk => xtrop(idxtrop:)
            else
              ! put result to intermediate buffer
              xtropblk => xbftrtmp
              xbftrtmp(1:lblk_troptmp) = 0d0
            end if

            ! loop over distributions of current Ms and IRREP 
            ! of AC and CC over ngastypes            
            first5 = .true.
            cac_loop: do
              if (.not.next_msgamdist_diag(first5,
     &              msc_dis_a,gmc_dis_a,
     &              nablk_cnt,
     &              cinfo_cnta,
     &              msc_a,igamc_a,nsym))
     &                 exit cac_loop
              first5 = .false.

              call ms2idxms(idxmsc_dis_a,msc_dis_a,
     &              cinfo_cnta,nablk_cnt)

              ! length of contraction
              call set_len_str(lstrcnt,ncblk_cnt,nablk_cnt,
     &              graphs,
     &              cinfo_cntc(1,2),idxmsc_dis_a,
     &                                  gmc_dis_a,cinfo_cntc(1,3),
     &              cinfo_cnta(1,2),idxmsc_dis_a,
     &                                  gmc_dis_a,cinfo_cnta(1,3),
     &              hpvxseq,.false.)

              if ( ncblk_cnt+nablk_cnt.gt.0 .and.
     &                   idxlist(0,lstrcnt,
     &                   ncblk_cnt+nablk_cnt,1).gt.0)
     &                   cycle cac_loop

              ! get Ms and IRREP distribution of op
              call merge_msgmdis(msopdis_c,gmopdis_c,
     &                                 ncblk_op,
     &                                 msc_dis_a,gmc_dis_a,
     &                                 mstropdis_c,gmtropdis_c,
     &                                 map_info_c)
              call merge_msgmdis(msopdis_a,gmopdis_a,
     &                                 nablk_op,
     &                                 msc_dis_a,gmc_dis_a,
     &                                 mstropdis_a,gmtropdis_a,
     &                                 map_info_a)

              call ms2idxms(idxmsopdis_c,msopdis_c,
     &              cinfo_opc,ncblk_op)
              call ms2idxms(idxmsopdis_a,msopdis_a,
     &              cinfo_opa,nablk_op)
                    
              call set_len_str(
     &                   lstrop,ncblk_op,nablk_op,
     &                   graphs,
     &                   cinfo_opc(1,2),idxmsopdis_c,
     &                                    gmopdis_c,cinfo_opc(1,3),
     &                   cinfo_opa(1,2),idxmsopdis_a,
     &                                    gmopdis_a,cinfo_opa(1,3),
     &                   hpvxseq,.false.)

              if ( ncblk_op+nablk_op.gt.0 .and.
     &                   idxlist(0,lstrop,
     &                             ncblk_op+nablk_op,1).gt.0)
     &                   cycle cac_loop

              ! get distribution index
              idxms =  msa2idxms4op(ms_op_tr_a(1),mstop,na_op,nc_op)
              if (ndis_op(igam_op_tr_a(1),idxms).gt.1) then
                idxdis =
     &                idx_msgmdst2(
     &                     iblkop,idxms,igam_op_tr_a(1),
     &                     cinfo_opc,idxmsopdis_c,
     &                              gmopdis_c,ncblk_op,
     &                     cinfo_opa,idxmsopdis_a,
     &                              gmopdis_a,nablk_op,
     &                     tra_op,me_op,nsym)
c     &                     .false.,me_op,nsym)
                idxop =   d_gam_ms_op(idxdis,igam_op_tr_a(1),idxms) + 1
     &                     - idxst_op+1
              end if

              xnrm = ddot(ielprd(lstrop,
     &                   ncblk_op+nablk_op),
     &                   xop(idxop),1,xop(idxop),1)
              if (xnrm.lt.1d-28) cycle cac_loop


              ! if we get here, trop will change:
              nonzero = .true.

              ! get igrphcnt,igrphext1->igrphop map
              ! for given ms and irreps
              lenmap = get_lenmap(lstrcnt,lstrtroptmp,
     &                   map_info_c,ncblk_op)
              call get_strmap_blk_c(map_tropcntc,
     &                   ncblk_cnt,ncblk_trop,ncblk_op,
     &                   cinfo_cntc,cinfo_tropc,lstrcnt,lstrtroptmp,
     &                   cinfo_cntc(1,2),cinfo_tropc(1,2),
     &                   idxmsc_dis_a,idxmstropdis_c,
     &                   gmc_dis_a,gmtropdis_c,map_info_c,
     &                   strmap_info,nsym,str_info%ngraph)

              lenmap = get_lenmap(lstrcnt(ncblk_cnt+1),
     &                                  lstrtroptmp(ncblk_trop+1),
     &                                  map_info_a,nablk_op)
              call get_strmap_blk_c(map_tropcnta,
     &                   nablk_cnt,nablk_trop,nablk_op,
     &                   cinfo_cnta,cinfo_tropa,
     &                     lstrcnt(ncblk_cnt+1),
     &                             lstrtroptmp(ncblk_trop+1),
     &                   cinfo_cnta(1,2),cinfo_tropa(1,2),
     &                   idxmsc_dis_a,idxmstropdis_a,
     &                   gmc_dis_a,gmtropdis_a,map_info_a,
     &                   strmap_info,nsym,str_info%ngraph)

              ! calculate the trace for this block
              call trace_op_blk(xfac*casign,
     &                   xtropblk,xop(idxop),
     &                   tra_op, tra_trop,
     &                   ncblk_op,nablk_op,ncblk_troptmp,nablk_troptmp,
     &                   ncblk_cnt,nablk_cnt,
     &                   cinfo_opc(1,3),cinfo_opa(1,3),
     &                   cinfo_troptmpc(1,3),cinfo_troptmpa(1,3),
     &                   lstrop,lstrtroptmp,
     &                   lstrcnt,
     &                   map_info_c, map_info_a,
     &                   map_tropcntc, map_tropcnta
     &                   )

            end do cac_loop
                  
            ! if necessary, reorder trop block:
            if (reo_trop.and.nonzero) then
              call reo_blk_wmaps_c(xtrop,xtropblk,
     &                  reo_info%sign_reo,
     &                  tra_trop,
     &                  ms_op_tr_c(2),ms_op_tr_a(2),
     &                                  igam_op_tr_c(2),igam_op_tr_a(2),
     &                  mstropdis_c,mstropdis_a,gmtropdis_c,gmtropdis_a,
     &                  ncblk_troptmp,nablk_troptmp,
     &                  cinfo_troptmpc,cinfo_troptmpa,
     &                  lstrtroptmp,
     &                  me_trop,iblktrop,
     &                  ncblk_trop,nablk_trop,
     &                  cinfo_tropc,cinfo_tropa,
     &                  reo_info%ncblk_reo,reo_info%nablk_reo,
     &                  reo_info%cinfo_reo_c,reo_info%cinfo_reo_a,
     &                  reo_info%ncblk_reo0,reo_info%nablk_reo0,
     &                  reo_info%cinfo_opreo0c,reo_info%cinfo_opreo0a,
     &                  reo_info%map_reo1c,reo_info%map_reo1a,
     &                  reo_info%map_reo2c,reo_info%map_reo2a,
     &                  nsym,str_info,strmap_info)
            end if

          end do caex_loop
              
        end do gam_loop
      end do ms_loop

      if (ntest.ge.1000) then
        if (iblktrop.gt.0
     &       ) then
          write(luout,*) 'operator(out) on exit (',trim(trop%name),')'
          call wrt_mel_buf(luout,5,xtrop,me_trop,
     &         iblktrop,iblktrop,str_info,orb_info)
        end if
      end if

      if (type_xret.eq.2) then
        xret(1) = xtrop(1)
      else if (type_xret.eq.1) then
        xret(1) = ddot(lentrop,xtrop,1,xtrop,1)
      end if

      ! put result to disc
      if (.not.buftrop) then
        call put_vec(fftrop,xtrop,idofftrop+idxst_trop,
     &                    idofftrop+idxst_trop-1+lentrop)
      end if

      deallocate(
     &     gmopdis_c, gmopdis_a,
     &     gmtropdis_c, gmtropdis_a,
     &     gmc_dis_c , gmc_dis_a ,
     &     msopdis_c, msopdis_a,
     &     mstropdis_c, mstropdis_a,
     &     msc_dis_c , msc_dis_a ,
     &     idxmsopdis_c, idxmsopdis_a,
     &     idxmstropdis_c, idxmstropdis_a,
     &     idxmsc_dis_c , idxmsc_dis_a ,
     &     lstrcnt,
     &     lstrop,lstrtroptmp
     &     )

      ifree = mem_flushmark()

      if (ntest.ge.100) then
        if (type_xret.ne.0)
     &       write(luout,*) 'xret on exit = ',xret(1)
      end if

      return
      end


