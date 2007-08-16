*----------------------------------------------------------------------*
      subroutine contr_op1op2_wmaps_c(xfac,casign,ffop1,ffop2,
     &     update,ffop1op2,xret,type_xret,
     &     op1,op2,op1op2,
     &     iblkop1,iblkop2,iblkop1op2,
     &     idoffop1,idoffop2,idoffop1op2,
     &       nca_blk,
     &       cinfo_op1c, cinfo_op1a, cinfo_op2c, cinfo_op2a,
     &       cinfo_op1op2c, cinfo_op1op2a,
     &       cinfo_ex1c, cinfo_ex1a, cinfo_ex2c, cinfo_ex2a,
     &       cinfo_cntc, cinfo_cnta,
     &       map_info_1c, map_info_1a,
     &       map_info_2c, map_info_2a,
     &       map_info_12c, map_info_12a,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     slightly improved version:
*       full blocks incore
*       addressing via precalculated maps
*       mainly written for testing maps
*
*      update: update matrix elements on ffop1op2
*              else overwrite
*
*     andreas, feb 2007
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
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_filinf.h'
      include 'def_strmapinf.h'
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
      integer, intent(in) ::
     &     type_xret,
     &     iblkop1, iblkop2, iblkop1op2,
     &     idoffop1,idoffop2,idoffop1op2,
     &     nca_blk(2,6),
     &     cinfo_op1c(nca_blk(1,1),3), cinfo_op1a(nca_blk(2,1),3),
     &     cinfo_op2c(nca_blk(1,2),3), cinfo_op2a(nca_blk(2,2),3),
     &     cinfo_op1op2c(nca_blk(1,3),3), cinfo_op1op2a(nca_blk(2,3),3),
     &     cinfo_ex1c(nca_blk(1,4),3), cinfo_ex1a(nca_blk(2,4),3),
     &     cinfo_ex2c(nca_blk(1,5),3), cinfo_ex2a(nca_blk(2,5),3),
     &     cinfo_cntc(nca_blk(1,6),3), cinfo_cnta(nca_blk(2,6),3),
     &     map_info_1c(*), map_info_1a(*),
     &     map_info_2c(*), map_info_2a(*),
     &     map_info_12c(*), map_info_12a(*),
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2
      type(filinf), intent(inout) ::
     &     ffop1,ffop2,ffop1op2
      type(operator), intent(in) ::
     &     op1, op2, op1op2
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      logical ::
     &     bufop1, bufop2, bufop1op2,
     &     first1, first2, first3, first4, first5
      integer ::
     &     nc_op1, na_op1, nc_op2, na_op2,
     &     nc_ex1, na_ex1, nc_ex2, na_ex2, 
     &     nc_op1op2, na_op1op2, nc_cnt, na_cnt,
     &     nsym, ifree,
     &     idxst_op1, idxst_op2, idxst_op1op2,
     &     idxop1, idxop2, idxop1op2,
     &     lenop1, lenop2, lenop1op2,
     &     mscmx_a, mscmx_c, msc_ac, msc_a, msc_c,
     &     msex1_a, msex1_c, msex2_a, msex2_c,
     &     igamc_ac, igamc_a, igamc_c,
     &     igamex1_a, igamex1_c, igamex2_a, igamex2_c,
     &     idxms, idxdis, lenmap
c     &     isignsum, , hpvx, hpvx2, ica
      real(8) ::
     &     xnrm
      real(8) ::
     &     cpu, sys, cpu0, sys0, cpu00, sys00

c      integer ::
c     &     iocc_op1(ngastp,2), iocc_op2(ngastp,2), iocc_op1op2(ngastp,2)

      real(8), pointer ::
     &     xop1(:), xop2(:), xop1op2(:)
      real(8), pointer ::
     &     xbf1(:), xbf2(:), xbf12(:)

      integer ::
c     &     irst_ext1(8*orb_info%ngas),irst_ext2(8*orb_info%ngas),
c     &     irst_cnt(8*orb_info%ngas),
     &     msbnd(2,3), igambnd(2,3),
     &     ms12i_a(3), ms12i_c(3), igam12i_a(3), igam12i_c(3),
     &     gmop1dis_c(nca_blk(1,1)), gmop1dis_a(nca_blk(2,1)),
     &     gmop2dis_c(nca_blk(1,2)), gmop2dis_a(nca_blk(2,2)),
     &     gmex1dis_c(nca_blk(1,4)), gmex1dis_a(nca_blk(2,4)),
     &     gmex2dis_c(nca_blk(1,5)), gmex2dis_a(nca_blk(2,5)),
     &     gmc_dis_c(nca_blk(1,6)), gmc_dis_a(nca_blk(2,6)),
     &     gmi_dis_c(nca_blk(1,3)), gmi_dis_a(nca_blk(2,3)),
     &     msop1dis_c(nca_blk(1,1)), msop1dis_a(nca_blk(2,1)),
     &     msop2dis_c(nca_blk(1,2)), msop2dis_a(nca_blk(2,2)),
     &     msex1dis_c(nca_blk(1,4)), msex1dis_a(nca_blk(2,4)),
     &     msex2dis_c(nca_blk(1,5)), msex2dis_a(nca_blk(2,5)),
     &     msc_dis_c(nca_blk(1,6)), msc_dis_a(nca_blk(2,6)),
     &     msi_dis_c(nca_blk(1,3)), msi_dis_a(nca_blk(2,3)),
     &     idxmsop1dis_c(nca_blk(1,1)), idxmsop1dis_a(nca_blk(2,1)),
     &     idxmsop2dis_c(nca_blk(1,2)), idxmsop2dis_a(nca_blk(2,2)),
     &     idxmsex1dis_c(nca_blk(1,4)), idxmsex1dis_a(nca_blk(2,4)),
     &     idxmsex2dis_c(nca_blk(1,5)), idxmsex2dis_a(nca_blk(2,5)),
     &     idxmsc_dis_c(nca_blk(1,6)), idxmsc_dis_a(nca_blk(2,6)),
     &     idxmsi_dis_c(nca_blk(1,3)), idxmsi_dis_a(nca_blk(2,3)),
     &     lstrex1(nca_blk(1,4)+nca_blk(2,4)),
     &     lstrex2(nca_blk(1,5)+nca_blk(2,5)),
     &     lstrcnt(nca_blk(1,6)+nca_blk(2,6)),
     &     lstrop1(nca_blk(1,1)+nca_blk(2,1)),
     &     lstrop2(nca_blk(1,2)+nca_blk(2,2)),
     &     lstrop1op2(nca_blk(1,3)+nca_blk(2,3))

      integer, pointer ::
     &     map_ex1ex2a(:), map_ex1ex2c(:),
     &     map_ex1cnta(:), map_ex1cntc(:),
     &     map_ex2cnta(:), map_ex2cntc(:)

      integer, pointer ::
     &     ndis_op1(:,:), d_gam_ms_op1(:,:,:), gam_ms_op1(:,:),
     &     ndis_op2(:,:), d_gam_ms_op2(:,:,:), gam_ms_op2(:,:),
     &     ndis_op1op2(:,:), d_gam_ms_op1op2(:,:,:), gam_ms_op1op2(:,:)

      type(graph), pointer ::
     &     graphs(:)

      integer, external ::
     &     ielsum, ielprd, idx_msgmdst2, get_lenmap, idxlist
      logical, external ::
     &     next_dist, next_msgamdist2
      real(8), external ::
     &     ddot

      if (ntest.gt.0) then
        call write_title(luout,wst_dbg_subr,
     &       'contr_op1op2_wmaps_c at work')
      end if
      if (ntest.ge.10) then
        write(luout,*) 'ffop1:   ',ffop1%name(1:len_trim(ffop1%name))
        write(luout,*) 'ffop2:   ',ffop2%name(1:len_trim(ffop2%name))
        write(luout,*) 'ffop1op2:',
     &       ffop1op2%name(1:len_trim(ffop1op2%name))
        write(luout,*) 'xfac = ',xfac
        write(luout,*) 'casign = ',casign
        if (type_xret.ne.0)
     &       write(luout,*) 'xret on entry = ',xret(1)
        write(luout,*) 'op1: ',op1%name(1:len_trim(op1%name)),
     &       ' block ',iblkop1
        write(luout,*) 'op2: ',op2%name(1:len_trim(op2%name)),
     &       ' block ',iblkop2
        if (iblkop1op2.gt.0) then
          write(luout,*) 'op1op2: ',
     &         op1op2%name(1:len_trim(op1op2%name)),
     &       ' block ',iblkop1op2
        else
          write(luout,*) 'op1op2: scalar'
        end if
      end if

      idxst_op1 = op1%off_op_occ(iblkop1) + 1
      lenop1    = op1%len_op_occ(iblkop1)
      idxst_op2 = op2%off_op_occ(iblkop2) + 1
      lenop2    = op2%len_op_occ(iblkop2)
      if (iblkop1op2.gt.0) then
        idxst_op1op2 = op1op2%off_op_occ(iblkop1op2) + 1
        lenop1op2    = op1op2%len_op_occ(iblkop1op2)
      else
        idxst_op1op2 = 1
        lenop1op2 = 1
      end if

      if (lenop1.le.0.or.lenop2.le.0.or.lenop1op2.le.0) then
        write(luout,*) 'lenop1, lenop2, lenop1op2: ',
     &                  lenop1, lenop2, lenop1op2
        call quit(1,'contr_op1op2_wmaps_c',
     &     'zero length for operator?')
      end if

      ifree = mem_setmark('contr1')

      if (ffop1%buffered.and.ffop1%incore(iblkop1).gt.0) then
        bufop1 = .true.
        xop1 => ffop1%buffer(idxst_op1:)
      else
        bufop1 = .false.
        ifree = mem_alloc_real(xbf1,lenop1,'xbf1')
        xop1 => xbf1
        call get_vec(ffop1,xop1,idoffop1+idxst_op1,
     &                          idoffop1+idxst_op1-1+lenop1)
      end if
      if (ffop2%buffered.and.ffop2%incore(iblkop2).gt.0) then
        bufop2 = .true.
        xop2 => ffop2%buffer(idxst_op2:)
      else
        bufop2 = .false.
        ifree = mem_alloc_real(xbf2,lenop2,'xbf2')
        xop2 => xbf2
        call get_vec(ffop2,xop2,idoffop2+idxst_op2,
     &                          idoffop2+idxst_op2-1+lenop2)
      end if

      if (ntest.ge.100) write(luout,*) ' bufop1/2: ',bufop1,bufop2

      ! get result vector as well (as we update)
      if (iblkop1op2.gt.0) then
        if (ffop1op2%buffered.and.ffop1op2%incore(iblkop1op2).gt.0) then
          bufop1op2 = .true.
          xop1op2 => ffop1op2%buffer(idxst_op1op2:)
        else
          bufop1op2 = .false.
          ifree = mem_alloc_real(xbf12,lenop1op2,'xbf12')
          xop1op2 => xbf12
          if (update) then
            ! read from disc
            call get_vec(ffop1op2,xop1op2,idoffop1op2+idxst_op1op2,
     &                             idoffop1op2+idxst_op1op2-1+lenop1op2)
          else
            ! init with zero
            xop1op2(1:lenop1op2) = 0d0
          end if
        end if
        if (ntest.ge.100) write(luout,*) ' bufop1op2: ',bufop1op2
      else
        bufop1op2 = .true.
        xop1op2 => xret
        if (ntest.ge.100) write(luout,*) ' result is scalar '
      end if

      if (ntest.ge.1000) then
        write(luout,*) 'operator 1'
        call wrt_op_buf(luout,5,xop1,op1,iblkop1,iblkop1,
     &                  str_info,orb_info)
        write(luout,*) 'operator 2'
        call wrt_op_buf(luout,5,xop2,op2,iblkop2,iblkop2,
     &                  str_info,orb_info)
        if (iblkop1op2.gt.0) then
          write(luout,*) 'operator 12 on entry'
          call wrt_op_buf(luout,5,xop1op2,op1op2,
     &                    iblkop1op2,iblkop1op2,
     &                    str_info,orb_info)
        end if
      end if

      graphs => str_info%g

      ndis_op1 => op1%off_op_gmox(iblkop1)%ndis
      gam_ms_op1 => op1%off_op_gmo(iblkop1)%gam_ms
      d_gam_ms_op1 => op1%off_op_gmox(iblkop1)%d_gam_ms
      ndis_op2 => op2%off_op_gmox(iblkop2)%ndis
      gam_ms_op2 => op2%off_op_gmo(iblkop2)%gam_ms
      d_gam_ms_op2 => op2%off_op_gmox(iblkop2)%d_gam_ms
      ndis_op1op2 => op1op2%off_op_gmox(iblkop1op2)%ndis
      gam_ms_op1op2 => op1op2%off_op_gmo(iblkop1op2)%gam_ms
      d_gam_ms_op1op2 => op1op2%off_op_gmox(iblkop1op2)%d_gam_ms

      call sum_occ(nc_op1,cinfo_op1c,nca_blk(1,1))
      call sum_occ(na_op1,cinfo_op1a,nca_blk(2,1))
      call sum_occ(nc_op2,cinfo_op2c,nca_blk(1,2))
      call sum_occ(na_op2,cinfo_op2a,nca_blk(2,2))
      call sum_occ(nc_op1op2,cinfo_op1op2c,nca_blk(1,3))
      call sum_occ(na_op1op2,cinfo_op1op2a,nca_blk(2,3))
      call sum_occ(nc_ex1,cinfo_ex1c,nca_blk(1,4))
      call sum_occ(na_ex1,cinfo_ex1a,nca_blk(2,4))
      call sum_occ(nc_ex2,cinfo_ex2c,nca_blk(1,5))
      call sum_occ(na_ex2,cinfo_ex2a,nca_blk(2,5))
      call sum_occ(nc_cnt,cinfo_cntc,nca_blk(1,6))
      call sum_occ(na_cnt,cinfo_cnta,nca_blk(2,6))

c dbg
c      print *,'casign = ',casign
c dbg

      ! set up maps (if necessary)
      call strmap_man_c(
     &     cinfo_ex1c(1,2),nca_blk(1,4),
     &     cinfo_ex2c(1,2),nca_blk(1,5),
     &     cinfo_op1op2c(1,2),nca_blk(1,3),map_info_12c,
     &     str_info,strmap_info,orb_info)
      call strmap_man_c(
     &     cinfo_ex2a(1,2),nca_blk(2,5),
     &     cinfo_ex1a(1,2),nca_blk(2,4),
     &     cinfo_op1op2a(1,2),nca_blk(2,3),map_info_12a,
     &     str_info,strmap_info,orb_info)
      ! get ex2,ex1 sequence, as well ???
      ! ...
      call strmap_man_c(
     &     cinfo_cntc(1,2),nca_blk(1,6),
     &     cinfo_ex1c(1,2),nca_blk(1,4),
     &     cinfo_op1c(1,2),nca_blk(1,1),map_info_1c,
     &     str_info,strmap_info,orb_info)
      call strmap_man_c(
     &     cinfo_cnta(1,2),nca_blk(2,6),
     &     cinfo_ex1a(1,2),nca_blk(2,4),
     &     cinfo_op1a(1,2),nca_blk(2,1),map_info_1a,
     &     str_info,strmap_info,orb_info)
      call strmap_man_c(
     &     cinfo_cnta(1,2),nca_blk(2,6),
     &     cinfo_ex2c(1,2),nca_blk(1,5),
     &     cinfo_op2c(1,2),nca_blk(1,2),map_info_2c,
     &     str_info,strmap_info,orb_info)
      call strmap_man_c(
     &     cinfo_cntc(1,2),nca_blk(1,6),
     &     cinfo_ex2a(1,2),nca_blk(2,5),
     &     cinfo_op2a(1,2),nca_blk(2,2),map_info_2a,
     &     str_info,strmap_info,orb_info)

      ! minimum Ms for ...
      msbnd(1,1) = -nc_op1 ! operator 1
      msbnd(1,2) = -nc_op2 ! operator 2        
      msbnd(1,3) = -nc_op1op2 ! product
      ! maximum Ms for ...
      msbnd(2,1) = -msbnd(1,1)
      msbnd(2,2) = -msbnd(1,2)
      msbnd(2,3) = -msbnd(1,3)
      ! max |Ms| for ...
      mscmx_a = na_cnt ! C(A)
      mscmx_c = nc_cnt ! C(C)
      ! minimum IRREP for operators
      igambnd(1,1) = 1
      igambnd(1,2) = 1
      igambnd(1,3) = 1
      ! maximum IRREP
      nsym = orb_info%nsym
      igambnd(2,1) = nsym
      igambnd(2,2) = nsym
      igambnd(2,3) = nsym
      ! loop Ms-cases of (Op1(A),Op2(A),Op1Op2(A))
      first1 = .true.
      ms_loop: do
        if (first1) then
          first1 = .false.
          ! initial Ms distribution
          ms12i_a(1:3) = msbnd(2,1:3)
        else
          ! next Ms distribution
          if (.not.next_dist(ms12i_a,3,msbnd,-2)) exit
        end if

        ms12i_c(1) = ms12i_a(1) + mstop1
        ms12i_c(2) = ms12i_a(2) + mstop2
        ms12i_c(3) = ms12i_a(3) + mstop1op2
        msc_ac = ms12i_a(1) + ms12i_a(2) - ms12i_a(3)

        if (mscmx_a+mscmx_c.lt.abs(msc_ac)) cycle ms_loop
          
        msc_loop: do msc_a = mscmx_a, -mscmx_a, -2
          msc_c = msc_ac - msc_a
          if (abs(msc_c).gt.mscmx_c) cycle
          ! Ms of string after lifting restrictions:
          !  Op1(C1,A1) -> Op1(C10,A10;CC,AC)
          msex1_c = ms12i_c(1) - msc_c
          if (abs(msex1_c).gt.nc_ex1)
     &         cycle msc_loop
          msex1_a = ms12i_a(1) - msc_a
          if (abs(msex1_a).gt.na_ex1)
     &         cycle msc_loop
          msex2_a = ms12i_a(2) - msc_c   ! other way 
          msex2_c = ms12i_c(2) - msc_a   ! round (!!)
          if (abs(msex2_c).gt.nc_ex2)
     &         cycle msc_loop
          if (abs(msex2_a).gt.na_ex2)
     &         cycle msc_loop
          if (ntest.ge.100) then
            write(luout,*) 'Current spin case:'
            write(luout,*) ' OP1/OP2/INT (C) ->',ms12i_c(1:3)
            write(luout,*) ' OP1/OP2/INT (A) ->',ms12i_a(1:3)
            write(luout,*) ' CNT(C)/CNT(A)   ->',msc_c,msc_a
          end if

          ! loop IRREP cases of (Op1(A),Op2(A),Interm)
          first2 = .true.
          gam_loop: do
            if (first2) then
              first2 = .false.
              ! initial IRREP distribution
              igam12i_a(1:3) = igambnd(1,1:3)
            else
              ! next IRREP distribution
              if (.not.next_dist(igam12i_a,3,igambnd,+1)) exit
            end if
            ! set up start addresses
            ! need to be modified, if more than one distribution
            ! exists, see below
            idxms = (na_op1-ms12i_a(1))/2 + 1
            idxop1 = gam_ms_op1(igam12i_a(1),idxms) + 1
     &             - idxst_op1+1
            idxms = (na_op2-ms12i_a(2))/2 + 1
            idxop2 = gam_ms_op2(igam12i_a(2),idxms) + 1
     &             - idxst_op2+1
            idxms = (na_op1op2-ms12i_a(3))/2 + 1
            if (iblkop1op2.gt.0)
     &           idxop1op2 = gam_ms_op1op2(igam12i_a(3),idxms) + 1
     &                - idxst_op1op2+1
            if (iblkop1op2.eq.0) idxop1op2 = 1

            igam12i_c(1) = multd2h(igam12i_a(1),igamtop1)
            igam12i_c(2) = multd2h(igam12i_a(2),igamtop2)
            igam12i_c(3) = multd2h(igam12i_a(3),igamtop1op2)

            igamc_ac = multd2h(igam12i_a(1),igam12i_a(2))
            igamc_ac = multd2h(igamc_ac,igam12i_a(3))

            gamc_loop: do igamc_a = 1, nsym
              igamc_c = multd2h(igamc_a,igamc_ac)

              ! IRREPs after lifting restrictions (cf. above)
              igamex1_a = multd2h(igam12i_a(1),igamc_a)
              igamex1_c = multd2h(igam12i_c(1),igamc_c)
              igamex2_a = multd2h(igam12i_a(2),igamc_c) !  !!
              igamex2_c = multd2h(igam12i_c(2),igamc_a) !  !!

              call atim_cs(cpu00,sys00)

              ! loop over distributions of current Ms and IRREP 
              ! of Aex1 and Cex1 over ngastypes
              first3 = .true.
              caex1_loop: do
                if (.not.next_msgamdist2(first3,
     &             msex1dis_c,msex1dis_a,gmex1dis_c,gmex1dis_a,
     &             nca_blk(1,4), nca_blk(2,4),
     &             cinfo_ex1c,cinfo_ex1a,
     &             msex1_c,msex1_a,igamex1_c,igamex1_a,nsym))
     &             exit caex1_loop
                first3 = .false.

                call ms2idxms(idxmsex1dis_c,msex1dis_c,
     &               cinfo_ex1c,nca_blk(1,4))
                call ms2idxms(idxmsex1dis_a,msex1dis_a,
     &               cinfo_ex1a,nca_blk(2,4))

                call set_len_str(lstrex1,nca_blk(1,4),nca_blk(2,4),
     &                  graphs,
     &                  cinfo_ex1c(1,2),idxmsex1dis_c,
     &                                 gmex1dis_c,cinfo_ex1c(1,3),
     &                  cinfo_ex1a(1,2),idxmsex1dis_a,
     &                                 gmex1dis_a,cinfo_ex1a(1,3),
     &                  hpvxseq,.false.)
                
                ! test C and A separately to avoid overflow
                if ( nca_blk(1,4)+nca_blk(2,4).gt.0 .and.
     &               idxlist(0,lstrex1,
     &                          nca_blk(1,4)+nca_blk(2,4),1).gt.0)
     &               cycle

                ! loop over distributions of current Ms and IRREP 
                ! of Aex2 and Cex2 over ngastypes            
                first4 = .true.
                caex2_loop: do
                  if (.not.next_msgamdist2(first4,
     &               msex2dis_c,msex2dis_a,gmex2dis_c,gmex2dis_a,
     &               nca_blk(1,5), nca_blk(2,5),
     &               cinfo_ex2c,cinfo_ex2a,
     &               msex2_c,msex2_a,igamex2_c,igamex2_a,nsym))
     &               exit caex2_loop
                  first4 = .false.

                  call ms2idxms(idxmsex2dis_c,msex2dis_c,
     &                 cinfo_ex2c,nca_blk(1,5))
                  call ms2idxms(idxmsex2dis_a,msex2dis_a,
     &                 cinfo_ex2a,nca_blk(2,5))

                  call set_len_str(lstrex2,nca_blk(1,5),nca_blk(2,5),
     &                 graphs,
     &                 cinfo_ex2c(1,2),idxmsex2dis_c,
     &                                gmex2dis_c,cinfo_ex2c(1,3),
     &                 cinfo_ex2a(1,2),idxmsex2dis_a,
     &                                gmex2dis_a,cinfo_ex2a(1,3),
     &                 hpvxseq,.false.)

                  if ( nca_blk(1,5)+nca_blk(2,5).gt.0.and.
     &                 idxlist(0,lstrex2,
     &                 nca_blk(1,5)+nca_blk(2,5),1).gt.0)
     &                 cycle

                  ! get Ms and IRREP distribution of intermediate
                  call merge_msgmdis(msi_dis_c,gmi_dis_c,
     &                                 nca_blk(1,3),
     &                                 msex1dis_c,gmex1dis_c,
     &                                 msex2dis_c,gmex2dis_c,
     &                                 map_info_12c)
                  call merge_msgmdis(msi_dis_a,gmi_dis_a,
     &                                 nca_blk(2,3),
     &                                 msex2dis_a,gmex2dis_a,
     &                                 msex1dis_a,gmex1dis_a,
     &                                 map_info_12a)

                  call ms2idxms(idxmsi_dis_c,msi_dis_c,
     &                   cinfo_op1op2c,nca_blk(1,3))
                  call ms2idxms(idxmsi_dis_a,msi_dis_a,
     &                   cinfo_op1op2a,nca_blk(2,3))

                  call set_len_str(
     &                   lstrop1op2,nca_blk(1,3),nca_blk(2,3),
     &                   graphs,
     &                   cinfo_op1op2c(1,2),idxmsi_dis_c,
     &                                    gmi_dis_c,cinfo_op1op2c(1,3),
     &                   cinfo_op1op2a(1,2),idxmsi_dis_a,
     &                                    gmi_dis_a,cinfo_op1op2a(1,3),
     &                   hpvxseq,.false.)
                  
                  ! get igrphext1,igrphext2->igrphop1op2 map
                  ! for given ms and irreps
                  ifree = mem_setmark('ex_str')
                  ! special product with map_info ...
                  lenmap = get_lenmap(lstrex1,lstrex2,
     &                   map_info_12c,nca_blk(1,3))
c dbg
c                  print *,'lenmap C: ',lenmap
c dbg
                  
                  ifree = mem_alloc_int(map_ex1ex2c,lenmap,'strmap_c')
                  ! for C: ex1,ex2 sequence !
                  call get_strmap_blk_c(map_ex1ex2c,
     &                 nca_blk(1,4),nca_blk(1,5),nca_blk(1,3),
     &                 cinfo_ex1c,cinfo_ex2c,lstrex1,lstrex2,
     &                 cinfo_ex1c(1,2),cinfo_ex2c(1,2),
     &                 idxmsex1dis_c,idxmsex2dis_c,
     &                 gmex1dis_c,gmex2dis_c,map_info_12c,
     &                 strmap_info,nsym,str_info%ngraph)

                  lenmap = get_lenmap(lstrex2(nca_blk(1,5)+1),
     &                                lstrex1(nca_blk(1,4)+1),
     &                   map_info_12a,nca_blk(2,3))
c dbg
c                  print *,'lenmap A: ',lenmap
c dbg
                  ifree = mem_alloc_int(map_ex1ex2a,lenmap,'strmap_a')
                  ! for A: ex2,ex1 sequence !
                  call get_strmap_blk_c(map_ex1ex2a,
     &                 nca_blk(2,5),nca_blk(2,4),nca_blk(2,3),
     &                 cinfo_ex2a,cinfo_ex1a,
     &                  lstrex2(nca_blk(1,5)+1),lstrex1(nca_blk(1,4)+1),
     &                 cinfo_ex2a(1,2),cinfo_ex1a(1,2),
     &                 idxmsex2dis_a,idxmsex1dis_a,
     &                 gmex2dis_a,gmex1dis_a,map_info_12a,
     &                 strmap_info,nsym,str_info%ngraph)

                  ! get distribution index
                  idxms = (na_op1op2-ms12i_a(3))/2 + 1
                  if (iblkop1op2.gt.0.and.
     &                 ndis_op1op2(igam12i_a(3),idxms).gt.1) then
                    idxdis =
     &                  idx_msgmdst2(
     &                   iblkop1op2,idxms,igam12i_a(3),
     &                   cinfo_op1op2c,idxmsi_dis_c,
     &                              gmi_dis_c,nca_blk(1,3),
     &                   cinfo_op1op2a,idxmsi_dis_a,
     &                              gmi_dis_a,nca_blk(2,3),
     &                   .false.,op1op2,nsym)
                         
                    idxop1op2 = 
     &                   d_gam_ms_op1op2(idxdis,igam12i_a(3),idxms) + 1
     &                   - idxst_op1op2+1
                  end if

                  ! loop over distributions of current Ms and IRREP 
                  ! of AC and CC over ngastypes            
                  first5 = .true.
                  cac_loop: do
                    if (.not.next_msgamdist2(first5,
     &                 msc_dis_c,msc_dis_a,gmc_dis_c,gmc_dis_a,
     &                 nca_blk(1,6), nca_blk(2,6),
     &                 cinfo_cntc,cinfo_cnta,
     &                 msc_c,msc_a,igamc_c,igamc_a,nsym))
     &                 exit cac_loop
                    first5 = .false.

                    ! length of contraction
                    call ms2idxms(idxmsc_dis_c,msc_dis_c,
     &                   cinfo_cntc,nca_blk(1,6))
                    call ms2idxms(idxmsc_dis_a,msc_dis_a,
     &                   cinfo_cnta,nca_blk(2,6))

                    ! length of contraction
                    call set_len_str(lstrcnt,nca_blk(1,6),nca_blk(2,6),
     &                  graphs,
     &                  cinfo_cntc(1,2),idxmsc_dis_c,
     &                                  gmc_dis_c,cinfo_cntc(1,3),
     &                  cinfo_cnta(1,2),idxmsc_dis_a,
     &                                  gmc_dis_a,cinfo_cnta(1,3),
     &                  hpvxseq,.false.)

                    if ( nca_blk(1,6)+nca_blk(2,6).gt.0 .and.
     &                   idxlist(0,lstrcnt,
     &                   nca_blk(1,6)+nca_blk(2,6),1).gt.0)
     &                   cycle cac_loop

                    ! get Ms and IRREP distribution of op1
                    call merge_msgmdis(msop1dis_c,gmop1dis_c,
     &                                 nca_blk(1,1),
     &                                 msc_dis_c,gmc_dis_c,
     &                                 msex1dis_c,gmex1dis_c,
     &                                 map_info_1c)
                    call merge_msgmdis(msop1dis_a,gmop1dis_a,
     &                                 nca_blk(2,1),
     &                                 msc_dis_a,gmc_dis_a,
     &                                 msex1dis_a,gmex1dis_a,
     &                                 map_info_1a)

                    call ms2idxms(idxmsop1dis_c,msop1dis_c,
     &                   cinfo_op1c,nca_blk(1,1))
                    call ms2idxms(idxmsop1dis_a,msop1dis_a,
     &                   cinfo_op1a,nca_blk(2,1))
                    
                    call set_len_str(
     &                   lstrop1,nca_blk(1,1),nca_blk(2,1),
     &                   graphs,
     &                   cinfo_op1c(1,2),idxmsop1dis_c,
     &                                    gmop1dis_c,cinfo_op1c(1,3),
     &                   cinfo_op1a(1,2),idxmsop1dis_a,
     &                                    gmop1dis_a,cinfo_op1a(1,3),
     &                   hpvxseq,.false.)

                    ! get distribution index
                    idxms = (na_op1-ms12i_a(1))/2 + 1
                    if (ndis_op1(igam12i_a(1),idxms).gt.1) then
                      idxdis =
     &                   idx_msgmdst2(
     &                     iblkop1,idxms,igam12i_a(1),
     &                     cinfo_op1c,idxmsop1dis_c,
     &                              gmop1dis_c,nca_blk(1,1),
     &                     cinfo_op1a,idxmsop1dis_a,
     &                              gmop1dis_a,nca_blk(2,1),
     &                     .false.,op1,nsym)
                      idxop1 = 
     &                     d_gam_ms_op1(idxdis,igam12i_a(1),idxms) + 1
     &                     - idxst_op1+1
                    end if

                    xnrm = ddot(ielprd(lstrop1,
     &                   nca_blk(1,1)+nca_blk(2,1)),
     &                   xop1(idxop1),1,xop1(idxop1),1)
                    if (xnrm.lt.1d-28) cycle cac_loop

                    ! get Ms and IRREP distribution of op2
                    ! remember: CNT^+ !
c dbg
c                    print *,'msc_dis_a:  ',msc_dis_a
c                    print *,'msex2dis_c: ',msex2dis_c
c dbg
                    call merge_msgmdis(msop2dis_c,gmop2dis_c,
     &                                 nca_blk(1,2),
     &                                 msc_dis_a,gmc_dis_a,
     &                                 msex2dis_c,gmex2dis_c,
     &                                 map_info_2c)
                    call merge_msgmdis(msop2dis_a,gmop2dis_a,
     &                                 nca_blk(2,2),
     &                                 msc_dis_c,gmc_dis_c,
     &                                 msex2dis_a,gmex2dis_a,
     &                                 map_info_2a)
c dbg
c                    print *,'msop2dis_c: ',msex2dis_c
c dbg

                    call ms2idxms(idxmsop2dis_c,msop2dis_c,
     &                   cinfo_op2c,nca_blk(1,2))
                    call ms2idxms(idxmsop2dis_a,msop2dis_a,
     &                   cinfo_op2a,nca_blk(2,2))

                    call set_len_str(
     &                   lstrop2,nca_blk(1,2),nca_blk(2,2),
     &                   graphs,
     &                   cinfo_op2c(1,2),idxmsop2dis_c,
     &                                    gmop2dis_c,cinfo_op2c(1,3),
     &                   cinfo_op2a(1,2),idxmsop2dis_a,
     &                                    gmop2dis_a,cinfo_op2a(1,3),
     &                   hpvxseq,.false.)

                    ! get distribution index
                    idxms = (na_op2-ms12i_a(2))/2 + 1
                    if (ndis_op2(igam12i_a(2),idxms).gt.1) then
                      idxdis =
     &                   idx_msgmdst2(
     &                     iblkop2,idxms,igam12i_a(2),
     &                     cinfo_op2c,idxmsop2dis_c,
     &                              gmop2dis_c,nca_blk(1,2),
     &                     cinfo_op2a,idxmsop2dis_a,
     &                              gmop2dis_a,nca_blk(2,2),
     &                     .false.,op2,nsym)
                      idxop2 = 
     &                     d_gam_ms_op2(idxdis,igam12i_a(2),idxms) + 1
     &                     - idxst_op2+1
                    end if

                    xnrm = ddot(ielprd(lstrop2,
     &                   nca_blk(1,2)+nca_blk(2,2)),
     &                   xop2(idxop2),1,xop2(idxop2),1)
                    if (xnrm.lt.1d-28) cycle cac_loop

                    ! get igrphcnt,igrphext1->igrphop1 map
                    ! for given ms and irreps
                    ifree = mem_setmark('cntstr')
                    lenmap = get_lenmap(lstrcnt,lstrex1,
     &                   map_info_1c,nca_blk(1,1))
                    ifree = mem_alloc_int(map_ex1cntc,lenmap,'strmap_c')
                    call get_strmap_blk_c(map_ex1cntc,
     &                   nca_blk(1,6),nca_blk(1,4),nca_blk(1,1),
     &                   cinfo_cntc,cinfo_ex1c,lstrcnt,lstrex1,
     &                   cinfo_cntc(1,2),cinfo_ex1c(1,2),
     &                   idxmsc_dis_c,idxmsex1dis_c,
     &                   gmc_dis_c,gmex1dis_c,map_info_1c,
     &                   strmap_info,nsym,str_info%ngraph)

                    lenmap = get_lenmap(lstrcnt(nca_blk(1,6)+1),
     &                                  lstrex1(nca_blk(1,4)+1),
     &                                  map_info_1a,nca_blk(2,1))
c dbg
c                    print *,'lstrcnt(A)',nca_blk(2,6),
c     &                 lstrcnt(nca_blk(1,6)+1:nca_blk(1,6)+nca_blk(2,6))
c                    print *,'lstrex1(A)',nca_blk(2,4),
c     &                 lstrex1(nca_blk(1,4)+1:nca_blk(1,4)+nca_blk(2,4))
c                    print *,'lenmap for map_ex1cnta: ',lenmap
c                    print *,'lstrop1(A)',nca_blk(2,1),
c     &                 lstrop1(nca_blk(1,1)+1:nca_blk(1,1)+nca_blk(2,1))
c dbg
                    ifree = mem_alloc_int(map_ex1cnta,lenmap,'strmap_a')
                    call get_strmap_blk_c(map_ex1cnta,
     &                   nca_blk(2,6),nca_blk(2,4),nca_blk(2,1),
     &                   cinfo_cnta,cinfo_ex1a,
     &                     lstrcnt(nca_blk(1,6)+1),
     &                             lstrex1(nca_blk(1,4)+1),
     &                   cinfo_cnta(1,2),cinfo_ex1a(1,2),
     &                   idxmsc_dis_a,idxmsex1dis_a,
     &                   gmc_dis_a,gmex1dis_a,map_info_1a,
     &                   strmap_info,nsym,str_info%ngraph)
c dbg
c                    print *,'map_ex1cnta(1): ',map_ex1cnta(1)
c dbg

                    ! get igrphcnt,igrphext2->igrphop2 map
                    ! for given ms and irreps
                    lenmap = get_lenmap(lstrcnt(nca_blk(1,6)+1),lstrex2,
     &                   map_info_2c,nca_blk(1,2))
                    ifree = mem_alloc_int(map_ex2cntc,lenmap,'strmap_c')
c dbg
c                  print *,'getting map_ex2cntc:'
c dbg
                    call get_strmap_blk_c(map_ex2cntc,
     &                   nca_blk(2,6),nca_blk(1,5),nca_blk(1,2),
     &                   cinfo_cnta,cinfo_ex2c,
     &                     lstrcnt(nca_blk(1,6)+1),lstrex2,
     &                   cinfo_cnta(1,2),cinfo_ex2c(1,2),
     &                   idxmsc_dis_a,idxmsex2dis_c,
     &                   gmc_dis_a,gmex2dis_c,map_info_2c,
     &                   strmap_info,nsym,str_info%ngraph)

                    lenmap = get_lenmap(lstrcnt,lstrex2(nca_blk(1,5)+1),
     &                   map_info_2a,nca_blk(2,2))
                    ifree = mem_alloc_int(map_ex2cnta,lenmap,'strmap_a')
                    call get_strmap_blk_c(map_ex2cnta,
     &                   nca_blk(1,6),nca_blk(2,5),nca_blk(2,2),
     &                   cinfo_cntc,cinfo_ex2a,
     &                       lstrcnt,lstrex2(nca_blk(1,5)+1),
     &                   cinfo_cntc(1,2),cinfo_ex2a(1,2),
     &                   idxmsc_dis_c,idxmsex2dis_a,
     &                   gmc_dis_c,gmex2dis_a,map_info_2a,
     &                   strmap_info,nsym,str_info%ngraph)                    

                    call atim_cs(cpu0,sys0)                    

                    ! make the contraction for this block
                    if (ntest.ge.100)
     &                   write(luout,*) 'calling blk1blk2',
     &                   lenop1,idxop1,
     &                   lenop2,idxop2,
     &                   lenop1op2,idxop1op2
c                    if (ntest.ge.1000) then
c                      write(luout,*) ' the maps:'
c                      write(luout,*) ' X1X2(A): '
c                      call prt_strmap_c(map_ex1ex2a,
c     &                     cinfo_ex2,cinfo_ex2,
c     &                     cinfo_ex2(1,3),cinfo_ex2(1,3),
c     &                     lstrex2(1,2),lstrex1(1,2),
c     &                     nca_blk(2,5),nca_blk(2,4))
c                      write(luout,*) ' X1X2(C): '
c                      call prt_strmap(map_ex1ex2c,
c     &                     iocc_ext1(1,1),iocc_ext2(1,1),
c     &                     lstrext1(1,1),lstrext2(1,1))
c                      write(luout,*) ' X1C(A): '
c                      call prt_strmap(map_ex1cnta,
c     &                     iocc_cnt(1,2),iocc_ext1(1,2),
c     &                     lstrcnt(1,2),lstrext1(1,2))
c                      write(luout,*) ' X1C(C): '
c                      call prt_strmap(map_ex1cntc,
c     &                     iocc_cnt(1,1),iocc_ext1(1,1),
c     &                     lstrcnt(1,1),lstrext1(1,1))
c                      write(luout,*) ' X2C(A): '
c                      call prt_strmap(map_ex2cnta,
c     &                     iocc_cnt(1,1),iocc_ext2(1,2),
c     &                     lstrcnt(1,1),lstrext2(1,2))
c                      write(luout,*) ' X2C(C): '
c                      call prt_strmap(map_ex2cntc,
c     &                     iocc_cnt(1,2),iocc_ext2(1,1),
c     &                     lstrcnt(1,2),lstrext2(1,1))
c                    end if
c dbg
c                    print *,'on call:'
c                    print *,'xop1: ',xop1(idxop1)
c                    print *,'xop2: ',xop2(idxop2)
c                    print *,'xop1op2: ',xop1op2(idxop1op2)
c dbg
                    call cntr_blk1blk2_wmaps_c(xfac*casign,
     &                   xop1op2(idxop1op2),
     &                                 xop1(idxop1),xop2(idxop2),
     &                   nca_blk,
     &                   cinfo_op1c(1,3),cinfo_op1a(1,3),
     &                   cinfo_op2c(1,3),cinfo_op2a(1,3),
     &                   cinfo_op1op2c(1,3),cinfo_op1op2a(1,3),
     &                   lstrop1,lstrop2,lstrop1op2,
     &                   lstrex1,lstrex2,lstrcnt,
     &                   map_info_12c, map_info_12a,
     &                   map_info_1c, map_info_1a,
     &                   map_info_2c, map_info_2a,
     &                   map_ex1ex2c, map_ex1ex2a,
     &                   map_ex1cntc, map_ex1cnta,
     &                   map_ex2cntc, map_ex2cnta
     &                   )                     
                    if (ntest.ge.100)
     &                   write(luout,*) 'after blk1blk2'
c dbg
c                    if (lenop1op2.eq.22)
c     &                   print *,'xop1op2: ',xop1op2(1:11)
c dbg

                    call atim_cs(cpu,sys)
                    cnt_kernel(1) = cnt_kernel(1)+cpu-cpu0
                    cnt_kernel(2) = cnt_kernel(2)+sys-sys0

                    ifree = mem_flushmark('cntstr')

                  end do cac_loop
                  
                  ifree = mem_flushmark('ex_str')
                end do caex2_loop
              end do caex1_loop

              call atim_cs(cpu,sys)
              cnt_dloop(1) = cnt_dloop(1)+cpu-cpu00
              cnt_dloop(2) = cnt_dloop(2)+sys-sys00

            end do gamc_loop
          end do gam_loop

        end do msc_loop
      end do ms_loop

      if (ntest.ge.1000) then
        if (iblkop1op2.gt.0
     &       ) then
          write(luout,*) 'operator 12 on exit'
          call wrt_op_buf(luout,5,xop1op2,op1op2,
     &         iblkop1op2,iblkop1op2,str_info,orb_info)
        end if
      end if
c dbg
c          write(luout,*) 'operator 12 on exit'
c          call wrt_op_buf(luout,2,xop1op2,op1op2,
c     &         iblkop1op2,iblkop1op2,str_info,orb_info)
c dbg

      if (type_xret.eq.2) then
        xret(1) = xop1op2(1)
      else if (type_xret.eq.1) then
        xret(1) = ddot(lenop1op2,xop1op2,1,xop1op2,1)
      end if

      ! put result to disc
      if (.not.bufop1op2) then
        call put_vec(ffop1op2,xop1op2,idoffop1op2+idxst_op1op2,
     &                    idoffop1op2+idxst_op1op2-1+lenop1op2)
      end if

      ifree = mem_flushmark()

      if (ntest.ge.100) then
        if (type_xret.ne.0)
     &       write(luout,*) 'xret on exit = ',xret(1)
      end if

      return
      end
*----------------------------------------------------------------------*
      subroutine cntr_blk1blk2_wmaps_c(xfac,            !prefactor   
     &     xop1op2,xop1,xop2,                     !buffers: res,op1,op2
     &     nca_blk,
     &     hpvxop1c, hpvxop1a,
     &     hpvxop2c, hpvxop2a,
     &     hpvxop1op2c, hpvxop1op2a,
     &     lstrop1, lstrop2, lstrop1op2,
     &     lstr_ex1, lstr_ex2, lstr_cnt,
     &     map_info_12c, map_info_12a,
     &     map_info_1c, map_info_1a,
     &     map_info_2c, map_info_2a,
     &     map_ex1ex2c, map_ex1ex2a,
     &     map_ex1cntc, map_ex1cnta,
     &     map_ex2cntc, map_ex2cnta
     &     )                     
*----------------------------------------------------------------------*
*     inner contraction routine, generation 2
*     - testing maps
*     - version c: interfaced to more general intermediates
*                  ("ngastp" can be larger than 4)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'hpvxseq.h'

      integer, parameter ::
     &     ntest = 00

      ! ngastp_op1op2a: number of non-zero entries in A occ of Op1Op2
      ! etc. for op1, op2
      ! lstrop1op2a(ngastp_op1op2a): # of strings
      ! etc.
      ! iocc_op1op2: occupation
      ! nstr_ex1/ex2/cnt,a/c: # strings
      ! map_ex1ex2a/c  mapping ex1*ex2->op1op2
      ! etc. for ex1cnt ex2cnt
      integer, intent(in) ::
     &     nca_blk(2,6),
     &     lstrop1(*), lstrop2(*), lstrop1op2(*),
     &     hpvxop1a(*), hpvxop2a(*), hpvxop1op2a(*),
     &     hpvxop1c(*), hpvxop2c(*), hpvxop1op2c(*),
     &     lstr_ex1(*), lstr_ex2(*), lstr_cnt(*),
     &     map_info_12c(*), map_info_12a(*),
     &     map_info_1c(*), map_info_1a(*),
     &     map_info_2c(*), map_info_2a(*),
     &     map_ex1ex2a(*), map_ex1ex2c(*),
     &     map_ex1cnta(*), map_ex1cntc(*),
     &     map_ex2cnta(*), map_ex2cntc(*)
      real(8), intent(in) ::
     &     xop1(*), xop2(*), xfac
      real(8), intent(inout) ::
     &     xop1op2(*)

      integer ::
     &     ielmap, ielmap2,
     &     idxop1op2a(nca_blk(2,3)), idxop1op2c(nca_blk(1,3)),
     &     nstr_ex1a1(nca_blk(2,1)), nstr_ex1c1(nca_blk(1,1)),
     &     nstr_ex1a12(nca_blk(2,3)), nstr_ex1c12(nca_blk(1,3)),
     &     nstr_ex2a2(nca_blk(2,2)), nstr_ex2c2(nca_blk(1,2)),
     &     nstr_ex2a12(nca_blk(2,3)), nstr_ex2c12(nca_blk(1,3)),
     &     nstr_cnta1(nca_blk(2,1)), nstr_cntc1(nca_blk(1,1)),
     &     nstr_cnta2(nca_blk(1,2)), nstr_cntc2(nca_blk(2,2)),
     &     nstr_ex1ex2a(nca_blk(2,3)), nstr_ex1ex2c(nca_blk(1,3)),
     &     nstr_ex1cnta(nca_blk(2,1)), nstr_ex1cntc(nca_blk(1,1)),
     &     nstr_ex2cnta(nca_blk(2,2)), nstr_ex2cntc(nca_blk(1,2)),
     &     ldim_op1a(nca_blk(2,1)), ldim_op1c(nca_blk(1,1)),
     &     ldim_op2a(nca_blk(2,2)), ldim_op2c(nca_blk(1,2)),
     &     ldim_op1op2a(nca_blk(2,3)), ldim_op1op2c(nca_blk(1,3))
      integer ::
     &     ioff, istr1, istr2, idx1, idx2, idx,
     &     ngastp_op1op2a, ngastp_op1op2c,
     &     ngastp_op1a,    ngastp_op1c,
     &     ngastp_op2a,    ngastp_op2c,
     &     nstr_ex1a_tot, nstr_ex1c_tot,
     &     nstr_ex2a_tot, nstr_ex2c_tot,
     &     nstr_cnta_tot, nstr_cntc_tot,
     &     isgnra, isgnrc, isgnr, isgnop1a, isgnop1c,
     &     isgnop2a, isgnop2c, isgn0, idxop1, idxop2, idxop1op2,
     &     idx0op1,idx0op2,
     &     istr_ex1a, istr_ex1c, istr_ex2a, istr_ex2c,
     &     istr_cnta, istr_cntc, icmp
      real(8) ::
     &     sgn

      integer, external ::
     &     idx_str_blk3, ielprd
c dbg
      integer lenop1, lenop2, lenop12

      lenop1  = ielprd(lstrop1,nca_blk(1,1)+nca_blk(2,1))
      lenop2  = ielprd(lstrop2,nca_blk(1,2)+nca_blk(2,2))
      lenop12 = ielprd(lstrop1op2,nca_blk(1,3)+nca_blk(2,3))
c dbg

      nstr_ex1c_tot = ielprd(lstr_ex1,nca_blk(1,4))
      nstr_ex1a_tot = ielprd(lstr_ex1(nca_blk(1,4)+1),nca_blk(2,4))
      nstr_ex2c_tot = ielprd(lstr_ex2,nca_blk(1,5))
      nstr_ex2a_tot = ielprd(lstr_ex2(nca_blk(1,5)+1),nca_blk(2,5))
      nstr_cntc_tot = ielprd(lstr_cnt,nca_blk(1,6))
      nstr_cnta_tot = ielprd(lstr_cnt(nca_blk(1,6)+1),nca_blk(2,6))

c dbg
c      print *,'nstr_ex1c_tot',nstr_ex1c_tot
c      print *,'nstr_ex1a_tot',nstr_ex1a_tot
c      print *,'nstr_ex2c_tot',nstr_ex2c_tot
c      print *,'nstr_ex2a_tot',nstr_ex2a_tot
c      print *,'nstr_cntc_tot',nstr_cntc_tot
c      print *,'nstr_cnta_tot',nstr_cnta_tot
c dbg

      ngastp_op1c = nca_blk(1,1)
      ngastp_op1a = nca_blk(2,1)
      ngastp_op2c = nca_blk(1,2)
      ngastp_op2a = nca_blk(2,2)
      ngastp_op1op2c = nca_blk(1,3)
      ngastp_op1op2a = nca_blk(2,3)

      ! C: ex1,ex2
      call set_strmapdim_c(nstr_ex1ex2c,nstr_ex1c12,nstr_ex2c12,
     &                   ngastp_op1op2c,
     &                   lstr_ex1,lstr_ex2,map_info_12c)
      ! A: ex2,ex1
      call set_strmapdim_c(nstr_ex1ex2a,nstr_ex2a12,nstr_ex1a12,
     &                   ngastp_op1op2a,
     &                   lstr_ex2(nca_blk(1,5)+1),
     &                    lstr_ex1(nca_blk(1,4)+1),map_info_12a)
      call set_strmapdim_c(nstr_ex1cntc,nstr_cntc1,nstr_ex1c1,
     &                   ngastp_op1c,
     &                   lstr_cnt,lstr_ex1,map_info_1c)
      call set_strmapdim_c(nstr_ex1cnta,nstr_cnta1,nstr_ex1a1,
     &                   ngastp_op1a,
     &                   lstr_cnt(nca_blk(1,6)+1),
     &                     lstr_ex1(nca_blk(1,4)+1),map_info_1a)
c dbg
c      print *,'jetzt kommt''s:'
c dbg
      call set_strmapdim_c(nstr_ex2cntc,nstr_cnta2,nstr_ex2c2,
     &                   ngastp_op2c,
     &                   lstr_cnt(nca_blk(1,6)+1),lstr_ex2,map_info_2c)
      call set_strmapdim_c(nstr_ex2cnta,nstr_cntc2,nstr_ex2a2,
     &                   ngastp_op2a,
     &                   lstr_cnt,lstr_ex2(nca_blk(1,5)+1),map_info_2a)
c dbg
c      print *,'nstr_cntc2: ',nstr_cntc2
c      print *,'nstr_cnta2: ',nstr_cnta2
c dbg
c dbg
c      print *,'ngastp_op1:',ngastp_op1c,ngastp_op1a
c      print *,'ngastp_op2:',ngastp_op2c,ngastp_op2a
c      print *,'ngastp_op1op2:',ngastp_op1op2c,ngastp_op1op2a
c dbg

      call set_op_ldim_c(ldim_op1c,ldim_op1a,
     &                   hpvxop1c,hpvxop1a,
     &                   lstrop1,nca_blk(1,1),nca_blk(2,1))
      call set_op_ldim_c(ldim_op2c,ldim_op2a,
     &                   hpvxop2c,hpvxop2a,
     &                   lstrop2,nca_blk(1,2),nca_blk(2,2))
      call set_op_ldim_c(ldim_op1op2c,ldim_op1op2a,
     &                   hpvxop1op2c,hpvxop1op2a,
     &                   lstrop1op2,nca_blk(1,3),nca_blk(2,3))
c dbg
c      print *,'ldim_op1c:',ldim_op1c(1:ngastp_op1c)
c      print *,'ldim_op1a:',ldim_op1a(1:ngastp_op1a)
c      print *,'ldim_op2c:',ldim_op2c(1:ngastp_op2c)
c      print *,'ldim_op2a:',ldim_op2a(1:ngastp_op2a)
c      print *,'ldim_op1op2c:',ldim_op1op2c(1:ngastp_op1op2c)
c      print *,'ldim_op1op2a:',ldim_op1op2a(1:ngastp_op1op2a)
c dbg

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &       'News from contr_op1op2_wmaps_c')
      end if

      ! loop over external A strings of operator 1
      ex1a: do istr_ex1a = 1, nstr_ex1a_tot

        ! loop over external A strings of operator 2
        ex2a: do istr_ex2a = 1, nstr_ex2a_tot

          ioff = 0
          isgnra = 1
          istr1 = istr_ex1a-1
          istr2 = istr_ex2a-1
          ! get sign and idxop1op2a
          do icmp = 1, ngastp_op1op2a
            idx1 = mod(istr1,nstr_ex1a12(icmp))+1
            idx2 = mod(istr2,nstr_ex2a12(icmp))+1
            idx  = (idx1-1)*nstr_ex2a12(icmp)+idx2
            ielmap = map_ex1ex2a(ioff+idx)
            if (ielmap.eq.0) cycle ex2a
            isgnra = isgnra*sign(1,ielmap)
            idxop1op2a(icmp) = abs(ielmap)-1
            ioff = ioff + nstr_ex1ex2a(icmp)
            istr1 = istr1/nstr_ex1a12(icmp)
            istr2 = istr2/nstr_ex2a12(icmp)
          end do

          ! loop over external C strings of operator 1
          ex1c: do istr_ex1c = 1, nstr_ex1c_tot

            ! loop over external C strings of operator 2
            ex2c: do istr_ex2c = 1, nstr_ex2c_tot
            
              ioff = 0
              isgnrc = 1
              istr1 = istr_ex2c-1
              istr2 = istr_ex1c-1
              ! get sign and idxop1op2c
              do icmp = 1, ngastp_op1op2c
                idx1 = mod(istr1,nstr_ex2c12(icmp))+1
                idx2 = mod(istr2,nstr_ex1c12(icmp))+1
                idx  = (idx1-1)*nstr_ex1c12(icmp)+idx2
                ielmap = map_ex1ex2c(ioff+idx)
                if (ielmap.eq.0) cycle ex2c
                isgnrc = isgnrc*sign(1,ielmap)
                idxop1op2c(icmp) = abs(ielmap)-1
                ioff = ioff + nstr_ex1ex2c(icmp)
                istr1 = istr1/nstr_ex2c12(icmp)
                istr2 = istr2/nstr_ex1c12(icmp)
              end do

              isgnr = isgnrc*isgnra

              idxop1op2 = idx_str_blk3(idxop1op2c,idxop1op2a,
     &                                 ldim_op1op2c,ldim_op1op2a,
     &                                 ngastp_op1op2c,ngastp_op1op2a)

              ! loop over contraction A string
              cnta: do istr_cnta = 1, nstr_cnta_tot
                ioff = 0
                isgnop1a = 1
                istr1 = istr_ex1a-1
                istr2 = istr_cnta-1
                idx0op1 = 1
                ! get sign and idxop1a
                do icmp = 1, ngastp_op1a
                  idx1 = mod(istr1,nstr_ex1a1(icmp))+1
                  idx2 = mod(istr2,nstr_cnta1(icmp))+1
                  idx  = (idx1-1)*nstr_cnta1(icmp)+idx2
                  ielmap = map_ex1cnta(ioff+idx)
                  if (ielmap.eq.0) cycle cnta
                  isgnop1a = isgnop1a*sign(1,ielmap)
                  idx0op1 = idx0op1+(abs(ielmap)-1)*ldim_op1a(icmp)
                  ioff = ioff + nstr_ex1cnta(icmp)
                  istr1 = istr1/nstr_ex1a1(icmp)
                  istr2 = istr2/nstr_cnta1(icmp)
                end do

                ioff = 0
                isgnop2c = 1
                istr1 = istr_ex2c-1
                istr2 = istr_cnta-1
                idx0op2 = 1
                ! get sign and idxop2c
                do icmp = 1, ngastp_op2c
                  idx1 = mod(istr1,nstr_ex2c2(icmp))+1
                  idx2 = mod(istr2,nstr_cnta2(icmp))+1
                  idx  = (idx1-1)*nstr_cnta2(icmp)+idx2
                  ielmap = map_ex2cntc(ioff+idx)
                  if (ielmap.eq.0) cycle cnta
                  isgnop2c = isgnop2c*sign(1,ielmap)
                  idx0op2 = idx0op2+(abs(ielmap)-1)*ldim_op2c(icmp)
                  ioff = ioff + nstr_ex2cntc(icmp)
                  istr1 = istr1/nstr_ex2c2(icmp)
                  istr2 = istr2/nstr_cnta2(icmp)
                end do

                isgn0 = isgnr*isgnop1a*isgnop2c

                ! loop over contraction C string
                cntc: do istr_cntc = 1, nstr_cntc_tot

                  ioff = 0
                  isgnop1c = 1
                  istr1 = istr_ex1c-1
                  istr2 = istr_cntc-1
                  idxop1 = idx0op1
                  ! get sign and idxop1c
                  do icmp = 1, ngastp_op1c
                    idx1 = mod(istr1,nstr_ex1c1(icmp))+1
                    idx2 = mod(istr2,nstr_cntc1(icmp))+1
                    idx  = (idx1-1)*nstr_cntc1(icmp)+idx2
                    ielmap = map_ex1cntc(ioff+idx)
                    if (ielmap.eq.0) cycle cntc
                    isgnop1c = isgnop1c*sign(1,ielmap)
                    idxop1 = idxop1+(abs(ielmap)-1)*ldim_op1c(icmp)
                    ioff = ioff + nstr_ex1cntc(icmp)
                    istr1 = istr1/nstr_ex1c1(icmp)
                    istr2 = istr2/nstr_cntc1(icmp)
                  end do

                  ioff = 0
                  isgnop2a = 1
                  istr1 = istr_ex2a-1
                  istr2 = istr_cntc-1
                  idxop2 = idx0op2
                  ! get sign and idxop2a
                  do icmp = 1, ngastp_op2a
                    idx1 = mod(istr1,nstr_ex2a2(icmp))+1
                    idx2 = mod(istr2,nstr_cntc2(icmp))+1
                    idx  = (idx1-1)*nstr_cntc2(icmp)+idx2
                    ielmap = map_ex2cnta(ioff+idx)
                    if (ielmap.eq.0) cycle cntc
                    isgnop2a = isgnop2a*sign(1,ielmap)
                    idxop2 = idxop2+(abs(ielmap)-1)*ldim_op2a(icmp)
                    ioff = ioff + nstr_ex2cnta(icmp)
                    istr1 = istr1/nstr_ex2a2(icmp)
                    istr2 = istr2/nstr_cntc2(icmp)
                  end do

                  sgn = xfac*dble(isgn0*isgnop1c*isgnop2a)

                  ! xfac contained in sgn
c dbg
                  if (idxop1.lt.1) stop 'range1'
                  if (idxop2.lt.1) stop 'range2'
                  if (idxop1op2.lt.1) stop 'range12'
                  if (idxop1.gt.lenop1) stop 'range1'
                  if (idxop2.gt.lenop2) stop 'range2'
                  if (idxop1op2.gt.lenop12) stop 'range12'
c dbg
                  xop1op2(idxop1op2) = xop1op2(idxop1op2)
     &                          + sgn * xop1(idxop1)
     &                                * xop2(idxop2)

                end do cntc           
              end do cnta

            end do ex2c
          end do ex1c

        end do ex2a
      end do ex1a
        
      return
      end

      
*----------------------------------------------------------------------*
      pure integer function idx_str_blk3(idxc,idxa,ldimc,ldima,
     &                                   ngastp_c,ngastp_a)
*----------------------------------------------------------------------*
*
*     version for compressed lists on idxc/a, lenc/a
*      
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     ngastp_c, ngastp_a,
     &     idxc(*), idxa(*), 
     &     ldimc(*), ldima(*)

      integer ::
     &     idx

      ! we have to add a one to get the index
      idx_str_blk3 = 1

      if ((ngastp_c*ngastp_a).eq.0) return

      idx_str_blk3 = idxc(1)*ldimc(1)+idxa(1)*ldima(1) + 1
      if ((ngastp_c*ngastp_a).eq.1) return

      do idx = 2, ngastp_c
        idx_str_blk3 = idx_str_blk3+idxc(idx)*ldimc(idx)
      end do
      do idx = 2, ngastp_a
        idx_str_blk3 = idx_str_blk3+idxa(idx)*ldima(idx)
      end do

      return
      end

      subroutine set_strmapdim_c(nstr1str2,nstr1,nstr2,
     &                         nblk12,
     &                         nstr1r,nstr2r,map12)

      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(out) ::
     &     nstr1str2(*), nstr1(*), nstr2(*)
      integer, intent(in) ::
     &     nblk12,
     &     nstr1r(*), nstr2r(*), map12(*)
            
      integer ::
     &     idxmap, iblk12, iblk, len1, len2, nblk

      idxmap = 0
      do iblk12 = 1, nblk12
c dbg
c        print *,'iblk12: ',iblk12
c        print *,'map12: ',map12(idxmap+1:idxmap+4)
c dbg
        idxmap = idxmap+1
        nblk = map12(idxmap) ! already adapted for arbitrary nblk>1
                             ! (not sure whether this works ...)
        len1 = 1
        do iblk = 1, nblk
          idxmap = idxmap+1
          len1 = len1*nstr1r(map12(idxmap))
        end do
        idxmap = idxmap+1
        len2 = 1
        nblk = map12(idxmap)
        do iblk = 1, nblk
          idxmap = idxmap+1
          len2 = len2*nstr2r(map12(idxmap))
        end do
        nstr1str2(iblk12) = len1*len2
        nstr1(iblk12) = len1
        nstr2(iblk12) = len2
c dbg
c        print *,'len1, len2: ',len1,len2
c dbg
      end do

      return
      end

      subroutine set_op_ldim_c(ldimc,ldima,hpvxc,hpvxa,nstr,nc,na)
      implicit none
      
      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(out) ::
     &     ldimc(nc), ldima(na)
      integer, intent(in) ::
     &     hpvxc(nc), hpvxa(na),
     &     nc, na, nstr(nc+na)

      integer ::
     &     idx, ldim, ioff, idx_hpvx, hpvx

      ldim = 1
      do idx_hpvx = 1, ngastp
        hpvx = hpvxseq(idx_hpvx)

        do idx = 1, nc
          if (hpvxc(idx).ne.hpvx) cycle
          ldimc(idx) = ldim
          ldim = ldim * nstr(idx)
        end do

        ioff = nc
        do idx = 1, na
          if (hpvxa(idx).ne.hpvx) cycle
          ldima(idx) = ldim
          ldim = ldim * nstr(ioff+idx)
        end do
      end do

      return
      end
      
