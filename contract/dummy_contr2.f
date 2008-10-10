*----------------------------------------------------------------------*
      subroutine dummy_contr2(flops,xmemtot,xmemblk,
     &     cnt_info,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     str_info,ngas,nsym)
*----------------------------------------------------------------------*
*     perform a dummy contraction for the operator blocks given
*     and estimate the flops and the memory consumption
*     CAVE: currently, restrictions are not handled correctly
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'multd2h.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_contraction_info.h'
      include 'hpvxseq.h'

      integer, parameter ::
     &     ntest = 00

      real(8), intent(out) ::
     &     flops, xmemtot, xmemblk
      type(contraction_info), target ::
     &     cnt_info
      integer, intent(in) ::
     &     mstop1,mstop2,mstop1op2,igamtop1,igamtop2,igamtop1op2
      integer, intent(in) ::
     &     ngas, nsym
      type(strinf), intent(in) ::
     &     str_info

      logical ::
     &     first1, first2, first3, first4, first5, fix_success
      integer ::
     &     nca(2,6),
     &     mscmx_a, mscmx_c, msc_ac, msc_a, msc_c,
     &     msex1_a, msex1_c, msex2_a, msex2_c,
     &     igamc_ac, igamc_a, igamc_c,
     &     igamex1_a, igamex1_c, igamex2_a, igamex2_c,
     &     icmp, maxidxms
      integer ::
     &     irst_ext(2,ngas,2,2,2),irst_cnt(2,ngas,2,2),
     &     ms12i_a(3), ms12i_c(3), igam12i_a(3), igam12i_c(3),
     &     msbnd(2,3), igambnd(2,3)
      integer, pointer ::
     &     gmex1dis_c(:), gmex1dis_a(:),
     &     gmex2dis_c(:), gmex2dis_a(:),
     &     gmc_dis_c(:), gmc_dis_a(:),
     &     gmi_dis_c(:), gmi_dis_a(:),
     &     msex1dis_c(:), msex1dis_a(:),
     &     msex2dis_c(:), msex2dis_a(:),
     &     msc_dis_c(:), msc_dis_a(:),
     &     msi_dis_c(:), msi_dis_a(:),
     &     idxmsex1dis_c(:), idxmsex1dis_a(:),
     &     idxmsex2dis_c(:), idxmsex2dis_a(:),
     &     idxmsc_dis_c(:), idxmsc_dis_a(:),
     &     idxmsi_dis_c(:), idxmsi_dis_a(:),
     &     lenex1(:), lenex2(:), lencnt(:),
     &     lenop1op2(:)
      integer, pointer ::
     &     ncblk_op1, nablk_op1, ncblk_ex1, nablk_ex1, 
     &     ncblk_op2, nablk_op2, ncblk_ex2, nablk_ex2, 
     &     ncblk_op1op2, nablk_op1op2, ncblk_op1op2tmp, nablk_op1op2tmp, 
     &     ncblk_cnt, nablk_cnt,
     &     cinfo_op1c(:,:),cinfo_op1a(:,:),
     &     cinfo_op2c(:,:),cinfo_op2a(:,:),
     &     cinfo_op1op2c(:,:),
     &     cinfo_op1op2a(:,:),
     &     cinfo_op1op2tmpc(:,:),
     &     cinfo_op1op2tmpa(:,:),
     &     cinfo_ex1c(:,:),cinfo_ex1a(:,:),
     &     cinfo_ex2c(:,:),cinfo_ex2a(:,:),
     &     cinfo_cntc(:,:),cinfo_cnta(:,:),
     &     map_info_1c(:),
     &     map_info_1a(:),
     &     map_info_2c(:),
     &     map_info_2a(:),
     &     map_info_12c(:),
     &     map_info_12a(:),
     &     lenstr_array(:,:,:)
      type(graph), pointer ::
     &     graphs(:)
      real(8) ::
     &     xlenex1, xlenex2, xlencnt, xlenop1op2

      logical, external ::
     &     next_dist, next_msgamdist2, idxlist
      logical, external ::
     &     zero_ivec

      if (ntest.ge.10) then
        write(luout,*) '========================'
        write(luout,*) ' here comes dummy_contr'
        write(luout,*) '========================'
      end if

      ! init
      flops = 0d0
      xmemtot = 0d0
      xmemblk = 0d0

      maxidxms =  str_info%max_idxms
      allocate(lenstr_array(nsym,str_info%max_idxms,str_info%ngraph))
      call set_lenstr_array(lenstr_array,nsym,
     &                      str_info%max_idxms,str_info)

      graphs => str_info%g

      ! set a few pointers:
      ncblk_op1 => cnt_info%ncblk_op1 ! 1,1
      nablk_op1 => cnt_info%nablk_op1 ! 2,1
      ncblk_ex1 => cnt_info%ncblk_ex1 ! 1,4
      nablk_ex1 => cnt_info%nablk_ex1 ! 2,4
      ncblk_op2 => cnt_info%ncblk_op2 ! 1,2 
      nablk_op2 => cnt_info%nablk_op2 ! 2,2
      ncblk_ex2 => cnt_info%ncblk_ex2 ! 1,5
      nablk_ex2 => cnt_info%nablk_ex2 ! 2,5
      ncblk_cnt => cnt_info%ncblk_cnt ! 1,3
      nablk_cnt => cnt_info%nablk_cnt ! 2,3
      ncblk_op1op2 => cnt_info%ncblk_op1op2 ! 1,6
      nablk_op1op2 => cnt_info%nablk_op1op2 ! 2,6
      ncblk_op1op2tmp => cnt_info%ncblk_op1op2tmp ! 1,7
      nablk_op1op2tmp => cnt_info%nablk_op1op2tmp ! 2,7

      cinfo_op1c => cnt_info%cinfo_op1c
      cinfo_op1a => cnt_info%cinfo_op1a
      cinfo_ex1c => cnt_info%cinfo_ex1c
      cinfo_ex1a => cnt_info%cinfo_ex1a
      cinfo_op2c => cnt_info%cinfo_op2c
      cinfo_op2a => cnt_info%cinfo_op2a
      cinfo_ex2c => cnt_info%cinfo_ex2c
      cinfo_ex2a => cnt_info%cinfo_ex2a
      cinfo_op1op2c => cnt_info%cinfo_op1op2c
      cinfo_op1op2a => cnt_info%cinfo_op1op2a
      cinfo_op1op2tmpc => cnt_info%cinfo_op1op2tmpc
      cinfo_op1op2tmpa => cnt_info%cinfo_op1op2tmpa
      cinfo_cntc => cnt_info%cinfo_cntc
      cinfo_cnta => cnt_info%cinfo_cnta

      map_info_1c => cnt_info%map_info_1c
      map_info_1a => cnt_info%map_info_1a
      map_info_2c => cnt_info%map_info_2c
      map_info_2a => cnt_info%map_info_2a
      map_info_12c => cnt_info%map_info_12c
      map_info_12a => cnt_info%map_info_12a

      call sum_occ(nca(1,1),cinfo_op1c,ncblk_op1)
      call sum_occ(nca(2,1),cinfo_op1a,nablk_op1)
      call sum_occ(nca(1,2),cinfo_op2c,ncblk_op2)
      call sum_occ(nca(2,2),cinfo_op2a,nablk_op2)
      call sum_occ(nca(1,3),cinfo_op1op2c,ncblk_op1op2)
      call sum_occ(nca(2,3),cinfo_op1op2a,nablk_op1op2)
      call sum_occ(nca(1,4),cinfo_ex1c, ncblk_ex1)
      call sum_occ(nca(2,4),cinfo_ex1a, nablk_ex1)
      call sum_occ(nca(1,5),cinfo_ex2c, ncblk_ex2)
      call sum_occ(nca(2,5),cinfo_ex2a, nablk_ex2)
      call sum_occ(nca(1,6),cinfo_cntc,ncblk_cnt)
      call sum_occ(nca(2,6),cinfo_cnta,nablk_cnt)

      allocate(
     &     gmex1dis_c(ncblk_ex1), gmex1dis_a(nablk_ex1),
     &     gmex2dis_c(ncblk_ex2), gmex2dis_a(nablk_ex2),
     &     gmc_dis_c(ncblk_cnt), gmc_dis_a(nablk_cnt),
     &     gmi_dis_c(ncblk_op1op2), gmi_dis_a(nablk_op1op2),
     &     msex1dis_c(ncblk_ex1), msex1dis_a(nablk_ex1),
     &     msex2dis_c(ncblk_ex2), msex2dis_a(nablk_ex2),
     &     msc_dis_c(ncblk_cnt), msc_dis_a(nablk_cnt),
     &     msi_dis_c(ncblk_op1op2), msi_dis_a(nablk_op1op2),
     &     idxmsex1dis_c(ncblk_ex1), idxmsex1dis_a(nablk_ex1),
     &     idxmsex2dis_c(ncblk_ex2), idxmsex2dis_a(nablk_ex2),
     &     idxmsc_dis_c(ncblk_cnt), idxmsc_dis_a(nablk_cnt),
     &     idxmsi_dis_c(ncblk_op1op2), idxmsi_dis_a(nablk_op1op2),
     &     lenex1(ncblk_ex1+nablk_ex1),
     &     lenex2(ncblk_ex2+nablk_ex2),
     &     lencnt(ncblk_cnt+nablk_cnt),
     &     lenop1op2(ncblk_op1op2+nablk_op1op2)
     &     )

      ! minimum Ms(A) for ...
      msbnd(1,1) = -nca(2,1) ! operator 1
      msbnd(1,2) = -nca(2,2) ! operator 2        
      msbnd(1,3) = -nca(2,3) ! intermediate
      ! maximum Ms(A) for ...
      msbnd(2,1) = -msbnd(1,1)
      msbnd(2,2) = -msbnd(1,2)
      msbnd(2,3) = -msbnd(1,3)
      ! max |Ms| for ...
      mscmx_a = nca(2,6) ! C(A)
      mscmx_c = nca(1,6) ! C(C)
      ! minimum IRREP for operators
      igambnd(1,1) = 1
      igambnd(1,2) = 1
      igambnd(1,3) = 1
      ! maximum IRREP 
      igambnd(2,1) = nsym
      igambnd(2,2) = nsym
      igambnd(2,3) = nsym
      ! loop Ms-cases of (Op1(A),Op2(A),Interm)
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
        if (abs(ms12i_c(1)).gt.nca(1,1)) cycle ms_loop
        if (abs(ms12i_c(2)).gt.nca(1,2)) cycle ms_loop
        if (abs(ms12i_c(3)).gt.nca(1,3)) cycle ms_loop

        msc_ac = ms12i_a(1) + ms12i_a(2) - ms12i_a(3)

        if (mscmx_a+mscmx_c.lt.abs(msc_ac)) cycle ms_loop
          
        msc_loop: do msc_a = mscmx_a, -mscmx_a, -2
          msc_c = msc_ac - msc_a
          if (abs(msc_c).gt.mscmx_c) cycle
          ! Ms of string after lifting restrictions:
          !  Op1(C1,A1) -> Op1(Cex1,Aex1;CC,AC)
          msex1_c = ms12i_c(1) - msc_c
          if (abs(msex1_c).gt.nca(1,4))
     &         cycle msc_loop
          msex1_a = ms12i_a(1) - msc_a
          if (abs(msex1_a).gt.nca(2,4))
     &         cycle msc_loop
          msex2_a = ms12i_a(2) - msc_c   ! other way 
          msex2_c = ms12i_c(2) - msc_a   ! round (!!)
          if (abs(msex2_c).gt.nca(1,5))
     &         cycle msc_loop
          if (abs(msex2_a).gt.nca(2,5))
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

              ! loop over distributions of current Ms and IRREP 
              ! of Aex1 and Cex1 over ngastypes
              first3 = .true.
              caex1_loop: do
                if (.not.next_msgamdist2(first3,
     &             msex1dis_c,msex1dis_a,gmex1dis_c,gmex1dis_a,
     &             ncblk_ex1, nablk_ex1,
     &             cinfo_ex1c,cinfo_ex1a,
     &             msex1_c,msex1_a,igamex1_c,igamex1_a,nsym,
     &             .false.,fix_success))
     &             exit caex1_loop
                first3 = .false.

                call ms2idxms(idxmsex1dis_c,msex1dis_c,
     &               cinfo_ex1c,ncblk_ex1)
                call ms2idxms(idxmsex1dis_a,msex1dis_a,
     &               cinfo_ex1a,nablk_ex1)

c                call set_len_str(lenex1,ncblk_ex1,nablk_ex1,
c     &                  graphs,
                call set_len_str2(lenex1,ncblk_ex1,nablk_ex1,
     &                  lenstr_array,nsym,maxidxms,
     &                  cinfo_ex1c(1,2),idxmsex1dis_c,
     &                                 gmex1dis_c,cinfo_ex1c(1,3),
     &                  cinfo_ex1a(1,2),idxmsex1dis_a,
     &                                 gmex1dis_a,cinfo_ex1a(1,3),
     &                  hpvxseq,.true.)
c dbg
c                print *,'nca_blk: ',ncblk_ex1,nablk_ex1
c                print *,'lenex1:   ',lenex1(1:ncblk_ex1+nablk_ex1)
c dbg

                if ( ncblk_ex1+nablk_ex1.gt.0 .and.
     &               idxlist(0,lenex1,ncblk_ex1+nablk_ex1,1).gt.0)
     &               cycle

                xlenex1 = 1d0
                do icmp = 1, ncblk_ex1+nablk_ex1
                  xlenex1 = xlenex1*dble(lenex1(icmp))
                end do

                ! loop over distributions of current Ms and IRREP 
                ! of Aex2 and Cex2 over ngastypes            
                first4 = .true.
                caex2_loop: do
                  if (.not.next_msgamdist2(first4,
     &               msex2dis_c,msex2dis_a,gmex2dis_c,gmex2dis_a,
     &               ncblk_ex2, nablk_ex2,
     &               cinfo_ex2c,cinfo_ex2a,
     &               msex2_c,msex2_a,igamex2_c,igamex2_a,nsym,
     &               .false.,fix_success))
     &               exit caex2_loop
                  first4 = .false.

                  call ms2idxms(idxmsex2dis_c,msex2dis_c,
     &                 cinfo_ex2c,ncblk_ex2)
                  call ms2idxms(idxmsex2dis_a,msex2dis_a,
     &                 cinfo_ex2a,nablk_ex2)

c                  call set_len_str(lenex2,ncblk_ex2,nablk_ex2,
c     &                 graphs,
                  call set_len_str2(lenex2,ncblk_ex2,nablk_ex2,
     &                 lenstr_array,nsym,maxidxms,
     &                 cinfo_ex2c(1,2),idxmsex2dis_c,
     &                                gmex2dis_c,cinfo_ex2c(1,3),
     &                 cinfo_ex2a(1,2),idxmsex2dis_a,
     &                                gmex2dis_a,cinfo_ex2a(1,3),
     &                 hpvxseq,.true.)

c dbg
c                  print *,'nca_blk: ',ncblk_ex2,nablk_ex2
c                  print *,'lenex2:   ',lenex2(1:ncblk_ex2+nablk_ex2)
c dbg
                  if ( ncblk_ex2+nablk_ex2.gt.0.and.
     &                 idxlist(0,lenex2,
     &                 ncblk_ex2+nablk_ex2,1).gt.0)
     &                 cycle

                  xlenex2 = 1d0
                  do icmp = 1, ncblk_ex2+nablk_ex2
                    xlenex2 = xlenex2*dble(lenex2(icmp))
                  end do

                  ! get Ms and IRREP distribution of intermediate
                  call merge_msgmdis(msi_dis_c,gmi_dis_c,
     &                                 ncblk_op1op2,
     &                                 msex1dis_c,gmex1dis_c,
     &                                 msex2dis_c,gmex2dis_c,
     &                                 map_info_12c)
                  call merge_msgmdis(msi_dis_a,gmi_dis_a,
     &                                 nablk_op1op2,
     &                                 msex2dis_a,gmex2dis_a,
     &                                 msex1dis_a,gmex1dis_a,
     &                                 map_info_12a)
c dbg
c                    print *,'gmi_dis_a: ',gmi_dis_a
c                    print *,'gmex1_dis_a:',gmex1dis_a
c                    print *,'gmex2_dis_a:',gmex2dis_a
c dbg

                  call ms2idxms(idxmsi_dis_c,msi_dis_c,
     &                   cinfo_op1op2c,ncblk_op1op2)
                  call ms2idxms(idxmsi_dis_a,msi_dis_a,
     &                   cinfo_op1op2a,nablk_op1op2)

                  ! length of intermediate
c                  call set_len_str(
c     &                   lenop1op2,ncblk_op1op2,nablk_op1op2,
c     &                   graphs,
                  call set_len_str2(
     &                   lenop1op2,ncblk_op1op2,nablk_op1op2,
     &                   lenstr_array,nsym,maxidxms,
     &                   cinfo_op1op2c(1,2),idxmsi_dis_c,
     &                                    gmi_dis_c,cinfo_op1op2c(1,3),
     &                   cinfo_op1op2a(1,2),idxmsi_dis_a,
     &                                    gmi_dis_a,cinfo_op1op2a(1,3),
     &                   hpvxseq,.true.)

                  xlenop1op2 = 1d0
                  do icmp = 1, ncblk_op1op2+nablk_op1op2
                    xlenop1op2 = xlenop1op2*dble(lenop1op2(icmp))
                  end do
c dbg
c                 print *,'xlenop1op2 ',xlenop1op2
c dbg
                    
                  xmemblk = max(xmemblk,xlenop1op2)

                  xmemtot = xmemtot+xlenop1op2

                  ! loop over distributions of current Ms and IRREP 
                  ! of AC and CC over ngastypes            
                  first5 = .true.
                  cac_loop: do
                    if (.not.next_msgamdist2(first5,
     &                 msc_dis_c,msc_dis_a,gmc_dis_c,gmc_dis_a,
     &                 ncblk_cnt, nablk_cnt,
     &                 cinfo_cntc,cinfo_cnta,
     &                 msc_c,msc_a,igamc_c,igamc_a,nsym,
     &                 .false.,fix_success))
     &                 exit cac_loop
                    first5 = .false.

                    call ms2idxms(idxmsc_dis_c,msc_dis_c,
     &                   cinfo_cntc,ncblk_cnt)
                    call ms2idxms(idxmsc_dis_a,msc_dis_a,
     &                   cinfo_cnta,nablk_cnt)

                    ! length of contraction
c                    call set_len_str(lencnt,ncblk_cnt,nablk_cnt,
c     &                  graphs,
                    call set_len_str2(lencnt,ncblk_cnt,nablk_cnt,
     &                  lenstr_array,nsym,maxidxms,
     &                  cinfo_cntc(1,2),idxmsc_dis_c,
     &                                  gmc_dis_c,cinfo_cntc(1,3),
     &                  cinfo_cnta(1,2),idxmsc_dis_a,
     &                                  gmc_dis_a,cinfo_cnta(1,3),
     &                  hpvxseq,.true.)

c dbg
c                    print *,'nca_blk: ',ncblk_cnt,nablk_cnt
c                    print *,'lencnt:   ',
c     &                   lencnt(1:ncblk_cnt+nablk_cnt)
c dbg
                    if ( ncblk_cnt+nablk_cnt.gt.0 .and.
     &                   idxlist(0,lencnt,
     &                   ncblk_cnt+nablk_cnt,1).gt.0)
     &                   cycle

                    xlencnt = 1d0
                    do icmp = 1, ncblk_cnt+nablk_cnt
                      xlencnt = xlencnt*dble(lencnt(icmp))
                    end do

                    ! Op1(C1,A1;Cc,Ac) Op2(C2,A2;Cc,Ac) -> Int(Ci,Ai)
                    ! so ....
                    flops = flops + xlenex1*xlenex2*xlencnt
c dbg
c                    print *,'xlenex1 = ',xlenex1
c                    print *,'xlenex2 = ',xlenex2
c                    print *,'xlencnt = ',xlencnt
c                    print *,'current blk:',xlenex1*xlenex2*xlencnt
c dbg

                  end do cac_loop
                end do caex2_loop
              end do caex1_loop

            end do gamc_loop
          end do gam_loop

        end do msc_loop
      end do ms_loop

      deallocate(
     &     gmex1dis_c, gmex1dis_a,
     &     gmex2dis_c, gmex2dis_a,
     &     gmc_dis_c, gmc_dis_a,
     &     gmi_dis_c, gmi_dis_a,
     &     msex1dis_c, msex1dis_a,
     &     msex2dis_c, msex2dis_a,
     &     msc_dis_c, msc_dis_a,
     &     msi_dis_c, msi_dis_a,
     &     idxmsex1dis_c, idxmsex1dis_a,
     &     idxmsex2dis_c, idxmsex2dis_a,
     &     idxmsc_dis_c, idxmsc_dis_a,
     &     idxmsi_dis_c, idxmsi_dis_a,
     &     lenex1, lenex2, lencnt, lenop1op2
     &     )

      deallocate(lenstr_array)

      return

      end
