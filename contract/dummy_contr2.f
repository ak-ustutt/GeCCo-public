*----------------------------------------------------------------------*
      subroutine dummy_contr2(flops,xmemtot,xmemblk,
     &     nca_blk,
     &     cinfo_op1c, cinfo_op1a, cinfo_op2c, cinfo_op2a,
     &     cinfo_op1op2c, cinfo_op1op2a,
     &     cinfo_10c, cinfo_10a, cinfo_20c, cinfo_20a,
     &     cinfo_cntc, cinfo_cnta,
     &     map_info_1c, map_info_1a,
     &     map_info_2c, map_info_2a,
     &     map_info_12c, map_info_12a,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     str_info,ngas,ihpvgas,nsym)
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
      include 'hpvxseq.h'

      integer, parameter ::
     &     ntest = 00

      real(8), intent(out) ::
     &     flops, xmemtot, xmemblk
      integer, intent(in) ::
     &     nca_blk(2,6),
     &     cinfo_op1c(nca(1,1),3), cinfo_op1a(nca_blk(2,1),3),
     &     cinfo_op2c(nca_blk(1,2),3), cinfo_op2a(nca_blk(2,2),3),
     &     cinfo_op1op2c(nca_blk(1,3),3), cinfo_op1op2a(nca_blk(2,3),3),
     &     cinfo_10c(nca_blk(1,4),3), cinfo_10a(nca_blk(2,4),3),
     &     cinfo_20c(nca_blk(1,5),3), cinfo_20a(nca_blk(2,5),3),
     &     cinfo_cntc(nca_blk(1,6),3), cinfo_cnta(nca_blk(2,6),3),
     &     map_info_1c(*), map_info_1a(*),
     &     map_info_2c(*), map_info_2a(*),
     &     map_info_12c(*), map_info_12a(*),
     &     mstop1,mstop2,mstop1op2,igamtop1,igamtop2,igamtop1op2
      integer, intent(in) ::
     &     ngas, nsym, ihpvgas(ngas)
      type(strinf), intent(in) ::
     &     str_info

      logical ::
     &     first1, first2, first3, first4, first5
      integer ::
     &     nca(2,6),
     &     mscmx_a, mscmx_c, msc_ac, msc_a, msc_c,
     &     ms10_a, ms10_c, ms20_a, ms20_c,
     &     igamc_ac, igamc_a, igamc_c,
     &     igam10_a, igam10_c, igam20_a, igam20_c,
     &     icmp
      integer ::
     &     irst_ext(2,ngas,2,2,2),irst_cnt(2,ngas,2,2),
     &     msbnd(2,3), igambnd(2,3),
     &     ms12i_a(3), ms12i_c(3), igam12i_a(3), igam12i_c(3),
     &     gm10dis_c(nca_blk(1,4)), gm10dis_a(nca_blk(2,4)),
     &     gm20dis_c(nca_blk(1,5)), gm20dis_a(nca_blk(2,5)),
     &     gmc_dis_c(nca_blk(1,6)), gmc_dis_a(nca_blk(2,6)),
     &     gmi_dis_c(nca_blk(1,3)), gmi_dis_a(nca_blk(2,3)),
     &     ms10dis_c(nca_blk(1,4)), ms10dis_a(nca_blk(2,4)),
     &     ms20dis_c(nca_blk(1,5)), ms20dis_a(nca_blk(2,5)),
     &     msc_dis_c(nca_blk(1,6)), msc_dis_a(nca_blk(2,6)),
     &     msi_dis_c(nca_blk(1,3)), msi_dis_a(nca_blk(2,3)),
     &     idxms10dis_c(nca_blk(1,4)), idxms10dis_a(nca_blk(2,4)),
     &     idxms20dis_c(nca_blk(1,5)), idxms20dis_a(nca_blk(2,5)),
     &     idxmsc_dis_c(nca_blk(1,6)), idxmsc_dis_a(nca_blk(2,6)),
     &     idxmsi_dis_c(nca_blk(1,3)), idxmsi_dis_a(nca_blk(2,3)),
     &     len10(nca_blk(1,4)+nca_blk(2,4)),
     &     len20(nca_blk(1,5)+nca_blk(2,5)),
     &     lencnt(nca_blk(1,6)+nca_blk(2,6)),
     &     lenop1op2(nca_blk(1,3)+nca_blk(2,3))
      type(graph), pointer ::
     &     graphs(:)
      real(8) ::
     &     xlen10, xlen20, xlencnt, xlenop1op2

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

      graphs => str_info%g
      
      call sum_occ(nca(1,1),cinfo_op1c,nca_blk(1,1))
      call sum_occ(nca(2,1),cinfo_op1a,nca_blk(2,1))
      call sum_occ(nca(1,2),cinfo_op2c,nca_blk(1,2))
      call sum_occ(nca(2,2),cinfo_op2a,nca_blk(2,2))
      call sum_occ(nca(1,3),cinfo_op1op2c,nca_blk(1,3))
      call sum_occ(nca(2,3),cinfo_op1op2a,nca_blk(2,3))
      call sum_occ(nca(1,4),cinfo_10c, nca_blk(1,4))
      call sum_occ(nca(2,4),cinfo_10a, nca_blk(2,4))
      call sum_occ(nca(1,5),cinfo_20c, nca_blk(1,5))
      call sum_occ(nca(2,5),cinfo_20a, nca_blk(2,5))
      call sum_occ(nca(1,6),cinfo_cntc,nca_blk(1,6))
      call sum_occ(nca(2,6),cinfo_cnta,nca_blk(2,6))

      ! minimum Ms for ...
      msbnd(1,1) = -nca(1,1) ! operator 1
      msbnd(1,2) = -nca(1,2) ! operator 2        
      msbnd(1,3) = -nca(1,3) ! intermediate
      ! maximum Ms for ...
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
        msc_ac = ms12i_a(1) + ms12i_a(2) - ms12i_a(3)

        if (mscmx_a+mscmx_c.lt.abs(msc_ac)) cycle ms_loop
          
        msc_loop: do msc_a = mscmx_a, -mscmx_a, -2
          msc_c = msc_ac - msc_a
          if (abs(msc_c).gt.mscmx_c) cycle
          ! Ms of string after lifting restrictions:
          !  Op1(C1,A1) -> Op1(C10,A10;CC,AC)
          ms10_c = ms12i_c(1) - msc_c
          if (abs(ms10_c).gt.nca(1,4))
     &         cycle msc_loop
          ms10_a = ms12i_a(1) - msc_a
          if (abs(ms10_a).gt.nca(2,4))
     &         cycle msc_loop
          ms20_a = ms12i_a(2) - msc_c   ! other way 
          ms20_c = ms12i_c(2) - msc_a   ! round (!!)
          if (abs(ms20_c).gt.nca(1,5))
     &         cycle msc_loop
          if (abs(ms20_a).gt.nca(2,5))
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
              igam10_a = multd2h(igam12i_a(1),igamc_a)
              igam10_c = multd2h(igam12i_c(1),igamc_c)
              igam20_a = multd2h(igam12i_a(2),igamc_c) !  !!
              igam20_c = multd2h(igam12i_c(2),igamc_a) !  !!

              ! loop over distributions of current Ms and IRREP 
              ! of A10 and C10 over ngastypes
              first3 = .true.
              ca10_loop: do
                if (.not.next_msgamdist2(first3,
     &             ms10dis_c,ms10dis_a,gm10dis_c,gm10dis_a,
     &             nca_blk(1,4), nca_blk(2,4),
     &             cinfo_10c,cinfo_10a,
     &             ms10_c,ms10_a,igam10_c,igam10_a,nsym)) exit ca10_loop
                first3 = .false.

                call ms2idxms(idxms10dis_c,ms10dis_c,
     &               cinfo_10c,nca_blk(1,4))
                call ms2idxms(idxms10dis_a,ms10dis_a,
     &               cinfo_10a,nca_blk(2,4))

                call set_len_str(len10,nca_blk(1,4),nca_blk(2,4),
     &                  graphs,
     &                  cinfo_10c(1,2),idxms10dis_c,
     &                                 gm10dis_c,cinfo_10c(1,3),
     &                  cinfo_10a(1,2),idxms10dis_a,
     &                                 gm10dis_a,cinfo_10a(1,3),
     &                  hpvxseq,.true.)
c dbg
c                print *,'nca_blk: ',nca_blk(1,4),nca_blk(2,4)
c                print *,'len10:   ',len10(1:nca_blk(1,4)+nca_blk(2,4))
c dbg

                if ( nca_blk(1,4)+nca_blk(2,4).gt.0 .and.
     &               idxlist(0,len10,nca_blk(1,4)+nca_blk(2,4),1).gt.0)
     &               cycle

                xlen10 = 1d0
                do icmp = 1, nca_blk(1,4)+nca_blk(2,4)
                  xlen10 = xlen10*dble(len10(icmp))
                end do

                ! loop over distributions of current Ms and IRREP 
                ! of A20 and C20 over ngastypes            
                first4 = .true.
                ca20_loop: do
                  if (.not.next_msgamdist2(first4,
     &               ms20dis_c,ms20dis_a,gm20dis_c,gm20dis_a,
     &               nca_blk(1,5), nca_blk(2,5),
     &               cinfo_20c,cinfo_20a,
     &               ms20_c,ms20_a,igam20_c,igam20_a,nsym))
     &               exit ca20_loop
                  first4 = .false.

                  call ms2idxms(idxms20dis_c,ms20dis_c,
     &                 cinfo_20c,nca_blk(1,5))
                  call ms2idxms(idxms20dis_a,ms20dis_a,
     &                 cinfo_20a,nca_blk(2,5))

                  call set_len_str(len20,nca_blk(1,5),nca_blk(2,5),
     &                 graphs,
     &                 cinfo_20c(1,2),idxms20dis_c,
     &                                gm20dis_c,cinfo_20c(1,3),
     &                 cinfo_20a(1,2),idxms20dis_a,
     &                                gm20dis_a,cinfo_20a(1,3),
     &                 hpvxseq,.true.)

c dbg
c                  print *,'nca_blk: ',nca_blk(1,5),nca_blk(2,5)
c                  print *,'len20:   ',len20(1:nca_blk(1,5)+nca_blk(2,5))
c dbg
                  if ( nca_blk(1,5)+nca_blk(2,5).gt.0.and.
     &                 idxlist(0,len20,
     &                 nca_blk(1,5)+nca_blk(2,5),1).gt.0)
     &                 cycle

                  xlen20 = 1d0
                  do icmp = 1, nca_blk(1,5)+nca_blk(2,5)
                    xlen20 = xlen20*dble(len20(icmp))
                  end do

                  ! get Ms and IRREP distribution of intermediate
                  call merge_msgmdis(msi_dis_c,gmi_dis_c,
     &                                 nca_blk(1,3),
     &                                 ms10dis_c,gm10dis_c,
     &                                 ms20dis_c,gm20dis_c,
     &                                 map_info_12c)
                  call merge_msgmdis(msi_dis_a,gmi_dis_a,
     &                                 nca_blk(2,3),
     &                                 ms20dis_a,gm20dis_a,
     &                                 ms10dis_a,gm10dis_a,
     &                                 map_info_12a)
c dbg
c                    print *,'gmi_dis_a: ',gmi_dis_a
c                    print *,'gm10_dis_a:',gm10dis_a
c                    print *,'gm20_dis_a:',gm20dis_a
c dbg

                  call ms2idxms(idxmsi_dis_c,msi_dis_c,
     &                   cinfo_op1op2c,nca_blk(1,3))
                  call ms2idxms(idxmsi_dis_a,msi_dis_a,
     &                   cinfo_op1op2a,nca_blk(2,3))

                  ! length of intermediate
                  call set_len_str(
     &                   lenop1op2,nca_blk(1,3),nca_blk(2,3),
     &                   graphs,
     &                   cinfo_op1op2c(1,2),idxmsi_dis_c,
     &                                    gmi_dis_c,cinfo_op1op2c(1,3),
     &                   cinfo_op1op2a(1,2),idxmsi_dis_a,
     &                                    gmi_dis_a,cinfo_op1op2a(1,3),
     &                   hpvxseq,.true.)

                  xlenop1op2 = 1d0
                  do icmp = 1, nca_blk(1,3)+nca_blk(2,3)
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
     &                 nca_blk(1,6), nca_blk(2,6),
     &                 cinfo_cntc,cinfo_cnta,
     &                 msc_c,msc_a,igamc_c,igamc_a,nsym))
     &                 exit cac_loop
                    first5 = .false.

                    call ms2idxms(idxmsc_dis_c,msc_dis_c,
     &                   cinfo_cntc,nca_blk(1,6))
                    call ms2idxms(idxmsc_dis_a,msc_dis_a,
     &                   cinfo_cnta,nca_blk(2,6))

                    ! length of contraction
                    call set_len_str(lencnt,nca_blk(1,6),nca_blk(2,6),
     &                  graphs,
     &                  cinfo_cntc(1,2),idxmsc_dis_c,
     &                                  gmc_dis_c,cinfo_cntc(1,3),
     &                  cinfo_cnta(1,2),idxmsc_dis_a,
     &                                  gmc_dis_a,cinfo_cnta(1,3),
     &                  hpvxseq,.true.)

c dbg
c                    print *,'nca_blk: ',nca_blk(1,6),nca_blk(2,6)
c                    print *,'lencnt:   ',
c     &                   lencnt(1:nca_blk(1,6)+nca_blk(2,6))
c dbg
                    if ( nca_blk(1,6)+nca_blk(2,6).gt.0 .and.
     &                   idxlist(0,lencnt,
     &                   nca_blk(1,6)+nca_blk(2,6),1).gt.0)
     &                   cycle

                    xlencnt = 1d0
                    do icmp = 1, nca_blk(1,6)+nca_blk(2,6)
                      xlencnt = xlencnt*dble(lencnt(icmp))
                    end do

                    ! Op1(C1,A1;Cc,Ac) Op2(C2,A2;Cc,Ac) -> Int(Ci,Ai)
                    ! so ....
                    flops = flops + xlen10*xlen20*xlencnt
c dbg
c                    print *,'xlen10 = ',xlen10
c                    print *,'xlen20 = ',xlen20
c                    print *,'xlencnt = ',xlencnt
c                    print *,'current blk:',xlen10*xlen20*xlencnt
c dbg

                  end do cac_loop
                end do ca20_loop
              end do ca10_loop

            end do gamc_loop
          end do gam_loop

        end do msc_loop
      end do ms_loop

      return

      end
