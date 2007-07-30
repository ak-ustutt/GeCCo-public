*----------------------------------------------------------------------*
      subroutine dummy_contr(flops,xmemtot,xmemblk,
     &                       iocc_op,iocc_ext,iocc_int,iocc_cnt,
     &                       irst_op,irst_int,
     &                       mstop,mstint,igamtop,igamtint,
     &                       str_info,ngas,ihpvgas,nsym)
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

      integer, parameter ::
     &     ntest = 00

      real(8), intent(out) ::
     &     flops, xmemtot, xmemblk
      integer, intent(in) ::
     &     iocc_op(ngastp,2,2), iocc_ext(ngastp,2,2),
     &     iocc_int(ngastp,2), iocc_cnt(ngastp,2),
     &     irst_op(2,ngas,2,2,2),irst_int(2,ngas,2,2), 
     &     mstop(2),mstint,igamtop(2),igamtint
      integer, intent(in) ::
     &     ngas, nsym, ihpvgas(ngas)
      type(strinf), intent(in) ::
     &     str_info

      logical ::
     &     first1, first2, first3, first4, first5
      integer ::
     &     mscmx_a, mscmx_c, msc_ac, msc_a, msc_c,
     &     ms10_a, ms10_c, ms20_a, ms20_c,
     &     igamc_ac, igamc_a, igamc_c,
     &     igam10_a, igam10_c, igam20_a, igam20_c,
     &     lenc1, lena1, lenc2, lena2, lencc, lenac, lenci, lenai
      integer ::
     &     irst_ext(2,ngas,2,2,2),irst_cnt(2,ngas,2,2),
     &     msbnd(2,3), igambnd(2,3),
     &     ms12i_a(3), ms12i_c(3), igam12i_a(3), igam12i_c(3),
     &     igam10dist(ngastp,2), igam20dist(ngastp,2),
     &     igamc_dist(ngastp,2), igami_dist(ngastp,2),
     &     ms10dist(ngastp,2), ms20dist(ngastp,2),
     &     msc_dist(ngastp,2), msi_dist(ngastp,2)

      integer, external ::
     &     ielsum
      logical, external ::
     &     next_dist, next_msgamdist

      if (ntest.ge.10) then
        write(luout,*) '========================'
        write(luout,*) ' here comes dummy_contr'
        write(luout,*) '========================'
        write(luout,*) ' ngas, nsym : ',ngas,nsym
        write(luout,*) 'OP 1 ',mstop(1),igamtop(1)
        call wrt_occ(luout,iocc_op(1,1,1))
        write(luout,*) 'OP 2 ',mstop(2),igamtop(2)
        call wrt_occ(luout,iocc_op(1,1,2))
        write(luout,*) 'CNT'
        call wrt_occ(luout,iocc_cnt)
        write(luout,*) 'RES ',mstint,igamtint
        call wrt_occ(luout,iocc_int)
      end if

      ! init
      flops = 0d0
      xmemtot = 0d0
      xmemblk = 0d0

      ! preliminary treatment of restrictions:
      call fit_restr(irst_ext(1,1,1,1,1),iocc_ext(1,1,1),
     &     irst_op(1,1,1,1,1),ihpvgas,ngas)
      call fit_restr(irst_ext(1,1,1,1,2),iocc_ext(1,1,2),
     &     irst_op(1,1,1,1,2),ihpvgas,ngas)
      call fit_restr(irst_cnt,iocc_cnt,
     &     irst_op(1,1,1,1,1),ihpvgas,ngas)
      ! end of preliminary code

      ! minimum Ms for ...
      msbnd(1,1) = -ielsum(iocc_op(1,1,1),ngastp) ! operator 1
      msbnd(1,2) = -ielsum(iocc_op(1,1,2),ngastp) ! operator 2        
      msbnd(1,3) = -ielsum(iocc_int,ngastp) ! intermediate
      ! maximum Ms for ...
      msbnd(2,1) = -msbnd(1,1)
      msbnd(2,2) = -msbnd(1,2)
      msbnd(2,3) = -msbnd(1,3)
      ! max |Ms| for ...
      mscmx_a = ielsum(iocc_cnt(1,2),ngastp) ! C(A)
      mscmx_c = ielsum(iocc_cnt(1,1),ngastp) ! C(C)
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

        ms12i_c(1) = ms12i_a(1) + mstop(1)
        ms12i_c(2) = ms12i_a(2) + mstop(2)
        ms12i_c(3) = ms12i_a(3) + mstint
        msc_ac = ms12i_a(1) + ms12i_a(2) - ms12i_a(3)

        if (mscmx_a+mscmx_c.lt.abs(msc_ac)) cycle ms_loop
          
        msc_loop: do msc_a = mscmx_a, -mscmx_a, -2
          msc_c = msc_ac - msc_a
          if (abs(msc_c).gt.mscmx_c) cycle
          ! Ms of string after lifting restrictions:
          !  Op1(C1,A1) -> Op1(C10,A10;CC,AC)
          ms10_c = ms12i_c(1) - msc_c
          if (abs(ms10_c).gt.ielsum(iocc_ext(1,1,1),ngastp))
     &         cycle msc_loop
          ms10_a = ms12i_a(1) - msc_a
          if (abs(ms10_a).gt.ielsum(iocc_ext(1,2,1),ngastp))
     &         cycle msc_loop
          ms20_a = ms12i_a(2) - msc_c   ! other way 
          ms20_c = ms12i_c(2) - msc_a   ! round (!!)
          if (abs(ms20_c).gt.ielsum(iocc_ext(1,1,2),ngastp))
     &         cycle msc_loop
          if (abs(ms20_a).gt.ielsum(iocc_ext(1,2,2),ngastp))
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

            igam12i_c(1) = multd2h(igam12i_a(1),igamtop(1))
            igam12i_c(2) = multd2h(igam12i_a(2),igamtop(2))
            igam12i_c(3) = multd2h(igam12i_a(3),igamtint)

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
                if (.not.next_msgamdist(first3,
     &             ms10_a,ms10_c,igam10_a,igam10_c,iocc_ext(1,1,1),nsym,
     &             ms10dist,igam10dist)) exit ca10_loop
                first3 = .false.

                call get_lenblk(lenc1,lena1,
     &               iocc_ext(1,1,1),irst_ext(1,1,1,1,1),
     &               ms10dist,igam10dist,
     &               str_info,ihpvgas,ngas)

                if (lenc1.eq.0.or.lena1.eq.0) cycle

                ! loop over distributions of current Ms and IRREP 
                ! of A20 and C20 over ngastypes            
                first4 = .true.
                ca20_loop: do
                  if (.not.next_msgamdist(first4,
     &               ms20_a,ms20_c,igam20_a,igam20_c,iocc_ext(1,1,2),
     &               nsym,
     &               ms20dist,igam20dist)) exit ca20_loop
                  first4 = .false.

                  call get_lenblk(lenc2,lena2,
     &                 iocc_ext(1,1,2),irst_ext(1,1,1,1,2),
     &                 ms20dist,igam20dist,
     &                 str_info,ihpvgas,ngas)

                  if (lenc2.eq.0.or.lena2.eq.0) cycle

                  ! loop over distributions of current Ms and IRREP 
                  ! of AC and CC over ngastypes            
                  first5 = .true.
                  cac_loop: do
                    if (.not.next_msgamdist(first5,
     &                 msc_a,msc_c,igamc_a,igamc_c,iocc_cnt,nsym,
     &                 msc_dist,igamc_dist)) exit cac_loop
                    first5 = .false.

                    ! length of contraction
                    call get_lenblk(lencc,lenac,
     &                   iocc_cnt,irst_cnt,msc_dist,igamc_dist,
     &                   str_info,ihpvgas,ngas)

                    if (lencc.eq.0.or.lenac.eq.0) cycle                    

                    ! Op1(C1,A1;Cc,Ac) Op2(C2,A2;Cc,Ac) -> Int(Ci,Ai)
                    ! so ....
                    flops = flops +
     &                   dble(lenc1)*dble(lena1)*
     &                   dble(lenc2)*dble(lena2)*
     &                   dble(lencc)*dble(lenac)

                    ! get Ms and IRREP distribution of intermediate
                    msi_dist(1:ngastp,1:2) =
     &                   ms10dist(1:ngastp,1:2) + ms20dist(1:ngastp,1:2)
                    call dirpr_gamdist(igami_dist,
     &                   igam10dist,igam20dist,ngastp*2)

                    call get_lenblk(lenci,lenai,
     &                   iocc_int,irst_int,msi_dist,igami_dist,
     &                   str_info,ihpvgas,ngas)
                    
                    xmemblk = max(xmemblk,dble(lenci)*dble(lenai))

                    xmemtot = xmemtot+dble(lenci)*dble(lenai)

                  end do cac_loop
                end do ca20_loop
              end do ca10_loop

            end do gamc_loop
          end do gam_loop

        end do msc_loop
      end do ms_loop

      return
      end
