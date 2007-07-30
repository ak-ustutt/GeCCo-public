*----------------------------------------------------------------------*
      subroutine contr_op1op2_simple(xfac,ffop1,ffop2,
     &     update,ffop1op2,xret,type_xret,
     &     op1,op2,op1op2,
     &     iblkop1,iblkop2,iblkop1op2,
     &     idoffop1,idoffop2,idoffop1op2,
     &     iocc_ext1,iocc_ext2,iocc_cnt,
     &     irst_op1,irst_op2,irst_op1op2,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     initial version:
*       full blocks incore
*       direct addressing
*       extremely transparent (I hope)
*       extremely slow (I fear)
*
*      update: update matrix elements on ffop1op2
*              else overwrite
*
*     andreas (not too proud), feb 2007
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
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'
      include 'ifc_operators.h'
      
      integer, parameter ::
     &     ntest = 000

      logical, intent(in) ::
     &     update
      real(8), intent(in) ::
     &     xfac
      real(8), intent(inout), target ::
     &     xret(1)
      integer, intent(in) ::
     &     type_xret,
     &     iblkop1, iblkop2, iblkop1op2,
     &     idoffop1,idoffop2,idoffop1op2,
     &     iocc_ext1(ngastp,2), iocc_ext2(ngastp,2), iocc_cnt(ngastp,2),
     &     irst_op1(*), irst_op2(*), irst_op1op2(*),
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2
      type(filinf), intent(inout) ::
     &     ffop1,ffop2,ffop1op2
      type(operator), intent(in) ::
     &     op1, op2, op1op2
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      logical ::
     &     bufop1, bufop2, bufop1op2,
     &     first1, first2, first3, first4, first5
      integer ::
     &     nsym, ifree,
     &     idxst_op1, idxst_op2, idxst_op1op2,
     &     idxop1, idxop2, idxop1op2,
     &     lenop1, lenop2, lenop1op2,
     &     na_op1, na_op2, na_op1op2,
     &     na_ext1, nc_ext1, na_ext2, nc_ext2, na_cnt, nc_cnt,
     &     mscmx_a, mscmx_c, msc_ac, msc_a, msc_c,
     &     ms10_a, ms10_c, ms20_a, ms20_c,
     &     igamc_ac, igamc_a, igamc_c,
     &     igam10_a, igam10_c, igam20_a, igam20_c,
     &     idxms, idxdis, maxac,
     &     isignsum
      real(8) ::
     &     xnrm, casign
      real(8) ::
     &     cpu, sys, cpu0, sys0, cpu00, sys00

      integer ::
     &     iocc_op1(ngastp,2), iocc_op2(ngastp,2), iocc_op1op2(ngastp,2)

      real(8), pointer ::
     &     xop1(:), xop2(:), xop1op2(:)
      real(8), pointer ::
     &     xbf1(:), xbf2(:), xbf12(:)

      integer ::
     &     irst_ext1(8*orb_info%ngas),irst_ext2(8*orb_info%ngas),
     &     irst_cnt(8*orb_info%ngas),
     &     msbnd(2,3), igambnd(2,3),
     &     ms12i_a(3), ms12i_c(3), igam12i_a(3), igam12i_c(3),
     &     igam10dist(ngastp,2), igam20dist(ngastp,2),
     &     igamc_dist(ngastp,2), igami_dist(ngastp,2),
     &     igamop1_dist(ngastp,2), igamop2_dist(ngastp,2),
     &     ms10dist(ngastp,2), ms20dist(ngastp,2),
     &     msc_dist(ngastp,2), msi_dist(ngastp,2),
     &     msop1_dist(ngastp,2), msop2_dist(ngastp,2),
     &     lstrext1(ngastp,2), lstrext2(ngastp,2), lstrcnt(ngastp,2),
     &     lstrop1(ngastp,2), lstrop2(ngastp,2), lstrop1op2(ngastp,2),
     &     igrphext1(ngastp,2), igrphext2(ngastp,2), igrphcnt(ngastp,2),
     &     igrphop1(ngastp,2), igrphop2(ngastp,2), igrphop1op2(ngastp,2)


      integer, external ::
     &     ielsum, ielprd, idx_msgmdst
      logical, external ::
     &     next_dist, next_msgamdist
      real(8), external ::
     &     ddot

      if (ntest.gt.0) then
        write(luout,*) '============================='
        write(luout,*) ' contr_op1op2_simple at work'
        write(luout,*) '============================='
      end if
      if (ntest.ge.10) then
        write(luout,*) 'ffop1:   ',trim(ffop1%name),' idoff: ',idoffop1
        write(luout,*) 'ffop2:   ',trim(ffop2%name),' idoff: ',idoffop2
        write(luout,*) 'ffop1op2:',trim(ffop1op2%name),
     &                                           ' idoff: ',idoffop1op2
        write(luout,*) 'xfac = ',xfac
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
      
      ! preliminary treatment of restrictions:
      call fit_restr(irst_ext1,iocc_ext1,
     &     irst_op1,orb_info%ihpvgas,orb_info%ngas)
      call fit_restr(irst_ext2,iocc_ext2,
     &     irst_op2,orb_info%ihpvgas,orb_info%ngas)
      call fit_restr(irst_cnt,iocc_cnt,
     &     irst_op1,orb_info%ihpvgas,orb_info%ngas)
      ! end of preliminary code

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

      if (.not.op1%dagger) then
        iocc_op1 = op1%ihpvca_occ(1:ngastp,1:2,iblkop1)
      else
        iocc_op1 = iocc_dagger(op1%ihpvca_occ(1:ngastp,1:2,iblkop1))
      end if
      if (.not.op2%dagger) then
        iocc_op2 = op2%ihpvca_occ(1:ngastp,1:2,iblkop2)
      else
        iocc_op2 = iocc_dagger(op2%ihpvca_occ(1:ngastp,1:2,iblkop2))
      end if
      if (iblkop1op2.gt.0) then
        if (.not.op1op2%dagger) then
          iocc_op1op2 = op1op2%ihpvca_occ(1:ngastp,1:2,iblkop1op2)
        else
          iocc_op1op2 = iocc_dagger(
     &         op1op2%ihpvca_occ(1:ngastp,1:2,iblkop1op2))
        end if
      else
        iocc_op1op2(1:ngastp,1:2) = 0
      end if

      ifree = mem_setmark('contr0')

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
     &                        idoffop1op2+idxst_op1op2-1+lenop1op2)
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

c dbg
c      xmerk = xop1op2(1)
c dbg

      na_op1 = ielsum(iocc_op1(1,2),ngastp)
      na_op2 = ielsum(iocc_op2(1,2),ngastp)
      na_op1op2 = ielsum(iocc_op1op2(1,2),ngastp)

      nc_ext1 = ielsum(iocc_ext1(1,1),ngastp)
      na_ext1 = ielsum(iocc_ext1(1,2),ngastp)
      nc_ext2 = ielsum(iocc_ext2(1,1),ngastp)
      na_ext2 = ielsum(iocc_ext2(1,2),ngastp)
      nc_cnt = ielsum(iocc_cnt(1,1),ngastp)
      na_cnt = ielsum(iocc_cnt(1,2),ngastp)
      maxac = max(na_ext1+na_ext2,nc_ext1+nc_ext2,
     &     na_ext1+na_cnt,na_ext2+nc_cnt,
     &     nc_ext1+nc_cnt,nc_ext2+na_cnt)

      ! sign from CA transpositions:
c      isignsum = nc_cnt*na_ext1+nc_cnt*nc_ext2+na_ext1*nc_ext2
c      casign = 1d0
c      if (mod(isignsum,2).eq.1) casign = -1d0
      isignsum = (na_ext1+nc_ext1+na_ext2+nc_ext2)*nc_cnt+
     &     na_ext1*nc_ext2
      casign = 1d0
      if (mod(isignsum,2).eq.1) casign = -1d0
c dbg
c      print *,'casign after CA: ',casign
c dbg

      ! minimum Ms for ...
      msbnd(1,1) = -ielsum(iocc_op1,ngastp) ! operator 1
      msbnd(1,2) = -ielsum(iocc_op2,ngastp) ! operator 2        
      msbnd(1,3) = -ielsum(iocc_op1op2,ngastp) ! product
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
          ms10_c = ms12i_c(1) - msc_c
          if (abs(ms10_c).gt.ielsum(iocc_ext1,ngastp))
     &         cycle msc_loop
          ms10_a = ms12i_a(1) - msc_a
          if (abs(ms10_a).gt.ielsum(iocc_ext1(1,2),ngastp))
     &         cycle msc_loop
          ms20_a = ms12i_a(2) - msc_c   ! other way 
          ms20_c = ms12i_c(2) - msc_a   ! round (!!)
          if (abs(ms20_c).gt.ielsum(iocc_ext2(1,1),ngastp))
     &         cycle msc_loop
          if (abs(ms20_a).gt.ielsum(iocc_ext2(1,2),ngastp))
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
            idxop1 = op1%off_op_gmo(iblkop1)%
     &                  gam_ms(igam12i_a(1),idxms) + 1
     &           - idxst_op1+1
            idxms = (na_op2-ms12i_a(2))/2 + 1
            idxop2 = op2%off_op_gmo(iblkop2)%
     &                  gam_ms(igam12i_a(2),idxms) + 1
     &           - idxst_op2+1
            idxms = (na_op1op2-ms12i_a(3))/2 + 1
            if (iblkop1op2.gt.0)
     &           idxop1op2 = op1op2%off_op_gmo(iblkop1op2)%
     &                  gam_ms(igam12i_a(3),idxms) + 1
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
              igam10_a = multd2h(igam12i_a(1),igamc_a)
              igam10_c = multd2h(igam12i_c(1),igamc_c)
              igam20_a = multd2h(igam12i_a(2),igamc_c) !  !!
              igam20_c = multd2h(igam12i_c(2),igamc_a) !  !!

              call atim_cs(cpu00,sys00)

              ! loop over distributions of current Ms and IRREP 
              ! of A10 and C10 over ngastypes
              first3 = .true.
              ca10_loop: do
                if (.not.next_msgamdist(first3,
     &             ms10_a,ms10_c,igam10_a,igam10_c,iocc_ext1(1,1),nsym,
     &             ms10dist,igam10dist)) exit ca10_loop
                first3 = .false.

                call get_lenblk_hpvca(lstrext1,igrphext1,
     &               iocc_ext1,irst_ext1,
     &               ms10dist,igam10dist,
     &               str_info,orb_info%ihpvgas,orb_info%ngas)
                
                ! test C and A separately to avoid overflow
                if (ielprd(lstrext1(1,1),ngastp).eq.0.or.
     &              ielprd(lstrext1(1,2),ngastp).eq.0) cycle

                ! loop over distributions of current Ms and IRREP 
                ! of A20 and C20 over ngastypes            
                first4 = .true.
                ca20_loop: do
                  if (.not.next_msgamdist(first4,
     &               ms20_a,ms20_c,igam20_a,igam20_c,iocc_ext2(1,1),
     &               nsym,
     &               ms20dist,igam20dist)) exit ca20_loop
                  first4 = .false.

                  call get_lenblk_hpvca(lstrext2,igrphext2,
     &                 iocc_ext2,irst_ext2,
     &                 ms20dist,igam20dist,
     &                 str_info,orb_info%ihpvgas,orb_info%ngas)

                  if (ielprd(lstrext2(1,1),ngastp).eq.0.or.
     &                ielprd(lstrext2(1,2),ngastp).eq.0) cycle

                  ! get Ms and IRREP distribution of intermediate
                  msi_dist(1:ngastp,1:2) =
     &                 ms10dist(1:ngastp,1:2) + ms20dist(1:ngastp,1:2)
                  call dirpr_gamdist(igami_dist,
     &                 igam10dist,igam20dist,ngastp*2)

                  call get_lenblk_hpvca(lstrop1op2,igrphop1op2,
     &                 iocc_op1op2,irst_op1op2,msi_dist,igami_dist,
     &                 str_info,orb_info%ihpvgas,orb_info%ngas)
                  
                  ! get distribution index
                  idxms = (na_op1op2-ms12i_a(3))/2 + 1
                  if (iblkop1op2.gt.0.and.
     &                 op1op2%off_op_gmox(iblkop1op2)%
     &                 ndis(igam12i_a(3),idxms).gt.1) then
                    idxdis =
     &                  idx_msgmdst(iblkop1op2,ms12i_a(3),igam12i_a(3),
     &                  msi_dist,igami_dist,.false.,op1op2,nsym)
                    idxop1op2 = op1op2%off_op_gmox(iblkop1op2)%
     &                   d_gam_ms(idxdis,igam12i_a(3),idxms) + 1
     &                   - idxst_op1op2+1
                  end if

                  ! loop over distributions of current Ms and IRREP 
                  ! of AC and CC over ngastypes            
                  first5 = .true.
                  cac_loop: do
                    if (.not.next_msgamdist(first5,
     &                 msc_a,msc_c,igamc_a,igamc_c,iocc_cnt,nsym,
     &                 msc_dist,igamc_dist)) exit cac_loop
                    first5 = .false.

                    ! length of contraction
                    call get_lenblk_hpvca(lstrcnt,igrphcnt,
     &                   iocc_cnt,irst_cnt,msc_dist,igamc_dist,
     &                   str_info,orb_info%ihpvgas,orb_info%ngas)

                    if (ielprd(lstrcnt(1,1),ngastp).eq.0.or.
     &                  ielprd(lstrcnt(2,1),ngastp).eq.0) cycle                    

                    ! get Ms and IRREP distribution of op1
                    msop1_dist(1:ngastp,1:2) =
     &                   ms10dist(1:ngastp,1:2) + msc_dist(1:ngastp,1:2)
                    call dirpr_gamdist(igamop1_dist,
     &                   igam10dist,igamc_dist,ngastp*2)
                    
                    call get_lenblk_hpvca(lstrop1,igrphop1,
     &                   iocc_op1,irst_op1,msop1_dist,igamop1_dist,
     &                   str_info,orb_info%ihpvgas,orb_info%ngas)

                  
                    ! get distribution index
                    idxms = (na_op1-ms12i_a(1))/2 + 1
                    if (op1%off_op_gmox(iblkop1)%
     &                   ndis(igam12i_a(1),idxms).gt.1) then
                      idxdis =
     &                   idx_msgmdst(iblkop1,ms12i_a(1),igam12i_a(1),
     &                   msop1_dist,igamop1_dist,.false.,op1,nsym)
                      idxop1 = op1%off_op_gmox(iblkop1)%
     &                     d_gam_ms(idxdis,igam12i_a(1),idxms) + 1
     &                     - idxst_op1+1
                    end if

                    xnrm = ddot(ielprd(lstrop1,2*ngastp),
     &                   xop1(idxop1),1,xop1(idxop1),1)
                    if (xnrm.lt.1d-28) cycle cac_loop

                    ! get Ms and IRREP distribution of op2
                    ! remember: CNT^+ !
                    msop2_dist(1:ngastp,1) =
     &                   ms20dist(1:ngastp,1) + msc_dist(1:ngastp,2)
                    msop2_dist(1:ngastp,2) =
     &                   ms20dist(1:ngastp,2) + msc_dist(1:ngastp,1)
                    call dirpr_gamdist(igamop2_dist(1,1),
     &                   igam20dist(1,1),igamc_dist(1,2),ngastp)
                    call dirpr_gamdist(igamop2_dist(1,2),
     &                   igam20dist(1,2),igamc_dist(1,1),ngastp)
                    
                    call get_lenblk_hpvca(lstrop2,igrphop2,
     &                   iocc_op2,irst_op2,msop2_dist,igamop2_dist,
     &                   str_info,orb_info%ihpvgas,orb_info%ngas)
                  
                    ! get distribution index
                    idxms = (na_op2-ms12i_a(2))/2 + 1
                    if (op2%off_op_gmox(iblkop2)%
     &                   ndis(igam12i_a(2),idxms).gt.1) then
                      idxdis =
     &                   idx_msgmdst(iblkop2,ms12i_a(2),igam12i_a(2),
     &                   msop2_dist,igamop2_dist,.false.,op2,nsym)
                      idxop2 = op2%off_op_gmox(iblkop2)%
     &                     d_gam_ms(idxdis,igam12i_a(2),idxms) + 1
     &                     - idxst_op2+1
                    end if

                    xnrm = ddot(ielprd(lstrop2,2*ngastp),
     &                   xop2(idxop2),1,xop2(idxop2),1)
                    if (xnrm.lt.1d-28) cycle cac_loop

                    call atim_cs(cpu0,sys0)

                    ! make the contraction for this block
c dbg:
c                    if (ntest.ge.100)
c     &                   write(luout,*) 'calling blk1blk2',
                        write(luout,*) 'calling blk1blk2',
     &                   lenop1,idxop1,
     &                   lenop2,idxop2,
     &                   lenop1op2,idxop1op2
                    call cntr_blk1blk2_0a
     &                   (xfac*casign,xop1op2(idxop1op2),
     &                                 xop1(idxop1),xop2(idxop2),
     &                   iocc_op1,iocc_op2,iocc_op1op2,
     &                   igrphop1,igrphop2,igrphop1op2,
     &                   lstrop1, lstrop2, lstrop1op2, 
     &                   iocc_ext1,iocc_ext2,iocc_cnt, 
     &                   irst_ext1,irst_ext2,irst_cnt, 
     &                   igrphext1,igrphext2,igrphcnt,
     &                   ms10dist,ms20dist,msc_dist, 
     &                   igam10dist,igam20dist,igamc_dist, 
     &                   na_ext1,nc_ext1,na_ext2,nc_ext2,na_cnt,nc_cnt,
     &                   maxac,
     &                   str_info,orb_info)
                    if (ntest.ge.100)
     &                   write(luout,*) 'after blk1blk2'

                    call atim_cs(cpu,sys)
                    cnt_kernel(1) = cnt_kernel(1)+cpu-cpu0
                    cnt_kernel(2) = cnt_kernel(2)+sys-sys0

                  end do cac_loop
                end do ca20_loop
              end do ca10_loop

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

      if (type_xret.eq.2) then
        xret(1) = xop1op2(1)
      else if (type_xret.eq.1) then
        xret(1) = ddot(lenop1op2,xop1op2,1,xop1op2,1)
      end if

      ! put result to disc
      if (.not.bufop1op2) then
        call put_vec(ffop1op2,xop1op2,idoffop1op2+idxst_op1op2,
     &                           idoffop1op2+idxst_op1op2-1+lenop1op2)
      end if

      ifree = mem_flushmark()

      if (ntest.ge.100) then
        if (type_xret.ne.0)
     &       write(luout,*) 'xret on exit = ',xret(1)
      end if

      return
      end
*----------------------------------------------------------------------*
      subroutine cntr_blk1blk2_0a(xfac,            !prefactor   
     &     xop1op2,xop1,xop2,                     !buffers: res,op1,op2
     &     iocc_op1,iocc_op2,iocc_op1op2,         ! occupations
     &     igrphop1,igrphop2,igrphop1op2,         ! which graphs needed
     &     lstrop1, lstrop2, lstrop1op2,          ! # of strings
     &     iocc_ex1,iocc_ex2,iocc_cnt,            ! occupations
     &     irst_ex1,irst_ex2,irst_cnt,            ! restrictions
     &     igrphex1,igrphex2,igrphcnt,            ! graphs
     &     msex1,   msex2,   mscnt,               ! Ms distrib
     &     igmex1, igmex2, igmcnt,             ! IRREP distrib
     &     na_ex1,nc_ex1,na_ex2,nc_ex2,na_cnt,nc_cnt,maxac,! # C/A
     &     str_info,orb_info)                     
*----------------------------------------------------------------------*
*     very first inner contraction routine
*     - with explicit addressing
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'hpvxseq.h'

      integer, parameter ::
     &     ntest = 00

      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     na_ex1, nc_ex1, na_ex2, nc_ex2, na_cnt, nc_cnt,maxac,
     &     iocc_ex1(ngastp,2), iocc_ex2(ngastp,2),
     &     iocc_cnt(ngastp,2), iocc_op1(ngastp,2),
     &     iocc_op2(ngastp,2), iocc_op1op2(ngastp,2),
     &     irst_ex1(2,orb_info%ngas,2,2), irst_ex2(2,orb_info%ngas,2,2),
     &     irst_cnt(2,orb_info%ngas,2,2),
     &     igrphex1(ngastp,2), igrphex2(ngastp,2), igrphcnt(ngastp,2),
     &     msex1(ngastp,2),    msex2(ngastp,2),    mscnt(ngastp,2),
     &     igmex1(ngastp,2),   igmex2(ngastp,2),   igmcnt(ngastp,2),
     &     igrphop1(ngastp,2), igrphop2(ngastp,2),igrphop1op2(ngastp,2),
     &     lstrop1(ngastp,2), lstrop2(ngastp,2), lstrop1op2(ngastp,2)
      real(8), intent(in) ::
     &     xop1(*), xop2(*), xfac
      real(8), intent(inout) ::
     &     xop1op2(*)

      logical ::
     &     firstex1a, firstex1c, firstex2a, firstex2c,
     &     firstcnta, firstcntc

      integer ::
     &     idxra(ngastp), idxrc(ngastp),
     &     idxop1a(ngastp), idxop2a(ngastp),
     &     idxop1c(ngastp), idxop2c(ngastp),
     &     idorbex1a(maxac), idspnex1a(maxac), idssex1a(maxac),
     &     idorbex2a(maxac), idspnex2a(maxac), idssex2a(maxac),
     &     idorbcnta(maxac), idspncnta(maxac), idsscnta(maxac),
     &     idorbex1c(maxac), idspnex1c(maxac), idssex1c(maxac),
     &     idorbex2c(maxac), idspnex2c(maxac), idssex2c(maxac),
     &     idorbcntc(maxac), idspncntc(maxac), idsscntc(maxac),
     &     idorbresa(maxac), idspnresa(maxac),
     &     idspcresa(maxac), idgamresa(maxac),
     &     idorbresc(maxac), idspnresc(maxac),
     &     idspcresc(maxac), idgamresc(maxac),
     &     idorbop1a(maxac), idspnop1a(maxac),
     &     idspcop1a(maxac), idgamop1a(maxac),
     &     idorbop1c(maxac), idspnop1c(maxac),
     &     idspcop1c(maxac), idgamop1c(maxac),
     &     idorbop2a(maxac), idspnop2a(maxac),
     &     idspcop2a(maxac), idgamop2a(maxac),
     &     idorbop2c(maxac), idspnop2c(maxac),
     &     idspcop2c(maxac), idgamop2c(maxac)
      integer ::
     &     lexlscrex1a(na_ex1,3),lexlscrex1c(nc_ex1,3),
     &     lexlscrex2a(na_ex2,3),lexlscrex2c(nc_ex2,3),
     &     lexlscrcnta(na_cnt,3),lexlscrcntc(nc_cnt,3)
      integer ::
     &     isgnra, isgnrc, isgnr, isgnop1a, isgnop1c,
     &     isgnop2a, isgnop2c, isgn0, idxop1, idxop2, idxr
      real(8) ::
     &     sgn

      logical, external ::
     &     next_tupel
      integer, external ::
     &     idx_str_blk, iordstr

c dbg
      integer ibop1, ibop2, ibop1op2, ielprd


c      print *,'maxac = ',maxac
c      print *,'lstrop1   ',lstrop1(1:ngastp,1)
c      print *,'          ',lstrop1(1:ngastp,2)
c      print *,'lstrop2   ',lstrop2(1:ngastp,1)
c      print *,'          ',lstrop2(1:ngastp,2)
c      print *,'lstrop1op2',lstrop1op2(1:ngastp,1)
c      print *,'          ',lstrop1op2(1:ngastp,2)
      ibop1 = ielprd(lstrop1,ngastp*2)
      ibop2 = ielprd(lstrop2,ngastp*2)
      ibop1op2 = ielprd(lstrop1op2,ngastp*2)
c      print *,'ibop1, ibop2, ibop1op2: ',ibop1, ibop2, ibop1op2
c      print *,'xop1:'
c      print '(x,5g12.6)',xop1(1:ibop1)
c      print *,'xop2:'
c      print '(x,5g12.6)',xop2(1:ibop2)
c      print *,'xop1op2: (entry)'
c      print '(x,5g12.6)',xop1op2(1:ibop1op2)
c dbg      
      firstex1a = .true.
c dbg
c      print *,'level ex1a'
c dbg
      ! loop over external A strings of operator 1
      do while(next_tupel(idorbex1a,idspnex1a,idssex1a,
     &     na_ex1,iocc_ex1(1,2),igrphex1(1,2),
     &     msex1(1,2),igmex1(1,2),firstex1a,
     &     str_info%igas_restr,
     &     orb_info%mostnd,orb_info%igamorb,
     &     orb_info%nsym,orb_info%ngas,orb_info%ngas_hpv,
     &                                         orb_info%idx_gas,
     &     hpvxseq,lexlscrex1a))
        firstex1a = .false.
c dbg
c        print *,'top A1'
c        if (na_ex1.gt.0) print *,'ex1a: ',idorbex1a(1:na_ex1)
c        if (na_ex1.gt.0) print *,'      ',idspnex1a(1:na_ex1)
c dbg

        ! loop over external A strings of operator 2
        firstex2a = .true.
c dbg
c        print *,'level ex2a'
c dbg
        do while(next_tupel(idorbex2a,idspnex2a,idssex2a,
     &       na_ex2,iocc_ex2(1,2),igrphex2(1,2),
     &       msex2(1,2),igmex2(1,2),firstex2a,
     &       str_info%igas_restr,
     &       orb_info%mostnd,orb_info%igamorb,
     &       orb_info%nsym,orb_info%ngas,orb_info%ngas_hpv,
     &                                         orb_info%idx_gas,
     &       hpvxseq,lexlscrex2a))
          firstex2a = .false.
c dbg
c          print *,'top A2'
c dbg

          ! join ext1 and ext2
          ! address on op1op2
          isgnra = iordstr(idorbresa,idspnresa,
     &         idorbex2a,idspnex2a,na_ex2,
     &         idorbex1a,idspnex1a,na_ex1)
c     &         idorbex1a,idspnex1a,na_ex1,
c     &         idorbex2a,idspnex2a,na_ex2)
          ! cycle, if non-existent
c dbg
c          if (isgnra.eq.0) print *,'cycling ra',isgnra
c dbg
          if (isgnra.eq.0) cycle

          call idx_tupel(idxra,
     &         iocc_op1op2(1,2),igrphop1op2(1,2),
     &         idspcresa,idorbresa,idspnresa,idgamresa,.true.,
     &         str_info,orb_info)

          firstex1c = .true.
          ! loop over external C strings of operator 1
c dbg
c          print *,'level ex1c'
c dbg
          do while(next_tupel(idorbex1c,idspnex1c,idssex1c,
     &         nc_ex1,iocc_ex1,igrphex1,
     &         msex1,igmex1,firstex1c,
     &         str_info%igas_restr,
     &         orb_info%mostnd,orb_info%igamorb,
     &         orb_info%nsym,orb_info%ngas,orb_info%ngas_hpv,
     &                                         orb_info%idx_gas,
     &         hpvxseq,lexlscrex1c))
            firstex1c = .false.
c dbg
c            print *,'top C1'
c dbg

            firstex2c = .true.
            ! loop over external C strings of operator 2
c dbg
c            print *,'level ex2c'
c dbg
            do while(next_tupel(idorbex2c,idspnex2c,idssex2c,
     &           nc_ex2,iocc_ex2,igrphex2,
     &           msex2,igmex2,firstex2c,
     &           str_info%igas_restr,
     &           orb_info%mostnd,orb_info%igamorb,
     &           orb_info%nsym,orb_info%ngas,orb_info%ngas_hpv,
     &                                         orb_info%idx_gas,
     &           hpvxseq,lexlscrex2c))
              firstex2c = .false.
c dbg
c              print *,'top C2'
c dbg

              isgnrc = iordstr(idorbresc,idspnresc,
     &             idorbex1c,idspnex1c,nc_ex1,
     &             idorbex2c,idspnex2c,nc_ex2)

              ! cycle, if non-existent
c dbg
c              if (isgnrc.eq.0) print *,'cycling rc',isgnrc
c dbg
              if (isgnrc.eq.0) cycle

              isgnr = isgnrc*isgnra
c dbg
c              print *,'nc_ex1,nc_ex2: ',nc_ex1,nc_ex2
c              print *,'na_ex1,na_ex2: ',na_ex1,na_ex2
c              print *,'bef idx_tupel',idorbresc(1:nc_ex1+nc_ex2)
c              print *,'bef idx_tupel',idspnresc(1:nc_ex1+nc_ex2)
c              call flush(6)
c dbg

              call idx_tupel(idxrc,
     &             iocc_op1op2,igrphop1op2,
     &             idspcresc,idorbresc,idspnresc,idgamresc,.true.,
     &             str_info,orb_info)

c dbg
c                print *,'bef idx_str_blk'
c                call flush(6)
c dbg
              idxr = idx_str_blk(idxrc,idxra,
     &             lstrop1op2,lstrop1op2(1,2),iocc_op1op2,hpvxseq)

c dbg
c                print *,'bef AC'
c                call flush(6)
c dbg

              firstcnta = .true.
              ! loop over contraction A string
c dbg
c              print *,'level cnta'
c dbg
              do while(next_tupel(idorbcnta,idspncnta,idsscnta,
     &             na_cnt,iocc_cnt(1,2),igrphcnt(1,2),
     &             mscnt(1,2),igmcnt(1,2),firstcnta,
     &             str_info%igas_restr,
     &             orb_info%mostnd,orb_info%igamorb,
     &             orb_info%nsym,orb_info%ngas,orb_info%ngas_hpv,
     &                                           orb_info%idx_gas,
     &             hpvxseq,lexlscrcnta))
                firstcnta = .false.
c dbg
c                print *,'top AC'
c                call flush(6)
c dbg
c dbg
c                print *,'1'
c                call flush(6)
c dbg
            
                ! join ext1 and cnt 
                ! address on op1
                isgnop1a = iordstr(idorbop1a,idspnop1a,
     &               idorbcnta,idspncnta,na_cnt,
     &               idorbex1a,idspnex1a,na_ex1)

c dbg
c                if (isgnop1a.eq.0) print *,'cycling op1a',isgnop1a
c dbg
                if (isgnop1a.eq.0) cycle

                ! join ext2 and cnt^+ 
                ! address on op2
                isgnop2c = iordstr(idorbop2c,idspnop2c,
     &               idorbcnta,idspncnta,na_cnt,
     &               idorbex2c,idspnex2c,nc_ex2)
c     &               idorbex2c,idspnex2c,nc_ex2,
c     &               idorbcnta,idspncnta,na_cnt)

c dbg
c                if (isgnop2c.eq.0) print *,'cycling op2c',isgnop2c
c dbg
                if (isgnop2c.eq.0) cycle

c dbg
c                print *,'this call?'
c                print *,iocc_op1(1:ngastp,2)
c dbg
                call idx_tupel(idxop1a,
     &               iocc_op1(1,2),igrphop1(1,2),
     &               idspcop1a,idorbop1a,idspnop1a,idgamop1a,.true.,
     &               str_info,orb_info)

c dbg
c                print *,'or this ?'
c                print *,iocc_op2(1:ngastp,1)
c dbg
                call idx_tupel(idxop2c,
     &               iocc_op2,igrphop2,
     &               idspcop2c,idorbop2c,idspnop2c,idgamop2c,.true.,
     &               str_info,orb_info)

                isgn0 = isgnr*isgnop1a*isgnop2c

                firstcntc = .true.
c dbg
c                print *,'level cntc'
c dbg
                do while(next_tupel(idorbcntc,idspncntc,idsscntc,
     &               nc_cnt,iocc_cnt,igrphcnt,
     &               mscnt,igmcnt,firstcntc,
     &               str_info%igas_restr,
     &               orb_info%mostnd,orb_info%igamorb,
     &               orb_info%nsym,orb_info%ngas,orb_info%ngas_hpv,
     &                                             orb_info%idx_gas,
     &               hpvxseq,lexlscrcntc))
                  firstcntc = .false.
c dbg
c              print *,'top CC'
c dbg

                  ! join ext1 and cnt 
                  ! address on op1
                  isgnop1c = iordstr(idorbop1c,idspnop1c,
     &                 idorbcntc,idspncntc,nc_cnt,
     &                 idorbex1c,idspnex1c,nc_ex1)

c dbg
c                  if (isgnop1c.eq.0) print *,'cycling op1c',isgnop1c
c dbg
                  if (isgnop1c.eq.0) cycle

                  ! join ext2 and cnt^+
                  ! address on op2
                  isgnop2a = iordstr(idorbop2a,idspnop2a,
     &                 idorbcntc,idspncntc,nc_cnt,
     &                 idorbex2a,idspnex2a,na_ex2)
c     &                 idorbex2a,idspnex2a,na_ex2,
c     &                 idorbcntc,idspncntc,nc_cnt)
                  
c dbg
c                  if (isgnop2a.eq.0) print *,'cycling op2a',isgnop2a
c dbg
                  if (isgnop2a.eq.0) cycle

                  sgn = xfac*dble(isgn0*isgnop1c*isgnop2a)

                  call idx_tupel(idxop1c,
     &                 iocc_op1,igrphop1,
     &                 idspcop1c,idorbop1c,idspnop1c,idgamop1c,.true.,
     &                 str_info,orb_info)
                  call idx_tupel(idxop2a,
     &                 iocc_op2(1,2),igrphop2(1,2),
     &                 idspcop2a,idorbop2a,idspnop2a,idgamop2a,.true.,
     &                 str_info,orb_info)

                  idxop1 = idx_str_blk(idxop1c,idxop1a,
     &                 lstrop1,lstrop1(1,2),iocc_op1,hpvxseq)
                  idxop2 = idx_str_blk(idxop2c,idxop2a,
     &                 lstrop2,lstrop2(1,2),iocc_op2,hpvxseq)
c dbg
c                  if (idxr.le.0) stop 'range 12'
c                  if (idxop1.le.0) stop 'range 1'
c                  if (idxop2.le.0) stop 'range 2'
c                  if (idxr.gt.ibop1op2) stop 'range 12e'
c                  if (idxop1.gt.ibop1) stop 'range 1e'
c                  if (idxop2.gt.ibop2) stop 'range 2e'
c dbg
c dbg
c                  if (idxr.eq.1)
c                  print '(a,3i5,4g15.8)','!!1',idxr,idxop1,idxop2,
c     &                 xop1op2(idxr),sgn,xop1(idxop1),xop2(idxop2)
c                  print '(a,6i5)',' ',isgnra,isgnrc,isgnop1a,
c     &                   isgnop1c,isgnop2a,isgnop2c
cc                  if (na_ex1+na_cnt.eq.2.and.
cc     &                na_ex2+nc_cnt.eq.2) then
c                    call prt_str_resort(luout,'OP1   :',
c     &                   nc_cnt,nc_ex1,na_cnt,na_ex1,
c     &                   idorbop1c,idorbop1a,idspnop1c,idspnop1a,
c     &                   idorbcntc,idorbcnta,idspncntc,idspncnta,
c     &                   idorbex1c,idorbex1a,idspnex1c,idspnex1a)
c                    print *,'signs: ',isgnop1c,isgnop1a
c                    call prt_str_resort(luout,'OP2   :',
c     &                   na_cnt,nc_ex2,nc_cnt,na_ex2,
c     &                   idorbop2c,idorbop2a,idspnop2c,idspnop2a,
c     &                   idorbcnta,idorbcntc,idspncnta,idspncntc,
c     &                   idorbex2c,idorbex2a,idspnex2c,idspnex2a)
c                    print *,'signs: ',isgnop2c,isgnop2a
c                    call prt_str_resort(luout,'OP1OP2:',
c     &                   nc_ex1,nc_ex2,na_ex2,na_ex1,
c     &                   idorbresc,idorbresa,idspnresc,idspnresa,
c     &                   idorbex1c,idorbex2a,idspnex1c,idspnex2a,
c     &                   idorbex2c,idorbex1a,idspnex2c,idspnex1a)
cc     &                   idorbex1c,idorbex1a,idspnex1c,idspnex1a,
cc     &                   idorbex2c,idorbex2a,idspnex2c,idspnex2a)
c                    print *,'signs: ',isgnrc,isgnra
cc                  end if
c dbg                  

                  ! xfac contained in sgn
                  xop1op2(idxr) = xop1op2(idxr)
     &                          + sgn * xop1(idxop1)
     &                                * xop2(idxop2)
c dbg
c                  if (idxr.eq.1)
c                  print *,'!!2',xop1op2(idxr)
c dbg                  
c dbg
c                  print *,'bottom CC'
c dbg
                end do                
              end do

            end do
          end do

        end do
      end do
        
c dbg                  
c      print *,'xop1op2: (exit)'
c      print '(x,5g12.6)',xop1op2(1:ibop1op2)
c dbg                  

      return
      end


* debug routines
      subroutine prt_str_resort(luout,label,
     &                   nc_ex,nc_cnt,na_ex,na_cnt,
     &                   idorbopc,idorbopa,idspnopc,idspnopa,
     &                   idorbexc,idorbexa,idspnexc,idspnexa,
     &                   idorbcntc,idorbcnta,idspncntc,idspncnta)
      
      implicit none

      integer, intent(in) ::
     &     luout,
     &     nc_ex,nc_cnt,na_ex,na_cnt,
     &     idorbopc(nc_ex+nc_cnt),idorbopa(na_ex+na_cnt),
     &     idspnopc(nc_ex+nc_cnt),idspnopa(na_ex+na_cnt),
     &     idorbexc(nc_ex),idorbexa(na_ex),
     &     idspnexc(nc_ex),idspnexa(na_ex),
     &     idorbcntc(nc_cnt),idorbcnta(na_cnt),
     &     idspncntc(nc_cnt),idspncnta(na_cnt)
      character, intent(in) ::
     &     label*(*)


      integer ::
     &     iel, ioff, ipos
      character ::
     &     spnstr*(nc_ex+nc_cnt+na_ex+na_cnt+1),
     &     spnstr_ex*(nc_ex+na_ex+1),
     &     spnstr_cnt*(nc_cnt+na_cnt+1),fmtstr*128

c      if (nc_ex.gt.0.and.na_ex.gt.0.and.nc_cnt.gt.0.and.na_cnt.gt.0)then
        write(fmtstr,'("(x,a,",i2,"i3,x,",i2,"i3,x,a,a,",'//
     &                     'i2,"i3,x,",i3,"i3,x,a,a,",'//
     &                     'i2,"i3,x,",i3,"i3,x,a)")')
     &     nc_ex+nc_cnt,na_ex+na_cnt,nc_ex,na_ex,nc_cnt,na_cnt

        ipos = 1
        do while(fmtstr(ipos:ipos).ne.")")
          if (fmtstr(ipos:ipos).eq."0") then
            fmtstr(ipos:ipos+4) = ' '
            ipos = ipos+4
          else
            ipos = ipos+1
          end if
        end do
c      else if (nc_ex.gt.0.and.na_ex.gt.0.and.nc_cnt.gt.0.and.na_cnt.gt.0)

      spnstr    (1:nc_ex+nc_cnt+na_ex+na_cnt+1) = ' '
      spnstr_ex (1:nc_ex+na_ex+1) = ' '
      spnstr_cnt(1:nc_cnt+na_cnt+1) = ' '

      do iel = 1, nc_ex+nc_cnt
        if (idspnopc(iel).eq.1)  spnstr(iel:iel)='+'
        if (idspnopc(iel).eq.-1) spnstr(iel:iel)='-'
        if (idspnopc(iel).eq.2)  spnstr(iel:iel)='2'
      end do
      ioff = nc_ex+nc_cnt+1
      do iel = 1, na_ex+na_cnt
        if (idspnopa(iel).eq.1)  spnstr(ioff+iel:ioff+iel)='+'
        if (idspnopa(iel).eq.-1) spnstr(ioff+iel:ioff+iel)='-'
        if (idspnopa(iel).eq.2)  spnstr(ioff+iel:ioff+iel)='2'
      end do

      do iel = 1, nc_ex
        if (idspnexc(iel).eq.1)  spnstr_ex(iel:iel)='+'
        if (idspnexc(iel).eq.-1) spnstr_ex(iel:iel)='-'
        if (idspnexc(iel).eq.2)  spnstr_ex(iel:iel)='2'
      end do
      ioff = nc_ex+1
      do iel = 1, na_ex
        if (idspnexa(iel).eq.1)  spnstr_ex(ioff+iel:ioff+iel)='+'
        if (idspnexa(iel).eq.-1) spnstr_ex(ioff+iel:ioff+iel)='-'
        if (idspnexa(iel).eq.2)  spnstr_ex(ioff+iel:ioff+iel)='2'
      end do

      do iel = 1,nc_cnt
        if (idspncntc(iel).eq.1)  spnstr_cnt(iel:iel)='+'
        if (idspncntc(iel).eq.-1) spnstr_cnt(iel:iel)='-'
        if (idspncntc(iel).eq.2)  spnstr_cnt(iel:iel)='2'
      end do
      ioff = nc_cnt+1
      do iel = 1, na_cnt
        if (idspncnta(iel).eq.1)  spnstr_cnt(ioff+iel:ioff+iel)='+'
        if (idspncnta(iel).eq.-1) spnstr_cnt(ioff+iel:ioff+iel)='-'
        if (idspncnta(iel).eq.2)  spnstr_cnt(ioff+iel:ioff+iel)='2'
      end do

c      print *,'fmt = ',fmtstr
      write(luout,fmtstr) label,
     &     idorbopc(1:nc_ex+nc_cnt),idorbopa(1:na_ex+na_cnt),spnstr,
     &     ' -> ',
     &     idorbexc(1:nc_ex),idorbexa(1:na_ex),spnstr_ex,
     &     ' * ',
     &     idorbcntc(1:nc_cnt),idorbcnta(1:na_cnt),spnstr_cnt
      
      return
      end
