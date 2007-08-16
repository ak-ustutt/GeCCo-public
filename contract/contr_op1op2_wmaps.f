*----------------------------------------------------------------------*
      subroutine contr_op1op2_wmaps(xfac,ffop1,ffop2,
     &     update,ffop1op2,xret,type_xret,
     &     op1,op2,op1op2,
     &     iblkop1,iblkop2,iblkop1op2,
     &     idoffop1,idoffop2,idoffop1op2,
     &     iocc_ext1,iocc_ext2,iocc_cnt,
     &     irst_op1,irst_op2,irst_op1op2,
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

      integer, parameter ::
     &     ntest = 00

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
      type(strmapinf), intent(inout) ::
     &     strmap_info
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
     &     idxms, idxdis,
     &     isignsum, lenmap, hpvx, hpvx2, ica
c     &     ngastp_op1op2a, ngastp_op1op2c,
c     &     ngastp_op1a,    ngastp_op1c,
c     &     ngastp_op2a,    ngastp_op2c,
c     &     nstrop1a, nstrop1c, nstrop2a, nstrop2c,
c     &     nstrop1op2a, nstrop1op2c,
c     &     nstrex1a, nstrex1c, nstrex2a, nstrex2c,
c     &     nstrcnta, nstrcntc
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

      integer, pointer ::
     &     map_ex1ex2a(:), map_ex1ex2c(:),
     &     map_ex1cnta(:), map_ex1cntc(:),
     &     map_ex2cnta(:), map_ex2cntc(:)


      integer, external ::
     &     ielsum, ielprd, idx_msgmdst
      logical, external ::
     &     next_dist, next_msgamdist
      real(8), external ::
     &     ddot

      if (ntest.gt.0) then
        write(luout,*) '============================'
        write(luout,*) ' contr_op1op2_wmaps at work'
        write(luout,*) '============================'
      end if
      if (ntest.ge.10) then
        write(luout,*) 'ffop1:   ',ffop1%name(1:len_trim(ffop1%name))
        write(luout,*) 'ffop2:   ',ffop2%name(1:len_trim(ffop2%name))
        write(luout,*) 'ffop1op2:',
     &       ffop1op2%name(1:len_trim(ffop1op2%name))
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

      na_op1 = ielsum(iocc_op1(1,2),ngastp)
      na_op2 = ielsum(iocc_op2(1,2),ngastp)
      na_op1op2 = ielsum(iocc_op1op2(1,2),ngastp)

      nc_ext1 = ielsum(iocc_ext1(1,1),ngastp)
      na_ext1 = ielsum(iocc_ext1(1,2),ngastp)
      nc_ext2 = ielsum(iocc_ext2(1,1),ngastp)
      na_ext2 = ielsum(iocc_ext2(1,2),ngastp)
      nc_cnt = ielsum(iocc_cnt(1,1),ngastp)
      na_cnt = ielsum(iocc_cnt(1,2),ngastp)

      casign = 1d0

      ! sign from CA transpositions:
      isignsum = (na_ext1+nc_ext1+na_ext2+nc_ext2)*nc_cnt+
     &     na_ext1*nc_ext2
      casign = 1d0
      if (mod(isignsum,2).eq.1) casign = -1d0
c dbg
c      print *,'ext1,ext2,cnt:'
c      call wrt_occ(luout,iocc_ext1)
c      call wrt_occ(luout,iocc_ext2)
c      call wrt_occ(luout,iocc_cnt)
c      print '(a,f6.2,a,3i4,a,2i4,a)',
c     &     'casign 1: ',casign,'(',na_ext1+nc_ext1,
c     &     na_ext2+nc_ext2,nc_cnt,'|',
c     &     na_ext1,nc_ext2,')'
c dbg
      ! sign for HPVX transpositions:
      isignsum = 0
      do ica = 1, 2
        do hpvx = 2, ngastp
          if (iocc_cnt(hpvx,ica).gt.0) then
            do hpvx2 = 1, hpvx-1
              isignsum = isignsum 
     &                 + iocc_cnt(hpvx,ica)*iocc_ext1(hpvx2,ica)
     &                 + iocc_cnt(hpvx,ica)*iocc_ext2(hpvx2,3-ica)
            end do
          end if
c          if (iocc_ext2(hpvx,ica).gt.0) then
c            do hpvx2 = 1, hpvx-1
c              isignsum = isignsum 
c     &                 + iocc_ext2(hpvx,ica)*iocc_ext1(hpvx2,ica)
c            end do
c          end if
c          if (iocc_ext1(hpvx,ica).gt.0) then
c            do hpvx2 = 1, hpvx-1
c              isignsum = isignsum 
c     &                 + iocc_ext1(hpvx,ica)*iocc_ext2(hpvx2,ica)
c            end do
c          end if
        end do
      end do
      do hpvx = 2, ngastp
        ! ordering: ex1c,ex2c, but ex2a,ex1a
        if (iocc_ext1(hpvx,1).gt.0) then
          do hpvx2 = 1, hpvx-1
            isignsum = isignsum 
     &           + iocc_ext1(hpvx,1)*iocc_ext2(hpvx2,1)
          end do
        end if
        if (iocc_ext2(hpvx,2).gt.0) then
          do hpvx2 = 1, hpvx-1
            isignsum = isignsum 
     &           + iocc_ext2(hpvx,2)*iocc_ext1(hpvx2,2)
          end do
        end if
      end do

      if (mod(isignsum,2).eq.1) casign = casign*(-1d0)
c dbg
c      print *,'casign = ',casign
c dbg

      ! get graph indices
      call get_grph4occ(igrphext1,iocc_ext1,irst_ext1,
     &     str_info,orb_info%ihpvgas,orb_info%ngas,1,.true.)
      call get_grph4occ(igrphext2,iocc_ext2,irst_ext2,
     &     str_info,orb_info%ihpvgas,orb_info%ngas,1,.true.)
      call get_grph4occ(igrphcnt,iocc_cnt,irst_cnt,
     &     str_info,orb_info%ihpvgas,orb_info%ngas,1,.true.)
      call get_grph4occ(igrphop1,iocc_op1,irst_op1,
     &     str_info,orb_info%ihpvgas,orb_info%ngas,1,.true.)
      call get_grph4occ(igrphop2,iocc_op2,irst_op2,
     &     str_info,orb_info%ihpvgas,orb_info%ngas,1,.true.)
      call get_grph4occ(igrphop1op2,iocc_op1op2,irst_op1op2,
     &     str_info,orb_info%ihpvgas,orb_info%ngas,1,.true.)
      ! set up maps (if necessary)
      call strmap_man(
     &     igrphext1,.false.,
     &     igrphext2,.false.,
     &     igrphop1op2,.false.,
     &     str_info,strmap_info,orb_info)
      ! get ex2,ex1 sequence, as well
      call strmap_man(
     &     igrphext2,.false.,
     &     igrphext1,.false.,
     &     igrphop1op2,.false.,
     &     str_info,strmap_info,orb_info)
      call strmap_man(
     &     igrphcnt,.false.,
     &     igrphext1,.false.,
     &     igrphop1,.false.,
     &     str_info,strmap_info,orb_info)
      call strmap_man(
     &     igrphcnt,.true.,
     &     igrphext2,.false.,
     &     igrphop2,.false.,
     &     str_info,strmap_info,orb_info)

c      ifree = mem_alloc_int(mapbuf_ex1ex2,maxblk_ex1ex2,'mapbuf12')
c      ifree = mem_alloc_int(mapbuf_ex1cnt,maxblk_ex1cnt,'mapbuf1c')
c      ifree = mem_alloc_int(mapbuf_ex2cnt,maxblk_ex2cnt,'mapbuf2c')

      ! number of GAS types that are non-zero in occupation
c      ngastp_op1op2c = ngastp - imltlist(0,iocc_op1op2(1,1),ngastp,1)
c      ngastp_op1op2a = ngastp - imltlist(0,iocc_op1op2(1,2),ngastp,1)
c      ngastp_op1c    = ngastp - imltlist(0,iocc_op1(1,1),ngastp,1)
c      ngastp_op1a    = ngastp - imltlist(0,iocc_op1(1,2),ngastp,1)
c      ngastp_op2c    = ngastp - imltlist(0,iocc_op2(1,1),ngastp,1)
c      ngastp_op2a    = ngastp - imltlist(0,iocc_op2(1,2),ngastp,1)

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
                  
                  ! get igrphext1,igrphext2->igrphop1op2 map
                  ! for given ms and irreps
                  ifree = mem_setmark('ex_str')
                  lenmap = 0
                  do hpvx = 1, ngastp
                    lenmap = lenmap + lstrext2(hpvx,1)*lstrext1(hpvx,1)
                  end do
                  ifree = mem_alloc_int(map_ex1ex2c,lenmap,'strmap_c')
                  ! for C: ex1,ex2 sequence !
                  call get_strmap_blk(map_ex1ex2c,
     &                 iocc_ext1,iocc_ext2,lstrext1(1,1),lstrext2(1,1),
     &                 igrphext1,igrphext2,
     &                 ms10dist,ms20dist,
     &                 igam10dist,igam20dist,
     &                 strmap_info,nsym,str_info%ngraph)

                  lenmap = 0
                  do hpvx = 1, ngastp
                    lenmap = lenmap + lstrext2(hpvx,2)*lstrext1(hpvx,2)
                  end do
                  ifree = mem_alloc_int(map_ex1ex2a,lenmap,'strmap_a')
                  ! for A: ex2,ex1 sequence !
                  call get_strmap_blk(map_ex1ex2a,
     &                 iocc_ext2(1,2),iocc_ext1(1,2),
     &                 lstrext2(1,2),lstrext1(1,2),
     &                 igrphext2(1,2),igrphext1(1,2),
     &                 ms20dist(1,2),ms10dist(1,2),
     &                 igam20dist(1,2),igam10dist(1,2),
     &                 strmap_info,nsym,str_info%ngraph)

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
     &                  ielprd(lstrcnt(1,2),ngastp).eq.0) cycle cac_loop                   

                    ! get Ms and IRREP distribution of op1
                    msop1_dist(1:ngastp,1:2) =
     &                   ms10dist(1:ngastp,1:2) + msc_dist(1:ngastp,1:2)
                    call dirpr_gamdist(igamop1_dist,
     &                   igam10dist,igamc_dist,ngastp*2)
                    
                    call get_lenblk_hpvca(lstrop1,igrphop1,
     &                   iocc_op1,irst_op1,msop1_dist,igamop1_dist,
     &                   str_info,orb_info%ihpvgas,orb_info%ngas)

                    ! get igrphext1,igrphcnt->igrphop1 map
                    ! for given ms and irreps
                  
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

                    ! get igrphext1,igrphcnt->igrphop2 map
                    ! for given ms and irreps
                  
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

                    ! get igrphcnt,igrphext1->igrphop1 map
                    ! for given ms and irreps
                    ifree = mem_setmark('cntstr')
                    lenmap = 0
                    do hpvx = 1, ngastp
                      lenmap = lenmap + lstrcnt(hpvx,1)*lstrext1(hpvx,1)
                    end do
                    ifree = mem_alloc_int(map_ex1cntc,lenmap,'strmap_c')
                    call get_strmap_blk(map_ex1cntc,
     &                   iocc_cnt,iocc_ext1,lstrcnt(1,1),lstrext1(1,1),
     &                   igrphcnt,igrphext1,
     &                   msc_dist,ms10dist,
     &                   igamc_dist,igam10dist,
     &                   strmap_info,nsym,str_info%ngraph)

                    lenmap = 0
                    do hpvx = 1, ngastp
                      lenmap = lenmap + lstrcnt(hpvx,2)*lstrext1(hpvx,2)
                    end do
                    ifree = mem_alloc_int(map_ex1cnta,lenmap,'strmap_a')
                    call get_strmap_blk(map_ex1cnta,
     &                   iocc_cnt(1,2),iocc_ext1(1,2),
     &                   lstrcnt(1,2),lstrext1(1,2),
     &                   igrphcnt(1,2),igrphext1(1,2),
     &                   msc_dist(1,2),ms10dist(1,2),
     &                   igamc_dist(1,2),igam10dist(1,2),
     &                   strmap_info,nsym,str_info%ngraph)

                    ! get igrphcnt,igrphext2->igrphop2 map
                    ! for given ms and irreps
                    lenmap = 0
                    do hpvx = 1, ngastp
                      lenmap = lenmap + lstrcnt(hpvx,2)*lstrext2(hpvx,1)
                    end do
                    ifree = mem_alloc_int(map_ex2cntc,lenmap,'strmap_c')
                    call get_strmap_blk(map_ex2cntc,
     &                   iocc_cnt(1,2),iocc_ext2,
     &                   lstrcnt(1,2),lstrext2,
     &                   igrphcnt(1,2),igrphext2,
     &                   msc_dist(1,2),ms20dist,
     &                   igamc_dist(1,2),igam20dist,
     &                   strmap_info,nsym,str_info%ngraph)

                    lenmap = 0
                    do hpvx = 1, ngastp
                      lenmap = lenmap + lstrcnt(hpvx,1)*lstrext2(hpvx,2)
                    end do                    
                    ifree = mem_alloc_int(map_ex2cnta,lenmap,'strmap_a')
                    call get_strmap_blk(map_ex2cnta,
     &                   iocc_cnt(1,1),iocc_ext2(1,2),
     &                   lstrcnt(1,1),lstrext2(1,2),
     &                   igrphcnt(1,1),igrphext2(1,2),
     &                   msc_dist(1,1),ms20dist(1,2),
     &                   igamc_dist(1,1),igam20dist(1,2),
     &                   strmap_info,nsym,str_info%ngraph)                    

                    call atim_cs(cpu0,sys0)                    

                    ! make the contraction for this block
                    if (ntest.ge.100)
     &                   write(luout,*) 'calling blk1blk2',
     &                   lenop1,idxop1,
     &                   lenop2,idxop2,
     &                   lenop1op2,idxop1op2
                    if (ntest.ge.1000) then
                      write(luout,*) ' the maps:'
                      write(luout,*) ' X1X2(A): '
                      call prt_strmap(map_ex1ex2a,
     &                     iocc_ext2(1,2),iocc_ext1(1,2),
     &                     lstrext2(1,2),lstrext1(1,2))
                      write(luout,*) ' X1X2(C): '
                      call prt_strmap(map_ex1ex2c,
     &                     iocc_ext1(1,1),iocc_ext2(1,1),
     &                     lstrext1(1,1),lstrext2(1,1))
                      write(luout,*) ' X1C(A): '
                      call prt_strmap(map_ex1cnta,
     &                     iocc_cnt(1,2),iocc_ext1(1,2),
     &                     lstrcnt(1,2),lstrext1(1,2))
                      write(luout,*) ' X1C(C): '
                      call prt_strmap(map_ex1cntc,
     &                     iocc_cnt(1,1),iocc_ext1(1,1),
     &                     lstrcnt(1,1),lstrext1(1,1))
                      write(luout,*) ' X2C(A): '
                      call prt_strmap(map_ex2cnta,
     &                     iocc_cnt(1,1),iocc_ext2(1,2),
     &                     lstrcnt(1,1),lstrext2(1,2))
                      write(luout,*) ' X2C(C): '
                      call prt_strmap(map_ex2cntc,
     &                     iocc_cnt(1,2),iocc_ext2(1,1),
     &                     lstrcnt(1,2),lstrext2(1,1))
                    end if
c dbg
c                    print *,'on call:'
c                    print *,'xop1: ',xop1(idxop1)
c                    print *,'xop2: ',xop2(idxop2)
c                    print *,'xop1op2: ',xop1op2(idxop1op2)
c dbg
                    call cntr_blk1blk2_wmaps(xfac*casign,
     &                   xop1op2(idxop1op2),
     &                                 xop1(idxop1),xop2(idxop2),
     &                   iocc_op1,iocc_op2,iocc_op1op2,        
     &                   lstrop1(1,2), lstrop1(1,1),
     &                       lstrop2(1,2), lstrop2(1,1),
     &                   lstrop1op2(1,2), lstrop1op2(1,1),
     &                   lstrext1(1,2), lstrext1(1,1),
     &                       lstrext2(1,2), lstrext2(1,1),
     &                   lstrcnt(1,2), lstrcnt(1,1),
     &                   map_ex1ex2a, map_ex1ex2c,
     &                   map_ex1cnta, map_ex1cntc,
     &                   map_ex2cnta, map_ex2cntc
     &                   )                     
                    if (ntest.ge.100)
     &                   write(luout,*) 'after blk1blk2'
c dbg
c                    print *,'xop1op2: ',xop1op2(idxop1op2)                    
c dbg

                    call atim_cs(cpu,sys)
                    cnt_kernel(1) = cnt_kernel(1)+cpu-cpu0
                    cnt_kernel(2) = cnt_kernel(2)+sys-sys0

                    ifree = mem_flushmark('cntstr')

                  end do cac_loop
                  
                  ifree = mem_flushmark('ex_str')
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
      subroutine cntr_blk1blk2_wmaps(xfac,            !prefactor   
     &     xop1op2,xop1,xop2,                     !buffers: res,op1,op2
     &     iocc_op1,iocc_op2,iocc_op1op2,         ! rest: see below
     &     lstrop1a, lstrop1c, lstrop2a, lstrop2c,
     &     lstrop1op2a, lstrop1op2c,
     &     lstr_ex1a, lstr_ex1c, lstr_ex2a, lstr_ex2c,
     &     lstr_cnta, lstr_cntc,
     &     map_ex1ex2a, map_ex1ex2c,
     &     map_ex1cnta, map_ex1cntc,
     &     map_ex2cnta, map_ex2cntc
     &     )                     
*----------------------------------------------------------------------*
*     inner contraction routine, generation 2
*     - testing maps
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
     &     iocc_op1(ngastp,2), iocc_op2(ngastp,2),iocc_op1op2(ngastp,2),
     &     lstrop1a(ngastp), lstrop1c(ngastp),
     &     lstrop2a(ngastp), lstrop2c(ngastp),
     &     lstrop1op2a(ngastp), lstrop1op2c(ngastp),
     &     lstr_ex1a(ngastp), lstr_ex1c(ngastp),
     &     lstr_ex2a(ngastp), lstr_ex2c(ngastp),
     &     lstr_cnta(ngastp), lstr_cntc(ngastp),
     &     map_ex1ex2a(*), map_ex1ex2c(*),
     &     map_ex1cnta(*), map_ex1cntc(*),
     &     map_ex2cnta(*), map_ex2cntc(*)
      real(8), intent(in) ::
     &     xop1(*), xop2(*), xfac
      real(8), intent(inout) ::
     &     xop1op2(*)

      integer ::
     &     ielmap, ielmap2,
     &     idxop1op2a(ngastp), idxop1op2c(ngastp),
c     &     idxop1a(ngastp), idxop2a(ngastp),
c     &     idxop1c(ngastp), idxop2c(ngastp),
     &     nstr_ex1a1(ngastp), nstr_ex1c1(ngastp),
     &     nstr_ex1a12(ngastp), nstr_ex1c12(ngastp),
     &     nstr_ex2a2(ngastp), nstr_ex2c2(ngastp),
     &     nstr_ex2a12(ngastp), nstr_ex2c12(ngastp),
     &     nstr_cnta1(ngastp), nstr_cntc1(ngastp),
     &     nstr_cnta2(ngastp), nstr_cntc2(ngastp),
     &     nstr_ex1ex2a(ngastp), nstr_ex1ex2c(ngastp),
     &     nstr_ex1cnta(ngastp), nstr_ex1cntc(ngastp),
     &     nstr_ex2cnta(ngastp), nstr_ex2cntc(ngastp),
     &     ldim_op1a(ngastp), ldim_op1c(ngastp),
     &     ldim_op2a(ngastp), ldim_op2c(ngastp),
     &     ldim_op1op2a(ngastp), ldim_op1op2c(ngastp)
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
     &     idx_str_blk2, ielprd

      nstr_ex1c_tot = ielprd(lstr_ex1c,ngastp)
      nstr_ex1a_tot = ielprd(lstr_ex1a,ngastp)
      nstr_ex2c_tot = ielprd(lstr_ex2c,ngastp)
      nstr_ex2a_tot = ielprd(lstr_ex2a,ngastp)
      nstr_cntc_tot = ielprd(lstr_cntc,ngastp)
      nstr_cnta_tot = ielprd(lstr_cnta,ngastp)

c dbg
c      print *,'nstr_ex1c_tot',nstr_ex1c_tot
c      print *,'nstr_ex1a_tot',nstr_ex1a_tot
c      print *,'nstr_ex2c_tot',nstr_ex2c_tot
c      print *,'nstr_ex2a_tot',nstr_ex2a_tot
c      print *,'nstr_cntc_tot',nstr_cntc_tot
c      print *,'nstr_cnta_tot',nstr_cnta_tot
c dbg

      ! C: ex1,ex2
      call set_strmapdim(nstr_ex1ex2c,nstr_ex1c12,nstr_ex2c12,
     &                   ngastp_op1op2c,
     &                   iocc_op1op2(1,1),
     &                   lstr_ex1c,lstr_ex2c)
      ! A: ex2,ex1
      call set_strmapdim(nstr_ex1ex2a,nstr_ex2a12,nstr_ex1a12,
     &                   ngastp_op1op2a,
     &                   iocc_op1op2(1,2),
     &                   lstr_ex2a,lstr_ex1a)
      call set_strmapdim(nstr_ex1cntc,nstr_cntc1,nstr_ex1c1,
     &                   ngastp_op1c,
     &                   iocc_op1(1,1),
     &                   lstr_cntc,lstr_ex1c)
      call set_strmapdim(nstr_ex1cnta,nstr_cnta1,nstr_ex1a1,
     &                   ngastp_op1a,
     &                   iocc_op1(1,2),
     &                   lstr_cnta,lstr_ex1a)
      call set_strmapdim(nstr_ex2cntc,nstr_cnta2,nstr_ex2c2,
     &                   ngastp_op2c,
     &                   iocc_op2(1,1),
     &                   lstr_cnta,lstr_ex2c)
      call set_strmapdim(nstr_ex2cnta,nstr_cntc2,nstr_ex2a2,
     &                   ngastp_op2a,
     &                   iocc_op2(1,2),
     &                   lstr_cntc,lstr_ex2a)
c dbg
c      print *,'ngastp_op1:',ngastp_op1c,ngastp_op1a
c      print *,'ngastp_op2:',ngastp_op2c,ngastp_op2a
c      print *,'ngastp_op1op2:',ngastp_op1op2c,ngastp_op1op2a
c dbg

      call set_op_ldim(ldim_op1c,ldim_op1a,
     &                 iocc_op1,iocc_op1(1,2),lstrop1c,lstrop1a)
      call set_op_ldim(ldim_op2c,ldim_op2a,
     &                 iocc_op2,iocc_op2(1,2),lstrop2c,lstrop2a)
      call set_op_ldim(ldim_op1op2c,ldim_op1op2a,
     &                 iocc_op1op2,iocc_op1op2(1,2),
     &                 lstrop1op2c,lstrop1op2a)
c dbg
c      print *,'ldim_op1c:',ldim_op1c(1:ngastp_op1c)
c      print *,'ldim_op1a:',ldim_op1a(1:ngastp_op1a)
c      print *,'ldim_op2c:',ldim_op2c(1:ngastp_op2c)
c      print *,'ldim_op2a:',ldim_op2a(1:ngastp_op2a)
c      print *,'ldim_op1op2c:',ldim_op1op2c(1:ngastp_op1op2c)
c      print *,'ldim_op1op2a:',ldim_op1op2a(1:ngastp_op1op2a)
c dbg

      if (ntest.ge.100) then
        write(luout,*) '============================'
        write(luout,*) 'News from contr_op1op2_wmaps'
        write(luout,*) '============================'
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

              idxop1op2 = idx_str_blk2(idxop1op2c,idxop1op2a,
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
      subroutine cntr_blk1blk2_wmaps_old(xfac,            !prefactor   
     &     xop1op2,xop1,xop2,                     !buffers: res,op1,op2
     &     iocc_op1,iocc_op2,iocc_op1op2,         ! rest: see below
     &     lstrop1a, lstrop1c, lstrop2a, lstrop2c,
     &     lstrop1op2a, lstrop1op2c,
     &     lstr_ex1a, lstr_ex1c, lstr_ex2a, lstr_ex2c,
     &     lstr_cnta, lstr_cntc,
     &     map_ex1ex2a, map_ex1ex2c,
     &     map_ex1cnta, map_ex1cntc,
     &     map_ex2cnta, map_ex2cntc
     &     )                     
*----------------------------------------------------------------------*
*     inner contraction routine, generation 2
*     - testing maps
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
     &     iocc_op1(ngastp,2), iocc_op2(ngastp,2),iocc_op1op2(ngastp,2),
     &     lstrop1a(ngastp), lstrop1c(ngastp),
     &     lstrop2a(ngastp), lstrop2c(ngastp),
     &     lstrop1op2a(ngastp), lstrop1op2c(ngastp),
     &     lstr_ex1a(ngastp), lstr_ex1c(ngastp),
     &     lstr_ex2a(ngastp), lstr_ex2c(ngastp),
     &     lstr_cnta(ngastp), lstr_cntc(ngastp),
     &     map_ex1ex2a(*), map_ex1ex2c(*),
     &     map_ex1cnta(*), map_ex1cntc(*),
     &     map_ex2cnta(*), map_ex2cntc(*)
      real(8), intent(in) ::
     &     xop1(*), xop2(*), xfac
      real(8), intent(inout) ::
     &     xop1op2(*)

      integer ::
     &     ielmap, ielmap2,
     &     idxop1op2a(ngastp), idxop1op2c(ngastp),
     &     idxop1a(ngastp), idxop2a(ngastp),
     &     idxop1c(ngastp), idxop2c(ngastp),
     &     nstr_ex1a1(ngastp), nstr_ex1c1(ngastp),
     &     nstr_ex1a12(ngastp), nstr_ex1c12(ngastp),
     &     nstr_ex2a2(ngastp), nstr_ex2c2(ngastp),
     &     nstr_ex2a12(ngastp), nstr_ex2c12(ngastp),
     &     nstr_cnta1(ngastp), nstr_cntc1(ngastp),
     &     nstr_cnta2(ngastp), nstr_cntc2(ngastp),
     &     nstr_ex1ex2a(ngastp), nstr_ex1ex2c(ngastp),
     &     nstr_ex1cnta(ngastp), nstr_ex1cntc(ngastp),
     &     nstr_ex2cnta(ngastp), nstr_ex2cntc(ngastp),
     &     ldim_op1a(ngastp), ldim_op1c(ngastp),
     &     ldim_op2a(ngastp), ldim_op2c(ngastp),
     &     ldim_op1op2a(ngastp), ldim_op1op2c(ngastp)
      integer ::
     &     ioff, istr1, istr2, idx1, idx2, idx,
     &     ngastp_op1op2a, ngastp_op1op2c,
     &     ngastp_op1a,    ngastp_op1c,
     &     ngastp_op2a,    ngastp_op2c,
     &     nstr_ex1a_tot, nstr_ex1c_tot,
     &     nstr_ex2a_tot, nstr_ex2c_tot,
     &     nstr_cnta_tot, nstr_cntc_tot,
     &     isgnra, isgnrc, isgnr, isgnop1a, isgnop1c,
     &     isgnop2a, isgnop2c, isgn0, idxop1, idxop2, idxop1op2
      integer ::
     &     idx_ex1ex2a, idx_ex1ex2c,
     &     idx_ex1cnta, idx_ex1cntc, idx_ex1cnta0, idx_ex1cntc0,
     &     idx_ex2cnta, idx_ex2cntc, idx_ex2cnta0, idx_ex2cntc0,
     &     istr_ex1a, istr_ex1c, istr_ex2a, istr_ex2c,
     &     istr_cnta, istr_cntc, icmp
      real(8) ::
     &     sgn

      integer, external ::
     &     idx_str_blk2, ielprd

      nstr_ex1c_tot = ielprd(lstr_ex1c,ngastp)
      nstr_ex1a_tot = ielprd(lstr_ex1a,ngastp)
      nstr_ex2c_tot = ielprd(lstr_ex2c,ngastp)
      nstr_ex2a_tot = ielprd(lstr_ex2a,ngastp)
      nstr_cntc_tot = ielprd(lstr_cntc,ngastp)
      nstr_cnta_tot = ielprd(lstr_cnta,ngastp)

      ! C: ex1,ex2
      call set_strmapdim(nstr_ex1ex2c,nstr_ex1c12,nstr_ex2c12,
     &                   ngastp_op1op2c,
     &                   iocc_op1op2(1,1),
     &                   lstr_ex1c,lstr_ex2c)
      ! A: ex2,ex1
      call set_strmapdim(nstr_ex1ex2a,nstr_ex2a12,nstr_ex1a12,
     &                   ngastp_op1op2a,
     &                   iocc_op1op2(1,2),
     &                   lstr_ex2a,lstr_ex1a)
      call set_strmapdim(nstr_ex1cntc,nstr_cntc1,nstr_ex1c1,
     &                   ngastp_op1c,
     &                   iocc_op1(1,1),
     &                   lstr_cntc,lstr_ex1c)
      call set_strmapdim(nstr_ex1cnta,nstr_cnta1,nstr_ex1a1,
     &                   ngastp_op1a,
     &                   iocc_op1(1,2),
     &                   lstr_cnta,lstr_ex1a)
      call set_strmapdim(nstr_ex2cntc,nstr_cnta2,nstr_ex2c2,
     &                   ngastp_op2c,
     &                   iocc_op2(1,1),
     &                   lstr_cnta,lstr_ex2c)
      call set_strmapdim(nstr_ex2cnta,nstr_cntc2,nstr_ex2a2,
     &                   ngastp_op2a,
     &                   iocc_op2(1,2),
     &                   lstr_cntc,lstr_ex2a)

      call set_op_ldim(ldim_op1c,ldim_op1a,
     &                 iocc_op1,iocc_op1(1,2),lstrop1c,lstrop1a)
      call set_op_ldim(ldim_op2c,ldim_op2a,
     &                 iocc_op2,iocc_op2(1,2),lstrop2c,lstrop2a)
      call set_op_ldim(ldim_op1op2c,ldim_op1op2a,
     &                 iocc_op1op2,iocc_op1op2(1,2),
     &                 lstrop1op2c,lstrop1op2a)

      if (ntest.ge.100) then
        write(luout,*) '============================'
        write(luout,*) 'News from contr_op1op2_wmaps'
        write(luout,*) '============================'
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

              idxop1op2 = idx_str_blk2(idxop1op2c,idxop1op2a,
     &                                 ldim_op1op2c,ldim_op1op2a,
     &                                 ngastp_op1op2c,ngastp_op1op2a)

              idx_ex1cnta = (istr_ex1a-1)*nstr_cnta_tot
              idx_ex2cntc = (istr_ex2c-1)*nstr_cnta_tot
              ! loop over contraction A string
              cnta: do istr_cnta = 1, nstr_cnta_tot
                idx_ex1cnta = idx_ex1cnta+1
                idx_ex2cntc = idx_ex2cntc+1
                ioff = 0
                isgnop1a = 1
                istr1 = istr_ex1a-1
                istr2 = istr_cnta-1
                ! get sign and idxop1a
                do icmp = 1, ngastp_op1a
                  idx1 = mod(istr1,nstr_ex1a1(icmp))+1
                  idx2 = mod(istr2,nstr_cnta1(icmp))+1
                  idx  = (idx1-1)*nstr_cnta1(icmp)+idx2
                  ielmap = map_ex1cnta(ioff+idx)
                  if (ielmap.eq.0) cycle cnta
                  isgnop1a = isgnop1a*sign(1,ielmap)
                  idxop1a(icmp) = abs(ielmap)-1
                  ioff = ioff + nstr_ex1cnta(icmp)
                  istr1 = istr1/nstr_ex1a1(icmp)
                  istr2 = istr2/nstr_cnta1(icmp)
                end do

                ioff = 0
                isgnop2c = 1
                istr1 = istr_ex2c-1
                istr2 = istr_cnta-1
                ! get sign and idxop2c
                do icmp = 1, ngastp_op2c
                  idx1 = mod(istr1,nstr_ex2c2(icmp))+1
                  idx2 = mod(istr2,nstr_cnta2(icmp))+1
                  idx  = (idx1-1)*nstr_cnta2(icmp)+idx2
                  ielmap = map_ex2cntc(ioff+idx)
                  if (ielmap.eq.0) cycle cnta
                  isgnop2c = isgnop2c*sign(1,ielmap)
                  idxop2c(icmp) = abs(ielmap)-1
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
                  ! get sign and idxop1c
                  do icmp = 1, ngastp_op1c
                    idx1 = mod(istr1,nstr_ex1c1(icmp))+1
                    idx2 = mod(istr2,nstr_cntc1(icmp))+1
                    idx  = (idx1-1)*nstr_cntc1(icmp)+idx2
                    ielmap = map_ex1cntc(ioff+idx)
                    if (ielmap.eq.0) cycle cntc
                    isgnop1c = isgnop1c*sign(1,ielmap)
                    idxop1c(icmp) = abs(ielmap)-1
                    ioff = ioff + nstr_ex1cntc(icmp)
                    istr1 = istr1/nstr_ex1c1(icmp)
                    istr2 = istr2/nstr_cntc1(icmp)
                  end do

                  ioff = 0
                  isgnop2a = 1
                  istr1 = istr_ex2a-1
                  istr2 = istr_cntc-1
                  ! get sign and idxop2a
                  do icmp = 1, ngastp_op2a
                    idx1 = mod(istr1,nstr_ex2a2(icmp))+1
                    idx2 = mod(istr2,nstr_cntc2(icmp))+1
                    idx  = (idx1-1)*nstr_cntc2(icmp)+idx2
                    ielmap = map_ex2cnta(ioff+idx)
                    if (ielmap.eq.0) cycle cntc
                    isgnop2a = isgnop2a*sign(1,ielmap)
                    idxop2a(icmp) = abs(ielmap)-1
                    ioff = ioff + nstr_ex2cnta(icmp)
                    istr1 = istr1/nstr_ex2a2(icmp)
                    istr2 = istr2/nstr_cntc2(icmp)
                  end do

                  sgn = xfac*dble(isgn0*isgnop1c*isgnop2a)

                  idxop1 = idx_str_blk2(idxop1c,idxop1a,
     &                                  ldim_op1c,ldim_op1a,
     &                                  ngastp_op1c,ngastp_op1a)
                  idxop2 = idx_str_blk2(idxop2c,idxop2a,
     &                                  ldim_op2c,ldim_op2a,
     &                                  ngastp_op2c,ngastp_op2a)

                  ! xfac contained in sgn
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
      pure integer function idx_str_blk2(idxc,idxa,ldimc,ldima,
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
      idx_str_blk2 = 1

      if ((ngastp_c*ngastp_a).eq.0) return

      idx_str_blk2 = idxc(1)*ldimc(1)+idxa(1)*ldima(1) + 1
      if ((ngastp_c*ngastp_a).eq.1) return

      do idx = 2, ngastp_c
        idx_str_blk2 = idx_str_blk2+idxc(idx)*ldimc(idx)
      end do
      do idx = 2, ngastp_a
        idx_str_blk2 = idx_str_blk2+idxa(idx)*ldima(idx)
      end do

      return
      end

      subroutine set_strmapdim(nstr1str2,nstr1,nstr2,
     &                         ngastp_nz,
     &                         iocc12,nstr1occ,nstr2occ)

      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(out) ::
     &     nstr1str2(*), nstr1(*), nstr2(*), ngastp_nz
      integer, intent(in) ::
     &     iocc12(ngastp), nstr1occ(ngastp), nstr2occ(ngastp)
            
      integer ::
     &     igastp, hpvx

      ngastp_nz = 0
      do igastp = 1, ngastp
        hpvx = hpvxseq(igastp)
        if (iocc12(hpvx).eq.0) cycle
        ngastp_nz = ngastp_nz+1
        nstr1(ngastp_nz) = nstr1occ(hpvx)
        nstr2(ngastp_nz) = nstr2occ(hpvx)
        nstr1str2(ngastp_nz) = nstr1occ(hpvx)*nstr2occ(hpvx)
      end do

      return
      end

      subroutine set_op_ldim(ldimc,ldima,iocc_c,iocc_a,nstrc,nstra)

      implicit none
      
      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(out) ::
     &     ldimc(*), ldima(*)
      integer, intent(in) ::
     &     iocc_c(ngastp), iocc_a(ngastp), nstrc(ngastp), nstra(ngastp)

      integer ::
     &     idxa, idxc, ldim, igastp, hpvx

      idxa = 0
      idxc = 0
      ldim = 1
      do igastp = 1, ngastp
        hpvx = hpvxseq(igastp)
        if (iocc_c(hpvx).gt.0) then
          idxc = idxc+1
          ldimc(idxc) = ldim
          ldim = ldim * nstrc(hpvx)
        end if
        if (iocc_a(hpvx).gt.0) then
          idxa = idxa+1
          ldima(idxa) = ldim
          ldim = ldim * nstra(hpvx)
        end if
      end do

      return
      end
      
