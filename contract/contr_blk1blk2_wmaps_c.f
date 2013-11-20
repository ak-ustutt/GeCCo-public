*----------------------------------------------------------------------*
      subroutine contr_blk1blk2_wmaps_c(xfac,            !prefactor   
     &     xop1op2,xop1,xop2,                     !buffers: res,op1,op2
     &     tra_op1,tra_op2,tra_op1op2,
     &     ncblk_op1,nablk_op1,ncblk_ex1,nablk_ex1,
     &     ncblk_op2,nablk_op2,ncblk_ex2,nablk_ex2,
     &     ncblk_cnt,nablk_cnt,ncblk_op1op2,nablk_op1op2,
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
      logical, intent(in) ::
     &     tra_op1, tra_op2, tra_op1op2
      integer, intent(in) ::
     &     ncblk_op1,nablk_op1,ncblk_ex1,nablk_ex1,
     &     ncblk_op2,nablk_op2,ncblk_ex2,nablk_ex2,
     &     ncblk_cnt,nablk_cnt,ncblk_op1op2,nablk_op1op2,
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
     &     idxop1op2a(nablk_op1op2), idxop1op2c(ncblk_op1op2)
      integer ::
     &     ireo_ex1a1(nablk_ex1),   ireo_ex1c1(ncblk_ex1),
     &     ireo_ex1a12(nablk_ex1),  ireo_ex1c12(ncblk_ex1),
     &     ireo_ex2a2(nablk_ex2),   ireo_ex2c2(ncblk_ex2),
     &     ireo_ex2a12(nablk_ex2),  ireo_ex2c12(ncblk_ex2),
     &     ireo_cnta1(nablk_cnt),   ireo_cntc1(ncblk_cnt),
     &     ireo_cnta2(nablk_cnt),   ireo_cntc2(ncblk_cnt)
      integer ::
     &     istr_ex1a1(nablk_op1),   istr_ex1c1(ncblk_op1),
     &     istr_ex1a12(nablk_op1op2),  istr_ex1c12(ncblk_op1op2),
     &     istr_ex2a2(nablk_op2),   istr_ex2c2(ncblk_op2),
     &     istr_ex2a12(nablk_op1op2),  istr_ex2c12(ncblk_op1op2),
     &     istr_cnta1(nablk_op1),   istr_cntc1(ncblk_op1),
     &     istr_cnta2(ncblk_op2),   istr_cntc2(nablk_op2)
      integer ::
     &     nstr_ex1a1(nablk_op1),   nstr_ex1c1(ncblk_op1),
     &     nstr_ex1a12(nablk_op1op2),  nstr_ex1c12(ncblk_op1op2),
     &     nstr_ex2a2(nablk_op2),   nstr_ex2c2(ncblk_op2),
     &     nstr_ex2a12(nablk_op1op2),  nstr_ex2c12(ncblk_op1op2),
     &     nstr_cnta1(nablk_op1),   nstr_cntc1(ncblk_op1),
     &     nstr_cnta2(ncblk_op2),   nstr_cntc2(nablk_op2),
     &     nstr_ex1ex2a(nablk_op1op2), nstr_ex1ex2c(ncblk_op1op2),
     &     nstr_ex1cnta(nablk_op1), nstr_ex1cntc(ncblk_op1),
     &     nstr_ex2cnta(nablk_op2), nstr_ex2cntc(ncblk_op2)
      integer ::
     &     ldim_op1a(nablk_op1), ldim_op1c(ncblk_op1),
     &     ldim_op2a(nablk_op2), ldim_op2c(ncblk_op2),
     &     ldim_op1op2a(nablk_op1op2), ldim_op1op2c(ncblk_op1op2)
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
c      integer lenop1, lenop2, lenop12
c
c      lenop1  = ielprd(lstrop1,ncblk_op1+nablk_op1)
c      lenop2  = ielprd(lstrop2,ncblk_op2+nablk_op2)
c      lenop12 = ielprd(lstrop1op2,ncblk_op1op2+nablk_op1op2)
c      print *,'lenop12: ',lenop12
c dbg

      nstr_ex1c_tot = ielprd(lstr_ex1,ncblk_ex1)
      nstr_ex1a_tot = ielprd(lstr_ex1(ncblk_ex1+1),nablk_ex1)
      nstr_ex2c_tot = ielprd(lstr_ex2,ncblk_ex2)
      nstr_ex2a_tot = ielprd(lstr_ex2(ncblk_ex2+1),nablk_ex2)
      nstr_cntc_tot = ielprd(lstr_cnt,ncblk_cnt)
      nstr_cnta_tot = ielprd(lstr_cnt(ncblk_cnt+1),nablk_cnt)

      ngastp_op1c = ncblk_op1
      ngastp_op1a = nablk_op1
      ngastp_op2c = ncblk_op2
      ngastp_op2a = nablk_op2
      ngastp_op1op2c = ncblk_op1op2
      ngastp_op1op2a = nablk_op1op2

      ! C: ex1,ex2
      call set_strmapdim_c2(nstr_ex1ex2c,nstr_ex1c12,nstr_ex2c12,
     &                   ireo_ex1c12,ireo_ex2c12,
     &                   ngastp_op1op2c,
     &                   lstr_ex1,lstr_ex2,map_info_12c)
      ! A: ex2,ex1
      call set_strmapdim_c2(nstr_ex1ex2a,nstr_ex2a12,nstr_ex1a12,
     &                   ireo_ex2a12,ireo_ex1a12,
     &                   ngastp_op1op2a,
     &                   lstr_ex2(ncblk_ex2+1),
     &                    lstr_ex1(ncblk_ex1+1),map_info_12a)
      call set_strmapdim_c2(nstr_ex1cntc,nstr_cntc1,nstr_ex1c1,
     &                   ireo_cntc1,ireo_ex1c1,
     &                   ngastp_op1c,
     &                   lstr_cnt,lstr_ex1,map_info_1c)
      call set_strmapdim_c2(nstr_ex1cnta,nstr_cnta1,nstr_ex1a1,
     &                   ireo_cnta1,ireo_ex1a1,
     &                   ngastp_op1a,
     &                   lstr_cnt(ncblk_cnt+1),
     &                     lstr_ex1(ncblk_ex1+1),map_info_1a)
      call set_strmapdim_c2(nstr_ex2cntc,nstr_cnta2,nstr_ex2c2,
     &                   ireo_cnta2,ireo_ex2c2,
     &                   ngastp_op2c,
     &                   lstr_cnt(ncblk_cnt+1),lstr_ex2,map_info_2c)
      call set_strmapdim_c2(nstr_ex2cnta,nstr_cntc2,nstr_ex2a2,
     &                   ireo_cntc2,ireo_ex2a2,
     &                   ngastp_op2a,
     &                   lstr_cnt,lstr_ex2(ncblk_ex2+1),map_info_2a)

      call set_op_ldim_c(ldim_op1c,ldim_op1a,
     &                   hpvxop1c,hpvxop1a,
     &                   lstrop1,ncblk_op1,nablk_op1,tra_op1)
      call set_op_ldim_c(ldim_op2c,ldim_op2a,
     &                   hpvxop2c,hpvxop2a,
     &                   lstrop2,ncblk_op2,nablk_op2,tra_op2)
      call set_op_ldim_c(ldim_op1op2c,ldim_op1op2a,
     &                   hpvxop1op2c,hpvxop1op2a,
     &                   lstrop1op2,ncblk_op1op2,nablk_op1op2,
     &                                                 tra_op1op2)

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'News from contr_op1op2_wmaps_c')
      end if

      ! loop over external A strings of operator 1
      ex1a: do istr_ex1a = 1, nstr_ex1a_tot

        ! break down into components
        istr_ex1a12(1:nablk_op1op2) = 1
        istr_ex1a1 (1:nablk_op1   ) = 1
        istr1 = istr_ex1a-1
        do icmp = 1, nablk_ex1
          idx1 = mod(istr1,lstr_ex1(ncblk_ex1+icmp))+1
          istr_ex1a12(ireo_ex1a12(icmp)) = idx1
          istr_ex1a1 (ireo_ex1a1 (icmp)) = idx1
          istr1 = istr1/lstr_ex1(ncblk_ex1+icmp)
        end do

        ! loop over external A strings of operator 2
        ex2a: do istr_ex2a = 1, nstr_ex2a_tot

          ! break down into components
          istr_ex2a12(1:nablk_op1op2) = 1
          istr_ex2a2 (1:nablk_op2   ) = 1
          istr1 = istr_ex2a-1
          do icmp = 1, nablk_ex2
            idx1 = mod(istr1,lstr_ex2(ncblk_ex2+icmp))+1
            istr_ex2a12(ireo_ex2a12(icmp)) = idx1
            istr_ex2a2 (ireo_ex2a2 (icmp)) = idx1
            istr1 = istr1/lstr_ex2(ncblk_ex2+icmp)
          end do

          ioff = 0
          isgnra = 1
          ! get sign and idxop1op2a
          do icmp = 1, ngastp_op1op2a
            idx1 = istr_ex1a12(icmp)
            idx2 = istr_ex2a12(icmp)
            idx  = (idx1-1)*nstr_ex2a12(icmp)+idx2
            ielmap = map_ex1ex2a(ioff+idx)
            if (ielmap.eq.0) cycle ex2a
            isgnra = isgnra*sign(1,ielmap)
            idxop1op2a(icmp) = abs(ielmap)-1
            ioff = ioff + nstr_ex1ex2a(icmp)
          end do

          ! loop over external C strings of operator 1
          ex1c: do istr_ex1c = 1, nstr_ex1c_tot

            ! break down into components
            istr_ex1c12(1:ncblk_op1op2) = 1
            istr_ex1c1 (1:ncblk_op1   ) = 1
            istr1 = istr_ex1c-1
            do icmp = 1, ncblk_ex1
              idx1 = mod(istr1,lstr_ex1(icmp))+1
              istr_ex1c12(ireo_ex1c12(icmp)) = idx1
              istr_ex1c1 (ireo_ex1c1 (icmp)) = idx1
              istr1 = istr1/lstr_ex1(icmp)
            end do

            ! loop over external C strings of operator 2
            ex2c: do istr_ex2c = 1, nstr_ex2c_tot

              ! break down into components
              istr_ex2c12(1:ncblk_op1op2) = 1
              istr_ex2c2 (1:ncblk_op2   ) = 1
              istr1 = istr_ex2c-1
              do icmp = 1, ncblk_ex2
                idx1 = mod(istr1,lstr_ex2(icmp))+1
                istr_ex2c12(ireo_ex2c12(icmp)) = idx1
                istr_ex2c2 (ireo_ex2c2 (icmp)) = idx1
                istr1 = istr1/lstr_ex2(icmp)
              end do

              ioff = 0
              isgnrc = 1
              ! get sign and idxop1op2c
              do icmp = 1, ngastp_op1op2c
                idx1 = istr_ex2c12(icmp)
                idx2 = istr_ex1c12(icmp)
                idx  = (idx1-1)*nstr_ex1c12(icmp)+idx2
                ielmap = map_ex1ex2c(ioff+idx)
                if (ielmap.eq.0) cycle ex2c
                isgnrc = isgnrc*sign(1,ielmap)
                idxop1op2c(icmp) = abs(ielmap)-1
                ioff = ioff + nstr_ex1ex2c(icmp)
              end do

              isgnr = isgnrc*isgnra

              idxop1op2 = idx_str_blk3(idxop1op2c,idxop1op2a,
     &                                 ldim_op1op2c,ldim_op1op2a,
     &                                 ngastp_op1op2c,ngastp_op1op2a)
c dbg
c                  if (idxop1op2.lt.1) stop 'range12-l'
c                  if (idxop1op2.gt.lenop12) then
c                  if (lenop12.eq.2) then
c                    print *,'-->',idxop1op2
c                    print *,'ngastp:',ngastp_op1op2c,ngastp_op1op2a
c                    print *,idxop1op2a
c                    print *,ldim_op1op2a
c                    print *,idxop1op2c
c                    print *,ldim_op1op2c
c                    stop 'range12-h'
c                  end if
c dbg

              ! loop over contraction A string
              cnta: do istr_cnta = 1, nstr_cnta_tot

                ! break down into components
                istr_cnta1 (1:nablk_op1) = 1
                istr_cnta2 (1:ncblk_op2) = 1
                istr1 = istr_cnta-1
                do icmp = 1, nablk_cnt
                  idx1 = mod(istr1,lstr_cnt(ncblk_cnt+icmp))+1
                  istr_cnta1(ireo_cnta1(icmp)) = idx1
                  istr_cnta2(ireo_cnta2(icmp)) = idx1
                  istr1 = istr1/lstr_cnt(ncblk_cnt+icmp)
                end do

                ioff = 0
                isgnop1a = 1
                idx0op1 = 1
                ! get sign and idxop1a
                do icmp = 1, ngastp_op1a
                  idx1 = istr_ex1a1(icmp)
                  idx2 = istr_cnta1(icmp)
                  idx  = (idx1-1)*nstr_cnta1(icmp)+idx2
                  ielmap = map_ex1cnta(ioff+idx)
                  if (ielmap.eq.0) cycle cnta
                  isgnop1a = isgnop1a*sign(1,ielmap)
                  idx0op1 = idx0op1+(abs(ielmap)-1)*ldim_op1a(icmp)
                  ioff = ioff + nstr_ex1cnta(icmp)
                end do

                ioff = 0
                isgnop2c = 1
                idx0op2 = 1
                ! get sign and idxop2c
                do icmp = 1, ngastp_op2c
                  idx1 = istr_ex2c2(icmp)
                  idx2 = istr_cnta2(icmp)
                  idx  = (idx1-1)*nstr_cnta2(icmp)+idx2
                  ielmap = map_ex2cntc(ioff+idx)
                  if (ielmap.eq.0) cycle cnta
                  isgnop2c = isgnop2c*sign(1,ielmap)
                  idx0op2 = idx0op2+(abs(ielmap)-1)*ldim_op2c(icmp)
                  ioff = ioff + nstr_ex2cntc(icmp)
                end do

                isgn0 = isgnr*isgnop1a*isgnop2c

                ! loop over contraction C string
                cntc: do istr_cntc = 1, nstr_cntc_tot

                  ! break down into components
                  istr_cntc1 (1:ncblk_op1) = 1
                  istr_cntc2 (1:nablk_op2) = 1
                  istr1 = istr_cntc-1
                  do icmp = 1, ncblk_cnt
                    idx1 = mod(istr1,lstr_cnt(icmp))+1
                    istr_cntc1(ireo_cntc1(icmp)) = idx1
                    istr_cntc2(ireo_cntc2(icmp)) = idx1
                    istr1 = istr1/lstr_cnt(icmp)
                  end do

                  ioff = 0
                  isgnop1c = 1
                  idxop1 = idx0op1
                  ! get sign and idxop1c
                  do icmp = 1, ngastp_op1c
                    idx1 = istr_ex1c1(icmp)
                    idx2 = istr_cntc1(icmp)
                    idx  = (idx1-1)*nstr_cntc1(icmp)+idx2
                    ielmap = map_ex1cntc(ioff+idx)
                    if (ielmap.eq.0) cycle cntc
                    isgnop1c = isgnop1c*sign(1,ielmap)
                    idxop1 = idxop1+(abs(ielmap)-1)*ldim_op1c(icmp)
                    ioff = ioff + nstr_ex1cntc(icmp)
                  end do

                  ioff = 0
                  isgnop2a = 1
                  idxop2 = idx0op2
                  ! get sign and idxop2a
                  do icmp = 1, ngastp_op2a
                    idx1 = istr_ex2a2(icmp)
                    idx2 = istr_cntc2(icmp)
                    idx  = (idx1-1)*nstr_cntc2(icmp)+idx2
                    ielmap = map_ex2cnta(ioff+idx)
                    if (ielmap.eq.0) cycle cntc
                    isgnop2a = isgnop2a*sign(1,ielmap)
                    idxop2 = idxop2+(abs(ielmap)-1)*ldim_op2a(icmp)
                    ioff = ioff + nstr_ex2cnta(icmp)
                  end do

                  sgn = xfac*dble(isgn0*isgnop1c*isgnop2a)

                  ! xfac contained in sgn
c dbg
c                  if (idxop1.lt.1) stop 'range1'
c                  if (idxop2.lt.1) stop 'range2'
c                  if (idxop1.gt.lenop1) stop 'range1'
c                  if (idxop2.gt.lenop2) stop 'range2'
c dbg
                  xop1op2(idxop1op2) = xop1op2(idxop1op2)
     &                          + sgn * xop1(idxop1)
     &                                * xop2(idxop2)
c dbg
c                  if (lenop12.eq.2) then
c                  print *,' sng ',isgnop1c,isgnop1a,
c     &                 isgnop2c,isgnop2a,isgnr
c                  print *,' +++ ',sgn,xop1(idxop1),xop2(idxop2),
c     &                                ' -> ',idxop1op2
c                  end if
c dbg

                end do cntc           
              end do cnta

            end do ex2c
          end do ex1c

        end do ex2a
      end do ex1a

c dbg
c      print *,'updated iop1op2: '
c      write(6,'(x,": ",4g20.12)') xop1op2(1:lenop12)
c dbg        
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

      if (ngastp_c.gt.0) idx_str_blk3 = idx_str_blk3 + idxc(1)*ldimc(1)
      if (ngastp_a.gt.0) idx_str_blk3 = idx_str_blk3 + idxa(1)*ldima(1)

c      idx_str_blk3 = idxc(1)*ldimc(1)+idxa(1)*ldima(1) + 1
cmh      if ((ngastp_c*ngastp_a).le.1) return

      do idx = 2, ngastp_c
        idx_str_blk3 = idx_str_blk3+idxc(idx)*ldimc(idx)
      end do
      do idx = 2, ngastp_a
        idx_str_blk3 = idx_str_blk3+idxa(idx)*ldima(idx)
      end do

      return
      end
      
