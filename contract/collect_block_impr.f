*----------------------------------------------------------------------*
      subroutine collect_block_impr(xblk,xop,
     &         kc_st0,kc_nd0,kc_len,ka_st0,ka_nd0,ka_len,k_len,kca_trp,
     &         lc_st0,lc_nd0,lc_len,la_st0,la_nd0,la_len,l_len,
     &         map_kl_c,map_kl_a,
     &         ncblk_op,nablk_op,
     &         ldim_opc,ldim_opa,
     &         ncblk_k,nablk_k,
     &         ncblk_l,nablk_l,
     &         lstr_kc,lstr_ka,
     &         lstr_lc,lstr_la,
     &         nstr_kc,nstr_ka,
     &         nstr_lc,nstr_la,
     &         ireo_kc,ireo_ka,
     &         ireo_lc,ireo_la)
*----------------------------------------------------------------------*
*     collect for a given batch of strings K and L the block Op(K,L)
*     from the full array Op(K*L)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      real(8), intent(in) ::
     &     xop(*)
      real(8), intent(out) ::
     &     xblk(*)
      logical, intent(in) ::
     &     kca_trp
      integer, intent(in) ::
     &     kc_st0,kc_nd0,kc_len,ka_st0,ka_nd0,ka_len,k_len,
     &     lc_st0,lc_nd0,lc_len,la_st0,la_nd0,la_len,l_len,
     &     ncblk_op, nablk_op,
     &     ncblk_k, nablk_k,
     &     ncblk_l, nablk_l
      integer, intent(in) ::
     &     map_kl_c(*), map_kl_a(*)
      integer, intent(in) ::
     &     ldim_opc(ncblk_op), ldim_opa(nablk_op),
     &     lstr_kc(ncblk_op),  lstr_ka(nablk_op),
     &     lstr_lc(ncblk_op),  lstr_la(nablk_op),
     &     nstr_lc(ncblk_op),  nstr_la(nablk_op),
     &     nstr_kc(ncblk_op),  nstr_ka(nablk_op),
     &     ireo_kc(ncblk_op),  ireo_ka(nablk_op),
     &     ireo_lc(ncblk_op),  ireo_la(nablk_op)

      real(8) ::
     &     xel
      integer ::
     &     ka, kc, kc_st, kc_nd, ka_st, ka_nd,
     &     la, lc, lc_st, lc_nd,
     &     ioff, istr1, istr2, icmp, idx1, idx2, idx, ielmap,
     &     isgn, isgna, idxop, idx0op,
     &     kdx0, kdx_st, ldx0, ldx, kldx, inc,
     &     lstr_kc1, ireo_kc1, lstr_lc1, ireo_lc1,
     &     lstr_ka1, ireo_ka1, lstr_la1, ireo_la1
      integer ::
     &     istr_kc(ncblk_op+1),  istr_ka(nablk_op+1),
     &     istr_lc(ncblk_op+1),  istr_la(nablk_op+1),
     &     istr_kc0(ncblk_op+1),  !istr_ka0(nablk_op+1),
     &     istr_lc0(ncblk_op+1),  istr_la0(nablk_op+1)
c dbg:
c     &     ,istr_kc1(ncblk_op+1)
c      logical ok

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'collect_block')
        write(lulog,*) 'k: from ',kc_st0,ka_st0,' to ',
     &                            kc_nd0,ka_nd0,' len: ',kc_len,ka_len
        write(lulog,*) 'l: from ',lc_st0,la_st0,' to ',
     &                            lc_nd0,la_nd0,' len: ',lc_len,la_len
        write(lulog,*) ' nstr_kc:  ',nstr_kc(1:ncblk_op)
        write(lulog,*) ' nstr_ka:  ',nstr_ka(1:nablk_op)
        write(lulog,*) ' nstr_lc:  ',nstr_lc(1:ncblk_op)
        write(lulog,*) ' nstr_la:  ',nstr_la(1:nablk_op)
        write(lulog,*) ' ireo_kc:  ',ireo_kc(1:ncblk_op)
        write(lulog,*) ' ireo_ka:  ',ireo_ka(1:nablk_op)
        write(lulog,*) ' ireo_lc:  ',ireo_lc(1:ncblk_op)
        write(lulog,*) ' ireo_la:  ',ireo_la(1:nablk_op)
        write(lulog,*) ' ncblk_op/nablk_op: ',ncblk_op,nablk_op
        write(lulog,*) ' ncblk_k/nablk_k:   ',ncblk_k, nablk_k
        write(lulog,*) ' ncblk_l/nablk_l:   ',ncblk_l, nablk_l
c dbg
c        print *,'map_kl_c: ',map_kl_c(1:4)
c        print *,'map_kl_a: ',map_kl_a(1:4)
c dbg
      end if

      ! init to zero
      xblk(1:k_len*l_len) = 0d0
      ! init strings (some will remain one)
      istr_ka(1:nablk_op) = 1
      istr_kc0(1:ncblk_op) = 1
      !istr_kc1(1:ncblk_op) = 1
      istr_la0(1:nablk_op) = 1
      istr_lc0(1:ncblk_op) = 1

      ireo_kc1=1
      lstr_kc1=1
      ireo_ka1=1
      lstr_ka1=1
      ireo_lc1=1
      lstr_lc1=1
      ireo_la1=1
      lstr_la1=1
      ! ireo_kc and lstr_kc contain sensible values only if 
      ! ncblk_k>0
      if (ncblk_k.gt.0) ireo_kc1=ireo_kc(1)
      if (ncblk_k.gt.0) lstr_kc1=lstr_kc(1)
      if (ncblk_l.gt.0) ireo_lc1=ireo_lc(1)
      if (ncblk_l.gt.0) lstr_lc1=lstr_lc(1)
      if (nablk_l.gt.0) ireo_la1=ireo_la(1)
      if (nablk_l.gt.0) lstr_la1=lstr_la(1)
      if (nablk_k.gt.0) ireo_ka1=ireo_ka(1)
      if (nablk_k.gt.0) lstr_ka1=lstr_ka(1)

      ! init LA start string (less one)
      istr1 = la_st0-1
      istr_la0(ireo_la1) = mod(la_st0-1,lstr_la1)
      do icmp = 2, nablk_l
        istr1 = istr1/lstr_la(icmp-1)
        istr_la0(ireo_la(icmp)) = mod(istr1,lstr_la(icmp))+1
      end do

      ! KA-strings
      if (.not.kca_trp) then
        ka_st = ka_st0
        ka_nd = ka_nd0
      else
        ka_st = 1
        ka_nd = ka_len
      end if

      ! init KA start string (less one)
      ! we can init istr_ka directly
      istr1 = ka_st-1
      istr_ka(ireo_ka1) = mod(ka_st-1,lstr_ka1)
      do icmp = 2, nablk_k
        istr1 = istr1/lstr_ka(icmp-1)
        istr_ka(ireo_ka(icmp)) = mod(istr1,lstr_ka(icmp))+1
      end do

c dbg
c      print *,'ka_st, ka_nd: ',ka_st, ka_nd
c dbg
      kdx0 = 0
      if (kca_trp) kdx0 = mod(ka_len-ka_st0+1,ka_len)
      ka_loop: do ka = ka_st, ka_nd
        if (.not.kca_trp) then
          kc_st = 1
          kc_nd = kc_len
          if (ka.eq.ka_st0) kc_st = kc_st0
          if (ka.eq.ka_nd0) kc_nd = kc_nd0
        else
          kc_st = kc_st0
          kc_nd = kc_nd0
          if (ka.lt.ka_st0) kc_st = kc_st0+1
          if (ka.gt.ka_nd0) kc_nd = kc_nd0-1
        end if

        ! init KC start string (less one)
        istr1 = kc_st-1
        istr_kc0(ireo_kc1) = mod(kc_st-1,lstr_kc1)
        do icmp = 2, ncblk_k
          istr1 = istr1/lstr_kc(icmp-1)
          istr_kc0(ireo_kc(icmp)) = mod(istr1,lstr_kc(icmp))+1
        end do

        if (kca_trp .and. ka.eq.ka_st0) kdx0 = 0

        ! the array which we collect to is addressed as
        ! (kdx,ldx), where these are the current relative indices
        ! in the string batches K, L
        ! get start point for K
        kdx_st = kdx0
        ! and increment for next round (we have to be cautious,
        ! as we might "cycle" inner loops!
        if (.not.kca_trp) kdx0 = kdx0+kc_nd-kc_st+1
        if (kca_trp) kdx0 = kdx0+1

        ! break down KA into components:
        istr_ka(ireo_ka1) = istr_ka(ireo_ka1)+1
        if (istr_ka(ireo_ka1).gt.lstr_ka1) then
          istr_ka(ireo_ka1) = 1
          inc_ka: do icmp = 2, nablk_k
            istr_ka(ireo_ka(icmp)) = istr_ka(ireo_ka(icmp))+1
            if (istr_ka(ireo_ka(icmp)).le.lstr_ka(icmp)) 
     &                                            exit inc_ka
            istr_ka(ireo_ka(icmp))=1
          end do inc_ka
        end if
        !istr1 = ka-1
        !do icmp = 1, nablk_k
        !  idx1 = mod(istr1,lstr_ka(icmp))+1
        !  istr_ka(ireo_ka(icmp)) = idx1
        !  istr1 = istr1/lstr_ka(icmp)
        !end do

        ! set first index (minus one) of first component
        istr_la(1:nablk_op) = istr_la0(1:nablk_op)

        ! LA-strings
c dbg
c        print *,'  kdx_st: ',kdx_st
c        print *,'  kc_st, kc_nd: ',kc_st, kc_nd
c        print *,'    la_st, la_nd: ',la_st0, la_nd0
c dbg
        ldx0 = 0
        la_loop: do la = la_st0, la_nd0
          lc_st = 1
          lc_nd = lc_len
          if (la.eq.la_st0) lc_st = lc_st0
          if (la.eq.la_nd0) lc_nd = lc_nd0
          ! offset for L batch
          ldx = ldx0
          ldx0 = ldx0+lc_nd-lc_st+1
c dbg
c          print *,'      ldx: ',ldx
c          print *,'      lc_st, lc_nd: ',lc_st, lc_nd
c dbg
          ! init LC start string (less one)
          istr1 = lc_st-1
          istr_lc0(ireo_lc1) = mod(lc_st-1,lstr_lc1)
          do icmp = 2, ncblk_l
            istr1 = istr1/lstr_lc(icmp-1)
            istr_lc0(ireo_lc(icmp)) = mod(istr1,lstr_lc(icmp))+1
          end do

          ! break down LA into components:
          istr_la(ireo_la1) = istr_la(ireo_la1)+1
          if (istr_la(ireo_la1).gt.lstr_la1) then
            istr_la(ireo_la1) = 1
            inc_la: do icmp = 2, nablk_l
              istr_la(ireo_la(icmp)) = istr_la(ireo_la(icmp))+1
              if (istr_la(ireo_la(icmp)).le.lstr_la(icmp)) 
     &                                                exit inc_la
              istr_la(ireo_la(icmp))=1
            end do inc_la
          end if
          ! break down LA into components:
          !istr1 = la-1
          !do icmp = 1, nablk_l
          !  idx1 = mod(istr1,lstr_la(icmp))+1
          !  istr_la(ireo_la(icmp)) = idx1
          !  istr1 = istr1/lstr_la(icmp)
          !end do
        
          ! get KA*LA
          ioff = 0
          isgna = 1
          idx0op = 1
          do icmp = 1, nablk_op
            idx1 = istr_la(icmp)
            idx2 = istr_ka(icmp)
            idx  = (idx1-1)*nstr_ka(icmp)+idx2
            ielmap = map_kl_a(ioff+idx)
            if (ielmap.eq.0) cycle la_loop
            isgna = isgna*sign(1,ielmap)
            idx0op = idx0op+(abs(ielmap)-1)*ldim_opa(icmp)
            ioff = ioff + nstr_ka(icmp)*nstr_la(icmp)
          end do

          ! set first index (minus one) of first component
          istr_lc(1:ncblk_op) = istr_lc0(1:ncblk_op)

          ! LC-strings
          lc_loop: do lc = lc_st, lc_nd
            ldx = ldx+1

            ! break down LC into components:
            istr_lc(ireo_lc1) = istr_lc(ireo_lc1)+1
            if (istr_lc(ireo_lc1).gt.lstr_lc1) then
              istr_lc(ireo_lc1) = 1
              inc_lc: do icmp = 2, ncblk_l
                istr_lc(ireo_lc(icmp)) = istr_lc(ireo_lc(icmp))+1
                if (istr_lc(ireo_lc(icmp)).le.lstr_lc(icmp)) 
     &                                                  exit inc_lc
                istr_lc(ireo_lc(icmp))=1
              end do inc_lc
            end if

            !istr1 = lc-1
            !do icmp = 1, ncblk_l
            !  idx1 = mod(istr1,lstr_lc(icmp))+1
            !  istr_lc(ireo_lc(icmp)) = idx1
            !  istr1 = istr1/lstr_lc(icmp)
            !end do

            ! KC-strings
            ! set index for xblk array:
c dbg
c            print *,'ldx,k_len,kdx_st: ',ldx,k_len,kdx_st
c dbg
            kldx = (ldx-1)*k_len + kdx_st
            inc  = 1
            if (kca_trp) kldx = (ldx-1)*k_len + kdx_st + 1 - ka_len
            if (kca_trp) inc = ka_len
c dbg
c            print *,'ldx,kldx,inc: ',ldx,kldx,inc,kca_trp
c dbg
            ! special loops
            if (ncblk_op.eq.0) then
              xel = dble(isgna)*xop(idx0op)
              kc_loop0: do kc = kc_st, kc_nd
                kldx = kldx+inc

                xblk(kldx) = xel

              end do kc_loop0

            else if (ncblk_op.eq.1) then
 
              kc_loop1: do kc = kc_st, kc_nd
                kldx = kldx+inc

                ! get LC*KC
                idx = (istr_lc(1)-1)*nstr_kc(1)+kc
                ielmap = map_kl_c(idx)
                if (ielmap.eq.0) cycle kc_loop1
                isgn = isgna*sign(1,ielmap)
                idxop = idx0op+(abs(ielmap)-1)*ldim_opc(1)

                xblk(kldx) = dble(isgn)*xop(idxop)

              end do kc_loop1

            else
            ! set first index (minus one) of first component
            istr_kc(1:ncblk_op) = istr_kc0(1:ncblk_op)
c dbg
c            print *,'kc_st',kc_st,' > ',istr_kc(1:ncblk_op)
c dbg

            kc_loop: do kc = kc_st, kc_nd
              kldx = kldx+inc

              ! break down KC into components:
              istr_kc(ireo_kc1) = istr_kc(ireo_kc1)+1
              if (istr_kc(ireo_kc1).gt.lstr_kc1) then
                istr_kc(ireo_kc1) = 1
                inc_kc: do icmp = 2, ncblk_k
                  istr_kc(ireo_kc(icmp)) = istr_kc(ireo_kc(icmp))+1
                  if (istr_kc(ireo_kc(icmp)).le.lstr_kc(icmp)) 
     &                                                  exit inc_kc
                  istr_kc(ireo_kc(icmp))=1
                end do inc_kc
              end if
c dbg
c              print *,'kc   ',kc,' > ',istr_kc(1:ncblk_op)
c dbg
c              istr1 = kc-1
c              do icmp = 1, ncblk_k
c                idx1 = mod(istr1,lstr_kc(icmp))+1
c                istr_kc1(ireo_kc(icmp)) = idx1
c                istr1 = istr1/lstr_kc(icmp)
c              end do
c              
c              ok = .true.
c              do icmp = 1, ncblk_op 
c                ok = ok.and.istr_kc(icmp).eq.istr_kc1(icmp)
c              end do
c              if (.not.ok) then
c                write(lulog,*) istr_kc0(1:ncblk_op)
c                write(lulog,*) istr_kc(1:ncblk_op)
c                write(lulog,*) istr_kc1(1:ncblk_op)
c          
c                write(lulog,*) 'kc_st, kc_nd, kc: ',kc_st,kc_nd, kc
c                write(lulog,*) 'lstr_kc: ',lstr_kc(1:ncblk_k)
c                write(lulog,*) 'ireo_kc: ',ireo_kc(1:ncblk_op)
c
c                call quit(1,'aus','ende - vorbei')
c              end if

              ! get LC*KC
              ioff = 0
              isgn = isgna
              idxop = idx0op
              do icmp = 1, ncblk_op
                idx1 = istr_lc(icmp)
                idx2 = istr_kc(icmp)
                idx = (idx1-1)*nstr_kc(icmp)+idx2
                ielmap = map_kl_c(ioff+idx)
                if (ielmap.eq.0) cycle kc_loop
                isgn = isgn*sign(1,ielmap)
                idxop = idxop+(abs(ielmap)-1)*ldim_opc(icmp)
                ioff = ioff+ nstr_kc(icmp)*nstr_lc(icmp)
              end do

c dbg
c              print *,'ka, kc, la, lc: ',ka, kc, la, lc
c              print *,'kldx, idxop, sign, val: ',
c     &             kldx,idxop,isgn,xop(idxop)
c dbg
              xblk(kldx) = dble(isgn)*xop(idxop)
              !if (isgn.eq.1) then
              !  xblk(kldx) =  xop(idxop)
              !else
              !  xblk(kldx) = -xop(idxop)
              !end if

            end do kc_loop

            end if

          end do lc_loop

        end do la_loop

      end do ka_loop

      return
      end
