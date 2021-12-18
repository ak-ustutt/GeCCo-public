*----------------------------------------------------------------------*
      subroutine trace_op_blk(xfac,            !prefactor   
     &     xtrop,xop,                          !buffers: res,op
     &     tra_op,tra_trop,
     &     ncblk_op,nablk_op,ncblk_trop,nablk_trop,
     &     ncblk_cnt,nablk_cnt,
     &     hpvxopc, hpvxopa,
     &     hpvxtropc, hpvxtropa,
     &     lstrop, lstr_trop,
     &     lstr_cnt,
     &     map_info_c, map_info_a,
     &     map_tropcntc, map_tropcnta
     &     )                     
*----------------------------------------------------------------------*
*     inner routine for taking the (partial) trace of an operator
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

      logical, intent(in) ::
     &     tra_op, tra_trop
      integer, intent(in) ::
     &     ncblk_op,nablk_op,ncblk_trop,nablk_trop,
     &     ncblk_cnt,nablk_cnt,
     &     lstrop(*), lstr_trop(*),
     &     hpvxopa(*), hpvxtropa(*),
     &     hpvxopc(*), hpvxtropc(*),
     &     lstr_cnt(*),
     &     map_info_c(*), map_info_a(*),
     &     map_tropcnta(*), map_tropcntc(*)
      real(8), intent(in) ::
     &     xop(*), xfac
      real(8), intent(inout) ::
     &     xtrop(*)

      integer ::
     &     ielmap, ielmap2,
     &     idxtropa(nablk_trop), idxtropc(ncblk_trop),
     &     nstr_tropa1(nablk_op),   nstr_tropc1(ncblk_op),
     &     nstr_cnta1(nablk_op),   nstr_cntc1(ncblk_op),
     &     nstr_tropcnta(nablk_op), nstr_tropcntc(ncblk_op),
     &     ldim_opa(nablk_op), ldim_opc(ncblk_op),
     &     ldim_tropa(nablk_trop), ldim_tropc(ncblk_trop)
      integer ::
     &     ioff, istr1, istr2, idx1, idx2, idx,
     &     nstr_tropa_tot, nstr_tropc_tot,
     &     nstr_cnta_tot, nstr_cntc_tot,
     &     isgnopa, isgnopc,
     &     idxop, idxtrop,
     &     idx0op,
     &     istr_tropa, istr_tropc, 
     &     istr_cnta, istr_cntc, icmp
      real(8) ::
     &     sgn

      integer, external ::
     &     idx_str_blk3, ielprd
c dbg
c      integer lenop, lentrop
c
c      lenop1  = ielprd(lstrop1,ncblk_op1+nablk_op1)
c      lenop2  = ielprd(lstrop2,ncblk_op2+nablk_op2)
c      lentrop = ielprd(lstr_trop,ncblk_trop+nablk_trop)
c      print *,'lenop12: ',lenop12
c dbg

      nstr_tropc_tot = ielprd(lstr_trop,ncblk_trop)
      nstr_tropa_tot = ielprd(lstr_trop(ncblk_trop+1),nablk_trop)
      nstr_cntc_tot = ielprd(lstr_cnt,ncblk_cnt)
      nstr_cnta_tot = ielprd(lstr_cnt(ncblk_cnt+1),nablk_cnt)

      call set_strmapdim_c(nstr_tropcntc,nstr_cntc1,nstr_tropc1,
     &                   ncblk_op,
     &                   lstr_cnt,lstr_trop,map_info_c)
      call set_strmapdim_c(nstr_tropcnta,nstr_cnta1,nstr_tropa1,
     &                   nablk_op,
     &                   lstr_cnt(ncblk_cnt+1),
     &                     lstr_trop(ncblk_trop+1),map_info_a)


      call set_op_ldim_c(ldim_opc,ldim_opa,
     &                   hpvxopc,hpvxopa,
     &                   lstrop,ncblk_op,nablk_op,tra_op)
      call set_op_ldim_c(ldim_tropc,ldim_tropa,
     &                   hpvxtropc,hpvxtropa,
     &                   lstr_trop,ncblk_trop,nablk_trop,
     &                                                 tra_trop)

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'News from trace_op_blk')
      end if

      ! loop over external A strings of operator 1
      tropa: do istr_tropa = 1, nstr_tropa_tot

        ! loop over contraction A string (= contraction C string!)
        cnta: do istr_cnta = 1, nstr_cnta_tot
          ioff = 0
          isgnopa = 1
          istr1 = istr_tropa-1
          istr2 = istr_cnta-1
          idx0op = 1
          ! get sign and idxop1a
          do icmp = 1, nablk_op
            idx1 = mod(istr1,nstr_tropa1(icmp))+1
            idx2 = mod(istr2,nstr_cnta1(icmp))+1
            idx  = (idx1-1)*nstr_cnta1(icmp)+idx2
            ielmap = map_tropcnta(ioff+idx)
            if (ielmap.eq.0) cycle cnta
            isgnopa = isgnopa*sign(1,ielmap)
            idx0op = idx0op+(abs(ielmap)-1)*ldim_opa(icmp)
            ioff = ioff + nstr_tropcnta(icmp)
            istr1 = istr1/nstr_tropa1(icmp)
            istr2 = istr2/nstr_cnta1(icmp)
          end do

          ! loop over external C strings of operator 1
          tropc: do istr_tropc = 1, nstr_tropc_tot

            idxtrop = idx_str_blk3(istr_tropc-1,istr_tropa-1,
     &                           ldim_tropc,ldim_tropa,
     &                           ncblk_trop,nablk_trop)


            ! loop over contraction C string
            istr_cntc = istr_cnta

            ioff = 0
            isgnopc = 1
            istr1 = istr_tropc-1
            istr2 = istr_cntc-1
            idxop = idx0op
            ! get sign and idxop1c
            do icmp = 1, ncblk_op
              idx1 = mod(istr1,nstr_tropc1(icmp))+1
              idx2 = mod(istr2,nstr_cntc1(icmp))+1
              idx  = (idx1-1)*nstr_cntc1(icmp)+idx2
              ielmap = map_tropcntc(ioff+idx)
              if (ielmap.eq.0) cycle tropc
              isgnopc = isgnopc*sign(1,ielmap)
              idxop = idxop+(abs(ielmap)-1)*ldim_opc(icmp)
              ioff = ioff + nstr_tropcntc(icmp)
              istr1 = istr1/nstr_tropc1(icmp)
              istr2 = istr2/nstr_cntc1(icmp)
            end do

            sgn = xfac*dble(isgnopa*isgnopc)

            xtrop(idxtrop) = xtrop(idxtrop)
     &                          + sgn * xop(idxop)
c dbg
c                  if (lenop12.eq.2) then
c            print *,' sng ',isgnopc,isgnopa
c            print *,' +++ ',idxop,sgn,xop(idxop),
c     &                                ' -> ',idxtrop
c                  end if
c dbg

          end do tropc

        end do cnta

      end do tropa

      return
      end

