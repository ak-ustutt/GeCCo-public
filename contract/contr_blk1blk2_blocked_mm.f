*----------------------------------------------------------------------*
      subroutine contr_blk1blk2_blocked_mm(xfac,         !prefactor   
     &     xop1op2,xop1,xop2,                     !buffers: res,op1,op2
     &     tra_op1,tra_op2,tra_op1op2,
     &     xscr,lenscr,lenblock,                   ! scratch
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
*     inner contraction routine, generation 3
*     - blocked reordering
*     - blocked dgemm
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'hpvxseq.h'

      integer, parameter ::
     &     ntest = 000

c      integer, parameter ::
c     &     lenblock = 1000
c     &     lenblock = 200

      ! lstrop1op2a(ncblk_op1op2a): # of strings
      ! etc.
      ! iocc_op1op2: occupation
      ! nstr_ex1/ex2/cnt,a/c: # strings
      ! map_ex1ex2a/c  mapping ex1*ex2->op1op2
      ! etc. for ex1cnt ex2cnt
      logical, intent(in) ::
     &     tra_op1, tra_op2, tra_op1op2
      integer, intent(in) ::
     &     lenblock,
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
     &     map_ex2cnta(*), map_ex2cntc(*),
     &     lenscr
      real(8), intent(in) ::
     &     xop1(*), xop2(*), xfac
      real(8), intent(inout) ::
     &     xop1op2(*),
     &     xscr(*)


      logical ::
     &     op1shorter
      integer ::
     &     ireo_ex1a1(nablk_ex1),   ireo_ex1c1(ncblk_ex1),
     &     ireo_ex1a12(nablk_ex1),  ireo_ex1c12(ncblk_ex1),
     &     ireo_ex2a2(nablk_ex2),   ireo_ex2c2(ncblk_ex2),
     &     ireo_ex2a12(nablk_ex2),  ireo_ex2c12(ncblk_ex2),
     &     ireo_cnta1(nablk_cnt),   ireo_cntc1(ncblk_cnt),
     &     ireo_cnta2(nablk_cnt),   ireo_cntc2(ncblk_cnt)
      integer ::
     &     nstr_ex1a1(nablk_op1),   nstr_ex1c1(ncblk_op1),
     &     nstr_ex1a12(nablk_op1op2),  nstr_ex1c12(ncblk_op1op2),
     &     nstr_ex2a2(nablk_op2),   nstr_ex2c2(ncblk_op2),
     &     nstr_ex2a12(nablk_op1op2),  nstr_ex2c12(ncblk_op1op2),
     &     nstr_cnta1(nablk_op1),   nstr_cntc1(ncblk_op1),
     &     nstr_cnta2(ncblk_op2),   nstr_cntc2(nablk_op2),
     &     nstr_ex1ex2a(nablk_op1op2), nstr_ex1ex2c(ncblk_op1op2),
     &     nstr_ex1cnta(nablk_op1), nstr_ex1cntc(ncblk_op1),
     &     nstr_ex2cnta(nablk_op2), nstr_ex2cntc(ncblk_op2),
     &     ldim_op1a(nablk_op1), ldim_op1c(ncblk_op1),
     &     ldim_op2a(nablk_op2), ldim_op2c(ncblk_op2),
     &     ldim_op1op2a(nablk_op1op2), ldim_op1op2c(ncblk_op1op2)
      integer ::
     &     lscr_op1op2, idxopsscr, idxoplscr, idxop1op2scr, idx, idxend,
     &     nstr_ex1a_tot, nstr_ex1c_tot,
     &     nstr_ex2a_tot, nstr_ex2c_tot,
     &     nstr_exsa_tot, nstr_exsc_tot,
     &     nstr_exla_tot, nstr_exlc_tot,
     &     nstr_cnta_tot, nstr_cntc_tot,
     &     lenop1, lenop2, lenop12,
     &     maxlen_cnt_batch, maxlen_exl_batch, maxlen_exs_batch,
     &     n_cnt_batch, n_exl_batch, n_exs_batch,
     &     cnt_batch, exl_batch, exs_batch, ncnt, nexl, nexs,
     &     istr_cntc_bst, istr_cntc_bnd, istr_cnta_bst, istr_cnta_bnd, 
     &     istr_exsc_bst, istr_exsc_bnd, istr_exsa_bst, istr_exsa_bnd, 
     &     istr_exlc_bst, istr_exlc_bnd, istr_exla_bst, istr_exla_bnd

      integer, external ::
     &     ielprd

      lenop1  = ielprd(lstrop1,ncblk_op1+nablk_op1)
      lenop2  = ielprd(lstrop2,ncblk_op2+nablk_op2)
      lenop12 = ielprd(lstrop1op2,ncblk_op1op2+nablk_op1op2)

      nstr_ex1c_tot = ielprd(lstr_ex1,ncblk_ex1)
      nstr_ex1a_tot = ielprd(lstr_ex1(ncblk_ex1+1),nablk_ex1)
      nstr_ex2c_tot = ielprd(lstr_ex2,ncblk_ex2)
      nstr_ex2a_tot = ielprd(lstr_ex2(ncblk_ex2+1),nablk_ex2)
      nstr_cntc_tot = ielprd(lstr_cnt,ncblk_cnt)
      nstr_cnta_tot = ielprd(lstr_cnt(ncblk_cnt+1),nablk_cnt)

      ! C: ex1,ex2
      call set_strmapdim_c2(nstr_ex1ex2c,nstr_ex1c12,nstr_ex2c12,
     &                   ireo_ex1c12,ireo_ex2c12,
     &                   ncblk_op1op2,
     &                   lstr_ex1,lstr_ex2,map_info_12c)
      ! A: ex2,ex1
      call set_strmapdim_c2(nstr_ex1ex2a,nstr_ex2a12,nstr_ex1a12,
     &                   ireo_ex2a12,ireo_ex1a12,
     &                   nablk_op1op2,
     &                   lstr_ex2(ncblk_ex2+1),
     &                    lstr_ex1(ncblk_ex1+1),map_info_12a)
      call set_strmapdim_c2(nstr_ex1cntc,nstr_cntc1,nstr_ex1c1,
     &                   ireo_cntc1,ireo_ex1c1,
     &                   ncblk_op1,
     &                   lstr_cnt,lstr_ex1,map_info_1c)
      call set_strmapdim_c2(nstr_ex1cnta,nstr_cnta1,nstr_ex1a1,
     &                   ireo_cnta1,ireo_ex1a1,
     &                   nablk_op1,
     &                   lstr_cnt(ncblk_cnt+1),
     &                     lstr_ex1(ncblk_ex1+1),map_info_1a)
      call set_strmapdim_c2(nstr_ex2cntc,nstr_cnta2,nstr_ex2c2,
     &                   ireo_cnta2,ireo_ex2c2,
     &                   ncblk_op2,
     &                   lstr_cnt(ncblk_cnt+1),lstr_ex2,map_info_2c)
      call set_strmapdim_c2(nstr_ex2cnta,nstr_cntc2,nstr_ex2a2,
     &                   ireo_cntc2,ireo_ex2a2,
     &                   nablk_op2,
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
     &                                               tra_op1op2)
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &       'News from contr_op1op2_blocked_mm')
      end if

      ! divide scratch into three blocks
      ! block for raw op1op2 batch needs at most lenblock^2
      ! (or #Ex1*#Ex2, if that is less)
      lscr_op1op2 = min(nstr_ex1c_tot*nstr_ex1a_tot*
     &                  nstr_ex2c_tot*nstr_ex2a_tot,lenblock*lenblock)
      ! blocks for resorted ops:
      ! is Op1 or Op2 list shorter? OpSHORT, OpLONG
      ! determine length of CNT blocks
      ! such that OpSH(CNT,EX) fits into scratch for all EX
      op1shorter = nstr_ex1c_tot*nstr_ex1a_tot.le.
     &             nstr_ex2c_tot*nstr_ex2a_tot
      if (op1shorter) then
        maxlen_exl_batch = min(nstr_ex2c_tot*nstr_ex2a_tot,lenblock)
        maxlen_exs_batch = min(nstr_ex1c_tot*nstr_ex1a_tot,lenblock)
        maxlen_cnt_batch = min(nstr_cntc_tot*nstr_cnta_tot,
     &       (lenscr-lscr_op1op2-1)/
     &       (    nstr_ex1c_tot*nstr_ex1a_tot + maxlen_exl_batch)+1)
        nstr_exsc_tot = nstr_ex1c_tot
        nstr_exsa_tot = nstr_ex1a_tot
        nstr_exlc_tot = nstr_ex2c_tot
        nstr_exla_tot = nstr_ex2a_tot
      else
        maxlen_exl_batch = min(nstr_ex1c_tot*nstr_ex1a_tot,lenblock)
        maxlen_exs_batch = min(nstr_ex2c_tot*nstr_ex2a_tot,lenblock)
        maxlen_cnt_batch = min(nstr_cntc_tot*nstr_cnta_tot,
     &       (lenscr-lscr_op1op2-1)/
     &       (    nstr_ex2c_tot*nstr_ex2a_tot + maxlen_exl_batch)+1)
        nstr_exsc_tot = nstr_ex2c_tot
        nstr_exsa_tot = nstr_ex2a_tot
        nstr_exlc_tot = nstr_ex1c_tot
        nstr_exla_tot = nstr_ex1a_tot
      end if

      if (maxlen_cnt_batch.lt.1)
     &     call quit(1,'contr_blk1blk2_blocked_mm',
     &     'insufficient scratch memory')

      idxopsscr = 1
      idxoplscr = idxopsscr +
     &     nstr_exsc_tot*nstr_exsa_tot*maxlen_cnt_batch
      idxop1op2scr = idxoplscr +
     &     maxlen_cnt_batch*maxlen_exl_batch
      idxend    = idxop1op2scr + maxlen_exl_batch*maxlen_exs_batch

      if (idxend.gt.lenscr)
     &     call quit(1,'contr_blk1blk2_blocked_mm',
     &                 'error in setting up batch lengthes')

      n_cnt_batch = (nstr_cntc_tot*nstr_cnta_tot-1)/maxlen_cnt_batch + 1

      n_exl_batch = (nstr_exlc_tot*nstr_exla_tot-1)/maxlen_exl_batch + 1

      n_exs_batch = (nstr_exsc_tot*nstr_exsa_tot-1)/maxlen_exs_batch + 1

      if (ntest.ge.100) then
        write(luout,*) 'batching info: '
        write(luout,*) 'lenop: ',lenop1,lenop2,lenop12
        write(luout,*) 'op1shorter: ',op1shorter
        write(luout,*) '# CNT    : ',nstr_cntc_tot*nstr_cnta_tot
        write(luout,*) '# EXlong : ',nstr_exlc_tot*nstr_exla_tot
        write(luout,*) '# EXshort: ',nstr_exsc_tot*nstr_exsa_tot
        write(luout,*) 'n_cnt_batch = ',n_cnt_batch
        write(luout,*) 'n_exl_batch = ',n_exl_batch
        write(luout,*) 'n_exs_batch = ',n_exs_batch
        write(luout,*) 'maxlen_cnt_batch = ',maxlen_cnt_batch
        write(luout,*) 'maxlen_exl_batch = ',maxlen_exl_batch
        write(luout,*) 'maxlen_exs_batch = ',maxlen_exs_batch
        write(luout,*) 'lenscr = ',lenscr
        write(luout,*) 'idxscr = ',idxopsscr, idxoplscr, idxop1op2scr
      end if

      ! loop over CNT blocks
      do cnt_batch = 1, n_cnt_batch

        call idxstnd_batch(ncnt,
     &       istr_cntc_bst,istr_cnta_bst,
     &       istr_cntc_bnd,istr_cnta_bnd,
     &       cnt_batch,nstr_cntc_tot,nstr_cnta_tot,
     &       maxlen_cnt_batch)

        if (ntest.ge.1000) then
          write(luout,*) 'cnt-batch: ',cnt_batch
          write(luout,*) 'maxlen_cnt_batch: ',maxlen_cnt_batch
          write(luout,*) 'ncnt = ',ncnt
          write(luout,*) '           ',istr_cntc_bst,istr_cntc_bnd
          write(luout,*) '           ',istr_cnta_bst,istr_cnta_bnd
        end if

        ! resort OpSHORT
        if (op1shorter) then
          call collect_block(xscr(idxopsscr),xop1,
     &         istr_cntc_bst,istr_cntc_bnd,nstr_cntc_tot,
     &         istr_cnta_bst,istr_cnta_bnd,nstr_cnta_tot,
     &         ncnt,.false.,
     &         1            ,nstr_ex1c_tot,nstr_ex1c_tot,
     &         1            ,nstr_ex1a_tot,nstr_ex1a_tot,
     &         nstr_ex1c_tot*nstr_ex1a_tot,
     &         map_ex1cntc,map_ex1cnta,
     &         ncblk_op1,nablk_op1,
     &         ldim_op1c,ldim_op1a,
     &         ncblk_cnt,nablk_cnt,
     &         ncblk_ex1,nablk_ex1,
     &         lstr_cnt,lstr_cnt(ncblk_cnt+1),
     &         lstr_ex1,lstr_ex1(ncblk_ex1+1),
     &         nstr_cntc1,nstr_cnta1,
     &         nstr_ex1c1,nstr_ex1a1,
     &         ireo_cntc1,ireo_cnta1,
     &         ireo_ex1c1,ireo_ex1a1)
        else
          call collect_block(xscr(idxopsscr),xop2,
     &         istr_cnta_bst,istr_cnta_bnd,nstr_cnta_tot, ! CA exchanged !
     &         istr_cntc_bst,istr_cntc_bnd,nstr_cntc_tot, ! CA exchanged !
     &         ncnt,.true.,
     &         1            ,nstr_ex2c_tot,nstr_ex2c_tot,
     &         1            ,nstr_ex2a_tot,nstr_ex2a_tot,
     &         nstr_ex2c_tot*nstr_ex2a_tot,
     &         map_ex2cntc,map_ex2cnta,
     &         ncblk_op2,nablk_op2,
     &         ldim_op2c,ldim_op2a,
     &         nablk_cnt,ncblk_cnt,      ! CA exchanged !
     &         ncblk_ex2,nablk_ex2,
     &         lstr_cnt(ncblk_cnt+1),lstr_cnt,    ! CA exchanged !
     &         lstr_ex2,lstr_ex2(ncblk_ex2+1),
     &         nstr_cnta2,nstr_cntc2,    ! CA exchanged !
     &         nstr_ex2c2,nstr_ex2a2,
     &         ireo_cnta2,ireo_cntc2,    ! CA exchanged !
     &         ireo_ex2c2,ireo_ex2a2)
        end if

        if (ntest.ge.1000) then
          write(luout,*) ncnt,nstr_exsc_tot*nstr_exsa_tot
          if (op1shorter)      write(luout,*) 'resorted operator 1'
          if (.not.op1shorter) write(luout,*) 'resorted operator 2'
          call wrtmat2(xscr(idxopsscr),ncnt,nstr_exsc_tot*nstr_exsa_tot,
     &                                 ncnt,nstr_exsc_tot*nstr_exsa_tot)
        end if

        ! loop over blocks of EX_LONG
        do exl_batch = 1, n_exl_batch

          call idxstnd_batch(nexl,
     &         istr_exlc_bst,istr_exla_bst,
     &         istr_exlc_bnd,istr_exla_bnd,
     &         exl_batch,nstr_exlc_tot,nstr_exla_tot,
     &         maxlen_exl_batch)

          if (ntest.ge.1000) then
            write(luout,*) '  exl_batch: ',exl_batch
            write(luout,*) '  maxlen_exl_batch: ',maxlen_exl_batch
            write(luout,*) '  nexl = ',nexl
            write(luout,*) '  istr_exlc_bst,istr_exla_bst:',
     &           istr_exlc_bst,istr_exla_bst
            write(luout,*) '  istr_exlc_bnd,istr_exla_bnd:',
     &           istr_exlc_bnd,istr_exla_bnd
          end if

          ! resort OpLONG -> (CNT,EX_LONG) for current block
          if (op1shorter) then
            call collect_block(xscr(idxoplscr),xop2,
     &         istr_cnta_bst,istr_cnta_bnd,nstr_cnta_tot, ! CA exchanged !
     &         istr_cntc_bst,istr_cntc_bnd,nstr_cntc_tot, ! CA exchanged !
     &         ncnt,.true.,
     &         istr_exlc_bst,istr_exlc_bnd,nstr_exlc_tot,
     &         istr_exla_bst,istr_exla_bnd,nstr_exla_tot,
     &         nexl,
     &         map_ex2cntc,map_ex2cnta,
     &         ncblk_op2,nablk_op2,
     &         ldim_op2c,ldim_op2a,
     &         nablk_cnt,ncblk_cnt,      ! CA exchanged !
     &         ncblk_ex2,nablk_ex2,
     &         lstr_cnt(ncblk_cnt+1),lstr_cnt,    ! CA exchanged !
     &         lstr_ex2,lstr_ex2(ncblk_ex2+1),
     &         nstr_cnta2,nstr_cntc2,    ! CA exchanged !
     &         nstr_ex2c2,nstr_ex2a2,
     &         ireo_cnta2,ireo_cntc2,    ! CA exchanged !
     &         ireo_ex2c2,ireo_ex2a2)
          else
            call collect_block(xscr(idxoplscr),xop1,
     &         istr_cntc_bst,istr_cntc_bnd,nstr_cntc_tot,
     &         istr_cnta_bst,istr_cnta_bnd,nstr_cnta_tot,
     &         ncnt,.false.,
     &         istr_exlc_bst,istr_exlc_bnd,nstr_exlc_tot,
     &         istr_exla_bst,istr_exla_bnd,nstr_exla_tot,
     &         nexl,
     &         map_ex1cntc,map_ex1cnta,
     &         ncblk_op1,nablk_op1,
     &         ldim_op1c,ldim_op1a,
     &         ncblk_cnt,nablk_cnt,
     &         ncblk_ex1,nablk_ex1,
     &         lstr_cnt,lstr_cnt(ncblk_cnt+1),
     &         lstr_ex1,lstr_ex1(ncblk_ex1+1),
     &         nstr_cntc1,nstr_cnta1,
     &         nstr_ex1c1,nstr_ex1a1,
     &         ireo_cntc1,ireo_cnta1,
     &         ireo_ex1c1,ireo_ex1a1)
          end if
          if (ntest.ge.1000) then
            if (op1shorter)      write(luout,*) 'resorted operator 2'
            if (.not.op1shorter) write(luout,*) 'resorted operator 1'
            call wrtmat2(xscr(idxoplscr),ncnt,nexl,ncnt,nexl)
          end if
          
          ! loop over blocks of EX_SHORT
          do exs_batch = 1, n_exs_batch

            call idxstnd_batch(nexs,
     &         istr_exsc_bst,istr_exsa_bst,
     &         istr_exsc_bnd,istr_exsa_bnd,
     &         exs_batch,nstr_exsc_tot,nstr_exsa_tot,
     &         maxlen_exs_batch)

            if (ntest.ge.1000) then
              write(luout,*) '    exs_batch: ',exs_batch
              write(luout,*) '    maxlen_exs_batch: ',maxlen_exs_batch
              write(luout,*) '    nexs = ',nexs
              write(luout,*) '    istr_exsc_bst,istr_exsa_bst:',
     &             istr_exsc_bst,istr_exsa_bst
              write(luout,*) '    istr_exsc_bnd,istr_exsa_bnd:',
     &             istr_exsc_bnd,istr_exsa_bnd
            end if

            if (op1shorter) then
              ! calculate Op1Op2(EX_SHORT,EX_LONG)
              idx = idxopsscr + (exs_batch-1)*maxlen_exs_batch*ncnt
              call dgemm('t','n',nexs,nexl,ncnt,
     &             xfac,xscr(idx),ncnt,
     &                  xscr(idxoplscr),ncnt,
     &              0d0,xscr(idxop1op2scr),nexs)

              if (ntest.ge.1000) then
                write(luout,*) 'result (fac was ',xfac,')'
                call wrtmat2(xscr(idxop1op2scr),nexs,nexl,nexs,nexl)
              end if

              ! resort Op1Op2
              call scatter_block(xop1op2,xscr(idxop1op2scr),
     &             istr_exsc_bst,istr_exsc_bnd,nstr_exsc_tot,
     &             istr_exsa_bst,istr_exsa_bnd,nstr_exsa_tot,
     &             nexs,
     &             istr_exlc_bst,istr_exlc_bnd,nstr_exlc_tot,
     &             istr_exla_bst,istr_exla_bnd,nstr_exla_tot,
     &             nexl,
     &             map_ex1ex2c,map_ex1ex2a,
     &             ncblk_op1op2,nablk_op1op2,
     &             ldim_op1op2c,ldim_op1op2a,
     &             ncblk_ex1,nablk_ex1,
     &             ncblk_ex2,nablk_ex2,
     &             lstr_ex1,lstr_ex1(ncblk_ex1+1),
     &             lstr_ex2,lstr_ex2(ncblk_ex2+1),
     &             nstr_ex1c12,nstr_ex1a12,
     &             nstr_ex2c12,nstr_ex2a12,
     &             ireo_ex1c12,ireo_ex1a12,
     &             ireo_ex2c12,ireo_ex2a12)

            else

              ! calculate Op1Op2(EX_LONG,EX_SHORT)
              idx = idxopsscr + (exs_batch-1)*maxlen_exs_batch*ncnt
              call dgemm('t','n',nexl,nexs,ncnt,
     &             xfac,xscr(idxoplscr),ncnt,
     &                  xscr(idx),ncnt,
     &              0d0,xscr(idxop1op2scr),nexl)

              if (ntest.ge.1000) then
                write(luout,*) 'result (fac was ',xfac,')'
                call wrtmat2(xscr(idxop1op2scr),nexl,nexs,nexl,nexs)
              end if

              ! resort Op1Op2
              call scatter_block(xop1op2,xscr(idxop1op2scr),
     &             istr_exlc_bst,istr_exlc_bnd,nstr_exlc_tot,
     &             istr_exla_bst,istr_exla_bnd,nstr_exla_tot,
     &             nexl,
     &             istr_exsc_bst,istr_exsc_bnd,nstr_exsc_tot,
     &             istr_exsa_bst,istr_exsa_bnd,nstr_exsa_tot,
     &             nexs,
     &             map_ex1ex2c,map_ex1ex2a,
     &             ncblk_op1op2,nablk_op1op2,
     &             ldim_op1op2c,ldim_op1op2a,
     &             ncblk_ex1,nablk_ex1,
     &             ncblk_ex2,nablk_ex2,
     &             lstr_ex1,lstr_ex1(ncblk_ex1+1),
     &             lstr_ex2,lstr_ex2(ncblk_ex2+1),
     &             nstr_ex1c12,nstr_ex1a12,
     &             nstr_ex2c12,nstr_ex2a12,
     &             ireo_ex1c12,ireo_ex1a12,
     &             ireo_ex2c12,ireo_ex2a12)

            end if

          end do

        end do

      end do

      return
      end

