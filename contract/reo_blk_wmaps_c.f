*----------------------------------------------------------------------*
      subroutine reo_blk_wmaps_c(xop_reo,xop_ori,
     &     len_reo,len_ori,
     &     sign_reo,
     &     tra_opreo,tra_opori,
     &     ms_op_c,ms_op_a,gm_op_c,gm_op_a,
     &     ms_i_dis_c,ms_i_dis_a,gm_i_dis_c,gm_i_dis_a,
     &     ncblk_opori,nablk_opori,
     &     cinfo_opori_c,cinfo_opori_a,
     &     lstr_opori,
     &     me_opreo, iblk_opreo,
     &     ncblk_opreo,nablk_opreo,
     &     cinfo_opreo_c,cinfo_opreo_a,
     &     ncblk_k,nablk_k,
     &     cinfo_k_c,cinfo_k_a,
     &     ncblk_i0,nablk_i0,
     &     cinfo_i0_c,cinfo_i0_a,
     &     map_info_to_ori_c, map_info_to_ori_a,
     &     map_info_to_reo_c, map_info_to_reo_a,
     &     nsym,str_info,strmap_info)
*----------------------------------------------------------------------*
*
*     reorder  OP([I1][I2]...) -> OP([I1'][I2']...)
*
*     initial version
*     andreas, oct 2007
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'multd2h.h'
      include 'hpvxseq.h'
      include 'def_filinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      ! buffer with originally ordered elements
      ! (single distribution, see below)
      real(8), intent(in) ::
     &     xop_ori(*)
      ! buffer with ALL distributions of current MS, IRREP
      real(8), intent(inout) ::
     &     xop_reo(*)
      ! global sign factor, and lengthes of buffers
      integer, intent(in) ::
     &     sign_reo, len_reo, len_ori
      ! transposed addressing?
      logical, intent(in) ::
     &     tra_opreo, tra_opori

      ! info about operator block in original order
      integer, intent(in) ::
     &     ncblk_opori, nablk_opori,
     &     ms_op_c, ms_op_a, gm_op_c, gm_op_a,
     &     ms_i_dis_c(ncblk_opori), ms_i_dis_a(nablk_opori), ! dis(MS)
     &     gm_i_dis_c(ncblk_opori), gm_i_dis_a(nablk_opori), ! dis(IRREP)
     &     cinfo_opori_c(ncblk_opori,3), cinfo_opori_a(nablk_opori,3),
     &     lstr_opori(ncblk_opori+nablk_opori)

      ! info about reordered operator block
      ! ms_opreo_c == ms_op_c etc.
      ! but we must loop over different distributions
      integer, intent(in) ::
     &     ncblk_opreo, nablk_opreo,
     &     cinfo_opreo_c(ncblk_opreo,3), cinfo_opreo_a(nablk_opreo,3)

      ! info about strings with resolved occupations [I0], [K]
      type(me_list), intent(in) ::
     &     me_opreo
      integer, intent(in) ::
     &     iblk_opreo,
     &     ncblk_k, nablk_k, ncblk_i0, nablk_i0,
     &     cinfo_k_c(ncblk_k,3), cinfo_k_a(nablk_k,3),
     &     cinfo_i0_c(ncblk_i0,3), cinfo_i0_a(nablk_i0,3)
      ! reorder mappings for (condensed) occupations:
      integer, intent(in) ::
     &     map_info_to_ori_c(*), map_info_to_ori_a(*),
     &     map_info_to_reo_c(*), map_info_to_reo_a(*)

      integer, intent(in) ::
     &     nsym
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info

      logical ::
     &     first, possible, fix_success
      integer ::
     &     ms_k_c_max, ms_k_a_max, ms_i0_c_max, ms_i0_a_max,
     &     ms_k_c,     ms_k_a,     ms_i0_c,     ms_i0_a,
     &     gm_k_c,     gm_k_a,     gm_i0_c,     gm_i0_a
      integer ::
     &     ngraph, lenmap, ifree
      integer ::
     &     istr_k_a, istr_k_c, istr_i0_a, istr_i0_c,
     &     nstr_k_c_tot, nstr_k_a_tot,
     &     nstr_i0_c_tot, nstr_i0_a_tot,
     &     idxms_op_a, idxdis,
     &     idx00opreo, idxst_opreo,
     &     idx0opreo, idx_opreo, ioff, istr1, istr2, idiv1, idiv2,
     &     isgnt, isgna,
     &     idx0opori, idx_opori, icmp, ielmap, idx, idx1, idx2

      integer ::
     &     ms_k_dis_c(ncblk_k), ms_k_dis_a(nablk_k),
     &     idxms_k_dis_c(ncblk_k), idxms_k_dis_a(nablk_k),
     &     gm_k_dis_c(ncblk_k), gm_k_dis_a(nablk_k),
     &     ms_i0_dis_c(ncblk_i0), ms_i0_dis_a(nablk_i0),
     &     idxms_i0_dis_c(ncblk_i0), idxms_i0_dis_a(nablk_i0),
     &     gm_i0_dis_c(ncblk_i0), gm_i0_dis_a(nablk_i0),
     &     idxms_ip_dis_c(ncblk_opreo), idxms_ip_dis_a(nablk_opreo),
     &     ms_ip_dis_c(ncblk_opreo), ms_ip_dis_a(nablk_opreo),
     &     gm_ip_dis_c(ncblk_opreo), gm_ip_dis_a(nablk_opreo),
     &     lstr_k(ncblk_k+nablk_k), lstr_i0(ncblk_i0+nablk_i0),
     &     lstr_opreo(ncblk_opreo+nablk_opreo),
     &     nstr_opreo_c(ncblk_opreo), nstr_opreo_a(nablk_opreo),
     &     nstr_opori_c(ncblk_opori), nstr_opori_a(nablk_opori),
     &     ldim_opreo_c(ncblk_opreo), ldim_opreo_a(nablk_opreo),
     &     ldim_opori_c(ncblk_opori), ldim_opori_a(nablk_opori),
     &     nstr_i0c1(ncblk_opori), nstr_i0a1(nablk_opori),
     &     nstr_i0c2(ncblk_opreo), nstr_i0a2(nablk_opreo),
     &     nstr_k_c1(ncblk_opori), nstr_k_a1(nablk_opori),
     &     nstr_k_c2(ncblk_opreo), nstr_k_a2(nablk_opreo),
     &     ireo_i0c1(ncblk_i0), ireo_i0a1(nablk_i0),
     &     ireo_i0c2(ncblk_i0), ireo_i0a2(nablk_i0),
     &     ireo_k_c1(ncblk_k), ireo_k_a1(nablk_k),
     &     ireo_k_c2(ncblk_k), ireo_k_a2(nablk_k),
     &     istr_i0c1(ncblk_opori), istr_i0a1(nablk_opori),
     &     istr_i0c2(ncblk_opreo), istr_i0a2(nablk_opreo),
     &     istr_k_c1(ncblk_opori), istr_k_a1(nablk_opori),
     &     istr_k_c2(ncblk_opreo), istr_k_a2(nablk_opreo)

      integer, pointer ::
     &     map_to_ori_c(:), map_to_ori_a(:),
     &     map_to_reo_c(:), map_to_reo_a(:)

      integer, pointer ::
     &     ndis_opreo(:,:), d_gam_ms_opreo(:,:,:),
     &     mapca(:), diag_idx(:), diag_ca(:)

      type(graph), pointer ::
     &     graphs(:)

      type(operator), pointer ::
     &     opreo
c dbg
c      logical ::
c     &     first_element
c dbg
      logical ::
     &     ms_fix

      logical, external ::
     &     next_msgamdist2, check_ms, check_gm, nondia_blk, nondia_distr
      integer, external ::
     &     get_lenmap, ielprd, idx_msgmdst2, idxlist,
     &     msa2idxms4op

c dbg
c      print *,'gm_op_c/a: ',gm_op_c,gm_op_a
c dbg

      if (ntest.ge.1000) then
        call write_title(luout,wst_dbg_subr,'reo_blk_wmaps_c')
        write(luout,*) 'sign_reo = ',sign_reo
      end if

      graphs => str_info%g
      ngraph =  str_info%ngraph
      opreo => me_opreo%op

      ms_fix = me_opreo%fix_vertex_ms

      if (me_opreo%diag_type.ne.0) then
        allocate(mapca(ncblk_opreo),diag_idx(ncblk_opreo),
     &           diag_ca(ncblk_opreo))
        if (nondia_blk(mapca,diag_idx,diag_ca,
     &                 opreo%ihpvca_occ(1,1,
     &                        (iblk_opreo-1)*opreo%njoined+1),
     &                 opreo%njoined,ncblk_opreo,
     &                 nablk_opreo,me_opreo%diag_type))
     &       call quit(1,'reo_blk_wmaps_c','non-diagonal block!')
      end if

      idxst_opreo = me_opreo%off_op_occ(iblk_opreo)
      ndis_opreo  => me_opreo%off_op_gmox(iblk_opreo)%ndis
      d_gam_ms_opreo  => me_opreo%off_op_gmox(iblk_opreo)%d_gam_ms

      call sum_occ(ms_k_c_max,cinfo_k_c,ncblk_k)
      call sum_occ(ms_k_a_max,cinfo_k_a,nablk_k)
      call sum_occ(ms_i0_c_max,cinfo_i0_c,ncblk_i0)
      call sum_occ(ms_i0_a_max,cinfo_i0_a,nablk_i0)

cmh      idxms_op_a = (ms_k_a_max+ms_i0_a_max - ms_op_a)/2 + 1
      ! also for msa_max.ne.msc_max:
      idxms_op_a = msa2idxms4op(ms_op_a,ms_op_c-ms_op_a,
     &             ms_k_a_max+ms_i0_a_max,ms_k_c_max+ms_i0_c_max)

      do ms_k_a = ms_k_a_max, -ms_k_a_max, -2
        ms_i0_a = ms_op_a - ms_k_a
        if (abs(ms_i0_a).gt.ms_i0_a_max) cycle

        do ms_k_c = ms_k_c_max, -ms_k_c_max, -2
          ms_i0_c = ms_op_c - ms_k_c
          if (abs(ms_i0_c).gt.ms_i0_c_max) cycle

          do gm_k_a = 1, nsym
            gm_i0_a = multd2h(gm_k_a,gm_op_a)
            do gm_k_c = 1, nsym
              gm_i0_c = multd2h(gm_k_c,gm_op_c)
c dbg
c              print *,'ms_k_a,ms_k_c: ',ms_k_a,gm_k_c
c              print *,'ms_i0_a,ms_i0_c: ',ms_i0_a,gm_i0_c
c              print *,'gm_k_a,gm_k_c: ',gm_k_a,gm_k_c
c              print *,'gm_i0_a,gm_i0_c: ',gm_i0_a,gm_i0_c
c dbg

               ! loop over symmetry and ms-distributions of shift string K
              first = .true.
              dis_loop: do
                if (.not.next_msgamdist2(first,
     &             ms_k_dis_c,ms_k_dis_a,gm_k_dis_c,gm_k_dis_a,
     &             ncblk_k,nablk_k,
     &             cinfo_k_c,cinfo_k_a,
     &             ms_k_c,ms_k_a,gm_k_c,gm_k_a,nsym,
     &             ms_fix,fix_success))
     &             exit dis_loop

                first = .false.
                if(ms_fix.and..not.fix_success)cycle dis_loop

                ! get symmetry and ms-distr. of I0
                call split_msgmdis(ms_i0_dis_c,gm_i0_dis_c,
     &                             possible,
     &                             ms_k_dis_c,gm_k_dis_c,
     &                             ms_i_dis_c,gm_i_dis_c,
     &                             ncblk_opori,
     &                             map_info_to_ori_c,.false.)
                if (.not.possible) cycle dis_loop

                call split_msgmdis(ms_i0_dis_a,gm_i0_dis_a,
     &                             possible,
     &                             ms_k_dis_a,gm_k_dis_a,
     &                             ms_i_dis_a,gm_i_dis_a,
     &                             nablk_opori,
     &                             map_info_to_ori_a,.true.)
                if (.not.possible) cycle dis_loop


c dbg                
c                  write(luout,*) '            I   GMD C',
c     &                 gm_i_dis_c(1:ncblk_opori)
c                  write(luout,*) '            K   GMD C',
c     &                 gm_k_dis_c(1:ncblk_k)
c                  write(luout,*) '            I0  GMD C',
c     &                 gm_i0_dis_c(1:ncblk_i0)
c                  write(luout,*) '            I   GMD A',
c     &                 gm_i_dis_a(1:nablk_opori)
c                  write(luout,*) '            K   GMD A',
c     &                 gm_k_dis_a(1:nablk_k)
c                  write(luout,*) '            I0  GMD A',
c     &                 gm_i0_dis_a(1:nablk_i0)
c                  write(luout,*) '            I   MSD C',
c     &                 ms_i_dis_c(1:ncblk_opori)
c                  write(luout,*) '            K   MSD C',
c     &                 ms_k_dis_c(1:ncblk_k)
c                  write(luout,*) '            I0  MSD C',
c     &                 ms_i0_dis_c(1:ncblk_i0)
c                  write(luout,*) '            I   MSD A',
c     &                 ms_i_dis_a(1:nablk_opori)
c                  write(luout,*) '            K   MSD A',
c     &                 ms_k_dis_a(1:nablk_k)
c                  write(luout,*) '            I0  MSD A',
c     &                 ms_i0_dis_a(1:nablk_i0)
c                  write(luout,*) ' possible: ',possible
c dbg
                if (ntest.ge.1500) then
                  write(luout,*) '               I  MS',ms_op_c,ms_op_a
                  write(luout,*) '               MSD C',
     &                 ms_i_dis_c(1:ncblk_opori)
                  write(luout,*) '               MSD A',
     &                 ms_i_dis_a(1:nablk_opori)
                  write(luout,*) 'current block: I0 MS',ms_i0_c,ms_i0_a
                  write(luout,*) '               IRREP',gm_i0_c,gm_i0_a
                  write(luout,*) '               MSD C',
     &                 ms_i0_dis_c(1:ncblk_i0)
                  write(luout,*) '               MSD A',
     &                 ms_i0_dis_a(1:nablk_i0)
                  write(luout,*) '               GMD C',
     &                 gm_i0_dis_c(1:ncblk_i0)
                  write(luout,*) '               GMD A',
     &                 gm_i0_dis_a(1:nablk_i0)
                  write(luout,*) '               K MS ',ms_k_c,ms_k_a
                  write(luout,*) '               IRREP',gm_k_c,gm_k_a
                  write(luout,*) '               MSD C',
     &                 ms_k_dis_c(1:ncblk_k)
                  write(luout,*) '               MSD A',
     &                 ms_k_dis_a(1:nablk_k)
                  write(luout,*) '               GMD C',
     &                 gm_k_dis_c(1:ncblk_k)
                  write(luout,*) '               GMD A',
     &                 gm_k_dis_a(1:nablk_k)
                end if

                if (.not.check_ms(ms_i0_dis_c,cinfo_i0_c,ncblk_i0))
     &               cycle dis_loop
                if (.not.check_ms(ms_i0_dis_a,cinfo_i0_a,nablk_i0))
     &               cycle dis_loop
                if (.not.check_gm(gm_i0_dis_c,cinfo_i0_c,ncblk_i0))
     &               cycle dis_loop
                if (.not.check_gm(gm_i0_dis_a,cinfo_i0_a,nablk_i0))
     &               cycle dis_loop
c dbg
c                print *,'ACCEPTED'
c                first_element = .true.
c dbg                
  
                ! get symmetry and ms-dist. of reordered operator I'
c dbg
c                print *,'merge_map1:',map_info_to_ori_c(1:8)
c                print *,'merge_map2:',map_info_to_reo_c(1:8)
c dbg
                call merge_msgmdis(ms_ip_dis_c,gm_ip_dis_c,
     &                             ncblk_opreo,
     &                             ms_i0_dis_c,gm_i0_dis_c,
     &                             ms_k_dis_c,gm_k_dis_c,
     &                             map_info_to_reo_c)
                call merge_msgmdis(ms_ip_dis_a,gm_ip_dis_a,
     &                             nablk_opreo,
     &                             ms_k_dis_a,gm_k_dis_a,
     &                             ms_i0_dis_a,gm_i0_dis_a,
     &                             map_info_to_reo_a)
c dbg
c                write(luout,*) '            I''  MSD C',
c     &                 ms_ip_dis_c(1:ncblk_opreo)
c                write(luout,*) '            I''  MSD A',
c     &                 ms_ip_dis_a(1:nablk_opreo)
c
c dbg

                ! reform MS-(times two)-values to idxms
                call ms2idxms(idxms_k_dis_c,ms_k_dis_c,
     &               cinfo_k_c,ncblk_k)
                call ms2idxms(idxms_k_dis_a,ms_k_dis_a,
     &               cinfo_k_a,nablk_k)

                call ms2idxms(idxms_i0_dis_c,ms_i0_dis_c,
     &               cinfo_i0_c,ncblk_i0)
                call ms2idxms(idxms_i0_dis_a,ms_i0_dis_a,
     &               cinfo_i0_a,nablk_i0)

c dbg
c                print *,'cinfo_opreo_c: ',
c     &               cinfo_opreo_c(1:ncblk_opreo,1)
c                print *,'cinfo_opreo_a: ',
c     &               cinfo_opreo_a(1:nablk_opreo,1)
c                print *,'gm_ip_dis_c: ',gm_ip_dis_c(1:ncblk_opreo)
c                print *,'gm_ip_dis_a: ',gm_ip_dis_a(1:nablk_opreo)
c                print *,'gm_i0_dis_c: ',gm_i0_dis_c(1:ncblk_i0)
c                print *,'gm_i0_dis_a: ',gm_i0_dis_a(1:nablk_i0)
c                print *,'gm_k_dis_c:  ',gm_k_dis_c(1:ncblk_k)
c                print *,'gm_k_dis_a:  ',gm_k_dis_a(1:nablk_k)
c dbg
                call ms2idxms(idxms_ip_dis_c,ms_ip_dis_c,
     &               cinfo_opreo_c,ncblk_opreo)
                call ms2idxms(idxms_ip_dis_a,ms_ip_dis_a,
     &               cinfo_opreo_a,nablk_opreo)

                call set_len_str(lstr_k,ncblk_k,nablk_k,
     &               graphs,
     &               cinfo_k_c(1,2),idxms_k_dis_c,
     &                               gm_k_dis_c,cinfo_k_c(1,3),
     &               cinfo_k_a(1,2),idxms_k_dis_a,
     &                               gm_k_dis_a,cinfo_k_a(1,3),
     &               hpvxseq,.false.)

                if (ncblk_k+nablk_k.gt.0 .and.
     &              idxlist(0,lstr_k,ncblk_k+nablk_k,1).gt.0) cycle

                call set_len_str(lstr_i0,ncblk_i0,nablk_i0,
     &               graphs,
     &               cinfo_i0_c(1,2),idxms_i0_dis_c,
     &                               gm_i0_dis_c,cinfo_i0_c(1,3),
     &               cinfo_i0_a(1,2),idxms_i0_dis_a,
     &                               gm_i0_dis_a,cinfo_i0_a(1,3),
     &               hpvxseq,.false.)

                if (ncblk_i0+nablk_i0.gt.0 .and.
     &              idxlist(0,lstr_i0,ncblk_i0+nablk_i0,1).gt.0) cycle
c dbg
c                print *,'cinfo_i0_c: ',cinfo_i0_c(1:ncblk_i0,1)
c                print *,'cinfo_i0_a: ',cinfo_i0_a(1:nablk_i0,1)
c                print *,'cinfo_i0_c: ',cinfo_i0_c(1:ncblk_i0,2)
c                print *,'cinfo_i0_a: ',cinfo_i0_a(1:nablk_i0,2)
c                print *,'cinfo_i0_c: ',cinfo_i0_c(1:ncblk_i0,3)
c                print *,'cinfo_i0_a: ',cinfo_i0_a(1:nablk_i0,3)
c                print *,'lstr_i0:',lstr_i0
c                print *,'lstr_k:',lstr_k
c dbg

                call set_len_str(lstr_opreo,ncblk_opreo,nablk_opreo,
     &               graphs,
     &               cinfo_opreo_c(1,2),idxms_ip_dis_c,
     &                               gm_ip_dis_c,cinfo_opreo_c(1,3),
     &               cinfo_opreo_a(1,2),idxms_ip_dis_a,
     &                               gm_ip_dis_a,cinfo_opreo_a(1,3),
     &               hpvxseq,.false.)
c dbg
c                print *,'cinfo_opreo_c:',cinfo_opreo_c(1:ncblk_opreo,2)
c                print *,'cinfo_opreo_a:',cinfo_opreo_a(1:nablk_opreo,2)
c                print *,'lstr_opreo:   ',lstr_opreo
c dbg

                if (ncblk_opreo+nablk_opreo.gt.0 .and.
     &              idxlist(0,lstr_opreo,
     &                        ncblk_opreo+nablk_opreo,1).gt.0) cycle

                call set_strmapdim_c2(
     &                 nstr_opori_c,nstr_i0c1,nstr_k_c1,
     &                 ireo_i0c1,ireo_k_c1,
     &                 ncblk_opori,
     &                 lstr_i0,lstr_k,map_info_to_ori_c)
                call set_strmapdim_c2(
     &                 nstr_opori_a,nstr_k_a1,nstr_i0a1,
     &                 ireo_k_a1,ireo_i0a1,
     &                 nablk_opori,
     &                 lstr_k(ncblk_k+1),
     &                       lstr_i0(ncblk_i0+1),map_info_to_ori_a)

                call set_strmapdim_c2(
     &                 nstr_opreo_c,nstr_i0c2,nstr_k_c2,
     &                 ireo_i0c2,ireo_k_c2,
     &                 ncblk_opreo,
     &                 lstr_i0,lstr_k,map_info_to_reo_c)
                call set_strmapdim_c2(
     &                 nstr_opreo_a,nstr_k_a2,nstr_i0a2,
     &                 ireo_k_a2,ireo_i0a2,
     &                 nablk_opreo,
     &                 lstr_k(ncblk_k+1),
     &                       lstr_i0(ncblk_i0+1),map_info_to_reo_a)
c dbg
c                write(*,'(a,10i4)') 'ireo_k_c1: ',ireo_k_c1
c                write(*,'(a,10i4)') 'ireo_k_c2: ',ireo_k_c2
c                write(*,'(a,10i4)') 'ireo_k_a1: ',ireo_k_a1
c                write(*,'(a,10i4)') 'ireo_k_a2: ',ireo_k_a2
c                write(*,'(a,10i4)') 'ireo_i0c1: ',ireo_i0c1
c                write(*,'(a,10i4)') 'ireo_i0c2: ',ireo_i0c2
c                write(*,'(a,10i4)') 'ireo_i0a1: ',ireo_i0a1
c                write(*,'(a,10i4)') 'ireo_i0a2: ',ireo_i0a2
c dbgend

                nstr_k_c_tot  = ielprd(lstr_k,ncblk_k)
                nstr_k_a_tot  = ielprd(lstr_k(ncblk_k+1),nablk_k)
                nstr_i0_c_tot = ielprd(lstr_i0,ncblk_i0)
                nstr_i0_a_tot = ielprd(lstr_i0(ncblk_i0+1),nablk_i0)
c dbg
c                print *,'nstr_k_c_tot:  ',nstr_k_c_tot
c                print *,'nstr_k_a_tot:  ',nstr_k_a_tot
c                print *,'nstr_i0_c_tot: ',nstr_i0_c_tot
c                print *,'nstr_i0_a_tot: ',nstr_i0_a_tot
c dbg
c                print *,'ncblk_k:  ',ncblk_k
c                print *,'nablk_k:  ',nablk_k
c                print *,'ncblk_i0: ',ncblk_i0
c                print *,'nablk_i0: ',nablk_i0
c                write(*,'(a,10i4)') 'lstr_i0: ',lstr_i0
c                write(*,'(a,10i4)') 'lstr_k:  ',lstr_k
c                write(*,'(a,10i4)')'map_info_to_ori_c:',
c     &            map_info_to_ori_c(1:10)
c                write(*,'(a,10i4)')'map_info_to_ori_a:',
c     &            map_info_to_ori_a(1:10)
c                write(*,'(a,10i4)')'map_info_to_reo_c:',
c     &            map_info_to_reo_c(1:10)
c                write(*,'(a,10i4)')'map_info_to_reo_a:',
c     &            map_info_to_reo_a(1:10)
c dbgend
c dbg

                call set_op_ldim_c(ldim_opori_c,ldim_opori_a,
     &                 cinfo_opori_c(1,3),cinfo_opori_a(1,3),
     &                 lstr_opori,ncblk_opori,nablk_opori,tra_opori)
                call set_op_ldim_c(ldim_opreo_c,ldim_opreo_a,
     &                 cinfo_opreo_c(1,3),cinfo_opreo_a(1,3),
     &                 lstr_opreo,ncblk_opreo,nablk_opreo,tra_opreo)

                ifree = mem_setmark('reostr')
                lenmap = get_lenmap(lstr_i0,lstr_k,
     &               map_info_to_ori_c,ncblk_opori)
                ifree = mem_alloc_int(map_to_ori_c,lenmap,'orimap_c')
                lenmap = get_lenmap(
     &               lstr_k(ncblk_k+1),lstr_i0(ncblk_i0+1),
     &               map_info_to_ori_a,nablk_opori)
                ifree = mem_alloc_int(map_to_ori_a,lenmap,'orimap_a')
                lenmap = get_lenmap(lstr_i0,lstr_k,
     &               map_info_to_reo_c,ncblk_opreo)
                ifree = mem_alloc_int(map_to_reo_c,lenmap,'reomap_c')
                lenmap = get_lenmap(
     &               lstr_k(ncblk_k+1),lstr_i0(ncblk_i0+1),
     &               map_info_to_reo_a,nablk_opreo)
                ifree = mem_alloc_int(map_to_reo_a,lenmap,'reomap_a')

                call get_strmap_blk_c(map_to_ori_c,
     &                 ncblk_i0,ncblk_k,ncblk_opori,
     &                 cinfo_i0_c,cinfo_k_c,lstr_i0,lstr_k,
     &                 cinfo_i0_c(1,2),cinfo_k_c(1,2),
     &                                 cinfo_opori_c(1,2),
     &                 idxms_i0_dis_c,idxms_k_dis_c,
     &                 gm_i0_dis_c,gm_k_dis_c,map_info_to_ori_c,
     &                 strmap_info,nsym,ngraph)
                call get_strmap_blk_c(map_to_ori_a,
     &                 nablk_k,nablk_i0,nablk_opori,
     &                 cinfo_k_a,cinfo_i0_a,
     &                    lstr_k(ncblk_k+1),lstr_i0(ncblk_i0+1),
     &                 cinfo_k_a(1,2),cinfo_i0_a(1,2),
     &                                cinfo_opori_a(1,2),
     &                 idxms_k_dis_a,idxms_i0_dis_a,
     &                 gm_k_dis_a,gm_i0_dis_a,map_info_to_ori_a,
     &                 strmap_info,nsym,ngraph)

                call get_strmap_blk_c(map_to_reo_c,
     &                 ncblk_i0,ncblk_k,ncblk_opreo,
     &                 cinfo_i0_c,cinfo_k_c,lstr_i0,lstr_k,
     &                 cinfo_i0_c(1,2),cinfo_k_c(1,2),
     &                                 cinfo_opreo_c(1,2),
     &                 idxms_i0_dis_c,idxms_k_dis_c,
     &                 gm_i0_dis_c,gm_k_dis_c,map_info_to_reo_c,
     &                 strmap_info,nsym,ngraph)
                call get_strmap_blk_c(map_to_reo_a,
     &                 nablk_k,nablk_i0,nablk_opreo,
     &                 cinfo_k_a,cinfo_i0_a,
     &                   lstr_k(ncblk_k+1),lstr_i0(ncblk_i0+1),
     &                 cinfo_k_a(1,2),cinfo_i0_a(1,2),
     &                                cinfo_opreo_a(1,2),
     &                 idxms_k_dis_a,idxms_i0_dis_a,
     &                 gm_k_dis_a,gm_i0_dis_a,map_info_to_reo_a,
     &                 strmap_info,nsym,ngraph)
c dbg
c                write(*,'(a,10i4)') 'map_to_ori_c: ',map_to_ori_c
c                write(*,'(a,10i4)') 'map_to_ori_a: ',map_to_ori_a
c                write(*,'(a,10i4)') 'map_to_reo_c: ',map_to_reo_c
c                write(*,'(a,10i4)') 'map_to_reo_a: ',map_to_reo_a
c dbgend


                if (me_opreo%diag_type.ne.0) then
                  ! skip non-diagonal distributions ...
                  if (nondia_distr(mapca,diag_idx,diag_ca,
     &                   ms_ip_dis_c,ms_ip_dis_a,
     &                   gm_ip_dis_c,gm_ip_dis_a,
     &                   ncblk_opreo,me_opreo%msdiag,
     &                   me_opreo%gamdiag)) then
                    ifree = mem_flushmark('reostr')
                    cycle
                  end if
                end if

                ! --> offset in xop_reo
                if (ndis_opreo(gm_op_a,idxms_op_a).gt.1) then
                  idxdis =
     &                 idx_msgmdst2(
     &                 iblk_opreo,idxms_op_a,gm_op_a,
     &                 cinfo_opreo_c,idxms_ip_dis_c,
     &                               gm_ip_dis_c,ncblk_opreo,
     &                 cinfo_opreo_a,idxms_ip_dis_a,
     &                               gm_ip_dis_a,nablk_opreo,
     &                 tra_opreo,-1,-1,me_opreo,nsym)
c     &                 .false.,me_opreo,nsym)
                  idx00opreo =
     &                 d_gam_ms_opreo(idxdis,gm_op_a,idxms_op_a) + 1
     &                                           - idxst_opreo
c dbg
c                  print *,'idxms_ip_dis_c: ',idxms_ip_dis_c
c                  print *,'idxms_ip_dis_a: ',idxms_ip_dis_a
c                  print *,'REO: ',idxdis,gm_op_a,idxms_op_a,'->',
c     &                 idx00opreo
c dbg
                else
                  idx00opreo =
     &                 d_gam_ms_opreo(1,gm_op_a,idxms_op_a) + 1
     &                                           - idxst_opreo
c dbg
c                  idxdis = 1
c                  print *,'REO: ',idxdis,gm_op_a,idxms_op_a,'->',
c     &                 idx00opreo
c dbg
                end if
c dbg
c                print *,'idx00opreo:',idx00opreo
c dbg

                ! loop over A strings
                k_a: do istr_k_a = 1, nstr_k_a_tot

                  ! break down into components
                  istr_k_a1(1:nablk_opori) = 1
                  istr_k_a2(1:nablk_opreo) = 1
                  istr1 = istr_k_a-1
                  do icmp = 1, nablk_k
                    idx1 = mod(istr1,lstr_k(ncblk_k+icmp))+1
                    istr_k_a1(ireo_k_a1(icmp)) = idx1
                    istr_k_a2(ireo_k_a2(icmp)) = idx1
                    istr1 = istr1/lstr_k(ncblk_k+icmp)
                  end do

                  i0_a: do istr_i0_a = 1, nstr_i0_a_tot

                    ! break down into components
                    istr_i0a1(1:nablk_opori) = 1
                    istr_i0a2(1:nablk_opreo) = 1
                    istr1 = istr_i0_a-1
                    do icmp = 1, nablk_i0
                      idx1 = mod(istr1,lstr_i0(ncblk_i0+icmp))+1
                      istr_i0a1(ireo_i0a1(icmp)) = idx1
                      istr_i0a2(ireo_i0a2(icmp)) = idx1
                      istr1 = istr1/lstr_i0(ncblk_i0+icmp)
                    end do

                    ! map K,I0 -> I
                    idx0opori = 1
                    ioff = 0
                    isgna = sign_reo
c dbg
c                    print *,'+-------------ORI-A----------------+'
c dbg
                    do icmp = 1, nablk_opori
                      idx1 = istr_i0a1(icmp)
                      idx2 = istr_k_a1(icmp)
                      idx  = (idx1-1)*nstr_k_a1(icmp)+idx2
c dbg
c                      print *,'icmp, ioff, idx:',icmp,ioff,idx
c dbgend
                      ielmap = map_to_ori_a(ioff+idx)
                      if (ielmap.eq.0) cycle i0_a
                      isgna = isgna*sign(1,ielmap)
                      idx0opori = idx0opori
     &                     + (abs(ielmap)-1)*ldim_opori_a(icmp)
c dbg
c                      print *,'icmp, idx_i0a, idx_k_a: ',icmp,idx1,idx2
c                      print *,'    ->idx_i_a: ',ielmap                      
c dbg
                      ioff = ioff + nstr_opori_a(icmp)
                    end do
                    
                    ! map K,I0 -> I'
                    idx0opreo = idx00opreo
                    ioff = 0
c                    do icmp = 1, nablk_opreo
c                      idiv2 = idiv2*nstr_k_a2(icmp)
c                    end do

c dbg
c                    print *,'+-------------REO-A----------------+'
c dbg
                    do icmp = 1, nablk_opreo
                      idx1 = istr_i0a2(icmp)
                      idx2 = istr_k_a2(icmp)
                      idx  = (idx1-1)*nstr_k_a2(icmp)+idx2
c dbg
c                      if (idx00opreo.eq.1) then
c                      print *,'istr2,idiv2,idx2:',istr2,idiv2,idx2
c                      print *,'idx1,idx2,idx,ioff:',idx1,idx2,idx,ioff
c                      print *,'map = ',map_to_reo_a(ioff+idx)
c                      end if
c dbg
c dbg
c                      print *,'icmp, ioff, idx:',icmp,ioff,idx
c dbgend
                      ielmap = map_to_reo_a(ioff+idx)
                      if (ielmap.eq.0) cycle i0_a
                      isgna = isgna*sign(1,ielmap)
                      idx0opreo = idx0opreo
     &                     + (abs(ielmap)-1)*ldim_opreo_a(icmp)
c dbg
c                      print *,'icmp, idx_i0a, idx_k_a: ',icmp,idx1,idx2
c                      print *,'    ->idx_ipa: ',ielmap
c dbg
                      ioff = ioff + nstr_opreo_a(icmp)
                    end do
c dbg
c                    print *,'+----------------------------------+'
c
c                    if (idx0opreo.eq.1.or.idx0opreo.eq.779)
c     &                   print *,'HIER HIER HIER: idx0opreo = ',
c     &                   idx0opreo
c dbg
c dbg
c                    print *,'idx0opreo: ',idx0opreo
c                    print *,'idx0opori: ',idx0opori
c dbg
                     
                    ! loop over C strings
                    k_c: do istr_k_c = 1, nstr_k_c_tot
c dbg
c                      print *,'istr_k_c = ',istr_k_c
c dbgend

                      ! break down into components
                      istr_k_c1(1:ncblk_opori) = 1
                      istr_k_c2(1:ncblk_opreo) = 1
                      istr1 = istr_k_c-1
                      do icmp = 1, ncblk_k
                        idx1 = mod(istr1,lstr_k(icmp))+1
                        istr_k_c1(ireo_k_c1(icmp)) = idx1
                        istr_k_c2(ireo_k_c2(icmp)) = idx1
                        istr1 = istr1/lstr_k(icmp)
c dbg
c                        print *,'icmp/idx1: ',icmp,idx1
c dbgend
                      end do

                      i0_c: do istr_i0_c = 1, nstr_i0_c_tot
c dbg
c                        print *,'istr_i0_c = ',istr_i0_c
c dbgend

                        ! break down into components
                        istr_i0c1(1:ncblk_opori) = 1
                        istr_i0c2(1:ncblk_opreo) = 1
                        istr1 = istr_i0_c-1
                        do icmp = 1, ncblk_i0
                          idx1 = mod(istr1,lstr_i0(icmp))+1
                          istr_i0c1(ireo_i0c1(icmp)) = idx1
                          istr_i0c2(ireo_i0c2(icmp)) = idx1
                          istr1 = istr1/lstr_i0(icmp)
                        end do

                        ! map K,I0 -> I
                        idx_opori = idx0opori
                        ioff = 0
                        isgnt = isgna
c dbg
c                        print *,'+-------------ORI-C----------------+'
c dbg
                        do icmp = 1, ncblk_opori
                          idx1 = istr_k_c1(icmp)
                          idx2 = istr_i0c1(icmp)
                          idx = (idx1-1)*nstr_i0c1(icmp)+idx2
c dbg
c                          print *,'icmp, ioff, idx:',icmp,ioff,idx
c dbgend
                          ielmap = map_to_ori_c(ioff+idx)
                          if (ielmap.eq.0) cycle i0_c
                          isgnt = isgnt*sign(1,ielmap)
                          idx_opori = idx_opori +
     &                         (abs(ielmap)-1)*ldim_opori_c(icmp)
c dbg
c                          print *,'    ->idx_i_c: ',ielmap
c dbgend
                          ioff = ioff + nstr_opori_c(icmp)
                        end do
                        ! map K,I0 -> I'
                        idx_opreo = idx0opreo
                        ioff = 0
c dbg
c                        print *,'+-------------REO-C----------------+'
c dbg
                        do icmp = 1, ncblk_opreo
                          idx1 = istr_k_c2(icmp)
                          idx2 = istr_i0c2(icmp)
                          idx = (idx1-1)*nstr_i0c2(icmp)+idx2
c dbg
c                          print *,'icmp, ioff, idx:',icmp,ioff,idx
c dbgend
                          ielmap = map_to_reo_c(ioff+idx)
                          if (ielmap.eq.0) cycle i0_c
                          isgnt = isgnt*sign(1,ielmap)
                          idx_opreo = idx_opreo +
     &                         (abs(ielmap)-1)*ldim_opreo_c(icmp)
c dbg
c                          print *,'    ->idx_ipc: ',ielmap              
c dbgend
                          ioff = ioff + nstr_opreo_c(icmp)
                        end do
c dbg
c                        if (!first_element ) then!.or.
c     &                       idx_opreo.eq.64 .or.
c     &                       idx_opreo.eq.253) then
c                          print *,'ori, reo, +/-A, +/-CA, val: ',
c     &                         idx_opori,idx_opreo,
c     &                         isgna,isgnt,xop_ori(idx_opori)
c                          first_element = .false.
c                        end if
c dbg
c dbg
c                        print *,'idx_opreo, idx_opori: ',
c     &                          idx_opreo,idx_opori
c dbgend

c dbg
c                        if (idx_opreo.lt.0.or.idx_opreo.gt.len_reo.or.
c     &                      idx_opori.lt.0.or.idx_opori.gt.len_ori) then
c                          print *,'RANGE dx_opreo/ori: ',
c     &                            idx_opreo,idx_opori
c                          print *,'max: ',len_reo,len_ori
c                        end if
c dbg
                        xop_reo(idx_opreo) = xop_reo(idx_opreo)
     &                       + dble(isgnt)*xop_ori(idx_opori)

                      end do i0_c
                    end do k_c

                  end do i0_a
                end do k_a

                ifree = mem_flushmark('reostr')

              end do dis_loop
              
            end do
          end do

        end do
      end do

      if (me_opreo%diag_type.ne.0) deallocate(mapca,diag_idx,diag_ca)

      return
      end
