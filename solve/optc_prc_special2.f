*----------------------------------------------------------------------*
      subroutine optc_prc_special2(me_grd,me_special,nspecial,
     &                            nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                            orb_info,str_info)
*----------------------------------------------------------------------*
*     experimental routine for testing special preconditioners
*     if nincore>1: xbuf1 contains the gradient vector on entry
*                   and should contain the preconditioned gradient
*                   on exit
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_orbinf.h'
      include 'hpvxseq.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'multd2h.h'
      include 'ifc_input.h'
      include 'par_opnames_gen.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nincore, lenbuf, nspecial
      type(me_list) ::
     &     me_grd
      type(me_list_array) ::
     &     me_special(nspecial)
      real(8), intent(inout) ::
     &     xbuf1(lenbuf), xbuf2(lenbuf), xbuf3(lenbuf)
      type(orbinf), intent(in), target ::
     &     orb_info
      type(strinf),intent(in), target ::
     &     str_info

      integer ::
     &     ngrd, nblk_grd, nbmat, nxmat, nfmat,
     &     idxmsa, msa, msc, gama, gamc,  ngam,
     &     msc_max, msa_max, idxdis,  idx,
     &     ncsub, nasub,
     &     nidx_cstr, nidx_astr, ms_cstr, ms_astr, gam_cstr, gam_astr,
     &     graph_cstr, graph_astr, len_cstr, len_astr,
     &     idx_grd, idx_b, idx_x,
     &     ngas, njoined, mstotal, gamtotal,
     &     iblk, iblk_off, idxms_bx, gam_bx, ld_bx, iblk_b,
     &     occ_b(ngastp,2)
      logical ::
     &     first, beyond_A

      type(me_list), pointer ::
     &     me_bmat, me_xmat, me_fmat
      type(graph), pointer ::
     &     graphs(:)

      integer ::
     &     f_hpvx(ngastp,2), msdst(ngastp,2), igamdst(ngastp,2),
     &     ioff_xsum(2*ngastp), nstr(2*ngastp), nincr(2*ngastp)

      integer, pointer ::
     &     hpvx_occ(:,:,:), idx_graph(:,:,:),
     &     ca_occ(:,:), occ_blk(:,:,:), igas_restr(:,:,:,:,:),
     &     off_grd_d_gam_ms(:,:,:), len_grd_d_gam_ms(:,:,:),
     &     off_bx_gam_ms(:,:),
     &     mostnd(:,:,:), igamorb(:), idx_gas(:), ngas_hpv(:),
     &     hpvx_csub(:), hpvx_asub(:), occ_csub(:), occ_asub(:),
     &     graph_csub(:), graph_asub(:), msdis_c(:), msdis_a(:),
     &     idxmsdis_c(:), idxmsdis_a(:), gamdis_c(:), gamdis_a(:),
     &     len_str(:)
      real(8), pointer ::
     &     f_dia(:), xmat(:), bmat(:)

      integer, external ::
     &     idxlist, iblk_occ
      logical,external ::
     &     next_msgamdist2

      if (nincore.ne.3)
     &     call quit(1,'optc_prc_special',
     &     'currently: special preconditioning for nincore==3, only')

      if (nspecial.lt.1)
     &     call quit(1,'optc_prc_special',
     &     'nspecial must be >= 1')

      me_bmat => me_special(1)%mel
      nbmat = me_bmat%len_op
      allocate(bmat(nbmat))
      ! should be open
      call vec_from_da(me_bmat%fhand,1,bmat,nbmat)

      beyond_A = me_bmat%op%name.ne.op_b_inv

      if (beyond_A) then
        if (nspecial.lt.3)
     &       call quit(1,'optc_prc_special','nspecial must be >= 3')

        me_xmat => me_special(2)%mel
        nxmat = me_xmat%len_op

        allocate(xmat(nxmat))
        ! should be open as well
        call vec_from_da(me_xmat%fhand,1,xmat,nxmat)

        me_fmat => me_special(3)%mel
        nfmat = me_fmat%len_op

        allocate(f_dia(2*orb_info%ntoob))

        ! xtract diagonal
        call onedia_from_op(f_dia,me_fmat,orb_info)
        
      else
        allocate(xmat(1),f_dia(1))
      end if
      
      ngam = orb_info%nsym
      ngas = orb_info%ngas
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      idx_gas => orb_info%idx_gas
      ngas_hpv => orb_info%ngas_hpv

      nblk_grd = me_grd%op%n_occ_cls
      njoined  = me_grd%op%njoined
      hpvx_occ => me_grd%op%ihpvca_occ
      idx_graph => me_grd%idx_graph
      ca_occ => me_grd%op%ica_occ
      mstotal = me_grd%mst
      gamtotal = me_grd%gamt
      graphs => str_info%g
      igas_restr => str_info%igas_restr

c      if (njoined.ne.2)
c     &     call quit(1,'optc_prc_special',
c     &     'strange -- expected njoined==2')

      ! loop over occupations of GRD
      do iblk = 1, nblk_grd
        iblk_off = (iblk-1)*njoined

        occ_blk =>
     &       hpvx_occ(1:ngastp,1:2,iblk_off+1:iblk_off+njoined)
        len_grd_d_gam_ms =>
     &       me_grd%len_op_gmox(iblk)%d_gam_ms
        off_grd_d_gam_ms =>
     &       me_grd%off_op_gmox(iblk)%d_gam_ms

        if (ntest.ge.100) then
          write(luout,*) 'Now caring for GRD block: '
          call wrt_occ_n(luout,occ_blk,njoined)
        end if

        call get_num_subblk(ncsub,nasub,
     &       hpvx_occ(1,1,iblk_off+1),njoined)

        ! find the relevant block of B
        if (njoined.eq.2) then
          iblk_b = 1
        else if (njoined.eq.1) then
          occ_b(1:ngastp,1) = occ_blk(1:ngastp,1,1)
          occ_b(1:ngastp,2) = occ_blk(1:ngastp,1,1)
          iblk_b = iblk_occ(occ_b,.false.,me_bmat%op)
          if (iblk_b.le.0)
     &         call quit(1,'optc_prc_special2',
     &         'did not find an appropriate block of B (precond)')
        else
          call quit(1,'optc_prc_special2',
     &         'gradient -- njoined > 2 ??')
        end if

        off_bx_gam_ms => me_bmat%off_op_gmo(iblk_b)%gam_ms

        ! special case: scalar:
        if (ncsub.eq.0.and.nasub.eq.0) then
          idx_grd = off_grd_d_gam_ms(1,1,1)+1
          idx_b   = off_bx_gam_ms(1,1)+1
          xbuf1(idx_grd) = xbuf1(idx_grd)/bmat(idx_b)
          cycle
        end if

        if (ncsub.ne.1.and.ncsub.ne.2.and.nasub.ne.1) then
          write(luout,*) 'ncsub, nasub: ',ncsub,nasub
          call quit(1,'optc_prc_special2','this is not what I expected')
        end if

        allocate(hpvx_csub(ncsub),hpvx_asub(nasub),
     &           occ_csub(ncsub), occ_asub(nasub),
     &           graph_csub(ncsub), graph_asub(nasub),
     &           msdis_c(ncsub),  msdis_a(nasub),
     &           idxmsdis_c(ncsub),  idxmsdis_a(nasub),
     &           gamdis_c(ncsub), gamdis_a(nasub),
     &           len_str(ncsub+nasub))

        msc_max = ca_occ(1,iblk)
        msa_max = ca_occ(2,iblk)

        ! set HPVX and OCC info
        call condense_occ(occ_csub, occ_asub,
     &                    hpvx_csub,hpvx_asub,
     &                    hpvx_occ(1,1,iblk_off+1),njoined,hpvxblkseq)
        ! do the same for the graph info
        call condense_occ(graph_csub, graph_asub,
     &                    hpvx_csub,hpvx_asub,
     &                    idx_graph(1,1,iblk_off+1),njoined,hpvxblkseq)

        ! loop over MS(A)
        idxmsa = 0
        msa_loop: do msa = msa_max, -msa_max, -2

          msc = msa + mstotal
          
          if (abs(msc).gt.msc_max) cycle msa_loop
          idxmsa = idxmsa+1

          gama_loop: do gama = 1, ngam
            
            gamc = multd2h(gama,gamtotal)

            if (ntest.ge.100)
     &           write(luout,*) 'current MS(A), MS(C), GAM(A), GAM(C):',
     &                           msa, msc, gama, gamc

            idxdis = 0
            first = .true.
            distr_loop: do

              if (.not.next_msgamdist2(first,
     &             msdis_c,msdis_a,gamdis_c,gamdis_a,
     &             ncsub,nasub,
     &             occ_csub,occ_asub,
     &             msc,msa,gamc,gama,ngam)) exit distr_loop

              first = .false.
              
              call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncsub)
              call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nasub)

              call set_len_str(len_str,ncsub,nasub,
     &                         graphs,
     &                         graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                         graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                         hpvxseq,.false.)

              if (idxlist(0,len_str,ncsub+nasub,1).gt.0)
     &             cycle distr_loop

              idxdis = idxdis+1

c test -- special insert 
              if (njoined.eq.1) then
                idx = idxlist(IPART,hpvx_csub,ncsub,1)
                if (idx.ne.1) stop '???'
                idxms_bx  = idxmsdis_c(idx)
                gam_bx = gamdis_c(idx)
                ld_bx  = len_str(idx)
                idx_b = off_bx_gam_ms(gam_bx,idxms_bx)+1
                len_astr = len_str(2)
                idx_grd = off_grd_d_gam_ms(idxdis,gama,idxmsa)+1
                call optc_prc_test(xbuf1(idx_grd),
     &               xbuf2,xbuf3,bmat(idx_b),xmat(idx_b),ld_bx,len_astr)

                cycle
              end if
              

              ! Gamma and Ms of B and X
              idx = idxlist(IHOLE,hpvx_csub,ncsub,1)
              if (idx.le.0)
     &             call quit(1,'optc_prc_special','no HOLE??')
              idxms_bx  = idxmsdis_c(idx)
              gam_bx = gamdis_c(idx)
              ld_bx  = len_str(idx)

              if (ld_bx.le.0) cycle

              ! occupation, Gamma and Ms of C-string
              idx = idxlist(IPART,hpvx_csub,ncsub,1)
              if (idx.le.0) then
                nidx_cstr = 0
                ms_cstr = 0
                gam_cstr = 1
                len_cstr = 1
                graph_cstr = 1 ! dummy
              else
                nidx_cstr = occ_csub(idx)
                ms_cstr   = msdis_c(idx)
                gam_cstr  = gamdis_c(idx)
                len_cstr = len_str(idx)
                graph_cstr = graph_csub(idx)
              end if

              if (len_cstr.le.0) cycle

              ! occupation, Gamma and Ms of A-string
              ! idx = 1,  we know that for sure
              nidx_astr = occ_asub(1)
              ms_astr   = msdis_a(1)
              gam_astr  = gamdis_a(1)
              len_astr = len_str(ncsub+1)
              graph_astr = graph_asub(1)

              if (len_astr.le.0) cycle

              ! offset of GRD
              idx_grd = off_grd_d_gam_ms(idxdis,gama,idxmsa)+1

              ! offset of appropriate block of B (and X)
              idx_b = off_bx_gam_ms(gam_bx,idxms_bx)+1
              idx_x = idx_b
              if (.not.beyond_A) idx_x = 1

              call optc_prc_special2_inner
     &             (xbuf1(idx_grd), beyond_A,
     &              xbuf2,xbuf3,bmat(idx_b),xmat(idx_x),f_dia,
     &              ld_bx,len_cstr, len_astr,
     &              nidx_cstr,ms_cstr,gam_cstr,
     &               igas_restr(1,1,1,1,graph_cstr),
     &               mostnd(1,1,idx_gas(IPART)),ngas_hpv(IPART),
     &              nidx_astr,ms_astr,gam_astr,
     &               igas_restr(1,1,1,1,graph_astr),
     &               mostnd(1,1,idx_gas(IHOLE)),ngas_hpv(IHOLE),
     &              igamorb,ngam,ngas)


            end do distr_loop

          end do gama_loop

        end do msa_loop


      end do

      deallocate(f_dia,xmat,bmat)

      return
      end


      subroutine optc_prc_test(grd,
     &     buf1,buf2,bmat,xmat,lenp,lenh)

      implicit none

      integer, intent(in) ::
     &     lenp, lenh
      real(8), intent(in) ::
     &     bmat(lenp,lenp), xmat(lenp,lenp)
      real(8), intent(inout) ::
     &     buf1(lenp,lenh), buf2(lenp,lenh)
      real(8), intent(inout) ::
     &     grd(lenp,lenh)

      real(8), parameter ::
     &     thrsh = 1d-12

      real(8), pointer ::
     &     scr(:,:), umat(:,:), vtmat(:,:), wrk(:), singval(:)

      integer ::
     &     idx, jdx, kdx, info, lwrk

      lwrk=max(1024,lenp*max(lenh,lenp))
      allocate(scr(lenp,lenp),umat(lenp,lenp),vtmat(lenp,lenp),
     &     wrk(lwrk),singval(lenp))

      scr = bmat + 10*xmat

      call dgesvd('A','A',lenp,lenp,
     &     scr,lenp,singval,
     &     umat,lenp,vtmat,lenp,
     &     wrk,lwrk,info)
c dbg
      print *,'singval: ',singval
c dbg

      do jdx = 1, lenh
        do idx = 1, lenp
          wrk(idx) = 0d0
          if (abs(singval(idx)).lt.thrsh) cycle
          do kdx = 1, lenp
            wrk(idx) = wrk(idx)+
     &           umat(kdx,idx)*grd(kdx,jdx)
          end do
          wrk(idx) = wrk(idx)/singval(idx)
        end do
        
        do idx = 1, lenp
          buf1(idx,jdx) = 0d0
          do kdx = 1, lenp
            buf1(idx,jdx) = buf1(idx,jdx) +
     &           vtmat(kdx,idx)*wrk(kdx)
          end do
        end do

      end do
      
c      call gaussj(scr,lenp,lenp)

c      do jdx = 1, lenh
c        do idx = 1, lenp
c          buf1(idx,jdx) = 0d0
c          do kdx = 1, lenp
c            buf1(idx,jdx) = buf1(idx,jdx)+
c     &           scr(kdx,idx)*grd(kdx,jdx)
c          end do
c        end do
c      end do

      grd = buf1

      deallocate(scr,wrk,umat,vtmat,singval)
      
      return
      end
