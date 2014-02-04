*----------------------------------------------------------------------*
      subroutine spin_prj_blk(buffer_out,
     &     buffer_in,fac,s2,
     &     mel,iblk,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     core routine for spin projector
*
*     buffer_in/buffer_out may be the same
*
*     matthias and andreas, april 2013
*     
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'
      include 'hpvxseq.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(me_list), intent(inout) ::
     &     mel
      integer, intent(in) ::
     &     iblk, s2
      real(8), intent(inout) ::
     &     buffer_in(*), buffer_out(*)
      real(8), intent(in) ::
     &     fac

      logical ::
     &     first, ms_fix, ms_fix_ok, ok
      integer ::
     &     nocc_cls, njoined, nsplc, isplc, 
     &     ifree, nblk, nbuff, idxmsa, idxmsc, ioff0, ioff,
     &     msa, msc, igama, igamc, idxa, idxc, ngam, lena, lenc,
     &     iblkoff, ncblk, nablk, msc_max, msa_max,
     &     istr, istr_csub0, istr_asub0, icmp, ngraph, maxbuf,
     &     imap, ld_spprj_c, ld_spprj_a, idx, nmap, idxmap,
     &     nel, nn, is2, ncoup, naa, nca, nomni,
     &     omnistrings(4)
      real(8) ::
     &     fac_off, fac_dia, relfac2

      type(graph), pointer ::
     &     graphs(:)
      type(operator), pointer ::
     &     op

      integer, pointer ::
     &     hpvx_csub(:),hpvx_asub(:),
     &     occ_csub(:), occ_asub(:),
     &     graph_csub(:), graph_asub(:),
     &     msdis_c(:),  msdis_a(:),
     &     msdis_c_tgt(:),  msdis_a_tgt(:),
     &     idxmsdis_c(:),  idxmsdis_a(:),
     &     gamdis_c(:), gamdis_a(:),
     &     len_str(:), gsign(:), csign(:), asign(:), 
     &     offsets(:), maps_c(:,:), maps_a(:,:), nmaps_c(:), nmaps_a(:),
     &     hpvx_occ(:,:,:), ca_occ(:,:), idx_graph(:,:,:),
     &     ldim_op_c(:,:), ldim_op_a(:,:),
     &     istr_csub(:,:), istr_asub(:,:), idxel(:),
     &     branch(:), psign(:)
      integer, pointer ::
     &     spprjmap_c(:), spprjmap_a(:)
      real(8), pointer ::
     &     value(:), coeff(:,:), overl(:)

      real(8), external ::
     &     ddot
      logical, external ::
     &     iocc_equal_n, next_msgamdist2, list_cmp
      integer, external ::
     &     ielprd, idx_msgmdst2, idx_str_blk3, std_spsign_msdis,
     &     msa2idxms4op, ibico, idxcount
      integer, parameter ::
     &     popcnt_par(16) = (/+1,-1,-1,+1,-1,+1,+1,-1,
     &                        -1,+1,+1,-1,+1,-1,-1,+1/)

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'spin_prj_blk')
      end if

      ifree = mem_setmark('spin_prj_blk')

      if (s2.ne.0) call quit(1,'spin_prj_blk','only singlets!')
      ! the following is redundant, but I don't care ...
      if (mel%mst.ne.0) then
        call quit(1,'spin_prj_blk',
     &       'ms.ne.0 not yet considered ')
      end if

      ms_fix = mel%fix_vertex_ms

      op  => mel%op

      ioff0 = mel%off_op_occ(iblk)

      njoined  = op%njoined
      iblkoff = (iblk-1)*njoined

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ngraph = str_info%ngraph

      hpvx_occ => op%ihpvca_occ
      idx_graph => mel%idx_graph
      ca_occ => op%ica_occ
      graphs => str_info%g

      call get_num_subblk(ncblk,nablk,
     &     op%ihpvca_occ(1,1,iblkoff+1),njoined)

      msc_max = ca_occ(1,iblk)
      msa_max = ca_occ(2,iblk)
      nel = msc_max + msa_max

      if (msc_max.ne.msa_max) 
     &      call quit(1,'spin_prj_blk',
     &                  'only particle-conserving operators!')

      ! just scalar? Delete if not singlet
      if (msc_max.eq.0.and.msa_max.eq.0) then
        if (s2.eq.0) then
          buffer_out(1) = buffer_in(1)
        else
          buffer_out(1) = 0d0
        end if
        ifree = mem_flushmark('spin_prj_blk')
        return
      end if

      ! walk through branching diagram to get number of coupling pathes
      allocate(branch(0:nel)) ! maximum 2*S
      branch(0) = 1
      branch(1:nel) = 0
      do nn = 1, nel
        if (mod(nn,2).eq.0) branch(0) = branch(1)
        branch(nn) = branch(nn-1) ! always = 1
        do is2 = mod(nn+1,2)+1, nn-2, 2
          branch(is2) = branch(is2-1) + branch(is2+1)
        end do
      end do
      ncoup = branch(s2)
      deallocate(branch)

      ! get number of possible alpha/beta permutations
      nsplc = 0
      if (mod(mel%mst+msc_max-msa_max,2).ne.0)
     &   call quit(1,'spin_prj_blk','impossible 2Ms value')
      !loop over distributions
      do naa = 0, msa_max
        nca = (mel%mst+2*naa+msc_max-msa_max)/2
        if (nca.lt.0.or.nca.gt.msc_max) cycle
        ! add number of permutations
        nsplc = nsplc + ibico(msc_max,nca) * ibico(msa_max,naa)
      end do

      allocate(hpvx_csub(ncblk),hpvx_asub(nablk),
     &         occ_csub(ncblk), occ_asub(nablk),
     &         graph_csub(ncblk), graph_asub(nablk),
     &         msdis_c(ncblk),  msdis_a(nablk),
     &         msdis_c_tgt(ncblk),  msdis_a_tgt(nablk),
     &         idxmsdis_c(ncblk),  idxmsdis_a(nablk),
     &         gamdis_c(ncblk), gamdis_a(nablk),
     &         len_str(ncblk+nablk), idxel(nsplc), value(nsplc),
     &         gsign(nsplc), csign(nsplc), asign(nsplc),
     &         offsets(nsplc), maps_c(nsplc,ncblk), maps_a(nsplc,nablk),
     &         nmaps_c(ncblk), nmaps_a(nablk),
     &         istr_csub(ncblk,nsplc),istr_asub(nablk,nsplc),
     &         ldim_op_c(ncblk,nsplc),ldim_op_a(nablk,nsplc),
     &         coeff(nsplc,ncoup),overl(ncoup),psign(nsplc))

      ! set HPVX and OCC info
      call condense_occ(occ_csub, occ_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  hpvx_occ(1,1,iblkoff+1),njoined,hpvxblkseq)
      ! do the same for the graph info
      call condense_occ(graph_csub, graph_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  idx_graph(1,1,iblkoff+1),njoined,hpvxblkseq)

      ! We will loop over the elements of the "sceleton part"
      ! of the ME list, i.e. either msa=msc=0 (for 2-body, 4-body, 
      ! ... operators) or msa=msc=1 (for 1-body, 3-body, ...
      ! operators). Needs to be generalized, if needed.
      idxmsa = 0

      msa = mod(msa_max,2)
      msc = msa + mel%mst ! <- should be zero (checked above)
      
      if (ntest.ge.100) then
        write(lulog,*) 'occupation summary for block: ', iblk
        write(lulog,*) 'occ_csub   = ',occ_csub(1:ncblk)
        write(lulog,*) 'occ_asub   = ',occ_asub(1:nablk)
        write(lulog,*) 'hpvx_csub  = ',hpvx_csub(1:ncblk)
        write(lulog,*) 'hpvx_asub  = ',hpvx_asub(1:nablk)
        write(lulog,*) 'graph_csub = ',graph_csub(1:ncblk)
        write(lulog,*) 'graph_asub = ',graph_asub(1:nablk)
        write(lulog,*) 'ncoup      = ',ncoup
        write(lulog,*) 'nsplc      = ',nsplc
      end if

      ! our target distributions
      call get_target_dis(ncblk,occ_csub,msdis_c_tgt)
      call get_target_dis(nablk,occ_asub,msdis_a_tgt)

      if (ntest.ge.100) then
        write(lulog,*) 'our target ms distribution: '
        write(lulog,*) ' C: ',msdis_c_tgt(1:ncblk)
        write(lulog,*) ' A: ',msdis_a_tgt(1:nablk)
      end if

      idxmsa = msa2idxms4op(msa,mel%mst,msa_max,msc_max)
      idxmsc = (msc_max-msc)/2 + 1

      ! set spin projection maps
      call strmap_man_spprj(
     &     maxbuf,nmaps_c,
     &     graph_csub,msdis_c_tgt,ncblk,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(spprjmap_c,maxbuf,'spprjmap_c')
      call strmap_man_spprj(
     &     maxbuf,nmaps_a,
     &     graph_asub,msdis_a_tgt,nablk,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(spprjmap_a,maxbuf,'spprjmap_a')
c dbg
c      print *,'maxbuf(A)=',maxbuf
c dbg

      ! Loop over Irrep of annihilator string.
      igama_loop: do igama = 1, ngam          

        igamc = multd2h(igama,mel%gamt)

        if (ntest.ge.100)
     &         write(lulog,*) 'MS(A), GAMMA(A): ',msa,igama,' len = ',
     &           mel%len_op_gmo(iblk)%gam_ms(igama,idxmsa)

        if (mel%len_op_gmo(iblk)%
     &         gam_ms(igama,idxmsa).le.0) cycle

        first = .true.
        distr_loop: do

          ! Loop over MS/Irrep-distributions of sceleton list
          if (.not.next_msgamdist2(first,
     &            msdis_c,msdis_a,gamdis_c,gamdis_a,
     &            ncblk, nablk,
     &            occ_csub,occ_asub,
     &            msc,msa,igamc,igama,ngam,
     &            ms_fix,ms_fix_ok)) exit
          first = .false.

c dbg
c          print *,'comparing present msdis:'
c          print *,'msdis_c:',msdis_c(1:ncblk)
c          print *,'msdis_a:',msdis_a(1:nablk)
c          print *,'to'
c          print *,'msdis_c_tgt:',msdis_c_tgt(1:ncblk)
c          print *,'msdis_a_tgt:',msdis_a_tgt(1:nablk)
c dbg          
          ! do we have a gamma distribution enumerator?
          ! for the moment we skip unless we get our target
          ! MS distribution
          if (.not.list_cmp(msdis_c,msdis_c_tgt,ncblk).or.
     &        .not.list_cmp(msdis_a,msdis_a_tgt,nablk)) cycle

c dbg
c          print *,'OK!'
c dbg
          if (ms_fix.and..not.ms_fix_ok) cycle

          call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncblk)
          call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nablk)

C          ! check for canon. sequence of MS-distributions
C          ok = .true.
C          do idx = 2, ncblk
C            ok = ok.and.(msdis_c(idx-1).le.msdis_c(idx))
C          end do
C          do idx = 2, nablk
C            ok = ok.and.(msdis_a(idx-1).le.msdis_a(idx))
C          end do
C
C          if (.not.ok) cycle

          call set_len_str(len_str,ncblk,nablk,
     &                     graphs,
     &                     graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                     graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                     hpvxseq,.false.)

          lenc = ielprd(len_str,ncblk)
          lena = ielprd(len_str(ncblk+1),nablk)

          if (lenc.eq.0.or.lena.eq.0) cycle

          ! generate info how to retrieve the elements that
          ! belong to a given spin-adapted element
          call set_spinprj_info(offsets,maps_c,maps_a,
     &            ldim_op_c,ldim_op_a,gsign,
     &            nsplc,mel,iblk,ioff0,igama,ngam,graphs,
     &            msdis_c,gamdis_c,hpvx_csub,occ_csub,graph_csub,ncblk,
     &            msdis_a,gamdis_a,hpvx_asub,occ_asub,graph_asub,nablk,
     &            s2,mel%mst,ncoup,coeff,nomni,omnistrings) 

          call get_spprjmap_blk(spprjmap_c,
     &         ncblk,occ_csub,len_str,graph_csub,idxmsdis_c,gamdis_c,
     &         strmap_info,ngam,ngraph)
          call get_spprjmap_blk(spprjmap_a,
     &         nablk,occ_asub,len_str(ncblk+1),
     &                                graph_asub,idxmsdis_a,gamdis_a,
     &         strmap_info,ngam,ngraph)

          idxa_loop: do idxa = 1, lena
            istr = idxa-1
            ioff = 0
            asign(1:nsplc) = gsign(1:nsplc)
            do icmp = 1, nablk
              nmap = nmaps_a(icmp)
              istr_asub0 = mod(istr,len_str(ncblk+icmp))
              do isplc = 1, nsplc
                idxmap = maps_a(isplc,icmp)
                if (idxmap.eq.0) then
                  istr_asub(icmp,isplc)=istr_asub0
                else 
                  imap = spprjmap_a(ioff+istr_asub0*nmap+idxmap)
                  if (imap.eq.0) asign(isplc)=0
                  istr_asub(icmp,isplc) = abs(imap)-1
                  asign(isplc) = asign(isplc)*sign(1,imap)
                end if
              end do
              istr = istr/len_str(ncblk+icmp)
              ioff = ioff + len_str(ncblk+icmp)*nmap
            end do
            idxc_loop: do idxc = 1, lenc
              istr = idxc-1
              ioff = 0
              csign(1:nsplc) = 1
              do icmp = 1, ncblk
                nmap = nmaps_c(icmp)
                istr_csub0 = mod(istr,len_str(icmp))
                do isplc = 1, nsplc
                  idxmap = maps_c(isplc,icmp)
                  if (idxmap.eq.0) then
                    istr_csub(icmp,isplc)=istr_csub0
                  else
                    imap = spprjmap_c(ioff+istr_csub0*nmap+idxmap)
                    if (imap.eq.0) csign(isplc) = 0
                    istr_csub(icmp,isplc) = abs(imap)-1
                    csign(isplc) = csign(isplc)*sign(1,imap)
                  end if
                end do
                istr = istr/len_str(icmp)
                ioff = ioff + len_str(icmp)*nmap
              end do

              do isplc = 1, nsplc

                idxel(isplc) = 0
                value(isplc) = 0d0
                if (csign(isplc)*asign(isplc).eq.0) cycle

                idxel(isplc) = offsets(isplc)
     &             + idx_str_blk3(istr_csub(1,isplc),istr_asub(1,isplc),
     &                            ldim_op_c(1,isplc),ldim_op_a(1,isplc),
     &                            ncblk,nablk)
                ! relative sign if index appears multiple times
                psign(isplc) = popcnt_par(idxcount(idxel(isplc),idxel,
     &                                             isplc,1))

                ! asign and csign are not needed because the
                ! alpha/beta strings are in the order in which they
                ! are considered in set_spinprj_info.
                ! exception is repeated indices, for which we assume
                ! that the first appearance has the correct order
                value(isplc) = buffer_in(idxel(isplc))*
     &               dble(psign(isplc))
c     &               dble(asign(isplc)*csign(isplc))

              end do

              ! check if element has been treated before:
              ! if any of the "omnipresent" elements is zero
              if (any(idxel(omnistrings(1:nomni)).eq.0)) cycle


              ! spin-adapt: project onto the desired spin components
              call dgemv('t',nsplc,ncoup,1d0,coeff,nsplc,
     &                   value,1,0d0,overl,1)
              call dgemv('n',nsplc,ncoup,1d0,coeff,nsplc,
     &                   overl,1,0d0,value,1)

              ! distribute values
              do isplc = 1, nsplc

                if (csign(isplc)*asign(isplc).eq.0) cycle

                buffer_out(idxel(isplc)) = value(isplc)*
     &               dble(psign(isplc))
c     &               dble(asign(isplc)*csign(isplc))

              end do


            end do idxc_loop
          end do idxa_loop

        end do distr_loop
          
      end do igama_loop
          
      deallocate(hpvx_csub,hpvx_asub,
     &         occ_csub, occ_asub,
     &         graph_csub, graph_asub,
     &         msdis_c,  msdis_a,
     &         msdis_c_tgt,  msdis_a_tgt,
     &         idxmsdis_c,  idxmsdis_a,
     &         gamdis_c, gamdis_a,
     &         len_str, idxel, value,
     &         gsign, csign, asign, offsets, 
     &         maps_c, maps_a,
     &         nmaps_c, nmaps_a,
     &         istr_csub,istr_asub,
     &         ldim_op_c,ldim_op_a,
     &         coeff,overl)

      ifree = mem_flushmark('spin_prj_blk')

      return

      contains

      subroutine get_target_dis(nocc,occ,tgtdis)

      implicit none

      integer, intent(in) ::
     &     nocc, occ(nocc)
      integer, intent(out) ::
     &     tgtdis(nocc)
      integer ::
     &     ii, jj, nodd, noddlarge, icnt
      logical ::
     &     assigned(nocc)

      ! assign dis=0 to even occupations, and
      ! count occurrences of odd numbers (and those larger than 1)
      assigned(1:nocc) = .false.
      nodd = 0
      noddlarge = 0
      do ii = 1, nocc
        if (mod(occ(ii),2).eq.1) then
          nodd = nodd + 1
          if (occ(ii).gt.1) noddlarge = noddlarge + 1
        else
          tgtdis(ii) = 0
          assigned(ii) = .true.
        end if
      end do

      if (all(assigned(1:nocc))) return

      ! odd occupations >1 must be assigned +1 distribution
      ! otherwise, handling of maps needs to be changed...
      if (noddlarge.gt.(nodd+1)/2) call warn('get_target_dis',
     &      'will not be able to assign enough +1 distributions')
      do jj = 1, min(noddlarge,(nodd+1)/2)
        do ii = 1, nocc
          if (assigned(ii)) cycle
          if (occ(ii).gt.1) then
            tgtdis(ii) = +1
            assigned(ii) = .true.
            exit
          end if
        end do
      end do

      ! now fill up with +1 distributions from the left
      ! and +1 distributions from the right
      icnt = min(noddlarge,(nodd+1)/2)
      do while (.not.all(assigned(1:nocc)))
        if (icnt.le.0) then ! +1 dist. on the left
          do ii = 1, nocc
            if (assigned(ii)) cycle
            tgtdis(ii) = +1
            assigned(ii) = .true.
            icnt = icnt+1
            exit
          end do
        else ! -1 dist. on the right
          do ii = nocc, 1, -1
            if (assigned(ii)) cycle
            tgtdis(ii) = -1
            assigned(ii) = .true.
            icnt = icnt-1
            exit
          end do
        end if
      end do

      return
      end subroutine

      end
