*----------------------------------------------------------------------*
      subroutine idx42str_nox(nstr,idxstr,
     &         idxprqs,igam,idss,igtp,
     &         orb_info,str_info,hlist,ihpvseq,error)
*----------------------------------------------------------------------*
*     a 4-tuple of indices (for a 2-electron integral, in type ordering) 
*     is given: (pq|rs), along with
*          igam: IRREP of orbital
*          igtp: H/P/V of orbital
*          idss: subspace of orbital (within H/P/V)
*     obtain all possible indices in string ordered operator array
*     to which this integral contributes:
*
*     Special version for non-antisymmetrized integrals
*
*     sign change indicated by negative index on idxstr
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'multd2h.h'
      include 'def_graph.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     idxprqs(4), igam(4), idss(4), igtp(4),
     &     ihpvseq(ngastp)
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(me_list), intent(in) ::
     &     hlist
      integer, intent(out) ::
     &     nstr, idxstr(*)
      logical, intent(out) ::
     &     error

      logical ::
     &     reo12, reo34, eqv12, eqv34, fail
      integer ::
     &     igamt, mst, iblk, iblk_ca, iblk_ac, ihpvdx, ihpv, ica,
     &     jdx, idxcnt, idstr_ca, idstr_ac, nel, msstr, lenlast,
     &     nmsd, imsd, igraph, idx1, idx2, idx3, idx4, ndup,
     &     idxms, ij, istr_ca, istr_ac

      integer ::
     &     iocc(ngastp,2,2), igmd(ngastp,2,2),
     &     msd(ngastp,2,2), ioff(ngastp,2,2),
     &     idx(2), len(2), ipos,!(2),
     &     idorb(4), idspn(4), idspc(4), idgam(4),
     &     occ_c(2), occ_a(2), idxms_c(2), idxms_a(2),
     &     gam_c(2), gam_a(2)

      type(graph), pointer ::
     &     curgraph
      type(operator), pointer ::
     &     hop

      integer, external ::
     &     iblk_occ, idx_msgmdst2, idx4sg
      logical, external ::
     &     allow_sbsp_dis

      if (ntest.ge.10) then
        write(lulog,*) '--------------------------'
        write(lulog,*) ' output from idx42str_nox'
        write(lulog,*) '--------------------------'
      end if
 
      nstr = 0

      hop => hlist%op

      ! Checks whether integral indices are the same, or if they need 
      ! to be reordered (want p.le.q and r.le.s)
      reo12 = idxprqs(1).gt.idxprqs(2)
      eqv12 = idxprqs(1).eq.idxprqs(2)
      reo34 = idxprqs(3).gt.idxprqs(4)
      eqv34 = idxprqs(3).eq.idxprqs(4)

      idx1 = 1
      idx2 = 2
      idx3 = 3
      idx4 = 4
      ! Reorder if necessary.
      !if (reo12) idx1 = 2
      !if (reo12) idx2 = 1
      !if (reo34) idx3 = 4
      !if (reo34) idx4 = 3

      ! Sets up the occupation for the integral i.e. the combination of
      ! creation and annihilation operators in each space.
      ! pr qs  ->  Type(p),Type(q) ; Type(r),Type(s)
      iocc(1:ngastp,1:2,1:2) = 0
      iocc(igtp(idx1),1,1) = 1
      iocc(igtp(idx3),2,1) = 1
      iocc(igtp(idx2),1,2) = 1
      iocc(igtp(idx4),2,2) = 1

      occ_c(1) = 1
      occ_c(2) = 1
      occ_a(1) = 1
      occ_a(2) = 1

      ! Sets up an offset matrix for the operator matrix.
      ioff(1:ngastp,1:2,1:2) = 5
      ioff(igtp(idx1),1,1) = 1
      ioff(igtp(idx3),2,1) = 3
      ioff(igtp(idx2),1,2) = 2
      ioff(igtp(idx4),2,2) = 4

      ! Deal with total symmetry of the integral.
      ! we actually mean: igama (=igamc)
      igamt = igam(1)
      igamt = multd2h(igamt,igam(2))

      ! Symmetries of each side of the integral.
      igmd(1:ngastp,1:2,1:2) = 1
      igmd(igtp(idx1),1,1) = igam(idx1)
      igmd(igtp(idx3),2,1) = igam(idx3)
      igmd(igtp(idx2),1,2) = igam(idx2)
      igmd(igtp(idx4),2,2) = igam(idx4)
 
      gam_c(1) = igam(idx1)
      gam_c(2) = igam(idx2)
      gam_a(1) = igam(idx3)
      gam_a(2) = igam(idx4)
 
      ! Which block of the operator does the passed integral represent?
      ! could be improved:
      iblk_ca = iblk_occ(iocc,.false.,hop,1)
      iblk_ac = iblk_occ(iocc,.true.,hop,1)

      ! no such block at all?
      if (iblk_ca.le.0.and.iblk_ac.le.0) return

      error=.false.

      if (ntest.ge.50) then
        write(lulog,*) 'input <pr|qs>: ',idxprqs(1:4)
        write(lulog,*) '        Gamma  ',igam(1:4)
        write(lulog,*) '        Subsp  ',idss(1:4)
        write(lulog,*) '        Types  ',igtp(1:4)
      end if
      if (ntest.ge.100) then
        write(lulog,*) ' igamt = ',igamt
        write(lulog,*) ' iocc:'
        call wrt_occ_n(lulog,iocc,2)
        write(lulog,*) ' igmd:'
        call wrt_occ_n(lulog,igmd,2)
        write(lulog,*) ' Blk(CA): ',iblk_ca
        write(lulog,*) ' Blk(AC): ',iblk_ac
        write(lulog,*) ' Flags: ',reo12,reo34,eqv12,eqv34
        write(lulog,*) ' ioff:'
        call wrt_occ_n(lulog,ioff,2)
      end if

      nstr = 0
      cnt_loop: do idxcnt = 1, 4
        ! ignore Pauli exclusion principle:
c        if ((eqv12.or.eqv34).and.(idxcnt.eq.1.or.idxcnt.eq.4))
c     &     cycle cnt_loop

        ! Associate spins with each orbital depending on the pass. 
        ! Alpha=1, beta=-1 => aaaa,abab, baba, bbbb.
        if (idxcnt.lt.4) idspn(1:4) = 1
        if (idxcnt.eq.2) idspn(idx2) = -1
        if (idxcnt.eq.2) idspn(idx4) = -1
        if (idxcnt.eq.3) idspn(idx1) = -1
        if (idxcnt.eq.3) idspn(idx3) = -1
        if (idxcnt.eq.4) idspn(1:4) = -1

        ! Get the total spin of the integral and associate an index.
        ! we actually mean: msa (=msc)
        if (idxcnt.eq.1) mst = 2
        if (idxcnt.eq.2) mst = 0
        if (idxcnt.eq.3) mst = 0
        if (idxcnt.eq.4) mst =-2
        idxms = (2-mst)/2+1

        ! Get orbital indices, symmetries, subspaces etc.
        ! set idspc etc....
        idorb(1) = idxprqs(idx1) 
        idorb(2) = idxprqs(idx2) 
        idorb(3) = idxprqs(idx3) 
        idorb(4) = idxprqs(idx4) 
        idgam(1) = igam(idx1)
        idgam(2) = igam(idx2)
        idgam(3) = igam(idx3)
        idgam(4) = igam(idx4)
        idspc(1) = idss(idx1)
        idspc(2) = idss(idx2)
        idspc(3) = idss(idx3)
        idspc(4) = idss(idx4)

        ! not needed for exchangeless integrals
        ! Loop to prevent the same string being counted twice.
        !if (idxcnt.eq.3.and.(eqv12.and.eqv34)) cycle cnt_loop

        ! Increment total number of strings, normal and adjoint.
        istr_ca = -1; istr_ac = -1
        if (iblk_ca.gt.0) then
          nstr = nstr+1
          istr_ca = nstr
        end if
        if (iblk_ac.gt.0) then
          nstr = nstr+1
          istr_ac = nstr
        end if

        ! make msd from idspn
        msd(1:ngastp,1:2,1:2) = 0
        msd(igtp(idx1),1,1) = idspn(1)
        msd(igtp(idx2),1,2) = idspn(2)
        msd(igtp(idx3),2,1) = idspn(3)
        msd(igtp(idx4),2,2) = idspn(4)

        idxms_c(1) = (2-idspn(1))/2+1
        idxms_c(2) = (2-idspn(2))/2+1
        idxms_a(1) = (2-idspn(3))/2+1
        idxms_a(2) = (2-idspn(4))/2+1

        if (ntest.ge.100) then
          write(lulog,*) ' idxcnt = ',idxcnt
          write(lulog,*) ' spins = ',idspn(1:4)
          call wrt_occ_n(lulog,msd,2)
        end if
        if (ntest.ge.100) then
          write(lulog,*) 'current psqr: ',idorb(1:4)
          write(lulog,*) '       Gamma  ',idgam(1:4)
          write(lulog,*) '       Subsp  ',idspc(1:4)
          write(lulog,*) '       spins  ',idspn(1:4)            
        end if

        ! Identify the index of the string in question and also its 
        ! adjoint.
        if (iblk_ca.gt.0) then
c          idstr_ca = idx_msgmdst(iblk_ca,mst,igamt,
c     &         msd,igmd,.false.,hlist,orb_info%nsym)
          fail = .true.
          idstr_ca = idx_msgmdst2(fail,iblk_ca,idxms,igamt,
     &         occ_c,idxms_c,gam_c,2,
     &         occ_a,idxms_a,gam_a,2,
     &         .false.,-1,-1,hlist,orb_info%nsym)
        else
          idstr_ca = -1
        end if
        if (iblk_ac.gt.0) then
c          idstr_ac = idx_msgmdst(iblk_ac,mst,igamt,
c     &         msd,igmd,.true.,hlist,orb_info%nsym)
          fail = .true.
          ! last argument line: No reordering upon transposition
          idstr_ac = idx_msgmdst2(fail,iblk_ac,idxms,igamt,
     &         occ_c,idxms_c,gam_c,2,
     &         occ_a,idxms_a,gam_a,2,
     &         .true.,(/1,2/),(/1,2/),hlist,orb_info%nsym)
        else
          idstr_ac = -1
        end if
        if (ntest.ge.100) then
          write(lulog,*) 'idstr_ca = ',idstr_ca
          write(lulog,*) 'idstr_ac = ',idstr_ac
        end if

        ! Locate offsets which were calculated earlier and store them
        ! in idxstr.
        if (istr_ca.gt.0)
     &      idxstr(istr_ca) =
     &         hlist%off_op_gmox(iblk_ca)%d_gam_ms(idstr_ca,igamt,idxms)
        if (istr_ac.gt.0)
     &      idxstr(istr_ac) =
     &         hlist%off_op_gmox(iblk_ac)%d_gam_ms(idstr_ac,igamt,idxms)

        if (ntest.ge.100) then
          write(lulog,*) 'idxstr(offsets) = '
          if (istr_ca.gt.0) write(lulog,*) idxstr(istr_ca)
          if (istr_ac.gt.0) write(lulog,*) idxstr(istr_ac)
        end if

        iblk = iblk_ca
        if (iblk.le.0) iblk = iblk_ac
        lenlast = 1
        ihpv_loop: do ihpvdx = 1, ngastp
         ihpv = ihpvseq(ihpvdx)
         idx(1:2) = 0
         len(1:2) = 1
         ij_loop: do ij = 1, 2
          ica_loop: do ica = 1,2
            nel = iocc(ihpv,ica,ij)
            if (iocc(ihpv,ica,ij).eq.0) cycle ica_loop
            ipos = ioff(ihpv,ica,ij)

            ! point to graph needed for current string
            if (iblk_ca.gt.0) then
              igraph = hlist%idx_graph(ihpv,ica,(iblk-1)*2+ij)
            else
              ! in this case we have to fetch it from the adj. block
              igraph = hlist%idx_graph(ihpv,3-ica,(iblk-1)*2+ij)
            end if
            curgraph => str_info%g(igraph)

            ! check for restrictions
            if (.not.allow_sbsp_dis(idspc(ipos),nel,
     &             orb_info%ngas_hpv(ihpv),
     &             str_info%igas_restr(1,1,1,1,igraph))) then
              !        ADAPT FOR OPEN-SHELL ^^^
              if (istr_ca.gt.0.and.istr_ac.gt.0) then
                nstr = nstr-2
              else
                nstr = nstr-1
              end if
              exit cnt_loop
            end if
c dbg
c            write(lulog,*) 'ica,ihpv,ij: ',ica,ihpv,ij
c            write(lulog,*) 'nel, ipos, igraph ',nel, ipos, igraph
c dbg

            ! get string index
            idx(ica) = idx4sg(nel,idspc(ipos),idorb(ipos),
     &             idspn(ipos),idgam(ipos),
     &             curgraph%y4sg,curgraph%yinf,
     &             curgraph%yssg,curgraph%wssg,
     &             curgraph%ioffstr_dgm,curgraph%ndis,
     &             orb_info%mostnd(1,1,orb_info%idx_gas(ihpv)),
     &             str_info%ispc_occ(igraph),orb_info%nsym,
     &             orb_info%ngas_hpv(ihpv))

            ! get total length within block
            msstr = (msd(ihpv,ica,ij)+iocc(ihpv,ica,ij))/2+1
            len(ica) =
     &             curgraph%lenstr_gm(igmd(ihpv,ica,ij),msstr)
      
          end do ica_loop
          ! string number is actual index - 1
          if (istr_ca.gt.0)
     &      idxstr(istr_ca) = idxstr(istr_ca) +
     &           ((idx(2))*len(1) + idx(1))*lenlast
          if (istr_ac.gt.0)
     &      idxstr(istr_ac) = idxstr(istr_ac) +
     &           ((idx(1))*len(2) + idx(2))*lenlast
          lenlast = lenlast*len(1)*len(2)
         end do ij_loop
        end do ihpv_loop

        if (istr_ca.gt.0) idxstr(istr_ca) = idxstr(istr_ca)+1
        if (istr_ac.gt.0) idxstr(istr_ac) = idxstr(istr_ac)+1

        ! ???
        !! Change sign if necessary.
        !if (idxcnt.eq.3.and.(eqv12.or.eqv34))
     &  !       idxstr(nstr-1:nstr) = -idxstr(nstr-1:nstr)

        if (ntest.ge.100) then
          write(lulog,*) 'idxstr(final) = '
          if (istr_ca.gt.0) write(lulog,*) idxstr(istr_ca)
          if (istr_ac.gt.0) write(lulog,*) idxstr(istr_ac)
        end if

        ! prelim: check, whether actually the same address resulted
        if (istr_ca.gt.0.and.istr_ac.gt.0) then
          ndup = 0
          check1: do jdx = 1, nstr-2
            if (idxstr(nstr-1).eq.idxstr(jdx)) ndup = ndup+1
            if (idxstr(nstr).eq.idxstr(jdx)) ndup = ndup+1
            if (ndup.eq.2) exit check1
          end do check1
          nstr = nstr-ndup
          if (idxstr(nstr-1).eq.idxstr(nstr)) nstr = nstr-1          
        end if

      end do cnt_loop

!      if ((reo12.and..not.reo34).or.(reo34.and..not.reo12))
!     &     idxstr(1:nstr) = -idxstr(1:nstr)

      if (ntest.ge.50) then
        write(lulog,*) 'on output:'
        write(lulog,*) ' nstr = ',nstr
        do jdx = 1, nstr
          write(lulog,*) '  idx: ',idxstr(jdx)
        end do
      end if

      return
      end

