*----------------------------------------------------------------------*
      subroutine idx42str(nstr,idxstr,
     &         idxprqs,igam,idss,igtp,
     &         orb_info,str_info,hop,ihpvseq,ierr)
*----------------------------------------------------------------------*
*     a 4-tuple of indices (for a 2-electron integral, in type ordering) 
*     is given: (pq|rs), along with
*          igam: IRREP of orbital
*          igtp: H/P/V of orbital
*          idss: subspace of orbital (within H/P/V)
*     obtain all possible indices in string ordered operator array
*     to which this integral contributes:
*         Ca   Aa:    <pr|qs> if p<r,q<s or with sign changing acc. to
*                             number of permutations neccessary to bring
*                             indices into this order
*                     <qs|pr> (as we do not exploit P/H symmetry)
*         CaCb AaAb:  <pr|qs> dto + all spin distrib.s
*                     <qs|pr> 
*         Cb   Ab:    <pr|qs> sign as CaAa case
*                     <qs|pr>
*
*     sign change indicated by negative index on idxstr
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'multd2h.h'
      include 'def_graph.h'
      include 'def_operator.h'
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
      type(operator), intent(in) ::
     &     hop
      integer, intent(out) ::
     &     nstr, idxstr(*)
      logical, intent(out) ::
     &     ierr

      logical ::
     &     reo12, reo34, eqv12, eqv34
      integer ::
     &     igamt, mst, iblk_ca, iblk_ac, ihpvdx, ihpv, ica,
     &     jdx, idxcnt, idstr_ca, idstr_ac, nel, msstr, lenlast,
     &     nmsd, imsd, igraph, idx1, idx2, idx3, idx4, ndup,
     &     idxms

      integer ::
     &     iocc(ngastp,2), igmd(ngastp,2), !igmdreo(ngastp,2),
     &     msd(ngastp,2), ioff(ngastp,2),
     &     idx(2), len(2), ipos,!(2),
     &     idorb(4), idspn(4), idspc(4), idgam(4)

      type(graph), pointer ::
     &     curgraph

      integer, external ::
     &     iblk_occ, idx_msgmdst, idx4sg
      logical, external ::
     &     allow_sbsp_dis


      if (ntest.ge.10) then
        write(luout,*) '----------------------'
        write(luout,*) ' output from idx42str'
        write(luout,*) '----------------------'
      end if

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
      if (reo12) idx1 = 2
      if (reo12) idx2 = 1
      if (reo34) idx3 = 4
      if (reo34) idx4 = 3

      ! Sets up the occupation for the integral i.e. the combination of
      ! creation and annihilation operators in each space.
      iocc(1:ngastp,1:2) = 0
      iocc(igtp(idx1),1) = 1
      iocc(igtp(idx2),1) = iocc(igtp(idx2),1)+1
      iocc(igtp(idx3),2) = 1
      iocc(igtp(idx4),2) = iocc(igtp(idx4),2)+1

      ! Sets up an offset matrix for the operator matrix.
      ioff(1:ngastp,1:2) = 5
      ioff(igtp(idx1),1) = 1
      ioff(igtp(idx2),1) = min(ioff(igtp(idx2),1),2)
      ioff(igtp(idx3),2) = 3
      ioff(igtp(idx4),2) = min(ioff(igtp(idx4),2),4)

      ! Deal with total symmetry of the integral.
      ! we actually mean: igama (=igamc)
      igamt = igam(1)
      igamt = multd2h(igamt,igam(2))

      ! Symmetries of each side of the integral.
      igmd(1:ngastp,1:2) = 1
      igmd(igtp(idx1),1) = igam(idx1)
      igmd(igtp(idx2),1) = multd2h(igmd(igtp(idx2),1),igam(idx2))
      igmd(igtp(idx3),2) = igam(idx3)
      igmd(igtp(idx4),2) = multd2h(igmd(igtp(idx4),2),igam(idx4))

      ! Which block of the operator does the passed integral represent?
      ! could be improved:
      iblk_ca = iblk_occ(iocc,.false.,hop)
      iblk_ac = iblk_occ(iocc,.true.,hop)
      if (iblk_ca.le.0.or.iblk_ac.lt.0) then
        write(luout,*) iblk_ca, iblk_ac
        call quit(1,'idx42str','something''s buggy!')
      end if

      ierr=.false.
      if(hop%formal_blk(iblk_ca))then
        ierr=.true.
        return
      endif

      if (ntest.ge.50) then
        write(luout,*) 'input <pr|qs>: ',idxprqs(1:4)
        write(luout,*) '        Gamma  ',igam(1:4)
        write(luout,*) '        Subsp  ',idss(1:4)
        write(luout,*) '        Types  ',igtp(1:4)
      end if
      if (ntest.ge.100) then
        write(luout,*) ' igamt = ',igamt
        write(luout,*) ' iocc:'
        call wrt_occ(luout,iocc)
        write(luout,*) ' igmd:'
        call wrt_occ(luout,igmd)
        write(luout,*) ' Blk(CA): ',iblk_ca
        write(luout,*) ' Blk(AC): ',iblk_ac
        write(luout,*) ' Flags: ',reo12,reo34,eqv12,eqv34
        write(luout,*) ' ioff:'
        call wrt_occ(luout,ioff)
      end if

      nstr = 0
      cnt_loop: do idxcnt = 1, 4
        ! If the spatial orbitals are equivalent for either electron 
        ! one or two and if the spins are the same (idxcnt=1 or 4) then
        ! there is no contribution due to Pauli exclusion principle:
        if ((eqv12.or.eqv34).and.(idxcnt.eq.1.or.idxcnt.eq.4))
     &     cycle cnt_loop

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

        ! Loop to prevent the same string being counted twice.
        if (idxcnt.eq.3.and.(eqv12.and.eqv34)) cycle cnt_loop

        ! Increment total number of strings, normal and adjoint.
        nstr = nstr+2

        ! make msd from idspn
        msd(1:ngastp,1:2) = 0
        msd(igtp(idx1),1) = idspn(1)
        msd(igtp(idx2),1) = msd(igtp(idx2),1)+idspn(2)
        msd(igtp(idx3),2) = idspn(3)
        msd(igtp(idx4),2) = msd(igtp(idx4),2)+idspn(4)

        ! correct entries for paired spins
        if (eqv12) idspn(1:2) = (/2,2/)
        if (eqv34) idspn(3:4) = (/2,2/)

        if (ntest.ge.100) then
          write(luout,*) ' idxcnt = ',idxcnt
          write(luout,*) ' spins = ',idspn(1:4)
          call wrt_occ(luout,msd)
        end if
        if (ntest.ge.100) then
          write(luout,*) 'current psqr: ',idorb(1:4)
          write(luout,*) '       Gamma  ',idgam(1:4)
          write(luout,*) '       Subsp  ',idspc(1:4)
          write(luout,*) '       spins  ',idspn(1:4)            
        end if

        ! Identify the index of the string in question and also its 
        ! adjoint.
        idstr_ca = idx_msgmdst(iblk_ca,mst,igamt,
     &         msd,igmd,.false.,hop,orb_info%nsym)
        idstr_ac = idx_msgmdst(iblk_ac,mst,igamt,
     &         msd,igmd,.true.,hop,orb_info%nsym)

        if (ntest.ge.100) then
          write(luout,*) 'idstr_ca = ',idstr_ca
          write(luout,*) 'idstr_ac = ',idstr_ac
          if (idstr_ca.lt.0.or.idstr_ac.lt.0)
     &         stop 'halt!'
        end if

        ! Locate offsets which were calculated earlier and store them
        ! in idxstr.
        idxstr(nstr-1) =
     &         hop%off_op_gmox(iblk_ca)%d_gam_ms(idstr_ca,igamt,idxms)
        idxstr(nstr)   =
     &         hop%off_op_gmox(iblk_ac)%d_gam_ms(idstr_ac,igamt,idxms)

        if (ntest.ge.100) then
          write(luout,*) 'idxstr(offsets) = ',idxstr(nstr-1:nstr)
        end if

        lenlast = 1
        ihpv_loop: do ihpvdx = 1, ngastp
          ihpv = ihpvseq(ihpvdx)
          idx(1:2) = 0
          len(1:2) = 1
          ica_loop: do ica = 1,2
            nel = iocc(ihpv,ica)
            if (iocc(ihpv,ica).eq.0) cycle ica_loop
            ipos = ioff(ihpv,ica)

            ! point to graph needed for current string
            igraph = hop%idx_graph(ihpv,ica,iblk_ca)
            curgraph => str_info%g(igraph)

            ! check for restrictions
            if (.not.allow_sbsp_dis(idspc(ipos),nel,
     &             orb_info%ngas_hpv(ihpv),
     &             str_info%igas_restr(1,1,1,igraph))) then
              nstr = nstr-2
              exit cnt_loop
            end if

            ! get string index
            idx(ica) = idx4sg(nel,idspc(ipos),idorb(ipos),
     &             idspn(ipos),idgam(ipos),
     &             curgraph%y4sg,curgraph%yinf,
     &             curgraph%yssg,curgraph%wssg,
     &             orb_info%mostnd(1,1,orb_info%idx_gas(ihpv)),
     &             str_info%ispc_occ(igraph),orb_info%nsym,
     &             orb_info%ngas_hpv(ihpv))

            ! get total length within block
            msstr = (msd(ihpv,ica)+iocc(ihpv,ica))/2+1
            len(ica) =
     &             curgraph%lenstr_gm(igmd(ihpv,ica),msstr)
      
          end do ica_loop
          ! string number is actual index - 1
          idxstr(nstr-1) = idxstr(nstr-1) +
     &           ((idx(2))*len(1) + idx(1))*lenlast
          idxstr(nstr)   = idxstr(nstr) +
     &           ((idx(1))*len(2) + idx(2))*lenlast
          lenlast = lenlast*len(1)*len(2)
        end do ihpv_loop

        idxstr(nstr-1:nstr) = idxstr(nstr-1:nstr)+1

        ! Change sign if necessary.
        if (idxcnt.eq.3.and.(eqv12.or.eqv34))
     &         idxstr(nstr-1:nstr) = -idxstr(nstr-1:nstr)

        if (ntest.ge.100) then
          write(luout,*) 'idxstr(final) = ',idxstr(nstr-1:nstr)
        end if

        ! prelim: check, whether actually the same address resulted
        ndup = 0
        check1: do jdx = 1, nstr-2
          if (idxstr(nstr-1).eq.idxstr(jdx)) ndup = ndup+1
          if (idxstr(nstr).eq.idxstr(jdx)) ndup = ndup+1
          if (ndup.eq.2) exit check1
        end do check1
        nstr = nstr-ndup
        if (idxstr(nstr-1).eq.idxstr(nstr)) nstr = nstr-1          
                    
      end do cnt_loop

      if ((reo12.and..not.reo34).or.(reo34.and..not.reo12))
     &     idxstr(1:nstr) = -idxstr(1:nstr)

      if (ntest.ge.50) then
        write(luout,*) 'on output:'
        write(luout,*) ' nstr = ',nstr
        do jdx = 1, nstr
          write(luout,*) '  idx: ',idxstr(jdx)
        end do
      end if

      return
      end

