*----------------------------------------------------------------------*
      subroutine read_mo_file_dalton_spc(buffer,ffinp,
     &       naux_min, naux_max, nh_min, nh_max,
     &       iy_int,typetab,ntypes,orb_info)
*----------------------------------------------------------------------*
*     read integrals from special DALTON environment and store
*     in a direct-access array addressed by iy_int and typetab
*
*     andreas, march 2008
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'par_dalton.h'

      integer, parameter::
     &     ntest= 00

      real(8), intent(out) ::
     &     buffer(*)
      integer, intent(in) ::
     &     ntypes, iy_int(*), typetab(24),
     &     naux_min, naux_max, nh_min, nh_max
      type(orbinf), intent(in), target ::
     &     orb_info
      type(filinf), intent(in) ::
     &     ffinp

      real(8) ::
     &     val, fac, abfac      ,cx,wx,c0,w0
      integer ::
     &     lu2in, iq, ip, ir, is, naux, ntoob, caborb, ngam, idx,
     &     idx_typ, idx_int, idxperm, iperm, ntotal, p,q, ld, idxpq,
     &     nocc, nh
      integer ::
     &     idxpqrs(4), idx_ord(4), nints(8), rank(4),
     &     rank1(4), rank2(4), perm(4), idx_ord1(4)
      integer, pointer ::
     &     reord(:), reost(:), igamorb(:)

      integer, external ::
     &     idx_int_graph, rank_ivec, ipermfac


      idxpq(p,q,ld) = (min(p,q)-1)*ld+max(p,q)

      abfac = 1d0
      if (ntypes.eq.6) abfac = -1d0

      ngam  = orb_info%nsym
      ntoob = orb_info%ntoob
      caborb = orb_info%caborb
      igamorb => orb_info%igamorb

      ntotal = ntoob+caborb
      nocc   = orb_info%norb_hpv(IHOLE,1)
      if (orb_info%nspin.eq.2)
     &     nocc = max(nocc,orb_info%norb_hpv(IHOLE,2))

      allocate(reord(ntoob+caborb))

      reost=>orb_info%ireost
      do idx = 1,ntoob
        reord(idx)=reost(idx)
      enddo
      do idx=ntoob+1,ntoob+caborb
        reord(idx)=idx
      enddo  

      lu2in = ffinp%unit

      ! rewind and over-read first two records
      rewind(lu2in)
      read(lu2in,*)
      read(lu2in,*)

      ! read from file to buffer
      do
          read(unit=lu2in,end=999,err=998,fmt=*) ip, iq, ir, is, val
          naux = 0
          if (ip.gt.ntoob) naux = naux+1
          if (iq.gt.ntoob) naux = naux+1
          if (ir.gt.ntoob) naux = naux+1
          if (is.gt.ntoob) naux = naux+1
          if (naux.lt.naux_min.or.naux.gt.naux_max) cycle

          idxpqrs(1) = reord(is)
          idxpqrs(2) = reord(ir)
          idxpqrs(3) = reord(iq)
          idxpqrs(4) = reord(ip)

          nh = 0
          if (idxpqrs(1).le.nocc) nh = nh+1
          if (idxpqrs(2).le.nocc) nh = nh+1
          if (idxpqrs(3).le.nocc) nh = nh+1
          if (idxpqrs(4).le.nocc) nh = nh+1
          if (nh.lt.nh_min.or.nh.gt.nh_max) cycle

          ! remove anti-hermiticity, if applicable
          fac = 1d0
c          if ((ip-1)*ntotal+ir.gt.(iq-1)*ntotal+is) fac = abfac
c          if (idxpq(ip,ir,ntotal).gt.idxpq(iq,is,ntotal)) fac = abfac
          if (idxpq(idxpqrs(2),idxpqrs(4),ntotal).gt.
     &        idxpq(idxpqrs(1),idxpqrs(3),ntotal)) fac = abfac

c dbg
          print *,'val = ',val, ' fac = ',fac
          print '(x,4i4,a,4i4)',ip,iq,ir,is,' --> ',idxpqrs(1:4)
c dbg

          idxperm = rank_ivec(rank,idxpqrs,4)+1

          idx_ord(rank(1)+1) = idxpqrs(1)
          idx_ord(rank(2)+1) = idxpqrs(2)
          idx_ord(rank(3)+1) = idxpqrs(3)
          idx_ord(rank(4)+1) = idxpqrs(4)

c dbg
          print '(x,4i4,a,4i4)',rank(1:4),' --> ',idx_ord(1:4)          
c dbg

          idx_int = idx_int_graph(idx_ord,4,iy_int,igamorb,ngam)
          idx_typ = abs(typetab(idxperm))
c dbg
          print *,'idx_int, idxtyp: ',idx_int, idx_typ
c          print *,'reord:'
c          print *,reord(1:10)
c          print *,'adr = ',(idx_int-1),ntypes,idx_typ
c          print *,'adr = ',(idx_int-1)*ntypes+idx_typ
c          print *,'------------------------------------------------'
c dbg
          buffer((idx_int-1)*ntypes+idx_typ) = fac*val

          if (ipermfac(idx_ord,4,.true.).eq.1) cycle

          perm = (/1,2,3,4/)

          ! check all permutations
          do iperm = 2, 24
            call next_perm(perm,4)
            call perm_mult(rank1,rank,perm,4)
            idx_ord1(rank1(1)+1) = idxpqrs(1)
            idx_ord1(rank1(2)+1) = idxpqrs(2)
            idx_ord1(rank1(3)+1) = idxpqrs(3)
            idx_ord1(rank1(4)+1) = idxpqrs(4)

            ! index quadruple invariant under this permutation?
            if (idx_ord1(1).ne.idx_ord(1) .or.
     &          idx_ord1(2).ne.idx_ord(2) .or.
     &          idx_ord1(3).ne.idx_ord(3) .or.
     &          idx_ord1(4).ne.idx_ord(4) ) cycle

            idxperm = rank_ivec(rank2,rank1,4)+1
            idx_typ = abs(typetab(idxperm))
c dbg
            print '(i4,a,4i4,a,i4,a,i4)',
     &           iperm,'perm: ',perm,' idx=',idxperm,' typ=',idx_typ
c          print *,'adr = ',(idx_int-1)*ntypes+idx_typ
c dbg

            buffer((idx_int-1)*ntypes+idx_typ) = fac*val
          end do

      end do

 999  write(lulog,*) 'read in finished'

      return
 998  call quit(0,'read_mo_file_dalton_spc','error reading integrals')
      end

