*----------------------------------------------------------------------*
      subroutine import_list_from_buffer(mel,buffer_in,
     &     nauxmin,nauxmax,
     &     fac_0,fac_ne0,scaling,
     &     iy_int,typetab,ntypes,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     set integral list on mel using the graph-indexed raw integral
*     buffer on buffer_in
*
*     andreas, march 2008
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_operators.h'
      include 'ifc_memman.h'


      integer, intent(in) ::
     &     nauxmin, nauxmax
      integer, intent(in) ::
     &     scaling
      type(me_list), intent(in) ::
     &     mel
      integer, intent(in) ::
     &     ntypes, iy_int(*), typetab(24)
      real(8), intent(in) ::
     &     fac_0, fac_ne0
      real(8), intent(in) ::
     &     buffer_in(*)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      
      logical ::
     &     first, close_again, blk_buf, scalar, dagger, xchange
      integer ::
     &     idoff, idxoff, idxoff_blk, iblk, lenblk, ifree, mmax,
     &     msamax, mscmax, idxms, ms, igam, idx_dis, ndis, did,
     &     idum, nel, mst, ngas, ngam, naux, lendis,
     &     njoined, idx_occ, ntotal
      integer ::
     &     msd(ngastp,2,mel%op%njoined), igamd(ngastp,2,mel%op%njoined),
     &     lexlscr(8,3), idorb(8), idorb2(8), idspn(8), idspc(8)
      integer, pointer ::
     &     occ(:,:,:), mostnd(:,:,:), idx_graph(:,:,:), igamorb(:),
     &     ngas_hpv(:), idx_gas(:), igas_restr(:,:,:,:,:)
      real(8) ::
     &     abfac, gfac
      real(8), pointer ::
     &     buffer_reo(:), curblk(:)

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      logical, external ::
     &     next_msgamdist
      real(8), external ::
     &     ddot

      ifree = mem_setmark('import_list')

      op => mel%op
      ffop => mel%fhand
      mst = mel%mst
      dagger = mel%op%dagger

      ngas = orb_info%ngas
      ngam = orb_info%nsym
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      ngas_hpv => orb_info%ngas_hpv
      idx_gas => orb_info%idx_gas

      ntotal = orb_info%ntoob+orb_info%caborb

      abfac = 1d0
      if (ntypes.eq.6) abfac = -1d0

      igas_restr => str_info%igas_restr

      njoined = op%njoined

      mmax = 0
      do iblk = 1, op%n_occ_cls
        if (op%formal_blk(iblk)) cycle
        idx_occ = (iblk-1)*njoined+1
        occ => op%ihpvca_occ(1:ngastp,1:2,idx_occ:idx_occ+njoined-1)
        naux = sum(occ(IEXTR,1:2,1:njoined))
        if (naux.lt.nauxmin .or. naux.gt.nauxmax) cycle

        mscmax = op%ica_occ(1,iblk)
        msamax = op%ica_occ(2,iblk)
        idxms = 0
        do ms = msamax, -msamax, -2
          if (abs(mst+ms).gt.mscmax) cycle
          idxms = idxms+1
          do igam = 1, orb_info%nsym
            mmax = max(mmax,mel%len_op_gmo(iblk)%gam_ms(igam,idxms))
          end do
        end do
      end do
      ifree = mem_alloc_real(buffer_reo,mmax,'buffer_reo')

      close_again = .false.
      if (ffop%unit.le.0) then
        close_again = .true.
        call file_open(ffop)
      end if

      do iblk = 1, op%n_occ_cls
        if (op%formal_blk(iblk)) cycle
        idx_occ = (iblk-1)*njoined+1

        occ => op%ihpvca_occ(1:ngastp,1:2,idx_occ:idx_occ+njoined-1)
        idx_graph =>
     &         mel%idx_graph(1:ngastp,1:2,idx_occ:idx_occ+njoined-1)
        naux = sum(occ(IEXTR,1:2,1:njoined))
c dbg
c        call wrt_occ_n(6,occ,njoined)
c        print *,'naux, (min, max): ',naux,'(',nauxmin,nauxmax,')'
c dbg
        if (naux.lt.nauxmin .or. naux.gt.nauxmax) cycle

        blk_buf = ffop%buffered
        if (blk_buf) blk_buf = blk_buf.and.ffop%incore(iblk).gt.0

        if (max(op%ica_occ(1,iblk),op%ica_occ(2,iblk)).ne.2) cycle

        scalar = max(op%ica_occ(1,iblk),op%ica_occ(2,iblk)).eq.0

        if (scalar) call quit(1,'import_list_from_buffer','scalar?')

        mscmax = op%ica_occ(1,iblk)
        msamax = op%ica_occ(2,iblk)
        nel = msamax+op%ica_occ(1,iblk)
        if (nel.ne.4) call quit(1,'import_list_from_buffer','nel!=4?')
        idxms = 0
        idxoff = 0
        idxoff_blk = 0
        do ms = msamax, -msamax, -2
          if (abs(ms+mst).gt.mscmax) cycle
          idxms = idxms+1

          xchange = ms.ne.0.or.scaling.gt.0

          gfac = fac_0
          if (ms.ne.0) gfac = fac_ne0
        
          do igam = 1, orb_info%nsym

            ! block offest and length:
            idxoff = mel%off_op_gmo(iblk)%gam_ms(igam,idxms)
            lenblk = mel%len_op_gmo(iblk)%gam_ms(igam,idxms)
            if (lenblk.eq.0) cycle

            ! point to current block
            if (.not.blk_buf) then
              curblk => buffer_reo
            else 
c              ioff = op%off_op_gmo(iblk)%gam_ms(igam,idxms)
              ! currently: idxoff should be valid here, as well
              curblk => ffop%buffer(idxoff+1:idxoff+lenblk)
            end if

c              write(luout,'(2x,a,i3,a,i2,a,i12)')
c     &           'Ms(A) = ',ms,'/2  IRREP(A) = ',igam,'  len = ',lenblk

            ! loop over distributions
            first = .true.
            ndis = mel%off_op_gmox(iblk)%ndis(igam,idxms)
            if (ndis.eq.0)
     &           call quit(1,'import_list_from_buffer','ndis=0?')
            distr_loop: do idx_dis = 1, ndis
                  
                lendis =
     &               mel%len_op_gmox(iblk)%d_gam_ms(idx_dis,igam,idxms)
                if (lendis.eq.0) cycle
                idxoff_blk =
     &               mel%off_op_gmox(iblk)%d_gam_ms(idx_dis,igam,idxms)
     &               - idxoff
                did = mel%off_op_gmox(iblk)%did(idx_dis,igam,idxms)
                call did2msgm(msd,igamd,did,occ,ngam,njoined)

c dbg
c                print *,'idxoff_blk, idxoff: ',idxoff_blk, idxoff
c dbg
                call inner_loop(curblk(idxoff_blk+1))
c dbg
c                print *,'current block:'
c                print '(5f16.8)',curblk(1:lenblk)
c dbg

            end do distr_loop

            ! save current block
            if (.not.blk_buf) then
              idoff = ffop%length_of_record*(ffop%current_record-1)
c dbg
c              print *,'putting to disk:',idoff,idxoff+1,idxoff+lenblk
c dbg
              call put_vec(ffop,curblk,
     &             idoff+idxoff+1,idoff+idxoff+lenblk)
            end if

          end do ! gam
        end do ! ms
      end do

      if (close_again) call file_close_keep(ffop)

      ifree = mem_flushmark()

      return

      contains

      subroutine inner_loop(curdisblk)

      implicit none

      include 'hpvxseq.h'

      real(8), intent(inout) ::
     &     curdisblk(*)

      real(8), parameter ::
     &     sp_fac = 1d0/3d0
      integer, parameter ::
     &     idx1(-1:2) = (/1,0,0,0/),
     &     idx2(-1:2) = (/1,0,0,1/)

      integer ::
     &     idxstr, idxperm, idx_typ, idx_int, ihelp,
     &     idx12, idx34, idx43, p,q, ld, idxpq
      integer ::
     &     idx_ord(4), rank(4)
      real(8) ::
     &     fac, facx

      integer, external ::
     &     idx_int_graph, rank_ivec
      logical, external ::
     &     next_tupel_ca

      idxpq(p,q,ld) = (min(p,q)-1)*ld+max(p,q)

      first = .true.
      idxstr = 0
      do while(next_tupel_ca(idorb,idspn,idspc,
     &     nel,njoined,occ,
     &     idx_graph,
     &     msd,igamd,first,
     &     igas_restr,
     &     mostnd,igamorb,
     &     ngam,ngas,
     &     ngas_hpv,idx_gas,
     &     hpvxseq,lexlscr))
        first = .false.
        idxstr = idxstr+1

        idx12 = idxpq(idorb(1),idorb(2), ntotal)
        idx34 = idxpq(idorb(3),idorb(4), ntotal)

        if (idspn(1).ge.idspn(2)) then
          idorb2(1) = idorb(1)
          idorb2(3) = idorb(2) 
          fac = gfac
        else
          idorb2(1) = idorb(2)
          idorb2(3) = idorb(1) 
          fac = -gfac
        end if
        if (idspn(3).ge.idspn(4)) then
          idorb2(2) = idorb(3) 
          idorb2(4) = idorb(4)
        else
          idorb2(2) = idorb(4) 
          idorb2(4) = idorb(3)
          fac = -fac
        end if

        if (idx12.gt.idx34    ) fac  = abfac * fac

        if (ms.eq.0.and.scaling.eq.2) then
          fac = (1d0+sp_fac*sp_fac)*fac
        end if
        
        idxperm = rank_ivec(rank,idorb2,4)+1
c dbg
c
c        print *,'idorb: ',idorb(1:4)
c        print *,'idspn: ',idspn(1:4)
c        print *,'idorb2: ',idorb2(1:4),'  fac = ',fac
c        print *,'rank :',rank(1:4),'  idxperm: ',idxperm
c dbg

        idx_ord(rank(1)+1) = idorb2(1)
        idx_ord(rank(2)+1) = idorb2(2)
        idx_ord(rank(3)+1) = idorb2(3)
        idx_ord(rank(4)+1) = idorb2(4)

        idx_int = idx_int_graph(idx_ord,4,iy_int,igamorb,ngam)
        idx_typ  = abs(typetab(idxperm))

c dbg
c        print *,'Xidorb2: ',idorb2(1:4)
c        print *,'rank :',rank(1:4),'  idxperm: ',idxperm
c dbg

c dbg
c        if (idxstr.le.10) then
c          if (idxstr.eq.1) print *,' -------e.g. -------------------'
c          print '(x,i4,a,4i4,a,2i4,2f16.10)',
c     &       idxstr,': (',idx_ord(1:4),') <- ',
c     &       idx_int,idx_typ, fac,
c     &       buffer_in((idx_int-1)*ntypes+idx_typ)
c          print *,'adr = ',(idx_int-1)*ntypes+idx_typ
c          if (idxstr.eq.10) print *,' ----------------- etc ---------'
c        end if
c dbg
        curdisblk(idxstr) = fac*buffer_in((idx_int-1)*ntypes+idx_typ)

        if (xchange) then

        if (idspn(1).ge.idspn(2)) then
          idorb2(1) = idorb(1)
          idorb2(3) = idorb(2) 
          fac = gfac
        else
          idorb2(1) = idorb(2)
          idorb2(3) = idorb(1) 
          fac = -gfac
        end if
c        if (idspn(4).ge.idspn(3)) then
        if (idspn(3).ge.idspn(4)) then
          idorb2(2) = idorb(4) 
          idorb2(4) = idorb(3)
        else
          idorb2(2) = idorb(3) 
          idorb2(4) = idorb(4)
          fac = -fac
        end if

        if (idx12.gt.idx34    ) fac  = abfac * fac

        if (ms.eq.0.and.scaling.eq.1) then
          fac = -sp_fac*fac
        else if (ms.eq.0.and.scaling.eq.2) then
          fac = -(2d0*sp_fac)*fac
        end if

        idxperm = rank_ivec(rank,idorb2,4)+1
c dbg
c        print *,'X: idorb: ',idorb(1:4)
c        print *,'X: idspn: ',idspn(1:4)
c        print *,'X: idorb2: ',idorb2(1:4),'  fac = ',fac
c        print *,'rank :',rank(1:4),'  idxperm: ',idxperm
c dbg

        idx_ord(rank(1)+1) = idorb2(1)
        idx_ord(rank(2)+1) = idorb2(2)
        idx_ord(rank(3)+1) = idorb2(3)
        idx_ord(rank(4)+1) = idorb2(4)

        idx_int = idx_int_graph(idx_ord,4,iy_int,igamorb,ngam)
        idx_typ  = abs(typetab(idxperm))

c dbg
c        print *,'Xidorb2: ',idorb2(1:4)
c        print *,'rank :',rank(1:4),'  idxperm: ',idxperm
c dbg

c dbg
c        if (idxstr.le.10) then
c          if (idxstr.eq.1) print *,' -------e.g. -------------------'
c          print '("X",i4,a,4i4,a,2i4,2f16.10)',
c     &       idxstr,': (',idx_ord(1:4),') <- ',
c     &       idx_int,idx_typ, -fac,
c     &       buffer_in((idx_int-1)*ntypes+idx_typ)
c          print *,'X adr = ',(idx_int-1)*ntypes+idx_typ
c          if (idxstr.eq.10) print *,' ----------------- etc ---------'
c        end if
c dbg

          curdisblk(idxstr) = curdisblk(idxstr) -
     &                   fac*buffer_in((idx_int-1)*ntypes+idx_typ)

        end if

      end do

      return
      end subroutine

      end
