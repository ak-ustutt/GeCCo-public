*------------------------------------------------------------------------*
      subroutine dia4op_ev(me_dia,ecore,xdia1,xdia2,str_info,orb_info)
*------------------------------------------------------------------------*
*     set up diagonal hamiltonian
*     ecore contains the core energy
*     xdia1 contains the diagonal one-electron Hamiltonian (Fock matrix)
*     xdia2 contains the diagonal two-electron Hamiltonian
*     matthias, fall 2009
*------------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'
      include 'ifc_baserout.h'
      include 'hpvxseq.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00,
     &     nlistmax = 120

      real(8), intent(in) ::
     &     ecore, xdia1(*), xdia2(*)
      type(me_list), intent(in) ::
     &     me_dia
      type(strinf), intent(in), target ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      logical ::
     &     first, first_str, open_close_ffdia, pair
      integer ::
     &     iblk, msamax, mscmax, ica, ihpvdx, ihpv, nsym, ngas,
     &     msa, msc, igama, igamc, igamstr, ms_str,
     &     idxms, igrph, ndis, idis, did, njoined,
     &     ncidx, naidx, maxidx, nidx, idx, ifree, jdx,
     &     nblk, nblkmax, idum, maxbuff,
     &     maxstrbuf, nstrbuf, nloop,
     &     istr, nouter, n_inner1, n_inner2, lenblk,
     &     iouter, ioffbuf, ioffbuf0, ioff_inner1, ioff_inner2,
     &     idxbuf, idx_inner1, x2_off, idxms2, ms2,
     &     ijoff, iincr, jincr, ni, nj, ilen, jlen, jhpv, iloop,
     %     jloop, jca, iorb, jorb, jdxms2, ims, jms, ioff, joff,
     &     kloop, jhpvdx, maxca, ii, jj, iofft

      integer ::
     &     msdst(ngastp,2), igamdst(ngastp,2),
     &     ioff_xsum(2*ngastp), nstr(2*ngastp), nincr(2*ngastp),
     &     iocc2(ngastp,2)
      real(8) ::
     &     xsum_outer, xsum_i1, fac, val,
     &     cpu, sys, wall, cpu0, sys0, wall0

      ! some pointers for easy typing and efficency
      integer, pointer ::
     &     iocc(:,:), idx_graph(:,:),
     &     igas_restr(:,:,:,:,:), mostnd(:,:,:),
     &     igamorb(:), ngas_hpv(:), idx_gas(:)

      ! allocatable stuff:
      integer, pointer ::
     &     idxorb(:),idxspn(:),idxdss(:),idxspn2(:)
      real(8), pointer ::
     &     buffer(:), xsum(:)
      integer, allocatable ::
     &     list(:,:,:,:)

      type(filinf), pointer ::
     &     ffdia
      type(operator), pointer ::
     &     op

      real(8), external ::
     &     dnrm2
      logical, external ::
     &     next_msgamdist, next_string

      ifree = mem_setmark('dia4op_ev')

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'this is dia4op_ev')
      end if

      op => me_dia%op
      ffdia => me_dia%fhand

      njoined = op%njoined

      if (njoined.ne.1)
     &     call quit(1,'dia4op_ev','njoined!=1 occurred!')

      call atim_csw(cpu0,sys0,wall0)
      
      igas_restr => str_info%igas_restr
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      ngas_hpv => orb_info%ngas_hpv
      idx_gas => orb_info%idx_gas
      nsym = orb_info%nsym
      ngas = orb_info%ngas

      open_close_ffdia = ffdia%unit.le.0
      if (open_close_ffdia) call file_open(ffdia)

      maxbuff=0
      do iblk = 1, op%n_occ_cls
        iocc => op%ihpvca_occ(1:ngastp,1:2,iblk)
c        ! currently only tested for pure excitations, so:
c        if (imltlist(0,iocc(1:,1),ngastp,1).lt.3.or.
c     &      imltlist(0,iocc(1:,2),ngastp,1).lt.3.or.
c     &      iocc(ihole,1).gt.0 .or.
c     &      iocc(ipart,2).gt.0 ) then
c          call wrt_occ(luout,op%ihpvca_occ(1,1,iblk))
c          call quit(1,'dia4op_ev',
c     &         'routine not tested for this kind of occupation')
c        end if
        maxbuff = max(maxbuff,me_dia%len_op_occ(iblk))
      end do

      ! we should of course use a better memory management
      ! when going to large calculations
      if (maxbuff.gt.ifree) then
        write(luout,*) 'maxbuff, ifree: ',maxbuff,ifree
        call mem_map(.true.)
        call quit(1,'dia4op_ev','inadequate memory management')
      end if

      if (ntest.ge.100)
     &     write(luout,*)'allocating result buffer of size: ',maxbuff
        
      ifree = mem_alloc_real(buffer,maxbuff,'buffer')
      allocate(list(orb_info%ntoob,2,0:nlistmax,2*ngastp))
      ! loop over operator elements
      occ_cls: do iblk = 1, op%n_occ_cls

        if(op%formal_blk(iblk))cycle

        ! buffer: start from new
        idxbuf = 0
        buffer(1:maxbuff) = 0d0

        iocc => op%ihpvca_occ(1:ngastp,1:2,iblk)
        iocc2 = iocc
        idx_graph => me_dia%idx_graph(1:ngastp,1:2,iblk)

        if (ntest.ge.100) then
          write(luout,*) 'current occupation class:'
          call wrt_occ(luout,op%ihpvca_occ(1,1,iblk))
        end if

        ncidx = op%ica_occ(1,iblk)
        naidx = op%ica_occ(2,iblk)
        maxidx = max(ncidx,naidx)
       
        mscmax = ncidx
        msamax = naidx

        ifree = mem_setmark('dia4op_ev_str')
        ifree = mem_alloc_int(idxorb,maxidx,'idxorb')
        ifree = mem_alloc_int(idxspn,maxidx,'idxspn')
        ifree = mem_alloc_int(idxspn2,maxidx,'idxspn2')
        ifree = mem_alloc_int(idxdss,maxidx,'idxdss')

        ioffbuf0 = 0
        ! inquire maximum number of strings:
        maxstrbuf = 0
        do msa = msamax, -msamax, -2
          msc = msa + me_dia%mst
          if (abs(msc).gt.mscmax) cycle
          do igama = 1, nsym
            igamc = multd2h(igama,me_dia%gamt)
            first = .true.
            distr_loop0: do
              if (.not.next_msgamdist(first,
     &           msa,msc,igama,igamc,iocc,nsym,
     &           msdst,igamdst)) exit distr_loop0
              first = .false.

              nstrbuf = 0
              do ica = 1, 2
                do ihpv = 1, ngastp
                  if (iocc(ihpv,ica).eq.0) cycle
                  igrph = idx_graph(ihpv,ica)
                  igamstr = igamdst(ihpv,ica)
                  idxms = (iocc(ihpv,ica)-msdst(ihpv,ica))/2+1
                  nstrbuf = nstrbuf +
     &                 str_info%g(igrph)%lenstr_gm(igamstr,idxms)
                end do
              end do
              maxstrbuf = max(maxstrbuf,nstrbuf)
            end do distr_loop0
          end do
        end do

        if (ntest.ge.100)
     &       write(luout,*) 'allocating xsum-buffer of size ',maxstrbuf
        ! buffer for orbital energy sums
        ifree = mem_alloc_real(xsum,maxstrbuf,'xsum')

        idxms = 0
        msa_loop: do msa = msamax, -msamax, -2

          ! C <-> A  means alpha <-> beta !!
          msc = msa + me_dia%mst          
          if (abs(msc).gt.mscmax) cycle msa_loop
          idxms = idxms+1

          igama_loop: do igama = 1, nsym
            igamc = multd2h(igama,me_dia%gamt)
            
            if (me_dia%len_op_gmo(iblk)%gam_ms(igama,idxms).eq.0)
     &           cycle igama_loop

c            idis = 0
c            first = .true.
            ! for pure P/H excitations we will not need this:
            ! separate loop / procedure ?
            ndis = me_dia%off_op_gmox(iblk)%ndis(igama,idxms)
            distr_loop: do idis = 1, ndis

              did = me_dia%off_op_gmox(iblk)%did(idis,igama,idxms)
            
              call did2msgm(msdst,igamdst,did,
     &               iocc,nsym,njoined)

c              if (.not.next_msgamdist(first,
c     &           msa,msc,igama,igamc,iocc,
c     &                   orb_info%nsym,
c     &           msdst,igamdst)) exit distr_loop
c              first = .false.
c              idis = idis+1

c              if (me_dia%len_op_gmox(iblk)%
c     &             d_gam_ms(idis,igama,idxms).eq.0)
c     &             cycle distr_loop
              
              if (ntest.ge.100) then
                write(luout,*) 'msc,msa,igamc,igama: ',
     &               msc,msa,igamc,igama
                write(luout,*) 'current distribution (MS,GAMMA):'
                call wrt_occ(luout,msdst)
                call wrt_occ(luout,igamdst)
              end if

              ! strings for HPV/CA
              nloop = 0
              istr = 0
              list(1:orb_info%ntoob,1:2,0:nlistmax,1:2*ngastp) = 0
              ! sequence: from outer to inner loop
              do ihpvdx = ngastp, 1, -1
                ihpv = hpvxseq(ihpvdx)
                do ica = 2, 1, -1 
                  if (iocc(ihpv,ica).eq.0) cycle
                  fac = 1d0
                  if (ica.eq.2) fac = -1d0
                  nloop = nloop+1
                  first_str = .true.
                  ioff_xsum(nloop) = istr
                  nidx = iocc(ihpv,ica)
                  igrph = idx_graph(ihpv,ica)
                  ms_str = msdst(ihpv,ica)
                  igamstr = igamdst(ihpv,ica)
c dbg
c                  print *,'ihpv,ica,iocc: ',ihpv,ica,nidx,igrph,ms_str,
c     &                 igamstr
c dbg
                  str_loop: do
                    if (.not.next_string(idxorb,idxspn,idxdss,
     &                 nidx,ms_str,igamstr,first_str,
     &                 igas_restr(1,1,1,1,igrph),
     &                 mostnd(1,1,idx_gas(ihpv)),igamorb,
     &                 nsym,ngas_hpv(ihpv))
     &                 ) exit str_loop

                    first_str = .false.
                                
                    istr = istr+1
c dbg
c          print *,'istr: ',istr
c          write(luout,'(x,a,8i3)') 'idxorb / idxspn: ',
c     &       idxorb(1:nidx), idxspn(1:nidx)
c dbgend
                    xsum(istr) = 0d0
                    idxspn2 = idxspn
                    pair = .false.
                    do idx = 1, nidx
                      ! use that orbitals are in ascending order
                      if (idxspn2(idx).eq.2.and.pair) then
                        idxspn2(idx) = -1
                        pair = .false.
                      else if (idxspn2(idx).eq.2) then
                        idxspn2(idx) = 1
                        pair = .true.
                      end if

                      ! update list orbital --> string indices
                      list(idxorb(idx),(idxspn2(idx)+3)/2,0,nloop)
     &                  = list(idxorb(idx),(idxspn2(idx)+3)/2,0,nloop)+1
                      if (list(idxorb(idx),(idxspn2(idx)+3)/2,0,nloop)
     &                    .gt.nlistmax) call quit(1,'dia4op_ev',
     &                    'Increase parameter nlistmax')
                      list(idxorb(idx),(idxspn2(idx)+3)/2,
     &                     list(idxorb(idx),(idxspn2(idx)+3)/2,0,nloop),
     &                     nloop) = istr - ioff_xsum(nloop)

                      ! need to patch for UHF:
                      xsum(istr) = xsum(istr) + fac*xdia1(idxorb(idx))
c dbg
c          print *,'added ',fac*xdia1(idxorb(idx))
c dbgend

                      do jdx = 1, idx-1
                        ms2 = idxspn2(idx) + idxspn2(jdx)
                        idxms2 = (2 - ms2)/2 + 1
                        x2_off = (idxms2 - 1) * orb_info%ntoob**2
                        xsum(istr) = xsum(istr) + xdia2(x2_off
     &                     + (idxorb(jdx)-1)*orb_info%ntoob+idxorb(idx))
c dbg
c                 write(luout,'(x,a,2i4,x,i4)') 'orb/idxms: ',
c     &               idxorb(jdx),idxorb(idx),idxms2
c                 print *,'added element ',x2_off
c     &                     + (idxorb(jdx)-1)*orb_info%ntoob+idxorb(idx),
c     &                   xdia2(x2_off
c     &                     + (idxorb(jdx)-1)*orb_info%ntoob+idxorb(idx))
c dbgend
                      end do
                    end do
                  end do str_loop
                  nstr(nloop) = istr-ioff_xsum(nloop)
                  ! do not consider empty strings
                  if (nstr(nloop).eq.0) then
                    nloop = nloop-1
                    iocc2(ihpv,ica) = 0
                  end if
     
                end do
              end do
c dbg
c              do idx = 1, nloop
c                print *,'idx, off, len: ',idx, ioff_xsum(idx), nstr(idx)
c                print *,'xsum = ',xsum(ioff_xsum(idx)+1:
c     &                                 ioff_xsum(idx)+nstr(idx))
c              end do
c dbg
              ! no contribution at all?
              if (nloop.eq.0) cycle

              ! only one loop?
              if (nloop.eq.1) then
                buffer(ioffbuf0+1:ioffbuf0+nstr(1)) =
     &               xsum(1:nstr(1))
                ioffbuf0 = ioffbuf0+nstr(1)
                cycle
              end if

              ! more than one loop:

              ! setup block-length information
              nouter = nloop-2
              n_inner2 = nstr(nloop)
              n_inner1 = nstr(nloop-1)
              lenblk = n_inner1*n_inner2
              do iouter = nouter, 1, -1
                nincr(iouter) = lenblk
                lenblk = lenblk*nstr(iouter)
              end do

              ! assemble the actual diagonal:
              ioffbuf = ioffbuf0
              ioff_inner1 = ioff_xsum(nouter+1)
              ioff_inner2 = ioff_xsum(nouter+2)
              do while(ioffbuf.lt.ioffbuf0+lenblk)
                ! collect contributions from formal outer loops 
                ! (over more than two strings):
                idxbuf = ioffbuf - ioffbuf0 !+1
                xsum_outer = 0d0
                do iouter = 1, nouter
                  idx = ioff_xsum(iouter) + idxbuf/nincr(iouter)+1
                  xsum_outer = xsum_outer + xsum(idx)
c dbg
c                  if (idx.gt.nstr(iouter).or.idx.lt.1) then
c                    write(6,*) 'range(O) iouter,idx,nstr: ',
c     &                 iouter,idx,nstr(iouter)
c                    write(6,*) 'range(O) '//
c     &                   'idxbuf,nincr(iouter),ioff_xsum(iouter): ',
c     &                   idxbuf,nincr(iouter),
c     &                 ioff_xsum(iouter)
c                    write(6,*) 'off,idx/nicr: ',
c     &                   ioff_xsum(iouter),idxbuf/nincr(iouter),
c     &                   ioff_xsum(iouter) + idxbuf/nincr(iouter)+1
c                  end if
c dbg
                  idxbuf = mod(idxbuf,nincr(iouter)) 
                end do
                ! explicit loops for the two innermost strings
                do idx_inner1 = ioff_inner1+1, ioff_inner1+n_inner1
                  xsum_i1 = xsum(idx_inner1) + xsum_outer
c dbg
c                  if (ioffbuf+1.lt.1) write(6,*)'range1:',ioffbuf+1
c                  if (ioffbuf+n_inner2.gt.maxbuff)
c     &                 write(6,*)'range2:',ioffbuf+n_inner2
c dbg
                  buffer(ioffbuf+1:ioffbuf+n_inner2) =
     &                xsum(ioff_inner2+1:ioff_inner2+n_inner2) + xsum_i1
                  ioffbuf = ioffbuf+n_inner2
                end do
              end do

              ! add two-electron integrals between separate strings:

              ! sequence: from inner to outer loop
              jloop = 0
              do jhpvdx = 1, ngastp
               jhpv = hpvxseq(jhpvdx)
               do jca = 1, 2
                if (iocc2(jhpv,jca).eq.0) cycle
                jloop = jloop + 1
                iloop = 0
                do ihpvdx = 1, jhpvdx
                 ihpv = hpvxseq(ihpvdx)
                 maxca = 2
                 if (ihpvdx.eq.jhpvdx) maxca = jca - 1
                 do ica = 1, maxca
                  fac = 1d0 - 2d0*abs(jca-ica)
                  if (iocc2(ihpv,ica).eq.0) cycle
                  iloop = iloop + 1
c dbg
c                  write(luout,'(a,6i4)') 'in loop:',
c     &                   iloop,jloop,ica,ihpv,jca,jhpv
c                  write(luout,'(a,2i4)') 'len:',
c     &                   nstr(nloop+1-iloop),nstr(nloop+1-jloop)
c dbgend
                  ! loop over all possible spin combinations ++/+-/-+/--
                  do jms = -1, 1, 2
                   if (abs(msdst(jhpv,jca)-jms).gt.iocc2(jhpv,jca))
     &                  cycle
                   do ims = -1, 1, 2
                    if (abs(msdst(ihpv,ica)-ims).gt.iocc2(ihpv,ica))
     &                   cycle
c dbg
c                    write(luout,'(a,2i4)') 'ims, jms: ',ims,jms
c dbgend
                    idxms2 = (2 - ims - jms)/2 + 1
                    x2_off = (idxms2 - 1) * orb_info%ntoob**2
                    ioff = sum(orb_info%igassh(1:nsym,
     &                          1:orb_info%ioff_gas(ihpv)))
                    joff = sum(orb_info%igassh(1:nsym,
     &                          1:orb_info%ioff_gas(jhpv)))
                    ilen = 1
                    nj = 1
                    ni = 1
                    do kloop = 1, iloop-1
                      ilen = ilen*nstr(nloop+1-kloop) !reverse loop order
                    end do
                    jlen = ilen
                    do kloop = iloop, jloop-1
                      jlen = jlen*nstr(nloop+1-kloop)
                    end do
                    do kloop = jloop+1, nloop
                      nj = nj*nstr(nloop+1-kloop)
                    end do
                    do kloop = iloop+1, jloop-1
                      ni = ni*nstr(nloop+1-kloop)
                    end do
                    iincr = nstr(nloop+1-iloop)*ilen
                    jincr = nstr(nloop+1-jloop)*jlen
                    !loop over diagonal elements of hamiltonian
                    do jorb = 1, orb_info%norb_hpv(jhpv,1)
                     do iorb = 1, orb_info%norb_hpv(ihpv,1)
                      val = fac*xdia2(x2_off
     &                     + (min(ioff+iorb,joff+jorb)-1)*orb_info%ntoob
     &                     + max(ioff+iorb,joff+jorb))
c dbg
c                      write(luout,'(2i4,E19.10)')
c     &                  ioff+iorb,joff+jorb,val
c                      write(luout,'(a,20i4)') 'list i: ',
c     &                  list(ioff+iorb,(ims+3)/2,
c     &                  1:list(ioff+iorb,(ims+3)/2,0,
c     &                  nloop+1-iloop),nloop+1-iloop)
c                      write(luout,'(a,20i4)') 'list j: ',
c     &                  list(joff+jorb,(jms+3)/2,
c     &                  1:list(joff+jorb,(jms+3)/2,0,
c     &                  nloop+1-jloop),nloop+1-jloop)
c dbgend
                      !loop over matching string indices
                      do idx = 1, list(ioff+iorb,(ims+3)/2,0,
     &                    nloop+1-iloop)
                       do jdx = 1, list(joff+jorb,(jms+3)/2,0,
     &                     nloop+1-jloop)
                        ijoff = (list(ioff+iorb,(ims+3)/2,idx,
     &                      nloop+1-iloop)-1)*ilen +
     &                      (list(joff+jorb,(jms+3)/2,jdx,
     &                      nloop+1-jloop)-1)*jlen
                        !loop over corresponding buffer elements
                        do ii = 0, ni-1
                         do jj = 0, nj-1
                          iofft = ioffbuf0 + ijoff + ii*iincr + jj*jincr
                          buffer(iofft+1:iofft+ilen) =
     &                        buffer(iofft+1:iofft+ilen) + val
                         end do
                        end do
                       end do
                      end do
                     end do
                    end do
                   end do
                  end do
                 end do
                end do
               end do
              end do

              ioffbuf0 = ioffbuf

            end do distr_loop

          end do igama_loop
        end do msa_loop

        ! add core energy
        buffer(1:me_dia%len_op_occ(iblk))
     &     = buffer(1:me_dia%len_op_occ(iblk)) + ecore
c dbg
c        if (ntest.ge.100)
c     &       print *,'final buffer: ',
c     &       buffer(1:me_dia%len_op_occ(iblk))
c dbg
        ! put buffer to disc
        call put_vec(ffdia,buffer,me_dia%off_op_occ(iblk)+1,
     &                            me_dia%off_op_occ(iblk)
     &                           +me_dia%len_op_occ(iblk))
        
        ifree = mem_flushmark()

      end do occ_cls

      if (open_close_ffdia) call file_close_keep(ffdia)

      ifree = mem_flushmark()
      deallocate(list)

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5)
     &     call prtim(luout,'time in dia4op_ev ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
