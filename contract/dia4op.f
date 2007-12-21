*------------------------------------------------------------------------*
      subroutine dia4op(me_dia,xdia,str_info,orb_info)
*------------------------------------------------------------------------*
*     set up diagonal hamiltonian
*     xdia contains the diagonal one-electron Hamiltonian (Fock matrix)
*     andreas, spring 2007
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
     &     ntest = 100

      real(8), intent(in) ::
     &     xdia(*)
      type(me_list), intent(in) ::
     &     me_dia
      type(strinf), intent(in), target ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      logical ::
     &     first, first_str, open_close_ffdia
      integer ::
     &     iocc_cls, msamax, mscmax, ica, ihpvdx, ihpv, nsym, ngas,
     &     msa, msc, igama, igamc, igamstr, ms_str,
     &     idxms, igrph, idis,
     &     ncidx, naidx, maxidx, nidx, idx, ifree,
     &     nblk, nblkmax, idum, maxbuff,
     &     maxstrbuf, nstrbuf, nloop,
     &     istr, nouter, n_inner1, n_inner2, lenblk,
     &     iouter, ioffbuf, ioffbuf0, ioff_inner1, ioff_inner2,
     &     idxbuf, idx_inner1

      integer ::
     &     msdst(ngastp,2), igamdst(ngastp,2),
     &     ioff_xsum(2*ngastp), nstr(2*ngastp), nincr(2*ngastp)
      real(8) ::
     &     xsum_outer, xsum_i1, fac,
     &     cpu, sys, wall, cpu0, sys0, wall0

      ! some pointers for easy typing and efficency
      integer, pointer ::
     &     iocc(:,:), idx_graph(:,:),
     &     igas_restr(:,:,:,:), mostnd(:,:,:),
     &     igamorb(:), ngas_hpv(:), idx_gas(:)

      ! allocatable stuff:
      integer, pointer ::
     &     idxorb(:),idxspn(:),idxdss(:)
      real(8), pointer ::
     &     buffer(:), xsum(:)
      type(filinf), pointer ::
     &     ffdia
      type(operator), pointer ::
     &     op

      real(8), external ::
     &     dnrm2
      logical, external ::
     &     next_msgamdist, next_string

      ifree = mem_setmark('dia4op')

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'this is dia4op')
      end if

      op => me_dia%op
      ffdia => me_dia%fhand

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
      do iocc_cls = 1, op%n_occ_cls
        iocc => op%ihpvca_occ(1:ngastp,1:2,iocc_cls)
        ! currently only tested for pure excitations, so:
        if (imltlist(0,iocc(1,1),ngastp,1).lt.3.or.
     &      imltlist(0,iocc(1,2),ngastp,1).lt.3.or.
     &      iocc(ihole,1).gt.0 .or.
     &      iocc(ipart,2).gt.0 ) then
          call wrt_occ(luout,op%ihpvca_occ(1,1,iocc_cls))
          call quit(1,'dia4op',
     &         'routine not tested for this kind of occupation')
        end if
        maxbuff = max(maxbuff,me_dia%len_op_occ(iocc_cls))
      end do

      ! we should of course use a better memory management
      ! when going to large calculations
      if (maxbuff.gt.ifree)
     &     call quit(1,'dia4op','inadequate memory management')

      if (ntest.ge.100)
     &     write(luout,*)'allocating result buffer of size: ',maxbuff
        
      ifree = mem_alloc_real(buffer,maxbuff,'buffer')
      ! loop over operator elements
      occ_cls: do iocc_cls = 1, op%n_occ_cls

        if(op%formal_blk(iocc_cls))cycle

        ! buffer: start from new
        idxbuf = 0

        iocc => op%ihpvca_occ(1:ngastp,1:2,iocc_cls)
        idx_graph => me_dia%idx_graph(1:ngastp,1:2,iocc_cls)

        if (ntest.ge.100) then
          write(luout,*) 'current occupation class:'
          call wrt_occ(luout,op%ihpvca_occ(1,1,iocc_cls))
        end if

        ncidx = op%ica_occ(1,iocc_cls)
        naidx = op%ica_occ(2,iocc_cls)
        maxidx = max(ncidx,naidx)
       
        mscmax = ncidx
        msamax = naidx

        ifree = mem_setmark('dia4op_str')
        ifree = mem_alloc_int(idxorb,maxidx,'idxorb')
        ifree = mem_alloc_int(idxspn,maxidx,'idxspn')
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

          igama_loop: do igama = 1, orb_info%nsym
            igamc = multd2h(igama,me_dia%gamt)
            
            if (me_dia%len_op_gmo(iocc_cls)%gam_ms(igama,idxms).eq.0)
     &           cycle igama_loop

            idis = 0
            first = .true.
            ! for pure P/H excitations we will not need this:
            ! separate loop / procedure ?
            distr_loop: do

              if (.not.next_msgamdist(first,
     &           msa,msc,igama,igamc,iocc,
     &                   orb_info%nsym,
     &           msdst,igamdst)) exit distr_loop
              first = .false.
              idis = idis+1

              if (me_dia%len_op_gmox(iocc_cls)%
     &             d_gam_ms(idis,igama,idxms).eq.0)
     &             cycle distr_loop
              
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
              ! sequence: from outer to inner loop
              do ihpvdx = ngastp, 1, -1
                ihpv = hpvxseq(ihpvdx)
                do ica = 1, 2                  
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
     &                 igas_restr(1,1,1,igrph),
     &                 mostnd(1,1,idx_gas(ihpv)),igamorb,
     &                 nsym,ngas_hpv(ihpv))
     &                 ) exit str_loop

                    first_str = .false.
                                
                    istr = istr+1
                    xsum(istr) = 0d0
                    do idx = 1, nidx
                      ! need to patch for UHF:
                      xsum(istr) = xsum(istr) + fac*xdia(idxorb(idx))
                    end do
                  end do str_loop
                  nstr(nloop) = istr-ioff_xsum(nloop)
                  ! do not consider empty strings
                  if (nstr(nloop).eq.0) nloop = nloop-1
                end do
              end do
c dbg
              do idx = 1, nloop
                print *,'idx, off, len: ',idx, ioff_xsum(idx), nstr(idx)
                print *,'xsum = ',xsum(ioff_xsum(idx)+1:
     &                                 ioff_xsum(idx)+nstr(idx))
              end do
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
                idxbuf = ioffbuf+1
                xsum_outer = 0d0
                do iouter = 1, nouter
                  idx = ioff_xsum(iouter) + idxbuf/nincr(iouter)+1
                  idxbuf = mod(idxbuf,nincr(iouter)) 
                  xsum_outer = xsum_outer + xsum(idx)
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
              ioffbuf0 = ioffbuf

            end do distr_loop

          end do igama_loop
        end do msa_loop

c dbg
        print *,'final buffer: ',buffer(1:me_dia%len_op_occ(iocc_cls))
c dbg
        ! put buffer to disc
        call put_vec(ffdia,buffer,me_dia%off_op_occ(iocc_cls)+1,
     &                            me_dia%off_op_occ(iocc_cls)
     &                           +me_dia%len_op_occ(iocc_cls))
        
        ifree = mem_flushmark()

      end do occ_cls

      if (open_close_ffdia) call file_close_keep(ffdia)

      ifree = mem_flushmark()

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5)
     &     call prtim(luout,'time in dia4op ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
