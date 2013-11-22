*------------------------------------------------------------------------*
      subroutine dia4op_r12(me_dia,h1dia,b2dia,b2off,x2dia,x2off,use_x,
     &     str_info,orb_info)
*------------------------------------------------------------------------*
*     set up diagonal hamiltonian with modifications for
*     R12-type of calculations
*     h1dia contains the diagonal one-electron Hamiltonian (Fock matrix)
*     b2dia contains the diagonal of the B matrix
*     x2dia contains the diagonal of the X matrix (if use_x is set true)
*
*     andreas, jan. 2008
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
     &     ntest = 00

      logical, intent(in) ::
     &     use_x
      real(8), intent(in) ::
     &     h1dia(*), b2dia(*), x2dia(*)
      integer, intent(in) ::
     &     b2off(*), x2off(*)
      type(me_list), intent(in) ::
     &     me_dia
      type(strinf), intent(in), target ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      logical ::
     &     first, first_str, open_close_ffdia, set_b
      integer ::
     &     iocc_cls, msamax, mscmax, ica, ihpvdx, ihpv, nsym, ngas,
     &     msa, msc, igama, igamc, igamstr, ms_str,
     &     idxms, igrph, idis,
     &     ncidx, naidx, maxidx, nidx, idx, jdx, kdx, ifree,
     &     nblk, nblkmax, idum, maxbuff, idxbuf,
     &     maxstrbuf, nstrbuf, nloop, loop_bx, loop_h, loop_p,
     &     istr, nouter, n_inner1, n_inner2, lenblk, ioff,
     &     ioff_h, ioff_p, ioff_bx, ilen,
     &     iouter, ioffbuf, ioffbuf0, ioff_inner1, ioff_inner2,
     &     idx_inner1

      integer, pointer ::
     &     lenstr_gm(:,:)

      integer ::
     &     msdst(ngastp,2), igamdst(ngastp,2),
     &     ioff_xsum(2*ngastp), nstr(2*ngastp), nincr(2*ngastp)
      real(8) ::
     &     xsum_outer, xsum_i1, fac, minel,
     &     cpu, sys, wall, cpu0, sys0, wall0

      ! some pointers for easy typing and efficency
      integer, pointer ::
     &     iocc(:,:), idx_graph(:,:),
     &     igas_restr(:,:,:,:,:), mostnd(:,:,:),
     &     igamorb(:), ngas_hpv(:), idx_gas(:)

      ! allocatable stuff:
      integer, pointer ::
     &     idxorb(:),idxspn(:),idxdss(:)
      real(8), pointer ::
     &     buffer(:), xsum(:), xmet(:)

      type(filinf), pointer ::
     &     ffdia
      type(operator), pointer ::
     &     op

      real(8), external ::
     &     dnrm2, fndmnx
      logical, external ::
     &     next_msgamdist, next_string

      ifree = mem_setmark('dia4op_r12')

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'this is dia4op_r12')
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
c        ! currently only tested for pure excitations, so:
c        if (imltlist(0,iocc(1,1),ngastp,1).lt.3.or.
c     &      imltlist(0,iocc(1,2),ngastp,1).lt.3.or.
c     &      iocc(ihole,1).gt.0 .or.
c     &      iocc(ipart,2).gt.0 ) then
c          call wrt_occ(lulog,op%ihpvca_occ(1,1,iocc_cls))
c          call quit(1,'dia4op',
c     &         'routine not tested for this kind of occupation')
c        end if
        maxbuff = max(maxbuff,me_dia%len_op_occ(iocc_cls))
      end do

      ! we should of course use a better memory management
      ! when going to large calculations
      if (maxbuff.gt.ifree)
     &     call quit(1,'dia4op_r12','inadequate memory management')

      if (ntest.ge.100)
     &     write(lulog,*)'allocating result buffer of size: ',maxbuff
        
      ifree = mem_alloc_real(buffer,maxbuff,'buffer')
      ! loop over operator elements
      occ_cls: do iocc_cls = 1, op%n_occ_cls

        if(op%formal_blk(iocc_cls))cycle

c        ! buffer: start from new
c        idxbuf = 0

        iocc => op%ihpvca_occ(1:ngastp,1:2,iocc_cls)
        idx_graph => me_dia%idx_graph(1:ngastp,1:2,iocc_cls)

        if (ntest.ge.100) then
          write(lulog,*) 'current occupation class:'
          call wrt_occ(lulog,op%ihpvca_occ(1,1,iocc_cls))
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
     &       write(lulog,*) 'allocating xsum-buffer of size ',maxstrbuf
        ! buffer for orbital energy sums
        ifree = mem_alloc_real(xsum,maxstrbuf,'xsum')
        ifree = mem_alloc_real(xmet,maxstrbuf,'xmet')

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
                write(lulog,*) 'msc,msa,igamc,igama: ',
     &               msc,msa,igamc,igama
                write(lulog,*) 'current distribution (MS,GAMMA):'
                call wrt_occ(lulog,msdst)
                call wrt_occ(lulog,igamdst)
              end if

              ! strings for HPV/CA
              nloop = 0
              istr = 0
              loop_bx = -1
              loop_h  = -1
              loop_p  = -1
              ! sequence: from outer to inner loop
              do ihpvdx = ngastp, 1, -1
                ihpv = hpvxseq(ihpvdx)
                do ica = 1, 2                  
                  if (iocc(ihpv,ica).eq.0) cycle
                  nloop = nloop+1
                  set_b = ihpv.eq.IHOLE.and.ica.eq.1
                  ! Just reusing the code from dia4op, not saying
                  ! that I now do the most general thing:
                  ! Set variables to address the loops associated
                  ! with r12-geminal, hole and particle indices
                  if (set_b) loop_bx = nloop
                  if (ihpv.eq.IHOLE.and.ica.eq.2)
     &                 loop_h = nloop
                  if (ihpv.eq.IPART.and.ica.eq.2)
     &                 loop_p = nloop
                  fac = 1d0
                  if (ica.eq.2) fac = -1d0
                  first_str = .true.
                  ioff_xsum(nloop) = istr
                  nidx = iocc(ihpv,ica)
                  if (set_b.and.nidx.ne.2)
     &                 call quit(1,'dia4op_r12','unexpected event!')
                  igrph = idx_graph(ihpv,ica)
                  ms_str = msdst(ihpv,ica)
                  igamstr = igamdst(ihpv,ica)
                  if (set_b)
     &                 lenstr_gm => str_info%g(igrph)%lenstr_gm
c dbg
c                  print *,'ihpv,ica,iocc: ',ihpv,ica,nidx,igrph,ms_str,
c     &                 igamstr
c dbg
                  if (.not.set_b) then
                    str_loop_h: do
                      if (.not.next_string(idxorb,idxspn,idxdss,
     &                   nidx,ms_str,igamstr,first_str,
     &                   igas_restr(1,1,1,1,igrph),
     &                   mostnd(1,1,idx_gas(ihpv)),igamorb,
     &                   nsym,ngas_hpv(ihpv))
     &                   ) exit str_loop_h

                      first_str = .false.
                                
                      istr = istr+1
                      xsum(istr) = 0d0
                      do idx = 1, nidx
                        ! need to patch for UHF:
                        xsum(istr) = xsum(istr) + fac*h1dia(idxorb(idx))
                      end do
                      xmet(istr) = 1d0
                    end do str_loop_h
                  else
                    ! exactly 2-membered strings, so we can 
                    ! a) calculate the index of the twodia block
                    idxms = (2 - ms_str)/2 + 1
                    idx = (idxms-1)*nsym + igamstr
                    ioff = b2off(idx)
                    ilen = lenstr_gm(igamstr,idxms) 
                    xsum(istr+1:istr+ilen) = !xsum(istr+1:istr+ilen) + 
     &                   b2dia(ioff+1:ioff+ilen)
                    if (use_x) then
                      ioff = x2off(idx)
                      xmet(istr+1:istr+ilen) =
     &                     x2dia(ioff+1:ioff+ilen)
                    else
                      xmet(istr+1:istr+ilen) = 1d0
                    end if
                    istr = istr+ilen
                  end if

                  nstr(nloop) = istr-ioff_xsum(nloop)
                  ! do not consider empty strings
                  if (nstr(nloop).eq.0) nloop = nloop-1
                end do
              end do

c dbg
c              do idx = 1, nloop
c                print *,'idx, off, len: ',idx, ioff_xsum(idx), nstr(idx)
c                print *,'xsum = ',xsum(ioff_xsum(idx)+1:
c     &                                 ioff_xsum(idx)+nstr(idx))
c                print *,'xmet = ',xmet(ioff_xsum(idx)+1:
c     &                                 ioff_xsum(idx)+nstr(idx))
c              end do
c              print *,'nloop',nloop
c              print *,'loop_h, loop_p, loop_bx: ',
c     &                 loop_h, loop_p, loop_bx
c dbg

              ! R12-doubles case:
              if (nloop.eq.2.and.loop_p.eq.-1) then

                if (loop_h.eq.-1.or.loop_bx.eq.-1)
     &               call quit(1,'dia4op_r12',
     &               'Something is wrong!')

                ! hard-wired implementation:
                ioff_h = ioff_xsum(loop_h)
                ioff_bx = ioff_xsum(loop_bx)
                idxbuf = 0
                do idx = 1, nstr(loop_h)
                  xsum_outer = xsum(ioff_h+idx)
                  do jdx = 1, nstr(loop_bx)
                    idxbuf = idxbuf+1
                    buffer(ioffbuf0+idxbuf) =
     &                   xsum_outer*xmet(ioff_bx+jdx)+xsum(ioff_bx+jdx)
                  end do
                end do
                ioffbuf0 = ioffbuf0+idxbuf
                cycle
              end if

              if (nloop.ne.3)then
c     &             call quit(1,'dia4op_r12',
c     &           'My assumptions are not fulfilled')
                write(lulog,*) 'Warning: R12 N-tuples Cycle 1'
                cycle
              endif

              ! R12 n-tuples case
              if (loop_h.eq.-1.or.loop_bx.eq.-1.or.loop_p.eq.-1)then
c     &             call quit(1,'dia4op_r12',
c     &               'Something is wrong (2)!')
                write(lulog,*) 'Warning: R12 N-tuples Cycle 2'
                cycle
              endif

              ! hard-wired implementation:
              ioff_h = ioff_xsum(loop_h)
              ioff_p = ioff_xsum(loop_p)
              ioff_bx = ioff_xsum(loop_bx)
              idxbuf = 0
              do idx = 1, nstr(loop_h)
                xsum_outer = xsum(ioff_h+idx)
                do jdx = 1, nstr(loop_bx)
                  xsum_i1 =
     &                 xsum_outer*xmet(ioff_bx+jdx)+xsum(ioff_bx+jdx)
                  do kdx = 1, nstr(loop_p)
                    idxbuf = idxbuf+1
                    buffer(ioffbuf0+idxbuf) =
     &                   xsum_i1+xsum(ioff_p+idx)
                  end do
                end do
              end do
              ioffbuf0 = ioffbuf0+idxbuf
              
            end do distr_loop

          end do igama_loop
        end do msa_loop

c dbg
c        if (ntest.ge.100)
c     &       print *,'final buffer: ',
c     &       buffer(1:me_dia%len_op_occ(iocc_cls))
c dbg
        ! shift
        ilen = me_dia%len_op_occ(iocc_cls)
        minel = fndmnx(buffer,ilen,-1)
        if (minel.lt.0d0) then
          write(lulog,*)
     &         'Negative lowest element of diagonal for block ',
     &         iocc_cls,' : ',minel
          write(lulog,*)
     &         'I shift by ',+0.5d0-minel          
          buffer(1:ilen) = buffer(1:ilen)+0.5d0-minel
        end if

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
     &     call prtim(lulog,'time in dia4op_r12 ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
