*----------------------------------------------------------------------*
      subroutine import_h2_molpro_ifc(lumoint,ifpos,hlist,fock,nfock,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     read through moints file of molpro and sort integral into list hlist;
*     no batch algorithm (consider writing a proper interface then)
*
*     andreas, january 2016
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'ifc_baserout.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     lumoint, nfock
      integer, intent(inout) ::
     &     ifpos
      real(8), intent(inout) ::
     &     fock(nfock)
      type(me_list), intent(in) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     int_disk, int_nonr, int_ordr, ioff, istr,
     &     idxst, idxnd, nblk, nblkmax, ifree, luerr, nbuff,
     &     len_op, idum, iblk, iblkst, idxfl
      logical ::
     &     error, closeit, avoid_dc

      real(8), pointer ::
     &     h2scr(:)

      real(8) ::
     &     val
      integer(2) ::
     &     idxrd(4)
      integer ::
     &     ip,iq,ir,is,nstr,ipq,irs,ioff_core,
     &     idxprqs(4),idss(4),igam(4),igtp(4), 
     &     lexlscr(4,3), idxstr(8)

      type(filinf), pointer ::
     &     ffham

c      type(operator), pointer ::
c     &     hop

      logical, pointer ::
     &     csorb(:)

      integer, pointer ::
     &     ireost(:), ihpvgas(:,:), igamorb(:),
     &     igasorb(:), idx_gas(:), iad_gas(:)

      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall, energy


      call atim_csw(cpu0,sys0,wall0)

      ffham => hlist%fhand

      if (ffham%unit.le.0) then
        call file_open(ffham)
        closeit = .true.
      else
        closeit = .false.
      end if

      ifree = mem_setmark('import_h2')
      
      iblkst = 1
      do iblk = hlist%op%n_occ_cls,1,-1
        if (hlist%op%ica_occ(1,iblk).eq.2) iblkst = iblk
      end do
      len_op = hlist%len_op
      nblk = min((len_op-hlist%off_op_occ(iblkst)-1)
     &                      / ffham%reclen+1,nblkmax)

      nbuff = len_op-hlist%off_op_occ(iblkst)

      write(lulog,'(x,a,f9.2,a)') 'size of buffer in geth2:   ',
     &     dble(nbuff)/128d0/1024d0, 'Mb'

      ifree = mem_alloc_real(h2scr,nbuff,'h2sort_buff')

      int_disk = 0
      int_nonr = 0
      int_ordr = 0

      ! note that the first ncore orbitals have been deleted
      ioff_core = orb_info%ncore_mpro

      ! dereference structure components for efficiency
      ireost => orb_info%ireost
      ihpvgas => orb_info%ihpvgas
      igamorb => orb_info%igamorb
      igasorb => orb_info%igasorb
      idx_gas => orb_info%idx_gas
      iad_gas => orb_info%iad_gas

      ! identify closed-shell orbitals (in molpro order)
      ! for setting up contributions to core Fock and core energy
      ! noting that "core" orbitals are absent and their contribution
      ! is already included in fock and reference energy
      allocate(csorb(orb_info%ntoob))
      call set_csorbs_molpro(csorb,orb_info)

      ! loop over file
      idxfl = ifpos

      idxst = hlist%off_op_occ(iblkst) + 1
      idxnd = idxst+nbuff-1
      ioff = -idxst+1
      h2scr(1:nbuff) = 0d0
      do 
        read(lumoint,pos=idxfl) val,idxrd(1:4)
        ip=idxrd(1); iq=idxrd(2); ir=idxrd(3); is=idxrd(4)

        ! if two-electron integrals are finished, we find a 0 for ir and is:
        if (ir.eq.0) exit

        idxfl = idxfl + 16 ! increment file position by length of entry

        ! the indices need an increment:
        if (ioff_core.ne.0) then
          ip=ip+ioff_core; iq=iq+ioff_core 
          ir=ir+ioff_core; is=is+ioff_core
        end if

        ! of interest for collecting Fock contributions?
        if (csorb(ip).or.csorb(iq).or.csorb(ir).or.csorb(is)) then
c dbg
c          if ((csorb(ip).and.(ip.eq.iq.or.ip.eq.is.or.ip.eq.ir)).or.
c     &        (csorb(ir).and.(ir.eq.is.or.ir.eq.iq)).or.
c     &        (csorb(iq).and.(iq.eq.is)))
c     &         print '(1x,a,4i4)','fock < ',ip,iq,ir,is
c dbg
          irs = (ir-1)*ir/2+is
          ipq = (ip-1)*ip/2+iq

          if (csorb(ip).and.ip.eq.iq)
     &         fock(irs) = fock(irs) + 2d0 * val
c dbg
c          if (csorb(ip).and.ip.eq.iq) print '(20x,a,2i4)','-> J1',ir,is
c dbg
          if (csorb(ip).and.ip.eq.ir) then
c dbg
c            if (iq.ge.is) print '(30x,a,2i4)','-> K1',iq,is
c            if (iq.lt.is) print '(30x,a,2i4)','-> K1',is,iq
c dbg
            if (iq.ge.is)
     &           fock((iq-1)*iq/2+is) = fock((iq-1)*iq/2+is) - val
            if (iq.lt.is)
     &           fock((is-1)*is/2+iq) = fock((is-1)*is/2+iq) - val
          end if
          avoid_dc = ip.eq.ir.and.iq.eq.is .or. ir.eq.is
          if (csorb(ip).and.ip.eq.is.and..not.avoid_dc) then
c dbg
c            if (iq.ge.ir) print '(30x,a,2i4)','-> K2',iq,ir
c            if (iq.lt.ir) print '(30x,a,2i4)','-> K2',ir,iq
c dbg
            if (iq.ge.ir)
     &           fock((iq-1)*iq/2+ir) = fock((iq-1)*iq/2+ir) - val
            if (iq.lt.ir)
     &           fock((ir-1)*ir/2+iq) = fock((ir-1)*ir/2+iq) - val
          end if

          if (csorb(ir).and.ipq.ne.irs.and.ir.eq.is)
     &         fock(ipq) = fock((ip-1)*ip/2+iq) + 2d0 * val
c dbg
c          if (csorb(ir).and.ipq.ne.irs.and.ir.eq.is)
c     &         print '(20x,a,2i4)','-> J2',ip,iq
c dbg

          avoid_dc = ip.eq.iq.and.ir.eq.is .or. ip.eq.iq
          if (csorb(ir).and.ir.eq.iq.and..not.avoid_dc) then
c dbg
c            if (ip.ge.is) print '(30x,a,2i4)','-> K3',ip,is
c            if (ip.lt.is) print '(30x,a,2i4)','-> K3',is,ip
c dbg
            if (ip.ge.is) 
     &           fock((ip-1)*ip/2+is) = fock((ip-1)*ip/2+is) - val
            if (ip.lt.is)
     &           fock((is-1)*is/2+ip) = fock((is-1)*is/2+ip) - val
          end if

          avoid_dc = ip.eq.iq.and.ir.eq.is .or. ir.eq.is
          if (csorb(is).and.is.eq.iq.and..not.avoid_dc) then
c dbg
c            if (ip.ge.ir) print '(30x,a,2i4)','-> K4',ip,ir
c            if (ip.lt.ir) print '(30x,a,2i4)','-> K4',ir,ip
c dbg
            if (ip.ge.ir) 
     &           fock((ip-1)*ip/2+ir) = fock((ip-1)*ip/2+ir) - val
            if (ip.lt.ir)
     &           fock((ir-1)*ir/2+ip) = fock((ir-1)*ir/2+ip) - val
          end if
        end if

        idxprqs(1) = ireost(ip)
        idxprqs(2) = ireost(ir)
        idxprqs(3) = ireost(iq)
        idxprqs(4) = ireost(is)
        igam(1) = igamorb(idxprqs(1))
        igam(2) = igamorb(idxprqs(2))
        igam(3) = igamorb(idxprqs(3))
        igam(4) = igamorb(idxprqs(4))
        idss(1) = igasorb(idxprqs(1))
        idss(2) = igasorb(idxprqs(2))
        idss(3) = igasorb(idxprqs(3))
        idss(4) = igasorb(idxprqs(4))
!! OPEN SHELL: NEEDS ADAPTATION
        igtp(1) = ihpvgas(idss(1),1)
        igtp(2) = ihpvgas(idss(2),1)
        igtp(3) = ihpvgas(idss(3),1)
        igtp(4) = ihpvgas(idss(4),1)
        if (iad_gas(idss(1)).ne.2.or.iad_gas(idss(3)).ne.2) cycle
        if (iad_gas(idss(2)).ne.2.or.iad_gas(idss(4)).ne.2) cycle
        idss(1) = idss(1)-idx_gas(igtp(1))+1
        idss(2) = idss(2)-idx_gas(igtp(2))+1
        idss(3) = idss(3)-idx_gas(igtp(3))+1
        idss(4) = idss(4)-idx_gas(igtp(4))+1

        ! generate string addresses of integrals 
        ! spin-orbital basis (w/o PH-symmetry) to which
        ! current (pq|rs) = <pr|qs> contributes
        call idx42str(nstr,idxstr,
     &           idxprqs,igam,idss,igtp,
     &           orb_info,str_info,hlist,hpvxseq,error)
        ! sign change if (CV and CP) xor (AV and AP)
        if ((igtp(1)*igtp(2).eq.6.and.igtp(3)*igtp(4).ne.6).or.
     &          (igtp(1)*igtp(2).ne.6.and.igtp(3)*igtp(4).eq.6))
     &           idxstr(1:nstr) = -idxstr(1:nstr)

        ! store integral in h2scr
        do istr = 1, nstr
          if (abs(idxstr(istr)).lt.idxst.or.
     &            abs(idxstr(istr)).gt.idxnd) cycle
          if (idxstr(istr).gt.0)
     &             h2scr(ioff+idxstr(istr)) =
     &             h2scr(ioff+idxstr(istr))+val
          if (idxstr(istr).lt.0)
     &             h2scr(ioff-idxstr(istr)) =
     &             h2scr(ioff-idxstr(istr))-val
        end do


        if (ip.eq.iq.or.ir.eq.is) cycle

        ! ip <-> iq.
        idxprqs(1) = ireost(iq)
        idxprqs(3) = ireost(ip)
        igam(1) = igamorb(idxprqs(1))
        igam(3) = igamorb(idxprqs(3))
        idss(1) = igasorb(idxprqs(1))
        idss(3) = igasorb(idxprqs(3))
        igtp(1) = ihpvgas(idss(1),1)
        igtp(3) = ihpvgas(idss(3),1)
        idss(1) = idss(1)-idx_gas(igtp(1))+1
        idss(3) = idss(3)-idx_gas(igtp(3))+1

        ! generate string addresses of integrals 
        ! spin-orbital basis (w/o PH-symmetry) to which
        ! current (pq|rs) = <pr|qs> contributes
        call idx42str(nstr,idxstr,
     &           idxprqs,igam,idss,igtp,
     &           orb_info,str_info,hlist,hpvxseq,error)
            ! sign change if (CV and CP) xor (AV and AP)
        if ((igtp(1)*igtp(2).eq.6.and.igtp(3)*igtp(4).ne.6).or.
     &          (igtp(1)*igtp(2).ne.6.and.igtp(3)*igtp(4).eq.6))
     &           idxstr(1:nstr) = -idxstr(1:nstr)

        ! store integrals in h2scr
        do istr = 1, nstr
c          if (abs(idxstr(istr)).lt.idxst.or.
c     &            abs(idxstr(istr)).gt.idxnd) cycle
          if (idxstr(istr).gt.0)
     &             h2scr(ioff+idxstr(istr)) =
     &             h2scr(ioff+idxstr(istr))+val
          if (idxstr(istr).lt.0)
     &             h2scr(ioff-idxstr(istr)) =
     &             h2scr(ioff-idxstr(istr))-val
        end do
          
      end do 

      ! remember last position
      ifpos = idxfl
          
      ! write reordered integral to disc
      call put_vec(ffham,h2scr,idxst,idxnd)
        
      write(lulog,*) 
      write(lulog,*) '2-el. integrals on disk: ',int_disk
      write(lulog,*) '   thereof nonredundant: ',int_nonr
      write(lulog,*) '    integrals reordered: ',!int_ordr,
     &     hlist%len_op-hlist%off_op_occ(iblkst)

      if (closeit)
     &     call file_close_keep(ffham)

      ! automatic deallocation:
      ifree = mem_flushmark('import_h2')
      deallocate(csorb)

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(lulog,'time in 2int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return

      end

