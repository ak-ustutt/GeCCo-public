*----------------------------------------------------------------------*
      subroutine import_h2_dalton64(hlist,str_info,orb_info)
*----------------------------------------------------------------------*
*     read through MOTWOINT file and sort integral into list hlist;
*     if memory is insufficient to keep the complete list,
*     we resort to a batch algorithm
*
*     patched version for 64 bit DALTON versions
*
*     andreas, january 2007
*
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
      include 'par_dalton.h'

      integer, parameter ::
     &     ntest = 00

      integer(8), parameter ::
     &     imsk16 = 65535,
     &     imsk08 =   255
      
      type(me_list), intent(in) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      ! dalton is i4:
      integer(8), pointer ::
     &     ibuf(:)
      real(8), pointer ::
     &     xbuf(:), h2scr(:)

      type(filinf) ::
     &     ffmo2
      logical ::
     &     closeit, error
c     &     first, first_str
      ! dalton is i4:
      integer(8) ::
     &     idxpq, idxrs, ip, iq, ir, is, len_
      integer(8) ::
     &     lbuf_, itrlevel
      integer ::
     &     lumo2, irat, lbuf,
     &     len, ii, nstr,
     &     int_disk, int_nonr, int_ordr, ioff, istr,
     &     idxst, idxnd, nblk, nblkmax, ifree, luerr, nbuff,
     &     len_op, idum, ipass, iblk, iblkst

      integer ::
     &     idxprqs(4),idss(4),igam(4),igtp(4), 
     &     lexlscr(4,3), idxstr(8)

      type(filinf), pointer ::
     &     ffham

c      type(operator), pointer ::
c     &     hop

      integer, pointer ::
     &     ireost(:), ihpvgas(:,:), igamorb(:),
     &     igasorb(:), idx_gas(:), iad_gas(:)

      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      ffham => hlist%fhand

      if (orb_info%ntoob.gt.255) then
        write(lulog,*) 'number of orbitals: ',orb_info%ntoob
        call quit(0,'import_h2_dalton',
     &       'not yet prepared for >255 orbitals')
      end if

      call atim_csw(cpu0,sys0,wall0)

      call file_init(ffmo2,motwoint,ftyp_sq_unf,0)
      call file_open(ffmo2)

      lumo2 = ffmo2%unit
      rewind lumo2

      read(lumo2)
      read(lumo2) lbuf_, itrlevel
      if (itrlevel.lt.10)
     &     call quit(0,'import_h2_dalton',
     &     'you must set DALTON transformation level to 10!')
      
      ifree = mem_setmark('import_h2')

      lbuf = lbuf_
c dbg
      print *,'lbuf ',lbuf,lbuf_
c dbg
      ifree = mem_alloc_real(xbuf,lbuf,'mo2_xbuff')
c      ifree = mem_alloc_int(ibuf,lbuf,'mo2_ibuff')
      ! i4 must be done by hand:
      irat = zirat()
      ifree = mem_register(lbuf/irat,'mo2_ibuff')
      allocate(ibuf(lbuf))

      if (ffham%unit.le.0) then
        call file_open(ffham)
        closeit = .true.
      else
        closeit = .false.
      end if
      
      nblkmax = ifree/ffham%reclen
      if (nblkmax.le.0) then
        write(lulog,*) 'free memory (words):  ',ifree
        write(lulog,*) 'block length (words): ',ffham%reclen
        call quit(1,'get_h2','not even 1 record fits into memory?')
      end if

cmh   determine number of first 2-el. block
      do iblk = hlist%op%n_occ_cls,1,-1
        if (hlist%op%ica_occ(1,iblk).eq.2) iblkst = iblk
      end do
      len_op = hlist%len_op
      nblk = min((len_op-hlist%off_op_occ(iblkst)-1)
     &                      / ffham%reclen+1,nblkmax)

      nbuff = min(len_op-hlist%off_op_occ(iblkst),nblk*ffham%reclen)

      write(lulog,*) 'number of incore-blocks in geth2: ',nblk
      write(lulog,'(x,a,f9.2,a)') 'size of buffer in geth2:   ',
     &     dble(nbuff)/128d0/1024d0, 'Mb'

      ifree = mem_alloc_real(h2scr,nbuff,'h2sort_buff')

      int_disk = 0
      int_nonr = 0
      int_ordr = 0

      ! dereference structure components for efficiency
      ireost => orb_info%ireost
      ihpvgas => orb_info%ihpvgas
      igamorb => orb_info%igamorb
      igasorb => orb_info%igasorb
      idx_gas => orb_info%idx_gas
      iad_gas => orb_info%iad_gas

      ! loop over batches of final integral file
      ipass = 0
      idxst = hlist%off_op_occ(iblkst) + 1
      do while(idxst.le.len_op)
        ipass = ipass+1
        idxnd = min(len_op,idxst-1+nbuff)
        ioff = -idxst+1
        h2scr(1:nbuff) = 0d0
        
        rewind lumo2
        luerr = lulog
        call mollab(motwolab,lumo2,luerr)
        do
          read(lumo2) xbuf(1:lbuf),ibuf(1:lbuf),len_
          len = len_
          int_disk = int_disk+len
          if (len.eq.0) cycle
          if (len.lt.0) exit
          ! get index (pq|rs)
          ! rs is always the same within record
          idxrs = iand(ishft(ibuf(1),-16),imsk16)
          ir = iand(ishft(idxrs,-8),imsk08)
          is = iand(idxrs,imsk08)
          idxprqs(2) = ireost(ir)
          idxprqs(4) = ireost(is)
          igam(2) = igamorb(idxprqs(2))
          igam(4) = igamorb(idxprqs(4))
          idss(2) = igasorb(idxprqs(2))
          idss(4) = igasorb(idxprqs(4))
!! OPEN SHELL: NEEDS ADAPTATION
          igtp(2) = ihpvgas(idss(2),1)
          igtp(4) = ihpvgas(idss(4),1)
          if (iad_gas(idss(2)).ne.2.or.iad_gas(idss(4)).ne.2) cycle
          idss(2) = idss(2)-idx_gas(igtp(2))+1
          idss(4) = idss(4)-idx_gas(igtp(4))+1
          do ii = 1, len
            idxpq = iand(ibuf(ii),imsk16)
            if (idxrs.gt.idxpq) cycle
            int_nonr = int_nonr+1
            ip = iand(ishft(idxpq,-8),imsk08)
            iq = iand(ibuf(ii),imsk08)
            
            idxprqs(1) = ireost(ip)
            idxprqs(3) = ireost(iq)
            igam(1) = igamorb(idxprqs(1))
            igam(3) = igamorb(idxprqs(3))
            idss(1) = igasorb(idxprqs(1))
            idss(3) = igasorb(idxprqs(3))
            igtp(1) = ihpvgas(idss(1),1)
            igtp(3) = ihpvgas(idss(3),1)
            if (iad_gas(idss(1)).ne.2.or.iad_gas(idss(3)).ne.2) cycle
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

            ! store integral in h2scr
            do istr = 1, nstr
              if (abs(idxstr(istr)).lt.idxst.or.
     &            abs(idxstr(istr)).gt.idxnd) cycle
              if (idxstr(istr).gt.0)
     &             h2scr(ioff+idxstr(istr)) =
     &             h2scr(ioff+idxstr(istr))+xbuf(ii)
              if (idxstr(istr).lt.0)
     &             h2scr(ioff-idxstr(istr)) =
     &             h2scr(ioff-idxstr(istr))-xbuf(ii)
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
              if (abs(idxstr(istr)).lt.idxst.or.
     &            abs(idxstr(istr)).gt.idxnd) cycle
              if (idxstr(istr).gt.0)
     &             h2scr(ioff+idxstr(istr)) =
     &             h2scr(ioff+idxstr(istr))+xbuf(ii)
              if (idxstr(istr).lt.0)
     &             h2scr(ioff-idxstr(istr)) =
     &             h2scr(ioff-idxstr(istr))-xbuf(ii)
            end do
          
          end do ! integrals in xbuf
          
        end do ! pass over DALTON integral file

        ! write reordered integral to disc
        call put_vec(ffham,h2scr,idxst,idxnd)
        idxst = idxnd+1
        
      end do ! pass over reordered integral file

      write(lulog,*) 'passes over integral file: ',ipass
      write(lulog,*) 
      write(lulog,*) '2-el. integrals on disk: ',int_disk
      write(lulog,*) '   thereof nonredundant: ',int_nonr
      write(lulog,*) '    integrals reordered: ',!int_ordr,
     &     hlist%len_op-hlist%off_op_occ(iblkst)

      if (closeit)
     &     call file_close_keep(ffham)
      call file_close_keep(ffmo2)

      ! deallocation by hand
      deallocate(ibuf)
      ifree = mem_dealloc('mo2_ibuff')
      ! automatic deallocation of rest:
      ifree = mem_flushmark('import_h2')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(lulog,'time in 2int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return

      end

