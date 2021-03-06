*----------------------------------------------------------------------*
      subroutine import_h2_cfour(hlist,str_info,orb_info)
*----------------------------------------------------------------------*
*     read through HF2 file and sort integral into list hlist;
*     if memory is insufficient to keep the complete list,
*     we resort to a batch algorithm
*
*     andreas, march 2016
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

      integer(8), pointer ::
     &     ibuf(:)
      real(8), pointer ::
     &     xbuf(:), h2scr(:)

      type(filinf) ::
     &     ffmo2
      logical ::
     &     closeit, error, found
c     &     first, first_str
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

      character(len=10) ::
     &     string

      type(filinf), pointer ::
     &     ffham

c      type(operator), pointer ::
c     &     hop

      integer, pointer ::
     &     ireost(:), ihpvgas(:,:), igamorb(:),
     &     igasorb(:), idx_gas(:), iad_gas(:), ireost_loc(:)

      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      ffham => hlist%fhand

      call atim_csw(cpu0,sys0,wall0)

      call file_init(ffmo2,trim(orb_info%name_intfile_ext),
     &               ftyp_sq_unf,0)
      call file_open(ffmo2)

      lumo2 = ffmo2%unit
      rewind lumo2

      found = .false.
      do while(.not.found)
        read(lumo2) string
        found=index(string,'**').ne.0
      end do
      
      ifree = mem_setmark('import_h2')

      lbuf = orb_info%bufflen
      if (orb_info%bufflen.ne.orb_info%ibufflen) 
     &       call quit(1,'import_h2_cfour','not yet prep.d')
      ifree = mem_alloc_real(xbuf,lbuf,'mo2_xbuff')
      ifree = mem_alloc_int(ibuf,lbuf,'mo2_ibuff')

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
      if (orb_info%ncore_ext.gt.0) then
        ! we have to modify ireost a bit
        ifree = mem_alloc_int(ireost_loc,orb_info%ntoob,'ireost_loc')
        ireost_loc = 0
        do ii = 1, orb_info%ntoob-orb_info%ncore_ext
          ireost_loc(ii)=ireost(orb_info%ext_fcreo(ii))
        end do
        ireost => ireost_loc
        if (ntest.ge.100) then
          write(lulog,'(x,"modified reorder array:")')
          write(lulog,'(x,5i5,x,5i5)') 
     &       ireost(1:orb_info%ntoob-orb_info%ncore_ext)
        end if
      end if
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

        if (ipass.gt.1) then        
          rewind lumo2

          found = .false.
          do while(.not.found)
            read(lumo2) string
            found=index(string,'**').ne.0
          end do

        end if

        do
          read(lumo2) xbuf(1:lbuf),ibuf(1:lbuf),len
          if (len.eq.0) cycle
          if (len.lt.0) exit
          int_disk = int_disk+len
          do ii = 1, len
            ip = iand(ibuf(ii),imsk16)
            iq = iand(ishft(ibuf(ii),-16),imsk16)
            ir = iand(ishft(ibuf(ii),-32),imsk16)
            is = iand(ishft(ibuf(ii),-48),imsk16)

            !NEEDED?? if (idxrs.gt.idxpq) cycle
            int_nonr = int_nonr+1
            
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
            igtp(1) = ihpvgas(idss(1),1)
            igtp(2) = ihpvgas(idss(2),1)
            igtp(3) = ihpvgas(idss(3),1)
            igtp(4) = ihpvgas(idss(4),1)
            if (iad_gas(idss(2)).ne.2.or.iad_gas(idss(4)).ne.2) cycle
            if (iad_gas(idss(1)).ne.2.or.iad_gas(idss(3)).ne.2) cycle
            idss(1) = idss(1)-idx_gas(igtp(1))+1
            idss(2) = idss(2)-idx_gas(igtp(2))+1
            idss(4) = idss(4)-idx_gas(igtp(4))+1
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

      ! automatic deallocation:
      ifree = mem_flushmark('import_h2')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(lulog,'time in 2int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return

      end

