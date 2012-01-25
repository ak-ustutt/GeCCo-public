*----------------------------------------------------------------------*
      subroutine import_h2_gamess(hlist,str_info,orb_info)
*----------------------------------------------------------------------*
*     read through MOINTS file and sort integral into list hlist;
*     if memory is insufficient to keep the complete list,
*     we resort to a batch algorithm
*
*     andreas, january 2007
*     matthias, spring 2010
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
      include 'par_gamess.h'

      integer, parameter ::
     &     ntest = 00

      integer(4), parameter ::
     &     imsk16 = 65535,
     &     imsk08 =   255
      
      type(me_list), intent(in) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      ! gamess is i8:
      integer(8), pointer ::
     &     ibuf(:)
      real(8), pointer ::
     &     xbuf(:), h2scr(:)

      type(filinf) ::
     &     ffmo2
      logical ::
     &     closeit, error, pack2
c     &     first, first_str
      ! gamess is i8:
      integer(8) ::
     &     idxpq, idxrs, ip, iq, ir, is, lbuf_, len_, itrlevel
      integer ::
     &     nlabmx, lumo2, irat, lbuf,
     &     len, ii, nstr, ilabel,
     &     int_disk, int_nonr, int_ordr, ioff, istr,
     &     idxst, idxnd, nblk, nblkmax, ifree, luerr, nbuff,
     &     len_op, idum, ipass, iblk, iblkst
      character*6 ::
     &     cmxao

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
     &     cpu0, sys0, wall0, cpu, sys, wall, energy

      ffham => hlist%fhand

      call atim_csw(cpu0,sys0,wall0)

      call file_init(ffmo2,moints,ftyp_sq_unf,0)
      call file_open(ffmo2)

      lumo2 = ffmo2%unit

      ! Dirty fix for determining whether LABSIZ was 1 or 2 in gamess:
      ! First let us assume that # of MOs was equal to # of AOs:
      pack2 = orb_info%ntoob.le.mxao
      if (pack2) then
        nlabmx = (nintmx+1)/2
      else
        nlabmx = nintmx
      end if
      ! Now we need to check the first orbital index: If zero: LABSIZ=2
      if (pack2) then
        rewind lumo2
        read (lumo2) energy
        read(lumo2) len_,idxpq
        if (len_.gt.0) then
          ip = ishft(idxpq,-56)
          if (ip.eq.0) then
            pack2 = .false.
            nlabmx = nintmx
            write(cmxao,'(i6)') mxao
            call warn('import_h2_gamess',
     &          'We assume that the number of AOs'//
     &          ' in the GAMESS run exceeded '//trim(cmxao))
          end if
        end if
      end if

      ifree = mem_setmark('import_h2')

      ifree = mem_alloc_real(xbuf,nintmx,'mo2_xbuff')
      ifree = mem_alloc_int(ibuf,nlabmx,'mo2_ibuff')

      if (ffham%unit.le.0) then
        call file_open(ffham)
        closeit = .true.
      else
        closeit = .false.
      end if
      
      nblkmax = ifree/ffham%reclen
      if (nblkmax.le.0) then
        write(luout,*) 'free memory (words):  ',ifree
        write(luout,*) 'block length (words): ',ffham%reclen
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

      write(luout,*) 'number of incore-blocks in geth2: ',nblk
      write(luout,'(x,a,f9.2,a)') 'size of buffer in geth2:   ',
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
        read (lumo2) energy
        do
c dbg
c          print *,'energy: ',energy
c          print *,'nlabmx,nintmx: ',nlabmx,nintmx
c dbgend
          read(lumo2) len_,ibuf(1:nlabmx),xbuf(1:nintmx)
c dbg
c          print *,'read ',len_,' integrals:'
c          do ii = 1, abs(len_)
c            if (.not.pack2) then
c              istr = ii
c            else
c              istr = (ii+1)/2
c            end if 
c            print *,ibuf(istr),xbuf(ii)
c          end do
c dbgend
          len = len_
          int_disk = int_disk+abs(len)
          if (len.eq.0) exit
          ! get index (pq|rs)

          do ii = 1, abs(len)
            if (pack2) then
              if (mod(ii,2).eq.0) then
                ilabel = ibuf(ii/2)
                ip = iand(ishft(ilabel,-24),imsk08)
                iq = iand(ishft(ilabel,-16),imsk08)
                ir = iand(ishft(ilabel, -8),imsk08)
                is = iand(ilabel,imsk08)
              else
                ilabel = ibuf(ii/2+1)
                ip = ishft(ilabel,-56)
                iq = iand(ishft(ilabel,-48),imsk08)
                ir = iand(ishft(ilabel,-40),imsk08)
                is = iand(ishft(ilabel,-32),imsk08)
              end if
            else
              ilabel = ibuf(ii)
              ip = ishft(ilabel,-48)
              iq = iand(ishft(ilabel,-32),imsk16)
              ir = iand(ishft(ilabel,-16),imsk16)
              is = iand(ilabel,imsk16)
            end if
c dbg
c            write(luout,'(a,4i3)') 'pqrs: ',ip,iq,ir,is
c dbgend
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
          if (len.lt.0) exit
          
        end do ! pass over DALTON integral file

        ! write reordered integral to disc
        call put_vec(ffham,h2scr,idxst,idxnd)
        idxst = idxnd+1
        
      end do ! pass over reordered integral file

      write(luout,*) 'passes over integral file: ',ipass
      write(luout,*) 
      write(luout,*) '2-el. integrals on disk: ',int_disk
      write(luout,*) '   thereof nonredundant: ',int_nonr
      write(luout,*) '    integrals reordered: ',!int_ordr,
     &     hlist%len_op-hlist%off_op_occ(iblkst)

      if (closeit)
     &     call file_close_keep(ffham)
      call file_close_keep(ffmo2)

      ! automatic deallocation:
      ifree = mem_flushmark('import_h2')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(luout,'time in 2int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return

      end

