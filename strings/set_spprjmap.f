*----------------------------------------------------------------------*
      subroutine set_spprjmap(init,grph,ityp,iocc,irestr,
     &                      idxspprj,ms_sc,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     set mapping for sceleton string vs. strings with all spin distr.
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'ifc_memman.h'
      include 'def_graph.h'
      include 'def_filinf.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'

      logical, intent(in) ::
     &     init
      type(graph), intent(in) ::
     &     grph
      integer, intent(in) ::
     &     ityp, iocc, irestr(*), idxspprj, ms_sc
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     nsym, ifree, ms,
     &     lenbuf, idxbuf, idxms, idxms_spprj, ioff_sym,
     &     idxmap, igam, nstr, nmaps, imaps,
     &     maxlen_blk, maxlenbuf
      integer, pointer ::
     &     mostnd(:,:,:), idx_gas(:), ngas_hpv(:), igamorb(:)

      integer, pointer ::
     &     buffer(:)

      if (ntest.ge.100) then
        write(luout,*) '==================='
        write(luout,*) 'this is set_spprjmap'
        write(luout,*) '==================='
        write(luout,*) 'idxspprj:   ',idxspprj
        write(luout,*) 'ms_sc:      ',ms_sc
        write(luout,*) 'irestr:     ',
     &       irestr(1:2*orb_info%ngas_hpv(ityp))
        write(luout,*) 'ityp, iocc: ',ityp, iocc
        write(luout,*) 'init:       ',init
      end if

      ifree = mem_setmark('set_spprjmap')

      nsym = orb_info%nsym
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      idx_gas => orb_info%idx_gas
      ngas_hpv => orb_info%ngas_hpv

c      idxms = (iocc - ms)/2 + 1
c      idxms_spprj = (iocc + ms)/2 + 1

      nmaps = 2**iocc - 1

      maxlen_blk = 0
      maxlenbuf = 0
      idxms = (iocc - ms_sc)/2 + 1
c      do ms = iocc, -iocc, -2
c      do ms = iocc, 0, -2
        lenbuf = 0
        do igam = 1, nsym
          nstr = grph%lenstr_gm(igam,idxms)
          if (nstr.ne.grph%lenstr_gm(igam,idxms)) then
            call quit(1,'set_spprjmap',
     &         'incompatible spin-spprjped string')
          end if
          maxlen_blk = max(maxlen_blk,nstr*nmaps)
          lenbuf = lenbuf+nstr*nmaps
        end do
        maxlenbuf = max(lenbuf,maxlenbuf)
c      end do
      
      lenbuf = maxlenbuf
      if (ntest.ge.100) then
        write(luout,*) 'allocating ',lenbuf,' integer words'
      end if
      ifree = mem_alloc_int(buffer,lenbuf,'buffer')

      if (init) then
        strmap_info%maxlen_blk_spprj(idxspprj) = maxlen_blk
        strmap_info%offsets_spprj(idxspprj)%ms(1:iocc+1) = -1
        strmap_info%offsets_spprj(idxspprj)%msgm(1:(iocc+1)*nsym) = -1
      end if

      idxmap = 0
c      do ms = iocc, -iocc, -2
      strmap_info%offsets_spprj(idxspprj)%ms(idxms)=strmap_info%idx_last
      idxbuf = 0
      do igam = 1, nsym
        nstr = grph%lenstr_gm(igam,idxms)
        strmap_info%offsets_spprj(idxspprj)%msgm(igam)
     &         = idxmap
c dbg
c          print *,'idxspprj,idxms,igam: ',idxspprj,idxms,igam
c          print *,'offset: ',idxmap
c          print *,'ityp,ngas_hpv(ityp)',ityp,ngas_hpv(ityp)
c dbg
        if (nstr.eq.0) cycle
        call set_spprjmap_kernel(buffer(idxbuf+1),
     &                 iocc,irestr,idxms,igam,
     &                 grph%y4sg,grph%yinf,
     &                 grph%yssg,grph%wssg,
     &                 grph%ioffstr_dgm,grph%ndis,
     &                 mostnd(1,1,idx_gas(ityp)),
     &                 nsym,ngas_hpv(ityp),igamorb)
        idxmap = idxmap + nstr*nmaps
        idxbuf = idxbuf + nstr*nmaps
      end do

c dbg
c        print *,'call to mem_iput in set_spprjmap: ',
c     &       strmap_info%idx_last, lenbuf
c dbg
      lenbuf = idxbuf
      call mem_iput(strmap_info%ffstrmap,buffer,
     &       strmap_info%idx_last+1,strmap_info%idx_last+lenbuf)

      strmap_info%idx_last = strmap_info%idx_last+lenbuf

      if (ntest.ge.150) then
        write(luout,*) 'maxlen:',strmap_info%maxlen_blk_spprj(idxspprj)
        write(luout,*) 'the MS offset array: '
        write(luout,'(1x,10i6)')
     &       strmap_info%offsets_spprj(idxspprj)%ms
        write(luout,*) 'the MS/IRREP offset array: '
        write(luout,'(1x,10i6)')
     &       strmap_info%offsets_spprj(idxspprj)%msgm
      end if

      ifree = mem_flushmark('set_spprjmap')

      return
      end
