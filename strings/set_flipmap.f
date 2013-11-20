*----------------------------------------------------------------------*
      subroutine set_flipmap(grph,ityp,iocc,irestr,
     &                      idxflip,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     set mapping for string vs. string with alpha/beta flipped
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'ifc_memman.h'
      include 'def_graph.h'
      include 'def_filinf.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'

      type(graph), intent(in) ::
     &     grph
      integer, intent(in) ::
     &     ityp, iocc, irestr(*), idxflip
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     nsym, ifree, ms,
     &     lenbuf, idxbuf, idxms, idxms_flip, ioff_sym,
     &     idxmap, igam, nstr,
     &     maxlen_blk, maxlenbuf
      integer, pointer ::
     &     mostnd(:,:,:), idx_gas(:), ngas_hpv(:), igamorb(:)

      integer, pointer ::
     &     buffer(:)

      if (ntest.ge.100) then
        write(lulog,*) '==================='
        write(lulog,*) 'this is set_flipmap'
        write(lulog,*) '==================='
        write(lulog,*) 'idxflip: ',idxflip
        write(lulog,*) 'irestr:  ',
     &       irestr(1:2*orb_info%nspin*orb_info%ngas)
        write(lulog,*) 'ityp, iocc: ',ityp, iocc
      end if

      ifree = mem_setmark('set_flipmap')

      nsym = orb_info%nsym
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      idx_gas => orb_info%idx_gas
      ngas_hpv => orb_info%ngas_hpv

      idxms = (iocc - ms)/2 + 1
      idxms_flip = (iocc + ms)/2 + 1

      maxlen_blk = 0
      maxlenbuf = 0
      idxms = 0
      do ms = iocc, -iocc, -2
c      do ms = iocc, 0, -2
        idxms = idxms+1
        lenbuf = 0
        do igam = 1, nsym
          nstr = grph%lenstr_gm(igam,idxms)
          if (nstr.ne.grph%lenstr_gm(igam,idxms)) then
            call quit(1,'set_flipmap',
     &         'incompatible spin-flipped string')
          end if
          maxlen_blk = max(maxlen_blk,nstr)
          lenbuf = lenbuf+nstr
        end do
        maxlenbuf = max(lenbuf,maxlenbuf)
      end do
      
      lenbuf = maxlenbuf
      if (ntest.ge.100) then
        write(lulog,*) 'allocating ',lenbuf,' integer words'
      end if
      ifree = mem_alloc_int(buffer,lenbuf,'buffer')
      strmap_info%maxlen_blk_flip(idxflip) = maxlen_blk
      strmap_info%offsets_flip(idxflip)%ms(1:iocc+1) = -1
      strmap_info%offsets_flip(idxflip)%msgm(1:(iocc+1)*nsym) = -1

      idxmap = 0
      idxms = 0
c      do ms = iocc, 0, -2
      do ms = iocc, -iocc, -2
        idxms = idxms+1
        strmap_info%offsets_flip(idxflip)%ms(idxms) = idxmap
        idxbuf = 0
        do igam = 1, nsym
          nstr = grph%lenstr_gm(igam,idxms)
          strmap_info%offsets_flip(idxflip)%msgm((idxms-1)*nsym+igam)
     &         = idxmap
c dbg
c          print *,'idxflip,idxms,igam: ',idxflip,idxms,igam
c          print *,'offset: ',idxmap
c          print *,'ityp,ngas_hpv(ityp)',ityp,ngas_hpv(ityp)
c dbg
          if (nstr.eq.0) cycle
          call set_flipmap_kernel(buffer(idxbuf+1),
     &                 iocc,irestr,idxms,igam,
     &                 grph%y4sg,grph%yinf,
     &                 grph%yssg,grph%wssg,
     &                 grph%ioffstr_dgm,grph%ndis,
     &                 mostnd(1,1,idx_gas(ityp)),
     &                 nsym,ngas_hpv(ityp),igamorb)
          idxmap = idxmap + nstr
          idxbuf = idxbuf + nstr
        end do

c dbg
c        print *,'call to mem_iput in set_flipmap: ',
c     &       strmap_info%idx_last, lenbuf
c dbg
        lenbuf = idxbuf
        call mem_iput(strmap_info%ffstrmap,buffer,
     &       strmap_info%idx_last+1,strmap_info%idx_last+lenbuf)

        strmap_info%idx_last = strmap_info%idx_last+lenbuf

      end do

      if (ntest.ge.150) then
        write(lulog,*) 'maxlen:',strmap_info%maxlen_blk_flip(idxflip)
        write(lulog,*) 'the MS offset array: '
        write(lulog,'(1x,10i6)')
     &       strmap_info%offsets_flip(idxflip)%ms
        write(lulog,*) 'the MS/IRREP offset array: '
        write(lulog,'(1x,10i6)')
     &       strmap_info%offsets_flip(idxflip)%msgm
      end if

      ifree = mem_flushmark('set_flipmap')

      return
      end
