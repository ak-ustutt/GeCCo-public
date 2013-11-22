*----------------------------------------------------------------------*
      subroutine set_fcmap(grph1,ityp1,iocc1,irestr1,
     &                     grph2,ityp2,iocc2,irestr2,
     &                     idxfc,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     set mapping for strings with different restrictions
*     I called this fc (frozen core) map, as this is the present case
*     but it will play a general role when dealing with restricted
*     strings
*
*     andreas, nov 2009
*
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
     &     grph1, grph2
      integer, intent(in) ::
     &     ityp1, iocc1, irestr1(*),
     &     ityp2, iocc2, irestr2(*),
     &     idxfc
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     nsym, ifree, ms,
     &     lenbuf, idxbuf, idxms, ioff_sym,
     &     idxmap, igam, nstr1,
     &     maxlen_blk, maxlenbuf
      integer, pointer ::
     &     mostnd(:,:,:), idx_gas(:), ngas_hpv(:), igamorb(:)

      integer, pointer ::
     &     buffer(:)

      if (ntest.ge.100) then
        write(lulog,*) '==================='
        write(lulog,*) 'this is set_fcmap'
        write(lulog,*) '==================='
        write(lulog,*) 'idxfc: ',idxfc
        write(lulog,*) 'ityp1, iocc1: ',ityp1, iocc1
        write(lulog,*) 'irestr1:  ',
     &       irestr1(1:2*orb_info%nspin*orb_info%ngas)
        write(lulog,*) 'ityp2, iocc2: ',ityp2, iocc2
        write(lulog,*) 'irestr2:  ',
     &       irestr2(1:2*orb_info%nspin*orb_info%ngas)
      end if

      ifree = mem_setmark('set_fcmap')

      if (ityp1.ne.ityp2)
     &     call quit(1,'set_fcmap','types are not equal')
      if (iocc1.ne.iocc2)
     &     call quit(1,'set_fcmap','occupations are not equal')

      nsym = orb_info%nsym
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      idx_gas => orb_info%idx_gas
      ngas_hpv => orb_info%ngas_hpv

      idxms = (iocc1 - ms)/2 + 1

      maxlen_blk = 0
      maxlenbuf = 0
      idxms = 0
      do ms = iocc1, -iocc1, -2
        idxms = idxms+1
        lenbuf = 0
        do igam = 1, nsym
          nstr1 = grph1%lenstr_gm(igam,idxms)
          ! we map idx_str2 = map(idx_str1) 
          maxlen_blk = max(maxlen_blk,nstr1) 
          lenbuf = lenbuf+nstr1
        end do
        maxlenbuf = max(lenbuf,maxlenbuf)
      end do
      
      lenbuf = maxlenbuf
      if (ntest.ge.100) then
        write(lulog,*) 'allocating ',lenbuf,' integer words'
      end if
      ifree = mem_alloc_int(buffer,lenbuf,'buffer')
      strmap_info%maxlen_blk_fc(idxfc) = maxlen_blk
      strmap_info%offsets_fc(idxfc)%ms(1:iocc1+1) = -1
      strmap_info%offsets_fc(idxfc)%msgm(1:(iocc1+1)*nsym) = -1

      idxmap = 0
      idxms = 0
      do ms = iocc1, -iocc1, -2
        idxms = idxms+1
        strmap_info%offsets_fc(idxfc)%ms(idxms) = idxmap
        idxbuf = 0
        do igam = 1, nsym
          nstr1 = grph1%lenstr_gm(igam,idxms)
          strmap_info%offsets_fc(idxfc)%msgm((idxms-1)*nsym+igam)
     &         = idxmap
c dbg
c          print *,'idxfc,idxms,igam: ',idxfc,idxms,igam
c          print *,'offset: ',idxmap
c          print *,'ityp,ngas_hpv(ityp)',ityp,ngas_hpv(ityp)
c dbg
          if (nstr1.eq.0) cycle
          call set_fcmap_kernel(buffer(idxbuf+1),
     &                 iocc1,irestr1,irestr2,
     &                 idxms,igam,
     &                 grph2%y4sg,grph2%yinf,
     &                 grph2%yssg,grph2%wssg,
     &                 grph2%ioffstr_dgm,grph2%ndis,
     &                 mostnd(1,1,idx_gas(ityp1)),
     &                 nsym,ngas_hpv(ityp1),igamorb)
          idxmap = idxmap + nstr1
          idxbuf = idxbuf + nstr1
        end do

        lenbuf = idxbuf
c dbg
c        print *,'call to mem_iput in set_fcmap: ',
c     &       strmap_info%idx_last, lenbuf
c dbg
        call mem_iput(strmap_info%ffstrmap,buffer,
     &       strmap_info%idx_last+1,strmap_info%idx_last+lenbuf)

        strmap_info%idx_last = strmap_info%idx_last+lenbuf

      end do

      if (ntest.ge.150) then
        write(lulog,*) 'maxlen:',strmap_info%maxlen_blk_fc(idxfc)
        write(lulog,*) 'the MS offset array: '
        write(lulog,'(1x,10i6)')
     &       strmap_info%offsets_fc(idxfc)%ms
        write(lulog,*) 'the MS/IRREP offset array: '
        write(lulog,'(1x,10i6)')
     &       strmap_info%offsets_fc(idxfc)%msgm
      end if

      ifree = mem_flushmark('set_fcmap')

      return
      end
