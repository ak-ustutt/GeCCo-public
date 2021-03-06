*----------------------------------------------------------------------*
      subroutine set_strmap(graph1,ityp1,iocc1,irestr1,
     &                      graph2,ityp2,iocc2,irestr2,
     &                      graph12,ityp12,iocc12,irestr12,
     &                      idxgrgr,strmap_info,orb_info)
*----------------------------------------------------------------------*
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
     &     graph1, graph2, graph12
      integer, intent(in) ::
     &     ityp1, iocc1, irestr1(*),
     &     ityp2, iocc2, irestr2(*),
     &     ityp12, iocc12, irestr12(*), idxgrgr
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     nsym, ifree, ipass,
     &     lenbuf, idxbuf, idxms1, idxms2, ioff_symtab,
     &     idxmap, igam2, igam1, nstr2, nstr1,
     &     maxlen_blk
      integer, pointer ::
     &     mostnd(:,:,:), idx_gas(:), ngas_hpv(:), igamorb(:)

      integer, pointer ::
     &     buffer(:)

      if (ntest.ge.100) then
        write(lulog,*) '=================='
        write(lulog,*) 'this is set_strmap'
        write(lulog,*) '=================='
        write(lulog,*) ' iocc1, iocc2, iocc12: ',
     &                   iocc1, iocc2, iocc12
        write(lulog,*) ' ityp1, ityp2, ityp12: ',
     &                   ityp1, ityp2, ityp12
        write(lulog,*) ' restriction on 1:'
        call wrt_srstr(lulog,irestr1,orb_info%ngas_hpv(ityp1),
     &                 orb_info%nspin)
        write(lulog,*) ' restriction on 2:'
        call wrt_srstr(lulog,irestr2,orb_info%ngas_hpv(ityp2),
     &                 orb_info%nspin)
        write(lulog,*) ' restriction on 12:'
        call wrt_srstr(lulog,irestr12,orb_info%ngas_hpv(ityp12),
     &                 orb_info%nspin)
      end if

      ! message on statistics file
      if (lustat.gt.0) 
     &write(lustat,'(">new map in space ",i1,": ",i2,",",i2," -> ",i2)')
     &       ityp1, iocc1, iocc2, iocc12

      ifree = mem_setmark('set_strmap')

      nsym = orb_info%nsym
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      idx_gas => orb_info%idx_gas
      ngas_hpv => orb_info%ngas_hpv

      ! consistency checks
      if (ityp1.ne.ityp2.or.ityp1.ne.ityp12.or.ityp2.ne.ityp12) then
        write(lulog,*) 'look at that: ',ityp1,ityp2,ityp12
        call quit(1,'set_strmap','incompatible graph types')
      end if
      if (iocc1+iocc2.ne.iocc12) then
        write(lulog,*) 'look at that: ',iocc1,' +',iocc2,' !=',iocc12
        call quit(1,'set_strmap','incompatible graph occupations')
      end if

      maxlen_blk = 0
      do ipass = 1, 2

        if (ipass.eq.2) then
          lenbuf = idxbuf
          if (ntest.ge.100) then
            write(lulog,*) 'allocating ',lenbuf,' integer words'
          end if
          ifree = mem_alloc_int(buffer,lenbuf,'buffer')
          strmap_info%maxlen_blk(idxgrgr) = maxlen_blk
          strmap_info%offsets(idxgrgr)%msms(1:(iocc1+1)*(iocc2+1)) = -1
        end if

        idxbuf = 0

        do idxms2 = 1, iocc2+1
          do idxms1 = 1, iocc1+1
            if (ipass.eq.2) then
              strmap_info%offsets(idxgrgr)%
     &             msms((idxms2-1)*(iocc1+1)+idxms1)
     &             = idxbuf
            end if
            idxmap = 0
            do igam2 = 1, nsym
              nstr2 = graph2%lenstr_gm(igam2,idxms2)
              do igam1 = 1, nsym
                nstr1 = graph1%lenstr_gm(igam1,idxms1)
                if (ntest.ge.100.and.ipass.eq.1)
     &            write(lulog,'(a,2i14,4i4)') ' nstr1,nstr2: ',
     &                    nstr1,nstr2,
     &                    igam1,igam2,idxms1,idxms2
                maxlen_blk = max(maxlen_blk,nstr1*nstr2)
                if (ipass.eq.2) then
                  strmap_info%offsets(idxgrgr)%
     &                 msmsgmgm(((idxms2-1)*(iocc1+1)+idxms1-1)
     &                 *nsym*nsym+(igam2-1)*nsym+igam1) = idxmap
                  if (nstr1.eq.0.or.nstr2.eq.0) cycle
                  call set_strmap_kernel(buffer(idxbuf+idxmap+1),
     &                 iocc1,irestr1,idxms1,igam1,
     &                 iocc2,irestr2,idxms2,igam2,
     &                 irestr12,
     &                 graph12%y4sg,graph12%yinf,
     &                 graph12%yssg,graph12%wssg,
     &                 graph12%ioffstr_dgm,graph12%ndis,
     &                 mostnd(1,1,idx_gas(ityp12)),
     &                 nsym,ngas_hpv(ityp12),igamorb)
                end if
                idxmap = idxmap + nstr1*nstr2
              end do
            end do
            idxbuf = idxbuf+idxmap
          end do
        end do
      end do

      if (ntest.ge.150) then
        write(lulog,*) 'the ms offset array'
        call wrtimat2(strmap_info%offsets(idxgrgr)%msms,
     &       iocc2+1,iocc1+1,iocc2+1,iocc1+1)
        write(lulog,*) 'the IRREP offset arrays:'
        do idxms2 = 1, iocc2+1
          do idxms1 = 1, iocc1+1
            write(lulog,*) 'ms = ',iocc2-(idxms2-1)*2,iocc1-(idxms1-1)*2
            ioff_symtab = buffer((idxms2-1)*iocc1+idxms1) +
     &           (iocc1+1)*(iocc2+1)
            call wrtimat2(strmap_info%offsets(idxgrgr)%
     &           msmsgmgm(((idxms2-1)*(iocc1+1)+idxms1-1)*nsym*nsym+1),
     &           nsym,nsym,nsym,nsym)
            if (ntest.ge.1000) then
              write(lulog,*) 'the maps'
              do igam2 = 1, nsym
                nstr2 = graph2%lenstr_gm(igam2,idxms2)
                do igam1 = 1, nsym
                  nstr1 = graph1%lenstr_gm(igam1,idxms1)
                  idxmap = strmap_info%offsets(idxgrgr)%
     &                 msms((idxms2-1)*(iocc1+1)+idxms1) +
     &                 strmap_info%offsets(idxgrgr)%
     &                 msmsgmgm(((idxms2-1)*(iocc1+1)
     &                 +idxms1-1)*nsym*nsym
     &                 +(igam2-1)*nsym+igam1)
                  write(lulog,*) 'size: ',nstr2,' x',nstr1
                  if (ntest.le.1000.and.nstr1*nstr2.gt.10000) then
                    write(lulog,*) 'too large for printing'
                    cycle
                  end if
                  call wrtimat2(buffer(idxmap+1),
     &                 nstr2,nstr1,nstr2,nstr1)
                end do
              end do
            end if
          end do
        end do
      end if

c dbg
c        print *,'call to mem_iput in set_strmap: ',
c     &       strmap_info%idx_last, lenbuf
c dbg
      call mem_iput(strmap_info%ffstrmap,buffer,
     &     strmap_info%idx_last+1,strmap_info%idx_last+lenbuf)

      strmap_info%idx_last = strmap_info%idx_last+lenbuf

      ifree = mem_flushmark()

      return
      end
