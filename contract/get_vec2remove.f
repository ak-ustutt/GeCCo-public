      subroutine get_vec2remove(pvecs,xmat,pdim,xdim,ms,igam,drank,
     &                          na,nc,gno,pass_spc,mel_spc,na1mx,nc1mx,
     &                          orb_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     given a (simple and purely active) high-rank operator (by na,nc),
*     a difference of ranks to a low-rank operator (drank) and
*     the transformation matrix xmat for orthogonalization of the set
*     of low-rank operators,
*     we determine the coefficients (pvecs) needed for expressing
*     the orthogonalized set of low-rank operators in terms of the
*     original high-rank operators (for a certain ms/igam block)
*
*     matthias, jun 2011
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'
      include 'hpvxseq.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     pdim,xdim,ms,igam,drank,na,nc,gno,na1mx,nc1mx
      real(8), intent(out) ::
     &     pvecs(pdim,xdim)
      real(8), intent(in) ::
     &     xmat(xdim,xdim)
      logical, intent(in) ::
     &     pass_spc
      type(me_list), intent(in) ::
     &     mel_spc

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info

      real(8) ::
     &     facinv, invmat(xdim,xdim), xnew(xdim,xdim)

      integer, target ::
     &     occ1(ngastp,2), occ2(ngastp,2), occ12(ngastp,2),
     &     idxg1(ngastp,2), idxg2(ngastp,2), idxg12(ngastp,2)
      integer ::
     &     map_info_c12(4), map_info_a12(4),
     &     occ_c1, occ_a1, graph_c1, graph_a1,
     &     hpvx_c1, hpvx_a1, idxmsdis_c1, idxmsdis_a1,
     &     occ_c12, occ_a12, graph_c12, graph_a12,
     &     hpvx_c12, hpvx_a12, idxmsdis_c12, idxmsdis_a12,
     &     len_str12(2), len_str1(2),
     &     ica, igraph, iop, maxbuf, ifree, ngam, nact,
     &     msa12, msc12, msa2, msc2, gama12, gamc12, gama2, gamc2,
     &     lenca12, lenca2, ioff12, ioff2, na2, nc2, nablk2, ncblk2,
     &     msc1, msa1, gama1, gamc1, nstra2, nstra12, istra2, istra12,
     &     nstrc1, nstrc2, nstrc12, istrc1, istrc2, istrc12,
     &     idx2, idx12, isgnc, isgn, ielmap, nfac, ii, nbuf, i_cls,
     &     idxms, ioff, jj

      integer, pointer ::
     &     pgraph, pocc,
     &     map_c12(:), map_a12(:),
     &     occ_c2(:), occ_a2(:), graph_c2(:), graph_a2(:),
     &     hpvx_c2(:), hpvx_a2(:), idxmsdis_c2(:), idxmsdis_a2(:),
     &     len_str2(:)

      real(8), pointer ::
     &     buffer(:)

      integer, external ::
     &     ielprd, ifac

      ifree = mem_setmark('vec2rem')
      ngam = orb_info%nsym
      nact = orb_info%nactel

      if (ntest.ge.100) then
        write(lulog,'(x,a)') '----------------------'
        write(lulog,'(x,a)') 'get_vec2remove at work'
        write(lulog,'(x,a)') '----------------------'
        write(lulog,*) 'input low-rank transformation matrix:'
        call wrtmat2(xmat,xdim,xdim,xdim,xdim)
      end if

      ! we assume the following operators:
      !  Op1: /0 0 drank 0\  Op2: /0 0 nc-drank 0\  Op12: /0 0 nc 0\
      !       \0 0 drank 0/       \0 0 na-drank 0/        \0 0 na 0/
      occ1 = 0
      occ12 = 0
      occ1(IVALE,1) = drank
      occ1(IVALE,2) = drank
      occ12(IVALE,1) = nc
      occ12(IVALE,2) = na
      occ2 = occ12 - occ1
      if (any(occ2(IVALE,1:2).lt.0)) call quit(1,'get_vec2remove',
     &           'partitioning not allowed')
      nc2 = occ2(IVALE,1)
      na2 = occ2(IVALE,2)
      if (nc2.gt.0) then
        ncblk2 = 1
      else
        ncblk2 = 0
      end if
      if (na2.gt.0) then
        nablk2 = 1
      else
        nablk2 = 0
      end if

      nfac = 1
      if (gno.eq.0) then
        ! determine factor Nfac = (Nact-na+drank)!/((Nact-na)!*drank!)
        ! only for particle-hole normal ordering
        do ii = nact-na+drank, nact-na+1, -1
          nfac = nfac*ii
        end do
        nfac = nfac/ifac(drank)
        xnew = xmat ! use plain trafo matrix
      else
        ! GNO: premultiply trafo matrix
        if (.not.pass_spc) call quit(1,'get_vec2remove',
     &       'seq. orth. with GNO requires special matrix')
        ! first read special matrix into buffer
        nbuf = 0
        do i_cls = 1, mel_spc%op%n_occ_cls
          nbuf = nbuf + mel_spc%len_op_occ(i_cls)
        end do
        ifree = mem_setmark('premultiplication')
        ifree = mem_alloc_real(buffer,nbuf,'buffer')
        call get_vec(mel_spc%fhand,buffer,1,nbuf) ! YAA: offset for records?
        ! which block do we need?
        if (na1mx.eq.1.and.nc1mx.eq.2) then
          i_cls = 1
        else if (na1mx.eq.2.and.nc1mx.eq.1) then
          i_cls = 2
        else if (na1mx.eq.1.and.nc1mx.eq.1) then
          ! do nothing
          xnew = xmat
          i_cls = 0
        else
          call quit(1,'get_vec2remove',
     &              'GNO with seq.orth. not prepared for this case yet')
        end if
        if (i_cls.gt.mel_spc%op%n_occ_cls)
     &     call quit(1,'get_vec2remove','mel_spc has too few blocks')
        if (i_cls.gt.0) then
          ! read in relevant block. we implicitly assume matrix symmetry
          if (ms.eq.1) then
            idxms = 1
          else if (ms.eq.-1) then
            idxms = 2
          else
            call quit(1,'get_vec2remove','not prepared for this ms')
          end if
          if (int(sqrt(dble(mel_spc%len_op_gmo(i_cls)%gam_ms(
     &                      igam,idxms)))).ne.xdim)
     &       call quit(1,'get_vec2remove','dimensions do not match')
          ioff = mel_spc%off_op_gmo(i_cls)%gam_ms(igam,idxms)
          do ii = 1, xdim
            do jj = 1, xdim
              invmat(jj,ii) = buffer(ioff+(ii-1)*xdim+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(lulog,*) 'read-in matrix for premultiplication:'
            call wrtmat2(invmat,xdim,xdim,xdim,xdim)
          end if
          ! premultiply the trafo matrix
          call dgemm('n','n',xdim,xdim,xdim,
     &               1d0,invmat,xdim,
     &               xmat,xdim,
     &               0d0,xnew,xdim)
        end if
        ifree = mem_flushmark()
        if (ntest.ge.100) then
          write(lulog,*) '(modified) transformation matrix:'
          call wrtmat2(xnew,xdim,xdim,xdim,xdim)
        end if
      end if
      facinv = 1d0/dble(nfac)

c dbg
c      print *,'nfac:',nfac
c      print *,'ncblk2,nablk2:',ncblk2,nablk2
c dbgend
      allocate(occ_c2(ncblk2),occ_a2(nablk2),
     &         graph_c2(ncblk2),graph_a2(nablk2),
     &         hpvx_c2(ncblk2),hpvx_a2(nablk2),
     &         idxmsdis_c2(ncblk2),idxmsdis_a2(nablk2),
     &         len_str2(ncblk2+nablk2))

      ! try to find the corresponding graphs
      idxg1 = 0
      idxg2 = 0
      idxg12 = 0
      do ica = 1, 2
        do iop = 1, 3
          select case(iop)
          case(1)
            pocc => occ1(IVALE,ica)
            pgraph => idxg1(IVALE,ica)
          case(2)
            pocc => occ2(IVALE,ica)
            pgraph => idxg2(IVALE,ica)
          case(3)
            pocc => occ12(IVALE,ica)
            pgraph => idxg12(IVALE,ica)
          end select
          if (pocc.eq.0) cycle
          do igraph = 1, str_info%ngraph
            if (IVALE.eq.str_info%ispc_typ(igraph).and.
     &          pocc.eq.str_info%ispc_occ(igraph)) then
              if (pgraph.ne.0) call quit(1,'get_vec2remove',
     &                       'more than one possible graph found!')
              pgraph = igraph
            end if
          end do
          if (pgraph.eq.0) call quit(1,'get_vec2remove',
     &                               'no graph found!')
        end do
      end do
c dbg
c      call wrt_occ(lulog,occ1)
c      call wrt_occ(lulog,occ2)
c      call wrt_occ(lulog,occ12)
c      call wrt_occ(lulog,idxg1)
c      call wrt_occ(lulog,idxg2)
c      call wrt_occ(lulog,idxg12)
c dbgend

      call condense_occ(occ_c1,occ_a1,
     &                  hpvx_c1,hpvx_a1,
     &                  occ1,1,hpvxblkseq)
      call condense_occ(occ_c2,occ_a2,
     &                  hpvx_c2,hpvx_a2,
     &                  occ2,1,hpvxblkseq)
      call condense_occ(occ_c12,occ_a12,
     &                  hpvx_c12,hpvx_a12,
     &                  occ12,1,hpvxblkseq)
      call condense_occ(graph_c1,graph_a1,
     &                  hpvx_c1,hpvx_a1,
     &                  idxg1,1,hpvxblkseq)
      call condense_occ(graph_c2,graph_a2,
     &                  hpvx_c2,hpvx_a2,
     &                  idxg2,1,hpvxblkseq)
      call condense_occ(graph_c12,graph_a12,
     &                  hpvx_c12,hpvx_a12,
     &                  idxg12,1,hpvxblkseq)

      call set_mapping_info(map_info_c12,map_info_a12,0,
     &                      occ1,1,.false.,
     &                      occ2,1,.false.,
     &                      occ12,(/1,1,1,1/),1,hpvxblkseq)

      call strmap_man_c(1,maxbuf,
     &     graph_c1,1,graph_c2,ncblk2,
     &     graph_c12,1,map_info_c12,
     &     str_info,strmap_info,orb_info)
c dbg
c      print *,'length of map_c12:',maxbuf
c dbgend
      ifree = mem_alloc_int(map_c12,maxbuf,'map_c12')
      call strmap_man_c(1,maxbuf,
     &     graph_a1,1,graph_a2,nablk2,
     &     graph_a12,1,map_info_a12,
     &     str_info,strmap_info,orb_info)
c dbg
c      print *,'length of map_a12:',maxbuf
c dbgend
      ifree = mem_alloc_int(map_a12,maxbuf,'map_a12')

      ! loop over Ms(A)/Gamma(A) for high-rank operator
      ioff12 = 0
      do msa12 = na, -na, -2
        msc12 = msa12 + ms
        if (abs(msc12).gt.nc) cycle
        do gama12 = 1, ngam
          gamc12 = multd2h(gama12,igam)

          call ms2idxms(idxmsdis_c12,msc12,occ_c12,1)
          call ms2idxms(idxmsdis_a12,msa12,occ_a12,1)
          call set_len_str(len_str12,1,1,str_info%g,
     &                     graph_c12,idxmsdis_c12,gamc12,hpvx_c12,
     &                     graph_a12,idxmsdis_a12,gama12,hpvx_a12,
     &                     hpvxseq,.false.)
          lenca12 = len_str12(1)*len_str12(2)
          if (lenca12.eq.0) cycle
c dbg
c          print *,'msa,gama,len_str:',msa12,gama12,len_str12(1:2)
c dbgend

          ! find all elements related to low-rank operator
          ioff2 = 0
          do msa2 = na2, -na2, -2
            msc2 = msa2 + ms
            if (abs(msc2).gt.nc2) cycle
            do gama2 = 1, ngam
              gamc2 = multd2h(gama2,igam)
              if (na2.eq.0.and.gama2.ne.1.or.
     &            nc2.eq.0.and.gamc2.ne.1) cycle

              call ms2idxms(idxmsdis_c2,msc2,occ_c2,ncblk2)
              call ms2idxms(idxmsdis_a2,msa2,occ_a2,nablk2)
              call set_len_str(len_str2,ncblk2,nablk2,str_info%g,
     &                     graph_c2,idxmsdis_c2,gamc2,hpvx_c2,
     &                     graph_a2,idxmsdis_a2,gama2,hpvx_a2,
     &                     hpvxseq,.false.)
              if (na2.eq.0.and.nc2.eq.0) then
                lenca2 = 1
              else if (na2.eq.0.or.nc2.eq.0) then
                lenca2 = len_str2(1)
              else
                lenca2 = len_str2(1)*len_str2(2)
              end if
              if (lenca2.eq.0) cycle
              ! can there be a matching element from unit tensor?
              msc1 = msc12 - msc2
              msa1 = msa12 - msa2
              gamc1 = multd2h(gamc12,gamc2)
              gama1 = multd2h(gama12,gama2)
              if (msc1.ne.msa1.or.gamc1.ne.gama1.or.
     &            abs(msc1).gt.drank.or.abs(msa1).gt.drank) then
                ioff2 = ioff2 + lenca2
                cycle
              end if
              call ms2idxms(idxmsdis_c1,msc1,occ_c1,1)
              call ms2idxms(idxmsdis_a1,msa1,occ_a1,1)
              call set_len_str(len_str1,1,1,str_info%g,
     &                     graph_c1,idxmsdis_c1,gamc1,hpvx_c1,
     &                     graph_a1,idxmsdis_a1,gama1,hpvx_a1,
     &                     hpvxseq,.false.)
c dbg
c          print *,'msa1,gama1,len_str1:',msa1,gama1,
c     &            len_str1(1:2)
c          print *,'msa2,gama2,len_str2:',msa2,gama2,
c     &            len_str2(1:ncblk2+nablk2)
c dbgend

              ! get maps
              call get_strmap_blk_c(map_c12,
     &             1,ncblk2,1,
     &             occ_c1,occ_c2,len_str1,len_str2,
     &             graph_c1,graph_c2,graph_c12,
     &             idxmsdis_c1,idxmsdis_c2,gamc1,gamc2,map_info_c12,
     &             strmap_info,ngam,str_info%ngraph)
              call get_strmap_blk_c(map_a12,
     &             1,nablk2,1,
     &             occ_a1,occ_a2,len_str1(2),len_str2(ncblk2+1),
     &             graph_a1,graph_a2,graph_a12,
     &             idxmsdis_a1,idxmsdis_a2,gama1,gama2,map_info_a12,
     &             strmap_info,ngam,str_info%ngraph)

              nstrc1 = ielprd(len_str1,1)
              nstrc2 = ielprd(len_str2,ncblk2)
              nstra2 = ielprd(len_str2(ncblk2+1),nablk2)
              nstrc12 = ielprd(len_str12,1)
              nstra12 = ielprd(len_str12(2),1)

              ! loop over C strings of low-rank operator
              do istrc2 = 1, nstrc2
                ! loop over C strings of unit tensor
                do istrc1 = 1, nstrc1
                  ielmap = map_c12((istrc2-1)*nstrc1+istrc1)
                  if (ielmap.eq.0) cycle
                  istrc12 = abs(ielmap)
                  isgnc = sign(1,ielmap)
                  ! loop over A strings of low-rank operator
                  do istra2 = 1, nstra2
                    ! unit tensor: A string = C string!
                    ielmap = map_a12((istra2-1)*nstrc1+istrc1)
                    if (ielmap.eq.0) cycle
                    istra12 = abs(ielmap)
                    isgn = isgnc * sign(1,ielmap)
                    idx2 = (istrc2-1)*nstra2+istra2
                    idx12 = (istrc12-1)*nstra12+istra12
c dbg
c                    print *,'idx2,idx12,istrc1:',idx2,idx12,istrc1
c dbgend
                    ! set elements of p vectors:
                    pvecs(ioff12+idx12,1:xdim)
     &                = pvecs(ioff12+idx12,1:xdim)
     &                + dble(isgn)*facinv*xnew(ioff2+idx2,1:xdim)
                  end do
                end do
              end do

              ioff2 = ioff2 + lenca2
            end do
          end do

          ioff12 = ioff12 + lenca12
        end do
      end do

      deallocate(occ_c2,occ_a2,graph_c2,graph_a2,hpvx_c2,hpvx_a2,
     &           idxmsdis_c2,idxmsdis_a2,len_str2)
      ifree = mem_flushmark()

      if (ntest.ge.100) then
        write(lulog,*) 'output vectors (to be projected out):'
        call wrtmat2(pvecs,pdim,xdim,pdim,xdim)
      end if

      return
      end
