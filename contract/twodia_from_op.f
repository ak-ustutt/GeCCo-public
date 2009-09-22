      subroutine twodia_from_op(x2dia,!offsets,nblk,
     &                          mel,
     &                          orb_info,str_info)
*-----------------------------------------------------------------------
*     Routine to extract the diagonal part of a two-electron operator.
*     GWR November 2007
*-----------------------------------------------------------------------

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'ifc_baserout.h'
      include 'hpvxseq.h'

      real(8), intent(out) ::
     &     x2dia(*)
c      integer, intent(out) ::
c     &     offsets(*), nblk
      type(me_list), intent(in) ::
     &     mel
      type(orbinf), intent(in), target ::
     &     orb_info
      type(strinf), intent(in), target ::
     &     str_info

      integer ::
     &     ifree, iocc_cls, nbuff, ioff_blk, ilen_blk, ms, idxms,
     &     imo_off, isym, idxcur, idx, jdx, igas, jgas, idxbuf,
     &     x2_off, len1, len2, idxdis, ioff, joff,
     &     njoined, ntoob, nsym, iblk_off, ncsub, nasub, ld_j
      integer ::
     &     ihpv(2), occ_temp(ngastp,2)
      integer, pointer ::
     &     ihpvgas(:,:), iad_gas(:), mostnd(:,:,:),
     &     hpvx_occ(:,:,:), idx_graph(:,:,:), idx_gas(:), 
     &     cinfo(:,:), ms_dis(:), gm_dis(:), idxms_dis(:), lstr(:),
     &     occ_blk_pnt(:,:,:), graph_blk_pnt(:,:,:)
      integer, allocatable ::
     &     map(:)
      type(graph), pointer ::
     &     graphs(:)

      logical ::
     &     loop(mel%op%n_occ_cls), first

      real(8), pointer ::
     &     buffer(:), curblk(:)

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop

      logical, external ::
     &     occ_is_diag_blk, next_msgamdist_diag
      integer, external ::
     &     idx_msgmdst2

      ifree = mem_setmark('twodia')

      op => mel%op
      ffop => mel%fhand

      ! Initialise some things.
      nbuff = 0
      loop(1:op%n_occ_cls) = .false.
      njoined = op%njoined
      hpvx_occ => op%ihpvca_occ
      graphs => str_info%g
      idx_graph => mel%idx_graph

      do iocc_cls = 1,op%n_occ_cls
        if(max(op%ica_occ(1,iocc_cls),
     &       op%ica_occ(2,iocc_cls)).ne.2  .or.
     &       op%formal_blk(iocc_cls)  .or. .not.
     &       occ_is_diag_blk(hpvx_occ(1,1,(iocc_cls-1)*njoined+1),
     &       njoined))
     &       loop(iocc_cls) = .true.          
      enddo

      ! Allocate space if the file, ffop, is not buffered.
      if(.not.ffop%buffered) then
        nbuff = 0
        do iocc_cls = 1,op%n_occ_cls
          if (loop(iocc_cls)) cycle

          nbuff = max(nbuff, mel%len_op_occ(iocc_cls))
        enddo

        ifree = mem_alloc_real(buffer,nbuff,'buffer')
      endif

      ! Initialise the array in which the diagonal elements will be held.
      ! This is the maximum length needed if no externals (I hope).
      x2dia(1:4*orb_info%ntoob**2) = 0d0

      mostnd => orb_info%mostnd
      ihpvgas => orb_info%ihpvgas
      iad_gas => orb_info%iad_gas
      idx_gas => orb_info%idx_gas
      nsym = orb_info%nsym
      ntoob = orb_info%ntoob

c      nblk = 0
      do iocc_cls = 1,op%n_occ_cls
        ! Two-electron operators of the appropriate type only.
        ! Formal blocks are ignored.
        if(loop(iocc_cls))cycle
        
        ioff_blk = mel%off_op_occ(iocc_cls)
        ilen_blk = mel%len_op_occ(iocc_cls)
        if(ffop%buffered)then
          idxbuf = ffop%idxrec(iocc_cls)+1
          if(idxbuf.lt.0)
     &         call quit(1,'twodia_from_op',
     &         'idxrec inconsistent for file '//trim(ffop%name))
          buffer => ffop%buffer(idxbuf:)
        else
          call get_vec(ffop,buffer,ioff_blk+1,ioff_blk+ilen_blk)
        endif

        if (ntest.ge.100) then
          write(luout,*) 'operator(in) (',trim(op%name),
     &         ',list=',trim(mel%label),')'
          call wrt_mel_buf(luout,5,buffer,mel,iocc_cls,iocc_cls,
     &                  str_info,orb_info)
        end if
          
        iblk_off = (iocc_cls-1)*njoined

        call get_num_subblk(ncsub,nasub,
     &       hpvx_occ(1,1,iblk_off+1),njoined)

        if (ncsub.ne.nasub)
     &       call quit(1,'twodia_from_op','haeehhh????')

        allocate(cinfo(ncsub,3),ms_dis(ncsub),gm_dis(ncsub),
     &       idxms_dis(ncsub),lstr(2*ncsub))

        occ_blk_pnt =>
     &         hpvx_occ(1:ngastp,1:2,iblk_off+1:iblk_off+njoined)
        graph_blk_pnt =>
     &        idx_graph(1:ngastp,1:2,iblk_off+1:iblk_off+njoined)

        ! set HPVX and OCC info
        call condense_occ(cinfo, cinfo,
     &                    cinfo(1,3),cinfo(1,3),
     &                    occ_blk_pnt,njoined,hpvxblkseq)
        ! do the same for the graph info
        call condense_occ(cinfo(1,2),cinfo(1,2),
     &                    cinfo(1,3),cinfo(1,3),
     &                    graph_blk_pnt,
     &                                njoined,hpvxblkseq)        

        ! Loop over the Ms of the annihilator string (equal to that of
        ! the creator as the operators should be symmetric).
        idxms = 0
        do ms = 2, -2, -2
          idxms = idxms + 1
c          if(ms.eq.-2) idxms = idxms+1
          x2_off = (idxms-1)*orb_info%ntoob**2
          ld_j   = orb_info%ntoob
c dbg
c          print *,'x2_off = ',x2_off
c          print *,'ld_j   = ',ld_j
c dbg

          ! Loop over irrep. of annihilator string.
          do isym = 1, nsym

c            nblk = nblk+1
c            offsets(nblk) = x2_off

            first = .true.
            distr_loop: do
              if (.not.next_msgamdist_diag(first,
     &              ms_dis,gm_dis,
     &              nasub,cinfo,
     &              ms,isym,nsym))
     &                 exit distr_loop
              first = .false.

              call ms2idxms(idxms_dis,ms_dis,
     &              cinfo,nasub)

              call set_len_str(lstr,nasub,nasub,
     &              graphs,
     &              cinfo(1,2),idxms_dis,gm_dis,cinfo(1,3),
     &              cinfo(1,2),idxms_dis,gm_dis,cinfo(1,3),
     &              hpvxseq,.false.)

              if (idxlist(0,lstr,ncsub+nasub,1).gt.0)
     &             cycle distr_loop
              
              idxdis = idx_msgmdst2(iocc_cls,idxms,isym,
     &                              cinfo,idxms_dis,gm_dis,nasub,
     &                              cinfo,idxms_dis,gm_dis,nasub,
     &                              .false.,mel,nsym)

              ! Offset for the buffer array.
              idxcur =
     &             mel%off_op_gmox(iocc_cls)%d_gam_ms(idxdis,isym,idxms)
     &             -ioff_blk+1
              curblk => buffer(idxcur:)

              ! Actual loops to extract diagonal elements from each block.
              if (nasub.eq.1) then
                ! Length of the block required in x2dia
                len1 = lstr(1)

                ! use map string index --> matrix element index
                allocate(map(len1))
                call set_map_diag(map,len1,ms,isym,mel,iocc_cls,
     &                            str_info,orb_info)
                do idx = 1, len1
                  x2dia(x2_off+map(idx)) = curblk((idx-1)*len1+idx)
                enddo
                deallocate(map)

              else if (nasub.eq.2) then
                len1 = lstr(1)
                len2 = lstr(2)
            
                igas  = idx_gas(cinfo(1,3))
                do while (iad_gas(igas).ne.2) 
                  igas = igas+1
                end do
                jgas  = idx_gas(cinfo(2,3))
                do while (iad_gas(jgas).ne.2) 
                  jgas = jgas+1
                end do

                joff = mostnd(1,gm_dis(2),jgas)-1
                ioff = mostnd(1,gm_dis(1),igas)-1
c dbg
c                print *,'igas,ioff: ',igas,ioff
c                print *,'jgas,joff: ',jgas,joff
c dbg    

                do jdx = 1, len2
                  do idx = 1, len1
c dbg
c                    print *,'idx, jdx : ',idx,jdx
c                    print *,'imo, jmo : ',idx+ioff,jdx+joff
c dbg
                    x2dia(x2_off+(joff+jdx-1)*ld_j+ioff+idx) =
     &                   curblk((jdx-1)*len1*len1*len2+
     &                          (jdx-1)*len1*len1+
     &                          (idx-1)*len1 + idx)
                  end do
                end do

              end if

            end do distr_loop

          enddo

        enddo

        deallocate(cinfo,idxms_dis,ms_dis,gm_dis,lstr)

      enddo

      if (ntest.ge.100) then
c        write(luout,*) 'blocks: ',nblk
c        write(luout,*) 'offsets:',offsets(1:nblk)
        write(luout,*) 'diagonal:'
        do idxms = 1, 3
          x2_off = (idxms-1)*ntoob**2
          write(luout,*) 'idxms = ',idxms
          call wrtmat2(x2dia(x2_off+1),ntoob,ntoob,ntoob,ntoob)
        end do
      end if

      ifree = mem_flushmark()

      return
      end
