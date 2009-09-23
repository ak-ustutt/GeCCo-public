*----------------------------------------------------------------------*
      subroutine import_2el_presort(buffer,ffpre,ffchain,
     &     incore,sbuffer,npass,len_bin,len_rec,
     &     ffinp,oplist,
     &     fac_aa,fac_ab,fac_abx,fac_ccaa,
     &     set_abx,set_ccaa,set_ca,set_pq_rs,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*
*     sort integral from file ffinp (special DALTON format)
*     - for incore==.true.  into buffer
*     - for incore==.false. into the chained direct access file ffpre
*                           (Yoshimine-type of sort)
*
*     andreas, may 2009
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'hpvxseq.h'
      include 'def_graph.h'
      include 'def_sort_buffer.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      
      integer, parameter ::
     &     ntest = 000

      type(me_list),intent(inout) ::
     &     oplist
      type(orbinf),intent(in),target ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(filinf), intent(inout) ::
     &     ffpre, ffchain, ffinp
      real(8), intent(in) ::
     &     fac_aa,fac_ab,fac_abx,fac_ccaa
      logical, intent(in) ::
     &     set_abx,set_ccaa,set_ca,set_pq_rs,incore

      integer, intent(in) ::
     &     npass, len_bin, len_rec
      real(8), intent(inout) ::
     &     buffer(*)
      type(sort_buffer), intent(inout) ::
     &     sbuffer

      real(8) ::
     &     fac, spfac(6)
      integer ::
     &     ngam, ntoob, caborb, ntotal, lu2in, lu_pre, lu_chain,
     &     ii, idx, idxst, n_so, max_so, idx_off, n_idx_cur, nx_max,
     &     ivl, ifc
      integer ::
     &     ipass,
     &     idxbin, idxbin_prev, first_bin, last_bin,
     &     bins_per_pass, bins_total, max_chain, len_per_rec,
     &     idxrec, jdx, kdx
      integer(8) ::
     &     maxlength, length, rec_wr, len_wr, idx_wr

      integer, parameter ::
     &     n_idx_sim = 16 ! more seems to give cache problems

      ! block look-up cache
      integer, parameter ::
     &     max_cache = 16
      integer ::
     &     len_cache,
     &     cache_tab(3,max_cache), cache_stat(3)

      integer, pointer ::
     &     igamorb(:), igasorb(:), hpvxgas(:,:),
     &     reost(:), reord(:)
      integer(2) ::
     &     idxprqs(4,n_idx_sim),
     &     igmprqs(4,n_idx_sim), ispprqs(4,n_idx_sim)
      integer(2), pointer ::
     &     idxbuf(:,:)
      integer, pointer ::
     &     idxlist(:), occ_hash(:),
     &     idxval(:), idxspin(:), idxperm(:),
     &     cur_idx(:), cur_chn(:), cur_len
      real(8), pointer ::
     &     valbuf(:), cur_val(:)

      type(sort_bin), pointer ::
     &     bin(:)

      integer, external ::
     &     max_rank_op

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'import_2el_presort')
        write(luout,*) 'incore = ',incore
      end if

      if (incore.and.npass.ne.1)
     &     call quit(luout,'import_2el_presort','npass>1 for incore?')

      ngam  = orb_info%nsym
      ntoob = orb_info%ntoob
      caborb = orb_info%caborb
      igamorb => orb_info%igamorb
      igasorb => orb_info%igasorb
      hpvxgas => orb_info%ihpvgas
      
      ntotal = ntoob+caborb

      allocate(reord(ntoob+caborb))

      reost=>orb_info%ireost
      do idx = 1,ntoob
        reord(idx)=reost(idx)
      enddo
      do idx=ntoob+1,ntoob+caborb
        reord(idx)=idx
      enddo  

      if (incore) then
        buffer(1:oplist%len_op) = 0d0
      else
        bins_per_pass = sbuffer%nbins
        len_per_rec   = sbuffer%max_length
        max_chain     = sbuffer%max_chain
        bins_total    = (oplist%len_op-1)/len_bin + 1
c dbg
c        print *,'bins_per_pass = ',bins_per_pass
c        print *,'bins_total    = ',bins_total
c        print *,'len_per_rec   = ',len_per_rec
c        print *,'len_bin       = ',len_bin
c dbg
        bin => sbuffer%bin
      end if

      lu2in = ffinp%unit
      if (.not.incore) then
        lu_pre = ffpre%unit
        lu_chain = ffchain%unit

        rewind(lu_chain)
        rec_wr = len_rec
        len_wr = len_per_rec
        idx_wr = max_chain
        write(lu_chain) rec_wr, len_wr, idx_wr
      end if

      ! rewind and get block length from first record
      rewind(lu2in)
      read(lu2in) maxlength

      max_so = 4*12*n_idx_sim
      allocate(idxbuf(4,maxlength),valbuf(maxlength),
     &         idxval(max_so),idxspin(max_so),idxperm(max_so),
     &         idxlist(max_so),occ_hash((ngastp*(ngastp+1)/2)**2))
      
      call set_occ_hash(occ_hash,oplist%op)

      nx_max = max_rank_op('XTOT',oplist%op,.true.)

      spfac(1:6) = (/fac_aa,fac_aa,fac_ab,fac_ab,fac_abx,fac_abx/)

      len_cache = 0
      cache_stat(1:3) = 0

      idxrec = 0
      last_bin = 0
      idxbin_prev = 0
      ! loop over passes over input file (hopefully 1)
      do ipass = 1, npass

        rewind(lu2in)
        read(lu2in) ! over-read first record

        if (.not.incore) then
          first_bin = last_bin+1
          last_bin  = min(last_bin+bins_per_pass,bins_total)

          ! reset buffer
          do idxbin = 1, bins_per_pass
            bin(idxbin)%chain(1) = 1  ! chain length + 1 
            bin(idxbin)%chain(2:max_chain) = 0
            bin(idxbin)%length = 0
          end do

        end if

        ! read from file to buffer
        do
          read(lu2in)
     &         length,idxbuf(1:4,1:length),valbuf(1:length)
          if (length.le.0) exit
      
          do idx_off = 0, length-1, n_idx_sim

            n_idx_cur = min(length-idx_off,n_idx_sim)

            ! convert to GeCCo ordering, and provide symmetry and space info
            call idx_reform(idxprqs,igmprqs,ispprqs,
     &                    idxbuf(1,idx_off+1),n_idx_cur,nx_max,
     &                                                fac_ccaa.lt.0d0,
     &                    reord,orb_info%igamorb,orb_info%igasorb,
     &                          orb_info%ihpvgas,ntotal) 
                                ! ^^^ used only for X -> spin-indep.

            ! provide expansion to spin-orbitals
            call so_expand(idxval,idxspin,idxperm,n_so,max_so,
     &                   idxprqs,n_idx_cur,
     &                   .true.,.true.,.true.,.true.,set_abx,
     &                   set_ccaa,set_ca)

            
            ! get indices on operator list
            call idx4so(idxlist,
c            call idx4so_ori(idxlist,
     &                idxprqs,igmprqs,ispprqs,n_idx_cur,
     &                idxval,idxspin,idxperm,n_so,
     &                fac_ccaa.lt.0d0,set_pq_rs,
     &                oplist,occ_hash,
     &                cache_tab,len_cache,max_cache,cache_stat,
     &                str_info,orb_info)

            if (incore) then
              ! direct distribution to output buffer
              do ii = 1, n_so
                ivl    = idxval(ii)
                ifc    = idxspin(ii)
                idx    = abs(idxlist(ii))
                fac    = dble(sign(1,idxlist(ii)))*spfac(ifc)
                if (idx.eq.0) cycle
                buffer(idx) = buffer(idx)+fac*valbuf(idx_off+ivl)
              end do
            else
c dbg
c              print *,'out-of-core part entered'
c dbg
              ! distribute into appropriate bins
              do ii = 1, n_so
                ivl    = idxval(ii)
                ifc    = idxspin(ii)
                idx    = abs(idxlist(ii))
                fac    = dble(sign(1,idxlist(ii)))*spfac(ifc)
                if (idx.eq.0) cycle
                idxbin = (idx-1)/len_bin + 1
                if (idxbin.lt.first_bin.or.idxbin.gt.last_bin) cycle
                ! correct for offset
                idxbin = idxbin - first_bin+1
c dbg
c                print *,'sorting to bin ',idxbin
c dbg
                ! point to present bin
                if (idxbin.ne.idxbin_prev) then
                  cur_len => bin(idxbin)%length
                  cur_val => bin(idxbin)%value
                  cur_idx => bin(idxbin)%index
                end if
                ! store value in bin
                jdx = cur_len+1
                cur_len = jdx
                cur_val(jdx) = fac*valbuf(idx_off+ivl)
                cur_idx(jdx) = idx

                ! empty full bin:
                if (jdx.eq.len_per_rec) then
c dbg
c                  print *,'empty bin # ',idxbin
c dbg
                  idxrec = idxrec + 1
                  len_wr = jdx  ! len_wr is always integer(8)

                  write(lu_pre,rec=idxrec) len_wr,
     &                                     cur_val(1:jdx),
     &                                     cur_idx(1:jdx)
                  cur_len = 0
                  
                  ! store chain information
                  cur_chn => bin(idxbin)%chain

                  kdx = cur_chn(1) + 1

                  if (kdx-1.gt.max_chain) then
                    call quit(1,'import_2el_presort',
     &                   'should not happen: chain too long!')
                  end if

                  cur_chn(1) = kdx
                  cur_chn(kdx) = idxrec
c dbg
c                  print *,'current chain: '
c                  print *,cur_chn(1:kdx)
c dbg
                end if

              end do ! loop over present batch of indices/values
    
            end if ! in-core / out-of-core switch
          
          end do ! loop over record

        end do ! loop over input file records
        ! empty buffers
        if (.not.incore) then

          do idxbin = 1, bins_per_pass
c dbg
c            print *,'final: empty bin # ',idxbin
c dbg
            cur_len => bin(idxbin)%length
            cur_val => bin(idxbin)%value
            cur_idx => bin(idxbin)%index

            if (cur_len.gt.0) then
              idxrec = idxrec + 1
              jdx    = cur_len
              len_wr = jdx      ! len_wr is always integer(8)

              write(lu_pre,rec=idxrec) len_wr,
     &                                 cur_val(1:jdx),
     &                                 cur_idx(1:jdx)
              cur_len = 0
                  
              ! store chain information
              cur_chn => bin(idxbin)%chain

              kdx = cur_chn(1) + 1

              if (kdx-1.gt.max_chain) then
                call quit(1,'import_2el_presort',
     &                   'should not happen: chain too long! (2)')
              end if

              cur_chn(1) = kdx
              cur_chn(kdx) = idxrec
c dbg
c              print *,'current chain: '
c              print *,cur_chn(1:kdx)
c dbg
            end if

          end do
          idxbin_prev = 0

          ! write chains to chain-file
          idx = 0
          do idxbin = first_bin, last_bin
            idx = idx+1
            idx_wr = idxbin
            len_wr = bin(idx)%chain(1)-1
            write(lu_chain) idx_wr,len_wr,bin(idx)%chain(2:len_wr+1)
          end do
        end if

      end do ! passes over input file

c dbg
      print *,'cache statistics:'
      print *,'max. length = ',max_cache
      print *,'used length = ',len_cache
      print *,'hits   = ',cache_stat(1)
      print *,'misses = ',cache_stat(2)
      print *,'misses (crit) = ',cache_stat(3)
c dbg
      if (.not.incore) then
        ! last record on chain file
        idx_wr = -1
        len_wr = 0
        write(lu_chain) idx_wr,len_wr
      end if

      deallocate(reord,idxbuf,valbuf,idxlist,occ_hash,
     &           idxval,idxspin,idxperm)

      if (ntest.ge.1000.and..not.incore) then
        write(luout,'(x,76("="))') 
        write(luout,*) 'dump of presort file:'
        write(luout,'(x,76("-"))') 
        call list_file_chda(ffpre,ffchain)
        write(luout,'(x,76("="))') 
      end if

      return

      end

