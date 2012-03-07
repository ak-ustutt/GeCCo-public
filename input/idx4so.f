*----------------------------------------------------------------------*
      subroutine idx4so(idxlist,
     &                  idxprqs,igmprqs,ispprqs,nidx,
     &                  idxval,idxspin,idxperm,n_so,
     &                  antiherm,pq_rs,
     &                  me,occ_hash,
     &                  cache_tab,len_cache,max_cache,cache_stat,
     &                  str_info,orb_info)
*----------------------------------------------------------------------*
*     input: list of spin-orbital indices
*     output: list of actual list-indices (if existent) + sign, if appr.
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'ioparam.h'
      include 'multd2h.h'
      include 'par_int_permute.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'hpvxseq.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf),intent(in),target ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(me_list),intent(inout) ::
     &     me
      integer, intent(in) ::
     &     n_so, nidx, occ_hash(100)
      integer(2), intent(inout) ::
     &     idxprqs(4,nidx), igmprqs(4,nidx), ispprqs(4,nidx)
      integer, intent(in) ::
     &     idxval(n_so), idxspin(n_so), idxperm(n_so)
      integer, intent(out) ::
     &     idxlist(n_so)
      logical, intent(in) ::
     &     antiherm,pq_rs

      logical ::
     &     reopr, reoqs, sgn_change, equal
      integer ::
     &     ii, i0, ipm, iind, nind, igraph, nstr, istr, jstr,
     &     gamt, mst, gama, idxmsa, idis,
     &     ip, iq, ir, is, ips, iqs, irs, iss, iscr, len, nj, nspin,
     &     igastp, jgastp, idx, isgn
c test
     &     ,s11,s21,s31,s12,s22,s32
      integer ::
     &     hpvx,
     &     hpvx_p, hpvx_q, hpvx_r, hpvx_s,
     &     gas_p, gas_q, gas_r, gas_s,
     &     gam_p, gam_q, gam_r, gam_s,
     &     ihash, ihash_last, ihash_blk, iblk, iblkoff, ioff
      integer ::
     &     ispin_raw(4),
     &     idxraw(6,n_so), string(3,2,4),
     &     idgam(4), idspc(4), strseq(4),
     &     string_list(3,2,4*n_so), idxlenstr(2,4*n_so),
     &     occ(ngastp,2,2), strl(2,2)
c      integer(2) :: string(3,2,4), string_list(3,2,4*n_so)

      integer ::
     &     nsym
      integer, pointer ::
     &     igasorb(:), igamorb(:), hpvxgas(:,:),
     &     idx_graph(:,:,:), gas_restr(:,:,:,:,:),
     &     ispc_occ(:), ngas_hpv(:), mostnd(:,:,:),
     &     ioff_gas(:), idx_gas(:),
     &     off_d_gam_ms(:,:,:), off_gam_ms(:,:)
      type(graph), pointer ::
     &     graphs(:), curgraph

      ! block look-up cache
c      integer, parameter ::
      integer ::
     &     max_cache
      integer ::
     &     len_cache, idx_cache,
c     &     cache_tab1(max_cache), cache_tab2(2,max_cache)
     &     cache_tab(3,max_cache), cache_stat(3)
      
      integer, external ::
     &     idx4sg, iblk_occ, search_list6
      logical, external ::
     &     allow_sbsp_dis, list_cmp


      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'idx4so')
        write(luout,*) 'n_so = ',n_so
      end if

      hpvxgas => orb_info%ihpvgas
      idx_graph => me%idx_graph

      nsym = orb_info%nsym
      nspin = orb_info%nspin
      nj = me%op%njoined

      nstr = 0

      ihash_last = -1
      idx_cache = len_cache
c      idx_cache = 0
c      len_cache = 0

      do ii = 1, n_so

        ! set spincase
        i0 = idxspin(ii)
        if (i0.eq.1) then
          ispin_raw = (/ 1, 1, 1, 1/)
        else if (i0.eq.2) then
          ispin_raw = (/-1,-1,-1,-1/)
        else if (i0.eq.3) then
          ispin_raw = (/ 1,-1, 1,-1/)
        else if (i0.eq.4) then
          ispin_raw = (/-1, 1,-1, 1/)
        else if (i0.eq.5) then
          ispin_raw = (/ 1,-1,-1, 1/)
        else if (i0.eq.6) then
          ispin_raw = (/-1, 1, 1,-1/)
        end if

        i0 = idxval(ii)
        ipm = idxperm(ii)

        sgn_change = antiherm.and.ipm.eq.2

        ! indices
        ip = idxprqs(int_permute(1,ipm),i0)
        ir = idxprqs(int_permute(2,ipm),i0)
        iq = idxprqs(int_permute(3,ipm),i0)
        is = idxprqs(int_permute(4,ipm),i0)

        ! irreps
        gam_p = igmprqs(int_permute(1,ipm),i0)
        gam_r = igmprqs(int_permute(2,ipm),i0)
        gam_q = igmprqs(int_permute(3,ipm),i0)
        gam_s = igmprqs(int_permute(4,ipm),i0)

        ! GAS
        gas_p = ispprqs(int_permute(1,ipm),i0)
        gas_r = ispprqs(int_permute(2,ipm),i0)
        gas_q = ispprqs(int_permute(3,ipm),i0)
        gas_s = ispprqs(int_permute(4,ipm),i0)          

        ips = ispin_raw(1)
        irs = ispin_raw(2)
        iqs = ispin_raw(3)
        iss = ispin_raw(4)

        if (ntest.ge.150) then
          write(luout,*) 'initial: ',ip,ir,iq,is
          write(luout,*) '         ',ips,irs,iqs,iss
        end if

        reopr = ip.gt.ir .or. (ip.eq.ir.and.ips.lt.irs)
        reoqs = iq.gt.is .or. (iq.eq.is.and.iqs.lt.iss)
        if (reopr) then
          iscr = ip
          ip = ir
          ir = iscr
          iscr = ips
          ips = irs
          irs = iscr
          iscr = gam_p
          gam_p = gam_r
          gam_r = iscr
          iscr = gas_p
          gas_p = gas_r
          gas_r = iscr
        end if
        if (reoqs) then
          iscr = iq
          iq = is
          is = iscr
          iscr = iqs
          iqs = iss
          iss = iscr
          iscr = gam_q
          gam_q = gam_s
          gam_s = iscr
          iscr = gas_q
          gas_q = gas_s
          gas_s = iscr
        end if

        if (ntest.ge.100) then
          write(luout,'(x,a,4i5,4x,4i3)') 'generated: ',
     &         ip,ir,iq,is,ips,irs,iqs,iss
        end if

        ! get HPVX
        if (nspin.eq.1) then
          hpvx_p = hpvxgas(gas_p,1)
          hpvx_r = hpvxgas(gas_r,1)
          hpvx_q = hpvxgas(gas_q,1)
          hpvx_s = hpvxgas(gas_s,1)
        else
          hpvx_p = hpvxgas(gas_p,(ips+1)/2+1)
          hpvx_r = hpvxgas(gas_r,(irs+1)/2+1)
          hpvx_q = hpvxgas(gas_q,(iqs+1)/2+1)
          hpvx_s = hpvxgas(gas_s,(iss+1)/2+1)
        end if

        ! get block of list
        ihash = (ngastp*(ngastp+1)/2)*
     &          (max(hpvx_p,hpvx_r)*(max(hpvx_p,hpvx_r)-1)/2+
     &           min(hpvx_p,hpvx_r) - 1 ) +
     &           max(hpvx_q,hpvx_s)*(max(hpvx_q,hpvx_s)-1)/2+
     &           min(hpvx_q,hpvx_s)

        ! turn into unique identifier for block-info cache
        !  gam_p, gam_r, gam_q
        !  ips, irs, iqs
        !  ihash
        ihash_blk = ihash*(8*8*8*8) + 
     &       (((gam_p-1)*8+gam_r-1)*8+gam_q-1)*8
        if (ips.gt.0) ihash_blk = ihash_blk+4
        if (irs.gt.0) ihash_blk = ihash_blk+2
        if (iqs.gt.0) ihash_blk = ihash_blk+1
Cc dbg
C        print '(a,6i3,i4,a,i10)',
C     &       'ihash_blk: ',gam_p,gam_r,gam_q,ips,irs,iqs,ihash,
C     &       '->',ihash_blk
Cc dbg
C
       ioff = -1
        do idx = 1, len_cache
c          if (ihash_blk.eq.cache_tab1(idx)) then
c            iblk   = cache_tab2(1,idx)
c            ioff   = cache_tab2(2,idx)
          if (ihash_blk.eq.cache_tab(1,idx)) then
            iblk   = cache_tab(2,idx)
            ioff   = cache_tab(3,idx)
            cache_stat(1) = cache_stat(1)+1 ! hit
Cc dbg
C            print *,'found: ',idx,'->iblk,ioff = ',iblk,ioff
Cc dbg
            exit
          end if
        end do

        if (ioff.lt.0) then
C          if (ihash.ne.ihash_last) then
c dbg
C            if (ioff.gt.0.and.iblk.ne.occ_hash(ihash)) then
C              print *,'iblk, iblk(hash): ',iblk,occ_hash(ihash)
C              stop 'err0'
C            end if
c dbg
            iblk = occ_hash(ihash)
c            ihash_last = ihash
C          end if
        end if

        idxraw(1:6,ii) = 0

        ! no success? skip that one
        if (iblk.le.0) cycle
        iblkoff = (iblk-1)*nj

        ! remember phase for reordered indices, if necessary
        idxraw(1,ii) = 1
        if (reopr) idxraw(1,ii) = -1
        if (reoqs) idxraw(1,ii) = -idxraw(1,ii)
C        if (reopr.xor.reoqs) idxraw(1,ii) = -1
        if (sgn_change)      idxraw(1,ii) = -idxraw(1,ii)

        if (ioff.ge.0) then
          idxraw(2,ii) = ioff
        else
          ! MS(A), GAM(A)
          idxmsa  = (2-iqs-iss)/2+1
          gama    = multd2h(gam_q,gam_s)

          ! more than one distribution?
          if (me%off_op_gmox(iblk)%ndis(gama,idxmsa).gt.1) then
            off_d_gam_ms => me%off_op_gmox(iblk)%d_gam_ms
            idis = get_dis()    ! local function
            idxraw(2,ii) = off_d_gam_ms(idis,gama,idxmsa)
            cache_stat(3) = cache_stat(3)+1 ! critical miss
          else
            off_gam_ms   => me%off_op_gmo(iblk)%gam_ms
            idxraw(2,ii) = off_gam_ms(gama,idxmsa)
          end if

          ! put to cache
          idx_cache = mod(idx_cache,max_cache)+1
          if (len_cache.lt.max_cache) len_cache = len_cache+1
c          cache_tab1(idx_cache) = ihash_blk
c          cache_tab2(1,idx_cache) = iblk
c          cache_tab2(2,idx_cache) = idxraw(2,ii)
          cache_stat(2) = cache_stat(2)+1 ! miss

          cache_tab(1,idx_cache) = ihash_blk
          cache_tab(2,idx_cache) = iblk
          cache_tab(3,idx_cache) = idxraw(2,ii)

        end if

c dbg   POINT 2
c        cycle
c dbg

        ! cut into individual strings
        if (.not.pq_rs) then
          if (hpvx_p.eq.hpvx_r) then
            nind = 1
            igraph = idx_graph(hpvx_p,1,iblkoff+1)
            string(1:3,1,1) = (/ip,ir,hpvx_p/)
            string(1:3,2,1) = (/ips,irs,igraph/)
          else
            nind = 2
            igraph = idx_graph(hpvx_p,1,iblkoff+1)
            string(1:3,1,1) = (/ip,0,hpvx_p/)
            string(1:3,2,1) = (/ips,0,igraph/)
            igraph = idx_graph(hpvx_r,1,iblkoff+1)
            string(1:3,1,2) = (/ir,0,hpvx_r/)
            string(1:3,2,2) = (/irs,0,igraph/)
          end if
          if (hpvx_q.eq.hpvx_s) then
            nind = nind+1
            igraph = idx_graph(hpvx_q,2,iblkoff+nj)
            string(1:3,1,nind) = (/iq,is,hpvx_q/)
            string(1:3,2,nind) = (/iqs,iss,igraph/)
          else
            nind = nind+2
            igraph = idx_graph(hpvx_q,2,iblkoff+nj)
            string(1:3,1,nind-1) = (/iq,0,hpvx_q/)
            string(1:3,2,nind-1) = (/iqs,0,igraph/)
            igraph = idx_graph(hpvx_s,2,iblkoff+nj)
            string(1:3,1,nind) = (/is,0,hpvx_s/)
            string(1:3,2,nind) = (/iss,0,igraph/)
          end if
        else
          nind = 4
          igraph = idx_graph(hpvx_p,1,iblkoff+1)
          string(1:3,1,1) = (/ip,0,hpvx_p/)
          string(1:3,2,1) = (/ips,0,igraph/)
          igraph = idx_graph(hpvx_r,1,iblkoff+nj)
          string(1:3,1,2) = (/ir,0,hpvx_r/)
          string(1:3,2,2) = (/irs,0,igraph/)
          igraph = idx_graph(hpvx_q,2,iblkoff+1)
          string(1:3,1,nind-1) = (/iq,0,hpvx_q/)
          string(1:3,2,nind-1) = (/iqs,0,igraph/)
          igraph = idx_graph(hpvx_s,2,iblkoff+nj)
          string(1:3,1,nind) = (/is,0,hpvx_s/)
          string(1:3,2,nind) = (/iss,0,igraph/)
        end if

c dbg   POINT 3
c        cycle
c dbg

        ! get ordering array for actual sequence of strings
        idx = 0
        do igastp = 1, ngastp
          jgastp = hpvxseq(igastp)
          do iind = 1, nind
            if (string(3,1,iind).ne.jgastp) cycle
            idx = idx+1
            strseq(iind) = idx
          end do
        end do

c dbg   POINT 4
c        cycle
c dbg
        if (ntest.ge.1000) then
          write(luout,*) 'strseq = ',strseq(1:4)
          write(luout,*) 'extracted strings:'
          write(luout,
     &         '(5x,i4,",",i4," (t=",i4,")",i4,",",i4," (g=",i4,")")')
     &         string(1:3,1:2,1:nind)          
        end if


        ! put string onto list for further processing
        do iind = 1, nind
          ! search string list for identical string
          jstr = -1 
c test
          s11=string(1,1,iind)
          s21=string(2,1,iind)
          s31=string(3,1,iind)
          s12=string(1,2,iind)
          s22=string(2,2,iind)
          s32=string(3,2,iind)
c test          
c          search: do istr = 1, nstr
c            if (string_list(1,1,istr).eq.s11.and.
c     &          string_list(2,1,istr).eq.s21.and.
c     &          string_list(3,1,istr).eq.s31.and.
c     &          string_list(1,2,istr).eq.s12.and.
c     &          string_list(2,2,istr).eq.s22.and.
c     &          string_list(3,2,istr).eq.s32
cc            if (string_list(1,1,istr).eq.string(1,1,iind).and.
cc     &          string_list(2,1,istr).eq.string(2,1,iind).and.
cc     &          string_list(3,1,istr).eq.string(3,1,iind).and.
cc     &          string_list(1,2,istr).eq.string(1,2,iind).and.
cc     &          string_list(2,2,istr).eq.string(2,2,iind).and.
cc     &          string_list(3,2,istr).eq.string(3,2,iind)
c     &         ) then
c              jstr = istr
c              exit search
c            end if
c          end do search
c          jstr = search_list6(string(1,1,iind),string_list,nstr)
c          search: do istr = 1, nstr
c            if (string_list(1,1,istr).ne.s11) cycle
c            if (string_list(2,1,istr).ne.s21) cycle
c            if (string_list(3,1,istr).ne.s31) cycle
c            if (string_list(1,2,istr).ne.s12) cycle
c            if (string_list(2,2,istr).ne.s22) cycle
c            if (string_list(3,2,istr).ne.s32) cycle
cc            if (string_list(1,1,istr).eq.string(1,1,iind).and.
cc     &          string_list(2,1,istr).eq.string(2,1,iind).and.
cc     &          string_list(3,1,istr).eq.string(3,1,iind).and.
cc     &          string_list(1,2,istr).eq.string(1,2,iind).and.
cc     &          string_list(2,2,istr).eq.string(2,2,iind).and.
cc     &          string_list(3,2,istr).eq.string(3,2,iind)
cc     &         ) then
c            jstr = istr
c            exit search
c          end do search
c          jstr = 1
          istr = 1
c          search: do while (istr .lt. nstr)
          search: do while (istr .le. nstr)
c          search: do istr = 1, nstr
c            istr = istr+1
            equal =           string_list(1,1,istr).eq.s11
            equal = equal.and.string_list(2,1,istr).eq.s21
            equal = equal.and.string_list(3,1,istr).eq.s31
            equal = equal.and.string_list(1,2,istr).eq.s12
            equal = equal.and.string_list(2,2,istr).eq.s22
            equal = equal.and.string_list(3,2,istr).eq.s32
            if (equal) exit search
            istr = istr+1
          end do search

          if (istr.gt.nstr) jstr=-1
          if (istr.le.nstr) jstr=istr

          if (jstr.gt.0) then
            ! string already on list
            idxraw(2+strseq(iind),ii) = jstr
          else
            ! add string to list
            nstr = nstr+1
            string_list(1:,1:,nstr) = string(1:,1:,iind)
            idxraw(2+strseq(iind),ii) = nstr
          end if
        end do

      end do

c dbg POINT 1
c      return
c dbg

      ! post-processing for ab strings
      do istr = 1, nstr
        if (string_list(1,1,istr).ne.string_list(2,1,istr)) cycle
        if (string_list(1,2,istr).eq.string_list(2,2,istr)) cycle
        string_list(1:2,2,istr) = (/2,2/)
      end do

      if (ntest.ge.100) then
        write(luout,*) 'I need ',nstr,' strings!'
        write(luout,
     &         '(5x,i4,",",i4," (t=",i4,")",i4,",",i4," (g=",i4,")")')
     &         string_list(1:3,1:2,1:nstr)          
        write(luout,*) 'idxraw:'
        write(luout,'(5x,i2,i10,x,4i4)') idxraw(1:6,1:n_so)
      end if

      ! set up a few pointers:
      graphs => str_info%g
      gas_restr => str_info%igas_restr
      ispc_occ  => str_info%ispc_occ
      igasorb => orb_info%igasorb
      igamorb => orb_info%igamorb
      
      ioff_gas => orb_info%ioff_gas
      idx_gas => orb_info%idx_gas
      ngas_hpv => orb_info%ngas_hpv
      mostnd => orb_info%mostnd

      ! process the strings
      do istr = 1, nstr
        len = 2
        if (string_list(2,1,istr).eq.0) len = 1
        hpvx = string_list(3,1,istr)
        igraph  = string_list(3,2,istr)
        curgraph => graphs(igraph)

        idspc(1) = igasorb(string_list(1,1,istr))-ioff_gas(hpvx)
        if (len.eq.2)
     &       idspc(2) = igasorb(string_list(2,1,istr))-ioff_gas(hpvx)

        if (.not.allow_sbsp_dis(idspc,len,
     &       ngas_hpv(hpvx),
     &       gas_restr(1,1,1,1,igraph))) then
      ! ADAPT FOR OPEN-SHELL^^^
          idxlenstr(1:2,istr) = -1
          cycle
        end if

        idgam(1) = igamorb(string_list(1,1,istr))
        gamt = idgam(1)
        if (len.eq.2) then
          idgam(2) = igamorb(string_list(2,1,istr))
          gamt = multd2h(gamt,idgam(2))
        end if

c test
        strl(1:2,1:2) = string_list(1:2,1:2,istr)
c test
        idxlenstr(1,istr) =
c     &       idx4sg(len,idspc,string_list(1:,1,istr),
c     &             string_list(1:,2,istr),idgam,
c test:
     &       idx4sg(len,idspc,strl(1:,1),
     &             strl(1:,2),idgam,
     &             curgraph%y4sg,curgraph%yinf,
     &             curgraph%yssg,curgraph%wssg,
     &             curgraph%ioffstr_dgm,curgraph%ndis,
     &             mostnd(1,1,idx_gas(hpvx)),
     &             ispc_occ(igraph),nsym,
     &             ngas_hpv(hpvx))

c dbg
c        if (abs(idxlenstr(1,istr)).gt.100000) then
c          print *,'fehler: ',idxlenstr(1,istr)
c          print *,'string: ',string_list(1:len,1,istr)
c          print *,'        ',string_list(1:len,2,istr)
c          print *,'    SP  ',idspc(1:len)
c          print *,'    GM  ',idgam(1:len)
c          print *,'hpvx,igraph: ',hpvx,igraph
c          print *,'ispc_occ: ',ispc_occ(igraph)
c          print *,'idx_gas(hpvx): ',idx_gas(hpvx)
c          print *,'mostnd(1:2,1,idx_gas(hpvx)): ',
c     &         mostnd(1:2,1,idx_gas(hpvx))
c          stop 'fehler'
c        end if
c dbg
        mst = string_list(1,2,istr)
        if (len.eq.2) mst = mst+string_list(2,2,istr)
        if (mst.eq.4) mst = 0 ! paired spin case

        idxlenstr(2,istr) = curgraph%lenstr_gm(gamt,(len-mst)/2+1)
        
      end do

      if (ntest.ge.100) then
        write(luout,*) 'obtained string-addresses (-1 = not allowed): '
        write(luout,'(2x,2i8,2x,2i8,2x,2i8,2x,2i8,2x,2i8)')
     &       idxlenstr(1:2,1:nstr)
      end if

      ! assemble addresses
      do ii = 1, n_so
        isgn = idxraw(1,ii)
        idxlist(ii) = 0
        if (isgn.eq.0) cycle
        idx = idxraw(2,ii)+1
        len = 1
        ! assemble full address (note: idxlenstr(1,ii) is index-1 )
        do iind = 3, 6
          if (idxraw(iind,ii).eq.0) exit
          if (idxlenstr(1,idxraw(iind,ii)).lt.0) isgn = 0
          idx = idx + idxlenstr(1,idxraw(iind,ii))*len
          len = len * idxlenstr(2,idxraw(iind,ii))
        end do
        idxlist(ii) = idx*isgn
      end do

      if (ntest.ge.100) then
        write(luout,*) 'final indices:'
        write(luout,'(x,i12)') idxlist(1:n_so)
      end if

      return

      contains

      integer function get_dis()

      integer ::
     &     nc, na, did, idx, idx_end
      integer ::
     &     cocc(2), cims(2), cgam(2),
     &     aocc(2), aims(2), agam(2)

      integer, pointer ::
     &     didarr(:,:,:)      

      integer, external ::
c     &     idx_msgmdst2
     &     msgmdid2

      if (hpvx_p.eq.hpvx_r) then
        nc = 1
        cocc(1) = 2
        cgam(1) = multd2h(gam_p,gam_r)
        cims(1) = (2-ips-irs)/2+1
      else
        nc = 2
        cocc(1:2) = (/1,1/)
        if (hpvxblkseq(hpvx_p).lt.hpvxblkseq(hpvx_r)) then
          cgam(1:2) = (/gam_p,gam_r/)
          cims(1:2) = (/(1-ips)/2+1,(1-irs)/2+1/)
        else
          cgam(1:2) = (/gam_r,gam_p/)
          cims(1:2) = (/(1-irs)/2+1,(1-ips)/2+1/)
        end if
      end if
      if (hpvx_q.eq.hpvx_s) then
        na = 1
        aocc(1) = 2
        agam(1) = gama
        aims(1) = idxmsa
      else
        na = 2
        aocc(1:2) = (/1,1/)
        if (hpvxblkseq(hpvx_q).lt.hpvxblkseq(hpvx_s)) then
          agam(1:2) = (/gam_q,gam_s/)
          aims(1:2) = (/(1-iqs)/2+1,(1-iss)/2+1/)
        else
          agam(1:2) = (/gam_s,gam_q/)
          aims(1:2) = (/(1-iss)/2+1,(1-iqs)/2+1/)
        end if
      end if

      did = msgmdid2(cocc,cims,cgam,nc,
     &               aocc,aims,agam,na,nsym)

      didarr => me%off_op_gmox(iblk)%did
      get_dis = -1
      idx_end = me%off_op_gmox(iblk)%ndis(gama,idxmsa)

      do idx = 1, idx_end
        if (didarr(idx,gama,idxmsa).eq.did) then
          get_dis = idx
          exit
        end if
      end do

      return

      end function

      end

      pure integer function search_list6(tgt,list,len)

      implicit none

      integer, intent(in) ::
     &     len, tgt(len), list(6*len)

      integer ::
     &     idx, ii
      logical ::
     &     equal

      idx = 0
      search_list6 = -1
      do ii = 1, len
        idx = idx+1
        equal = tgt(1).eq.list(idx)
        idx = idx+1
        equal = equal.and.tgt(2).eq.list(idx)
        idx = idx+1
        equal = equal.and.tgt(3).eq.list(idx)
        idx = idx+1
        equal = equal.and.tgt(4).eq.list(idx)
        idx = idx+1
        equal = equal.and.tgt(5).eq.list(idx)
        idx = idx+1
        equal = equal.and.tgt(6).eq.list(idx)
c        if (equal) search_list6 = ii
        if (equal) exit
      end do

      return
      end
