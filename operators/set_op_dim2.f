*----------------------------------------------------------------------*
      subroutine set_op_dim2(ipass,mel,str_info,ngam)
*----------------------------------------------------------------------*
*
*     set up the dimension arrays for ME-list mel operator mel%op
*
*     new version: general intermediates added
*
*     the operator has a total Ms and IRREP, as given in "mel".
*     it is stored as
*
*          Op(Cx,Ax,Cp,Ap,Ch,Ah,Cv,Av)
*
*     where C are indices associated with creation strings and
*           A are indices associated with annihilation strings
*     in (p)article, (h)ole, and (v)alence spaces
*
*     for intermediates: Cx,Ax,etc. may break down in several blocks
*     
*     in this routine we define the sequence in which the different blocks
*     of the operator are to be stored; with "different blocks" we mean
*     all the different cases of Ms and symmetry labels distributed among
*     the strings (Cp,Ch,Cv,Ap,Ah,Av) such that they add up to the 
*     required total Ms and IRREP
*
*     ipass = 1   set len_op_gmo, len_op_occ, 
*                     off_op_gmo, off_op_occ, len_op
*                 and get dimension info for off_op_gmox
*
*     ipass = 2   set off_op_gmox
*
*     the len_op_gmox is only needed for operators, where more than
*     one class of spaces (I mean H/P/V) is occupied for C or for A
*     usual H->P excitation operators (e.g. T-operators of SR-CC) do
*     not fall into this group and ipass=2 may be skipped
*
*     andreas, end of 2006
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'hpvxseq.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 000
      
      integer, intent(in) ::
     &     ipass, ngam
      type(strinf), intent(in) ::
     &     str_info
      type(me_list), intent(inout) ::
     &     mel

      logical ::
     &     first, ms_fix, fix_success
      integer ::
     &     idxstr, idxstr_tot, idxdis, iblk, iblkoff, nexc,
     &     msa_max, msc_max, njoined, nblk,
     &     msa, msc, idxmsa, idxmsa2, igama, igamc,
     &     nasub, ncsub, icmp,
     &     did, iexc, igam, len_blk, ld_blk, idx, jdx, tot_c, tot_a

      integer, pointer ::
     &     hpvx_csub(:), hpvx_asub(:),
     &     occ_csub(:), occ_asub(:),
     &     graph_csub(:), graph_asub(:),
     &     msdis_c(:),  msdis_a(:),
     &     idxmsdis_c(:),  idxmsdis_a(:),
     &     gamdis_c(:), gamdis_a(:),
     &     len_str(:)
      integer, pointer ::
     &     hpvx_occ(:,:,:), idx_graph(:,:,:),
     &     ca_occ(:,:)
      type(graph), pointer ::
     &     graphs(:)
      type(operator), pointer ::
     &     op

      integer ::
c     &     msd(ngastp,2), igamd(ngastp,2)
     &     msd(ngastp,2,mel%op%njoined), igamd(ngastp,2,mel%op%njoined)

      logical, external ::
     &     next_msgamdist2
      integer, external ::
     &     msgmdid2, msgmdid, msa2idxms4op
      
      op => mel%op

      if (ntest.gt.5) then
        call write_title(luout,wst_dbg_subr,'set_op_dim')
        write(luout,*) ' ipass = ',ipass
        write(luout,*) ' ME-list  = ',trim(mel%label)
        write(luout,*) ' operator = ',trim(op%name)
        write(luout,*) ' IRREP    = ',mel%gamt
        write(luout,*) ' Ms       = ',mel%mst
      end if

      idxstr = 0
      idxstr_tot = 0
      njoined = op%njoined
      nblk = op%n_occ_cls
      ms_fix = mel%fix_vertex_ms
      hpvx_occ => op%ihpvca_occ
      idx_graph => mel%idx_graph
      ca_occ => op%ica_occ
      graphs => str_info%g

c dbg
c      print *,'set dim, fix =',ms_fix
c dbg

      ! we better initialize some of the key arrays
      if (ipass.eq.1) then
        mel%len_op_occ(1:nblk) = 0
        mel%off_op_occ(1:nblk) = 0
        do iblk = 1, nblk
          mel%off_op_gmox(iblk)%maxd = 0
          mel%len_op_gmo(iblk)%gam_ms = 0
          mel%off_op_gmo(iblk)%gam_ms = 0
        end do
      end if
      if (ipass.eq.2) then
        do iblk = 1, nblk
          mel%len_op_gmox(iblk)%d_gam_ms = 0
          mel%off_op_gmox(iblk)%d_gam_ms = 0
          mel%ld_op_gmox(iblk)%d_gam_ms = 0
        end do
      end if

      ! loop over occupations (= blocks)
      occ_cls: do iblk = 1, nblk

        if (ntest.ge.100.and.op%formal_blk(iblk))
     &     write(luout,*) 'skipping formal block (#',iblk,')'
        if (op%formal_blk(iblk)) cycle

        iblkoff = (iblk-1)*njoined
        if (ntest.ge.150) then
          write(luout,*) 'class: ',iblk
          call wrt_occ_n(luout,hpvx_occ(1,1,iblkoff+1),njoined)
        end if
        
        ! find the number sub-blocks for C and A
        call get_num_subblk(ncsub,nasub,hpvx_occ(1,1,iblkoff+1),njoined)

        ! HPVX-info, OCC-info, 
        ! MS/IRREP-distributions, and string lengths
        allocate(hpvx_csub(ncsub),hpvx_asub(nasub),
     &           occ_csub(ncsub), occ_asub(nasub),
     &           graph_csub(ncsub), graph_asub(nasub),
     &           msdis_c(ncsub),  msdis_a(nasub),
     &           idxmsdis_c(ncsub),  idxmsdis_a(nasub),
     &           gamdis_c(ncsub), gamdis_a(nasub),
     &           len_str(ncsub+nasub))

        msc_max = ca_occ(1,iblk)
        msa_max = ca_occ(2,iblk)

        ! set HPVX and OCC info
        call condense_occ(occ_csub, occ_asub,
     &                    hpvx_csub,hpvx_asub,
     &                    hpvx_occ(1,1,iblkoff+1),njoined,hpvxblkseq)
        ! do the same for the graph info
        call condense_occ(graph_csub, graph_asub,
     &                    hpvx_csub,hpvx_asub,
     &                    idx_graph(1,1,iblkoff+1),njoined,hpvxblkseq)

c dbg
c        print *,'graph_csub: ',graph_csub(1:ncsub)
c        print *,'graph_asub: ',graph_asub(1:nasub)
c        print *,'hpvx_csub: ',hpvx_csub(1:ncsub)
c        print *,'hpvx_asub: ',hpvx_asub(1:nasub)
c dbg
        ! set offsets
        mel%off_op_occ(iblk) = idxstr

        if (ipass.eq.1) then
          mel%off_op_gmox(iblk)%maxd = 0
        else
          mel%off_op_gmox(iblk)%
     &       d_gam_ms(1:mel%off_op_gmox(iblk)%maxd,1:ngam,1:msa_max)=-1
          mel%off_op_gmox(iblk)%
     &       did(1:mel%off_op_gmox(iblk)%maxd,1:ngam,1:msa_max) = 0
          mel%off_op_gmox(iblk)%
     &       ndis(1:ngam,1:msa_max) = 0
        end if

        ! loop over Ms of A-string (fixes Ms of C-string)
        idxmsa = 0
        msa_loop: do msa = msa_max, -msa_max, -2

          ! C <-> A  means alpha <-> beta !!
          msc = msa + mel%mst
          ! to make this point clear:
          ! usually, we have mst=0 operators, which implies
          ! Ms(C) == Ms(A)

          if (abs(msc).gt.msc_max) cycle msa_loop
          idxmsa = idxmsa+1

          ! test indexing routine
          idxmsa2 = msa2idxms4op(msa,mel%mst,msa_max,msc_max)

          if (idxmsa.ne.idxmsa2)
     &         call quit(1,'set_op_dim2','bug in msa2idxms4op!')
c dbg
c          print *,'msa, msc: ',msa,msc
c dbg

          ! loop over IRREP of A-string (fixes IRREP of C-string)
          igama_loop: do igama = 1, ngam
            igamc = multd2h(igama,mel%gamt)
            
            ! store the current position in offset array
            mel%off_op_gmo(iblk)%gam_ms(igama,idxmsa) = idxstr

            ! now, to be general, we have to loop over all
            ! possible MS and IRREP distributions over X/H/P/V spaces
            ! for both C and A strings 
            ! (for intermediates, njoined>1, these are further 
            ! subdivided into substrings, so even more work to come)

            idxdis = 0            

            first = .true.
            distr_loop: do
                  
              if (.not.next_msgamdist2(first,
     &            msdis_c,msdis_a,gamdis_c,gamdis_a,
     &            ncsub, nasub,
     &            occ_csub,occ_asub,
     &            msc,msa,igamc,igama,ngam,
     &            ms_fix,fix_success))
     &           exit distr_loop
c              if(ms_fix.and..not.fix_success)cycle distr_loop
c dbg
c              print *,'top of dist_loop'
c dbg
              first = .false.

              call ms2idxms(idxmsdis_c,msdis_c,occ_csub,ncsub)
              call ms2idxms(idxmsdis_a,msdis_a,occ_asub,nasub)

c dbg fix by mh ???
              if (associated(graphs)) then
              print *,'assograph!'
c dbg original
              call set_len_str(len_str,ncsub,nasub,
     &                         graphs,
     &                         graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                         graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                         hpvxseq,.false.)
c dbg resume fix
              end if
c dbg end fix

c dbg
c              print *,'current dis:'
c              write(luout,*) idxmsdis_c(1:ncsub)
c              write(luout,*) gamdis_c(1:ncsub)
c              write(luout,*) idxmsdis_a(1:nasub)
c              write(luout,*) gamdis_a(1:nasub)
c              print *,'graphs c:',graph_csub(1:ncsub)
c              print *,'graphs a:',graph_asub(1:nasub)
c              print *,'len_str: ',len_str(1:ncsub+nasub)
c dbg
              ld_blk = 1
              do icmp = 1, ncsub
                ld_blk = ld_blk*len_str(icmp)
c dbg
c                print *,'icmp, len_str: ',icmp,len_str(icmp)
c dbg
              end do
              len_blk = ld_blk
              do icmp = ncsub+1, ncsub+nasub
                len_blk = len_blk*len_str(icmp)
c dbg
c                print *,'icmp, len_str: ',icmp,len_str(icmp)
c dbg
              end do

              if (len_blk.le.0) cycle distr_loop
              
              ! increment distribution index
              idxdis = idxdis+1
              
              if (njoined.eq.1.and.ntest.ge.150.or.
     &             (ms_fix.and.ipass.eq.2)) then
                if(.not.ms_fix)
     &               write(luout,*) 'current MS and IRREP distr:'
                call expand_occ(msd,idx_graph(1,1,iblkoff+1),
     &                    ncsub,nasub,
     &                    msdis_c,msdis_a,
     &                    hpvx_csub,hpvx_asub,
     &                    njoined)
                if(.not.ms_fix)
     &               call wrt_occ_n(luout,msd,njoined)
                call expand_occ(igamd,idx_graph(1,1,iblkoff+1),
     &                    ncsub,nasub,
     &                    gamdis_c,gamdis_a,
     &                    hpvx_csub,hpvx_asub,
     &                    njoined)
                if(.not.ms_fix)
     &               call wrt_occ_n(luout,igamd,njoined)
                if (njoined.eq.1) then
                  did = msgmdid(hpvx_occ(1,1,iblkoff+1),
     &                        msd,igamd,ngam)
                  write(luout,*) 'DID old: ',did
                end if
                write(luout,*) 'current idxdis = ',idxdis
              end if
c dbg
              if(ms_fix)then
                do idx = 1, njoined
                  tot_c = 0
                  tot_a = 0
                  do jdx = 1, ngastp
                    tot_c = tot_c + msd(jdx,1,idx)
                    tot_a = tot_a + msd(jdx,2,idx)
                  enddo
                  if(tot_c.ne.tot_a) cycle distr_loop
                enddo
              endif
c dbg

              ! save current offset
              if (ipass.eq.2) then
                mel%off_op_gmox(iblk)%
     &               d_gam_ms(idxdis,igama,idxmsa)=idxstr
                ! get ID of current distr
                did = msgmdid2(occ_csub,idxmsdis_c,gamdis_c,ncsub,
     &                         occ_asub,idxmsdis_a,gamdis_a,nasub,ngam)

c dbg
c                if(trim(mel%op%name).eq.'G_Z')
c     &               print *,'did',did
c dbg
                ! save ID of current distr
                mel%off_op_gmox(iblk)%
     &               did(idxdis,igama,idxmsa) = did
                if (ntest.ge.150) then
                  write(luout,*) 'current did = ',did
                  write(luout,*) idxmsdis_c(1:ncsub)
                  write(luout,*) gamdis_c(1:ncsub)
                  write(luout,*) idxmsdis_a(1:nasub)
                  write(luout,*) gamdis_a(1:nasub)
                end if
              end if 

              ! increment string element index
              idxstr = idxstr+len_blk
              idxstr_tot = idxstr_tot+len_blk

              if (ipass.eq.2) then
                mel%len_op_gmox(iblk)%
     &               d_gam_ms(idxdis,igama,idxmsa) = len_blk
                mel%ld_op_gmox(iblk)%
     &               d_gam_ms(idxdis,igama,idxmsa) = ld_blk
              end if

              if (ntest.ge.150) then
                write(luout,*) 'current block length: ',len_blk
              end if
              
            end do distr_loop

            mel%len_op_gmo(iblk)%gam_ms(igama,idxmsa) = idxstr -
     &           mel%off_op_gmo(iblk)%gam_ms(igama,idxmsa)

            mel%off_op_gmox(iblk)%maxd =
     &           max(mel%off_op_gmox(iblk)%maxd,idxdis)

            if (ipass.eq.2) then
              mel%off_op_gmox(iblk)%ndis(igama,idxmsa) = idxdis
            end if

          end do igama_loop

        end do msa_loop

        mel%len_op_occ(iblk) = idxstr - mel%off_op_occ(iblk)

        deallocate(hpvx_csub,hpvx_asub,
     &           occ_csub, occ_asub,
     &           graph_csub, graph_asub,
     &           msdis_c,  msdis_a,
     &           idxmsdis_c,  idxmsdis_a,
     &           gamdis_c, gamdis_a,
     &           len_str)

      end do occ_cls

      mel%len_op = idxstr_tot

      if (ntest.ge.100) then
        if (ipass.eq.1) then
          write(luout,*) 'total number of operator elements: ',
     &         mel%len_op
          write(luout,*) 'length per occupation class:'
          call iwrtma(mel%len_op_occ,nblk,1,nblk,1)
          write(luout,*) 'offsets per occupation class:'
          call iwrtma(mel%off_op_occ,nblk,1,nblk,1)
          write(luout,*) 'info per occupation class, IRREP, MS:'
          do iblk = 1, nblk
            if (op%formal_blk(iblk)) cycle
            nexc = min(ca_occ(1,iblk),
     &                 ca_occ(2,iblk))
            write(luout,*) 'occ-class: ',iblk
            write(luout,*) 'lengths:'
            call iwrtma(mel%len_op_gmo(iblk)%gam_ms,
     &           ngam,nexc+1,ngam,nexc+1)
            write(luout,*) 'offsets:'
            call iwrtma(mel%off_op_gmo(iblk)%gam_ms,
     &           ngam,nexc+1,ngam,nexc+1)
          end do
        else
          write(luout,*) 'info per occupation class, DISTR, IRREP, MS:'
          write(luout,*) 'offsets:'
          do iblk = 1, nblk
            if (op%formal_blk(iblk)) cycle
            nexc = min(ca_occ(1,iblk),
     &                 ca_occ(2,iblk))
            write(luout,*) 'occ-class: ',iblk
            iblkoff = (iblk-1)*njoined
            call wrt_occ_n(luout,op%ihpvca_occ(1,1,iblkoff+1),njoined)
            do iexc = 1, nexc+1
              do igam = 1, ngam
                if (mel%off_op_gmox(iblk)%ndis(igam,iexc).eq.0) cycle
                write(luout,*) iexc,igam,' -> ',
     &               mel%off_op_gmox(iblk)%
     &               d_gam_ms(1:mel%off_op_gmox(iblk)%
     &               ndis(igam,iexc),igam,iexc)
              end do
            end do            
          end do
          write(luout,*) 'distribution IDs:'
          do iblk = 1, nblk
            if (op%formal_blk(iblk)) cycle
            nexc = min(ca_occ(1,iblk),
     &                 ca_occ(2,iblk))
            write(luout,*) 'occ-class: ',iblk
            iblkoff = (iblk-1)*njoined
            call wrt_occ_n(luout,op%ihpvca_occ(1,1,iblkoff+1),njoined)
            do iexc = 1, nexc+1
              do igam = 1, ngam
                if (mel%off_op_gmox(iblk)%ndis(igam,iexc).eq.0) cycle
                write(luout,*) iexc,igam,' -> ',
     &               mel%off_op_gmox(iblk)%
     &               did(1:mel%off_op_gmox(iblk)%
     &               ndis(igam,iexc),igam,iexc)
              end do
            end do
          end do

        end if
        
      end if

      return
      end
*----------------------------------------------------------------------*
