*----------------------------------------------------------------------*
      subroutine tex_contr(lutex,lencnt,first,newline,contr,op_info)
*----------------------------------------------------------------------*
*     write info on contraction in TeX style onto unit lutex
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, intent(in) ::
     &     lutex
      integer, intent(out) ::
     &     lencnt
      logical, intent(in) ::
     &     first, newline
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      type(operator), pointer ::
     &     op
      
      real(8) ::
     &     fac
      logical ::
     &     dagger
      integer ::
     &     njres, maxnj, maxset, isvtx, idxoff, nset, nj, narc, nxarc,
     &     ipos, nocc, ioff, ivtx, iarc, iarc0, ixarc, idx, ifac, p, q
      integer, pointer ::
     &     occ(:,:,:), occset(:,:,:,:), typset(:), idx0set(:,:,:,:),
     &     idxoff_cnt(:,:,:), idxoff_opn(:,:,:), occ_scr(:,:),
     &     arc_seq(:)
      character ::
     &     str*1024

      integer, external ::
     &     sign_global

c dbg
c      integer ::
c     &     nr12, nttld, nt1x, idx_t, idx_tb, idx_r, idx_tt, idx_ttb,
c     &     idx_op, iblkop
c dbg

      allocate(arc_seq(contr%narc+1))

      ! get a nice sequence of contractions (as to get "nice-looking diagrams")
      call get_nice_arc_seq(arc_seq,contr)

      fac = contr%fac*dble(sign_global(arc_seq,contr,op_info))

      op => op_info%op_arr(contr%idx_res)%op
      njres = op%njoined
      idxoff = (contr%iblk_res-1)*njres 
      occ => op%ihpvca_occ(1:ngastp,1:2,idxoff+1:idxoff+njres)

      ipos = 1
c dbg
c      nr12 = 0
c      nttld = 0
c      nt1x = 0
c      idx_t = 2
c      idx_tb = 8
c      idx_r = 15
c      idx_tt = 4
c      idx_ttb = 5
c      do ivtx = 1, contr%nvtx
c        idx_op = contr%vertex(ivtx)%idx_op
c        iblkop = contr%vertex(ivtx)%iblk_op
c        if (idx_op.eq.idx_t.and.iblkop.eq.2) nt1x = nt1x+1
c        if (idx_op.eq.idx_tb.and.iblkop.eq.2) nt1x = nt1x+1
c        if (idx_op.eq.idx_r) nr12 = nr12+1
c        if (idx_op.eq.idx_tt) nttld = nttld+1
c        if (idx_op.eq.idx_ttb) nttld = nttld+1        
c      end do
c      if (nt1x.gt.0) then
c        write(str,'("D")')
c      else if (nttld.gt.0) then
c        write(str,'("C")')
c      else if (nr12.gt.0) then
c        write(str,'("B")')
c      else
c        write(str,'("A")')
c      end if
c      ipos = 2
c dbg

      maxnj = njres
      do isvtx = 1, contr%nsupvtx
        maxnj = max(maxnj,contr%joined(0,isvtx))
      end do
      maxset = max(1,contr%narc + contr%nxarc)

      allocate(occset(ngastp,2,maxnj,maxset),
     &        idx0set(ngastp,2,maxnj,maxset),
     &         typset(maxset),occ_scr(ngastp,2))

      lencnt = 0

      if (first) then
        nset = 1
        typset(1) = 1
        occset(1:ngastp,1:2,1:njres,1) = occ(1:ngastp,1:2,1:njres)
        lencnt = lencnt + max(sum(occ(1:ngastp,1,1:njres)),
     &                        sum(occ(1:ngastp,2,1:njres))) + 3
        idx0set(1:ngastp,1:2,1:njres,1) = 0
        call tex_op(str,op%name,contr%dagger,
     &       1,typset,occset,idx0set,njres)
        ipos = len_trim(str)+1
        write(str(ipos:),'("&=")')
      else if (newline) then
        if (fac.ge.0d0) write(str(ipos:),'("\\&+")')
        if (fac.lt.0d0) write(str(ipos:),'("\\&-")')
        lencnt = lencnt+1
      else
        if (fac.ge.0d0) write(str(ipos:),'("+")')
        if (fac.lt.0d0) write(str(ipos:),'("-")')
        lencnt = lencnt+1
      end if

      if (abs(contr%fac-1d0).gt.1d-12) then
        call real2rat(p,q,contr%fac)
        ipos = len_trim(str)+1
        write(str(ipos:),'("\frac{",i6,"}{",i6,"}")') p,q
        lencnt = lencnt+2
      end if

      narc = contr%narc
      nxarc = contr%nxarc

      ! offsets for indices
      allocate(idxoff_cnt(ngastp,2,max(1,narc)),
     &         idxoff_opn(ngastp,2,max(1,nxarc)))
      
      idxoff_cnt(1:ngastp,1:2,1) = 0      
      occ_scr = 0
      do iarc0 = 1, narc-1
        iarc = arc_seq(iarc0)
        occ_scr = occ_scr + contr%arc(iarc)%occ_cnt
        idxoff_cnt(1:ngastp,1:2,iarc0+1) = occ_scr
      end do

      idxoff_opn(1:ngastp,1:2,1) = 0
      occ_scr = 0
      do ixarc = 1, nxarc-1
        occ_scr = occ_scr + contr%xarc(ixarc)%occ_cnt
        idxoff_opn(1:ngastp,1:2,ixarc+1) = occ_scr
      end do

      do isvtx = 1, contr%nsupvtx

        ! get operator info
        nj = contr%joined(0,isvtx)
        op => op_info%op_arr( contr%vertex(
     &                               contr%joined(1,isvtx) )%idx_op
     &                      ) %op
        dagger = contr%vertex(contr%joined(1,isvtx))%dagger

        ! set up index info
        do idx = 1, nj
          ivtx = contr%joined(idx,isvtx)
          ! scan arcs for contributions
          occset = 0
          idx0set = 0
          typset = 0
          nset = 0
          do iarc0 = 1, narc
            iarc = arc_seq(iarc0)
            if (contr%arc(iarc)%link(1).eq.ivtx) then
              nset = nset+1
              typset(nset) = 2
              occset(1:ngastp,1:2,idx,nset) = contr%arc(iarc)%occ_cnt
              idx0set(1:ngastp,1:2,idx,nset) =
     &                                   idxoff_cnt(1:ngastp,1:2,iarc0)
            end if
            if (contr%arc(iarc)%link(2).eq.ivtx) then
              nset = nset+1
              typset(nset) = 3
              occset(1:ngastp,1:2,idx,nset) =
     &                           iocc_dagger(contr%arc(iarc)%occ_cnt)
              idx0set(1:ngastp,1:2,idx,nset) =
     &                       iocc_dagger(idxoff_cnt(1:ngastp,1:2,iarc0))
            end if
          end do

          ! scan xarcs for contributions
          do ixarc = 1, nxarc
            if (contr%xarc(ixarc)%link(1).eq.ivtx) then
              nset = nset+1
              typset(nset) = 1
              occset(1:ngastp,1:2,idx,nset) = contr%xarc(ixarc)%occ_cnt
              idx0set(1:ngastp,1:2,idx,nset) =
     &                                   idxoff_opn(1:ngastp,1:2,ixarc)
            end if
          end do
          
        end do

        ! generate operator symbol and indices
        ipos = len_trim(str)+1
        call tex_op(str(ipos:),op%name,dagger,
     &       nset,typset,occset,idx0set,nj)

        lencnt = lencnt + 1 + 
     &       max( sum(occset(1:ngastp,1,1:nj,1:nset)),
     &            sum(occset(1:ngastp,2,1:nj,1:nset)) )

      end do

      write(lutex,'(a)') trim(str)

      deallocate(occset,idx0set,typset,occ_scr,
     &     idxoff_cnt,idxoff_opn,arc_seq)

      return
      end
