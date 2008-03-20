*----------------------------------------------------------------------*
      subroutine topomap4contr(mode,base_tm,base_em,
     &                         topomap,eqv_map,neqv,idx_eqv,
     &                         contr,occ_vtx)
*----------------------------------------------------------------------*
*     generate a few maps needed for topology analysis
*
*     mode==1: topomap only   mode==2: all quantities
*
*     base_tm: base for setting up the packed occ_cnt; should be greater
*           than the maximum occupation appearing in the present 
*           contraction;       if -1 is provided
*           we will use an appropriate one (which might, however, not
*           be comparable to other contractions for which another
*           base might be chose automatically)
*
*     base_em: dto. for equiv. map (see below) must be larger than
*           max. block number in present contraction
*
*     topomap(nvtx,nvtx) :  who is connected to who?
*                           the index is a packed occupation, which
*                           allows to identify contractions of same
*                           type and to establish an ordering scheme
*
*     eqv_map(nvtx) :       a packed operator + block index which
*                           allows to identify similar vertices and
*                           to establish an ordering scheme
*
*     neqv(nvtx) :          number of equivalent vertices (i.e. same
*                           op+blk and commuting) for first occurrence
*                           -1 for second etc occurrence of vertex
*
*     idx_eqv(nvtx,nvtx):   idx_eqv(1:neqv(ivtx),ivtx) contains
*                           the equivalent vertices (ivtx is first
*                           occurrence of vertex)
*
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     mode, base_tm, base_em
      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,contr%nvtx)
      integer, intent(out) ::
     &     topomap(contr%nvtx,contr%nvtx),eqv_map(contr%nvtx),
     &     neqv(contr%nvtx),idx_eqv(contr%nvtx,contr%nvtx)

      logical ::
     &     ok
      integer ::
     &     nvtx, narc, ivtx, jvtx, kvtx, iarc, jarc, idx,
     &     icpack, ibase, maxblk, idx_op, iblk_op, jdx_op, jblk_op
      integer, pointer ::
     &     occ_cnt(:,:)

      integer, external ::
     &     maxblk_in_contr, ifndmax, int_pack

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'this is topomap4contr')
      end if
      
      if (mode.ne.1.and.mode.ne.2.and.mode.ne.3)
     &     call quit(1,'topomap4contr','illegal mode parameter')

      nvtx = contr%nvtx
      narc = contr%narc

      topomap = 0
      do iarc = 1, narc
        ivtx = contr%arc(iarc)%link(1)
        jvtx = contr%arc(iarc)%link(2)
        idx = iarc
        ! look wether connection is similar to a previous one
        occ_cnt => contr%arc(iarc)%occ_cnt
        ! ignore zero contractions
        if (.not.iocc_nonzero(occ_cnt)) cycle
        do jarc = 1, iarc-1
          if (iocc_equal(occ_cnt,.false.,
     &                   contr%arc(jarc)%occ_cnt,.false.)) then
            idx = jarc
            exit
          end if
        end do
        if (base_tm.lt.0) then
          ! old code, not completely understood any more :-(
          ibase = ifndmax(occ_vtx(1,1,jvtx),1,ngastp*2,1)+1
          icpack = int_pack(contr%arc(idx)%occ_cnt,ngastp*2,ibase)
          topomap(ivtx,jvtx) = icpack
          ibase = ifndmax(occ_vtx(1,1,ivtx),1,ngastp*2,1)+1
          icpack = int_pack(contr%arc(idx)%occ_cnt,ngastp*2,ibase)
          topomap(jvtx,ivtx) = icpack
        else
          icpack = int_pack(contr%arc(idx)%occ_cnt,ngastp*2,base_tm)
          topomap(ivtx,jvtx) = icpack
          topomap(jvtx,ivtx) = icpack
        end if
      end do

      if (mode.eq.1) return

      if (mode.eq.2) then
        if (base_em.lt.0) then
          maxblk = maxblk_in_contr(contr)
        else
          maxblk = base_em-1
        end if
        do ivtx = 1, nvtx
          idx_op = contr%vertex(ivtx)%idx_op
          iblk_op = contr%vertex(ivtx)%iblk_op
          idx = idx_op*(maxblk+1) + iblk_op
          eqv_map(ivtx) = idx
        end do
        return
      end if

      neqv(1:nvtx) = 1
      do ivtx = 1, nvtx
        idx_eqv(1,ivtx) = ivtx
        idx_eqv(2:nvtx,ivtx) = 0
      end do
      
      if (base_em.lt.0) then
        maxblk = maxblk_in_contr(contr)
      else
        maxblk = base_em-1
      end if
c      maxblk = maxblk_in_contr(contr)
c      maxop  = maxop_in_contr(contr)
      do ivtx = 1, nvtx
        idx_op = contr%vertex(ivtx)%idx_op
        iblk_op = contr%vertex(ivtx)%iblk_op
        idx = idx_op*(maxblk+1) + iblk_op
        eqv_map(ivtx) = idx
        do jvtx = 1, ivtx-1
          if (neqv(jvtx).lt.0) cycle
          jdx_op = contr%vertex(jvtx)%idx_op
          jblk_op = contr%vertex(jvtx)%iblk_op
          ! vertices identical ...
          if (idx_op.eq.jdx_op.and.
     &        iblk_op.eq.jblk_op) then
            ! ... and commuting? ...
            ok = .true.
            do kvtx = ivtx-1,jvtx,-1
              ok = ok.and.topomap(kvtx,ivtx).eq.0
            end do
            ! If so: equivalent!
            if (ok) then
              neqv(jvtx) = neqv(jvtx)+1
              neqv(ivtx) = -1
              idx_eqv(neqv(jvtx),jvtx) = ivtx
            end if          
          end if          
        end do
      end do

      if (ntest.ge.100) then
        write(luout,*) 'topo-map:'
        call iwrtma(topomap,nvtx,nvtx,nvtx,nvtx)
        write(luout,*) 'equiv-map:'
        call iwrtma(eqv_map,1,nvtx,1,nvtx)
        write(luout,*) 'neqv:'
        call iwrtma(neqv,1,nvtx,1,nvtx)
        write(luout,*) 'idx_eqv:'
        call iwrtma(idx_eqv,nvtx,nvtx,nvtx,nvtx)
      end if

      return
      end
