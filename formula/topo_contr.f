*----------------------------------------------------------------------*
      subroutine topo_contr(ieqvfac,resort,vtx_reo,
     &     contr,occ_vtx,fix_vtx)
*----------------------------------------------------------------------*
*     analyze the topology of the contraction contr
*     ieqvfac is the number of identical permutations of the vertex
*     sequence
*     fix_vtx means: do not consider permutations of this node for
*       evaluation of ieqvfac; note that this does not affect the
*       resort array, i.e. fixed nodes will be sorted to their canonical
*       position as well
*     resort is set true if the sequence is not in canonical order
*     vtx_reo gives the reodering information:
*       vtx_reo(idx_reo) = idx_ori
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     ieqvfac, vtx_reo(*)
      logical, intent(out) ::
     &     resort
      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,*)
      logical, intent(in) ::
     &     fix_vtx(*)

      logical ::
     &     ok, new, idag, jdag
      integer ::
     &     nvtx, ivtx, jvtx, kvtx, lvtx, narc, iarc, jarc, idx, nsame,
     &     idx_op, iblk_op, jdx_op, jblk_op, maxblk, ibase, icpack,
     &     itopo
     &     ,ieqvfac2
      
      integer, pointer ::
     &     topomap(:,:), eqv_map(:), occ_cnt(:,:), scr(:), svmap(:),
     &     neqv(:), idx_eqv(:,:), svertex(:)

      integer, external ::
     &     ifac, maxblk_in_contr, ifndmax, int_pack


c dbg
c      print *,'on entry'
c      call wrt_occ_n(luout,occ_vtx,contr%nvtx)
c dbg

      nvtx = contr%nvtx
      narc = contr%narc
      
      allocate(topomap(nvtx,nvtx),eqv_map(nvtx),scr(nvtx),
     &     neqv(nvtx),idx_eqv(nvtx,nvtx),svmap(nvtx),svertex(nvtx))

      call svmap4contr2(svmap,contr)
      svertex = contr%svertex

      topomap = 0
      do iarc = 1, narc
        ivtx = contr%arc(iarc)%link(1)
        jvtx = contr%arc(iarc)%link(2)
        idx = iarc
        ! look whether connection is similar to a previous one
        occ_cnt => contr%arc(iarc)%occ_cnt
        ! ignore zero contractions
        if (.not.iocc_nonzero(occ_cnt)) cycle
        do jarc = 1, iarc-1
c dbg
c          print *,'> ',iarc,jarc,iocc_equal(occ_cnt,.false.,
c     &                   contr%arc(jarc)%occ_cnt,.false.)
c dbg
          if (iocc_equal(occ_cnt,.false.,
     &                   contr%arc(jarc)%occ_cnt,.false.)) then
            idx = jarc
            exit
          end if
        end do
        ibase = ifndmax(occ_vtx(1,1,jvtx),1,ngastp*2,1)+1
        icpack = int_pack(contr%arc(idx)%occ_cnt,ngastp*2,ibase)
        topomap(ivtx,jvtx) = icpack
        ibase = ifndmax(occ_vtx(1,1,ivtx),1,ngastp*2,1)+1
        icpack = int_pack(contr%arc(idx)%occ_cnt,ngastp*2,ibase)
        topomap(jvtx,ivtx) = icpack
      end do

      neqv(1:nvtx) = 1
      do ivtx = 1, nvtx
        idx_eqv(1,ivtx) = ivtx
        idx_eqv(2:nvtx,ivtx) = 0
      end do
      
      maxblk = maxblk_in_contr(contr)
c      maxop  = maxop_in_contr(contr)
      do ivtx = 1, nvtx
        idx_op = contr%vertex(ivtx)%idx_op
        iblk_op = contr%vertex(ivtx)%iblk_op
        idag = contr%vertex(ivtx)%dagger
        idx = idx_op*(maxblk+1) + iblk_op
        eqv_map(ivtx) = idx
cmh     do we need to distinguish Op^+ from Op in eqv_map???
c        if (idag) eqv_map(ivtx)=eqv_map(ivtx)+1000
cmhend
        do jvtx = 1, ivtx-1
          if (neqv(jvtx).lt.0) cycle
          jdx_op = contr%vertex(jvtx)%idx_op
          jblk_op = contr%vertex(jvtx)%iblk_op
          jdag = contr%vertex(jvtx)%dagger
          ! vertices identical ...
          if (idx_op.eq.jdx_op.and.
     &        idag.eq.jdag.and.
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

      ieqvfac = 1
      ivtx_loop: do ivtx = 1, nvtx
        if (fix_vtx(ivtx).or.neqv(ivtx).lt.0) cycle
        nsame = 1
        kvtx = ivtx
        do jvtx = 2, neqv(ivtx)
          lvtx = idx_eqv(jvtx,ivtx)
          itopo = topo_cmp2(topomap(1:nvtx,kvtx),topomap(1:nvtx,lvtx),
     &         eqv_map,nvtx)
          kvtx = lvtx
c dbg
c          print *,'comparing'
c          print *,'1: ',topomap(1:nvtx,ivtx)
c          print *,'2: ',topomap(1:nvtx,kvtx)
c          print *,'result -> ',itopo
c dbg
          if (itopo.lt.0) then
            ieqvfac = ieqvfac*ifac(nsame)
            nsame = 1
          else if (itopo.eq.0) then
            nsame = nsame+1
          else if (itopo.gt.0) then
            ieqvfac = -1
            exit ivtx_loop
          end if
        end do
        ieqvfac = ieqvfac*ifac(nsame)
      end do ivtx_loop

      if (ntest.ge.100) then
        write(luout,*) '-> fac = ',ieqvfac
      end if

c      ! look for equivalent tupels of lines in the topo-map:
c      ieqvfac2 = 1
c      do ivtx = 1, nvtx
c        if (fix_vtx(ivtx).or.topomap(1,ivtx).lt.0) cycle
c        nsame = 1
c        do jvtx = ivtx+1, nvtx
c          if (eqv_map(ivtx).eq.eqv_map(jvtx)) then
c            ! compare lines including permutations that leave
c            ! eqv_map invariant ...
c            if (topo_cmp(topomap(1:nvtx,ivtx),topomap(1:nvtx,jvtx),
c     &          eqv_map,nvtx)) then
c              nsame = nsame+1
c              topomap(1,jvtx) = -1 ! ignore next time
c            end if
c          end if
c        end do
c        ieqvfac2 = ieqvfac2*ifac(nsame)
c      end do
c
c      if (ntest.ge.100) then
c        write(luout,*) '-> fac = ',ieqvfac2
c        if (ieqvfac2.ne.ieqvfac) write(luout,*) 'they differ!'
c      end if

      ! suggest a reordering array for operator sequence:
      ! the sequence of operators may be changed if a zero
      ! occurs on the first off-diagonal band of topomap
      ! ensure that the numbering on eqv_map is ascending for
      ! these operators
      do ivtx = 1, nvtx
        vtx_reo(ivtx) = ivtx
      end do
      resort = .false.
      new = .true.
      if (new) then
      ! loop over vertices
      do jvtx = 2, nvtx
c dbg
c        print *,'jvtx = ',jvtx
c        print *,'present sequence: ',eqv_map(1:nvtx)
c        print *,'present reo     : ',vtx_reo(1:nvtx)
c        print *,'present svmap   : ',svmap(1:nvtx)
c        print *,'present svertex : ',svertex(1:nvtx)
c dbg
        ivtx = jvtx-1
        ! maximum "upper" position to that we may shift vertex(jvtx):
        do while (ivtx.gt.0.and.topomap(ivtx,jvtx).eq.0.and.! must commute
     &            (svmap(ivtx).eq.0.or.svmap(jvtx).eq.0.or. ! and not interchange
     &             svmap(ivtx).eq.svmap(jvtx)).and.         !  external lines
     &             svertex(ivtx).ne.svertex(jvtx))
          ivtx = ivtx-1
c dbg fix by mh
          if (ivtx.eq.0) exit
c dbg end of fix


        end do
        ivtx = ivtx+1 ! actual "upper" position is ivtx+1
c dbg
c        print *,'upper = ',ivtx
c dbg
        ! now find the uppermost postion such that vertex(jvtx) has lower
        ! rank than the succeding vertex
        do while(ivtx.lt.jvtx .and. 
     &           ( eqv_map(jvtx).gt.eqv_map(ivtx) .or.
     &            (eqv_map(jvtx).eq.eqv_map(ivtx).and.
     &             topo_cmp2(topomap(1:nvtx,jvtx),
     &                       topomap(1:nvtx,ivtx),
     &                       eqv_map,nvtx).ge.0
     &            )
     &           )
     &          )
          ivtx = ivtx+1
        end do
c dbg
c        print *,'insert at ',ivtx
c dbg

        ! insert jvtx at position ivtx
        if (ivtx.ne.jvtx) then
          resort = .true.
          call shift_imat(topomap,jvtx,ivtx,nvtx)
          call shift_ivec(eqv_map,jvtx,ivtx,nvtx)
          call shift_ivec(vtx_reo,jvtx,ivtx,nvtx)
          call shift_ivec(svmap,jvtx,ivtx,nvtx)
          call shift_ivec(svertex,jvtx,ivtx,nvtx)
        end if

      end do
c dbg
c      print *,'2nd round'
c dbg
      do jvtx = nvtx-1, 1, -1
c dbg
c        print *,'jvtx = ',jvtx
c        print *,'present sequence: ',eqv_map(1:nvtx)
c        print *,'present reo     : ',vtx_reo(1:nvtx)
c        print *,'present svmap   : ',svmap(1:nvtx)
c        print *,'present svertex : ',svertex(1:nvtx)
c dbg
        ivtx = jvtx+1
        ! maximum "upper" position to that we may shift vertex(jvtx):
        do while (ivtx.le.nvtx.and.topomap(ivtx,jvtx).eq.0.and.! must commute
     &            (svmap(ivtx).eq.0.or.svmap(jvtx).eq.0.or. ! and not interchange
     &             svmap(ivtx).eq.svmap(jvtx)).and.         !  external lines
     &             svertex(ivtx).ne.svertex(jvtx))
          ivtx = ivtx+1
c dbg fix by mh
          if (ivtx.gt.nvtx) exit
c dbg end fix
        end do
        ivtx = ivtx-1 ! actual "lower" position is ivtx-1
c dbg
c        print *,'lower = ',ivtx
c dbg
        ! now find the lowermost postion such that vertex(jvtx) has higher
        ! rank than the preceding vertex
        do while(ivtx.gt.jvtx .and. 
     &           ( eqv_map(jvtx).lt.eqv_map(ivtx) .or.
     &            (eqv_map(jvtx).eq.eqv_map(ivtx).and.
     &             topo_cmp2(topomap(1:nvtx,jvtx),
     &                       topomap(1:nvtx,ivtx),
     &                       eqv_map,nvtx).lt.0
     &            )
     &           )
     &          )
          ivtx = ivtx-1
        end do
c dbg
c        print *,'insert at ',ivtx
c dbg

        ! insert jvtx at position ivtx
        if (ivtx.ne.jvtx) then
          resort = .true.
          call shift_imat(topomap,jvtx,ivtx,nvtx)
          call shift_ivec(eqv_map,jvtx,ivtx,nvtx)
          call shift_ivec(vtx_reo,jvtx,ivtx,nvtx)
          call shift_ivec(svmap,jvtx,ivtx,nvtx)
          call shift_ivec(svertex,jvtx,ivtx,nvtx)
        end if

      end do
      ok = .true.

      else
      ! looks like bubble sort:
      ! the sort is restricted to certain transpositions, so this
      ! is the most error-save variant
      resort = .false.
      do jvtx = 1, nvtx
        ok = .true.
        do ivtx = 1, nvtx-1
          if (topomap(ivtx,ivtx+1).eq.0.and.
     &        (eqv_map(ivtx).gt.eqv_map(ivtx+1).or.
     &         eqv_map(ivtx).eq.eqv_map(ivtx+1).and.
     &         topo_cmp2(topomap(1:nvtx,ivtx),topomap(1:nvtx,ivtx+1),
     &         eqv_map,nvtx).gt.0) .and.
     &        (svmap(ivtx).eq.0.or.svmap(ivtx+1).eq.0.or.
     &         svmap(ivtx).eq.svmap(ivtx+1)) .and.
     &         svertex(ivtx).ne.svertex(ivtx+1)) then
            ok = .false.
            resort = .true.
            ! interchange the two indices
            ! (a) in topomap
            scr(1:nvtx) = topomap(1:nvtx,ivtx)
            topomap(1:nvtx,ivtx) = topomap(1:nvtx,ivtx+1)
            topomap(1:nvtx,ivtx+1) = scr(1:nvtx)
            scr(1:nvtx) = topomap(ivtx,1:nvtx)
            topomap(ivtx,1:nvtx) = topomap(ivtx+1,1:nvtx)
            topomap(ivtx+1,1:nvtx) = scr(1:nvtx)
            ! (b) in eqv_map
            scr(1) = eqv_map(ivtx)
            eqv_map(ivtx) = eqv_map(ivtx+1)
            eqv_map(ivtx+1) = scr(1)
            ! (c) in vtx_reo
            scr(1) = vtx_reo(ivtx)
            vtx_reo(ivtx) = vtx_reo(ivtx+1)
            vtx_reo(ivtx+1) = scr(1)
            ! (d) in svmap
            scr(1) = svmap(ivtx)
            svmap(ivtx) = svmap(ivtx+1)
            svmap(ivtx+1) = scr(1)
          end if
        end do
        if (ok) exit
      end do
      end if

      if (.not.ok)
     &     call quit(1,'topo_contr','restricted sort in problems')

      if (ntest.ge.100) then
        write(luout,*) 'suggested reordering of operators'
        write(luout,*) vtx_reo(1:nvtx)
      end if

      deallocate(topomap,eqv_map,scr,neqv,idx_eqv)

      return

      contains

      logical function topo_cmp2(top1,top2,eqv,nel)

      implicit none

      integer, intent(in) ::
     &     nel, top1(nel), top2(nel), eqv(nel)
c      type(cntr_arc), intent(in) ::
c     &     arc(*)

      integer ::
     &     iel, jel, ieqv, npick, nconn1, nconn2
      integer ::
     &     pick1(nel), pick2(nel), eqv_cp(nel)

      integer, external ::
     &     imltlist

      ! number of vertices connected to
      nconn1 = nvtx-imltlist(0,top1,nvtx,1)
      nconn2 = nvtx-imltlist(0,top2,nvtx,1)

      if (nconn1.gt.nconn2) then
        topo_cmp2 = -1
        return
      else if (nconn1.lt.nconn2) then
        topo_cmp2 = +1
        return
      end if

      eqv_cp = eqv ! get a copy
      do iel = 1, nel
        ! already processed this one?
        if (eqv_cp(iel).le.0) cycle
        ! current eqv number is:
        ieqv = eqv_cp(iel)
        ! pick all elements from eqv places:
        npick = 1
        pick1(1) = top1(iel)
        pick2(1) = top2(iel)
        do jel = iel+1, nel
          if (eqv_cp(jel).eq.ieqv) then
            eqv_cp(jel) = -1 ! mark this element as "already processed"
            npick = npick+1
            pick1(npick) = top1(jel)
            pick2(npick) = top2(jel)
          end if
        end do
        ! compare
        if (npick.eq.1) then
          if (pick1(1).lt.pick2(1)) then
            topo_cmp2 = -1
          else if (pick1(1).eq.pick2(1)) then
            topo_cmp2 = 0
          else
            topo_cmp2 = +1
          end if
        else
          ! sort the two lists
          call isort(pick1,npick,+1)
          call isort(pick2,npick,+1)
          topo_cmp2 = 0
          do jel = 1, npick
            if (pick1(jel).lt.pick2(jel)) then
              topo_cmp2 = -1
              exit
            else if (pick1(jel).gt.pick2(jel)) then
              topo_cmp2 = +1
              exit
            end if
          end do
        end if
        ! if different, we can leave the loop
        if (topo_cmp2.ne.0) exit
      end do

      return
      end function

      subroutine shift_ivec(vec,idx1,idx2,len)

      implicit none

      integer, intent(in) ::
     &     len, idx1, idx2
      integer, intent(inout) ::
     &     vec(len)

      integer ::
     &     ihelp, idx, inc

      inc = +1
      if (idx1.gt.idx2) inc = -1

      ihelp = vec(idx1)
      do idx = idx1, idx2-inc, inc
        vec(idx) = vec(idx+inc)
      end do
      vec(idx2) = ihelp

      return
      end subroutine

      subroutine shift_imat(mat,idx1,idx2,len)

      implicit none

      integer, intent(in) ::
     &     len, idx1, idx2
      integer, intent(inout) ::
     &     mat(len,len)

      integer ::
     &     ihelp(len), idx, inc

      inc = +1
      if (idx1.gt.idx2) inc = -1

      ! shift rows
      ihelp(1:len) = mat(idx1,1:len)
      do idx = idx1, idx2-inc, inc
        mat(idx,1:len) = mat(idx+inc,1:len)
      end do
      mat(idx2,1:len) = ihelp(1:len)
      ! shift columns
      ihelp(1:len) = mat(1:len,idx1)
      do idx = idx1, idx2-inc, inc
        mat(1:len,idx) = mat(1:len,idx+inc)
      end do
      mat(1:len,idx2) = ihelp(1:len)

      return
      end subroutine

c      logical function topo_cmp(top1,top2,eqv,nel)
c
c      implicit none
c
c      integer, intent(in) ::
c     &     nel, top1(nel), top2(nel), eqv(nel)
c
c      integer ::
c     &     iel, jel, ieqv, npick
c      integer ::
c     &     pick1(nel), pick2(nel), eqv_cp(nel)
c
c      topo_cmp = .true.
c      eqv_cp = eqv ! get a copy
c      do iel = 1, nel
c        ! already processed this one?
c        if (eqv_cp(iel).le.0) cycle
c        ! current eqv number is:
c        ieqv = eqv_cp(iel)
c        ! pick all elements from eqv places:
c        npick = 1
c        pick1(1) = top1(iel)
c        pick2(1) = top2(iel)
c        do jel = iel+1, nel
c          if (eqv_cp(jel).eq.ieqv) then
c            eqv_cp(jel) = -1 ! mark this element as "already processed"
c            npick = npick+1
c            pick1(npick) = top1(jel)
c            pick2(npick) = top2(jel)
c          end if
c        end do
c        ! compare
c        if (npick.eq.1) then
c          topo_cmp = topo_cmp.and.pick1(1).eq.pick2(1)
c        else
c          ! sort the two lists
c          call isort(pick1,npick,+1)
c          call isort(pick2,npick,+1)
c          do jel = 1, npick
c            topo_cmp = topo_cmp.and.pick1(jel).eq.pick2(jel)
c          end do
c        end if
c        ! if different, we can leave the loop
c        if (.not.topo_cmp) exit
c      end do
c      
c      return
c      end function

      end
