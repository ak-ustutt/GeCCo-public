*----------------------------------------------------------------------*
      subroutine reorder_supvtx_x(possible,
     &     modify_contr,set_reord_list,reo_info,
     &     contr,occ_vtx,idxop12)
*----------------------------------------------------------------------*
*     analoguous to reorder_supvtx_x but now for external arcs
*     (awaiting a more elegant general routine)
*
*     check whether an external vertex is linked to more than one 
*     component of a super vertex
*     if this is the case, we may shift the involved occupations to
*     the super-vertex closest to the vertex they are commonly 
*     contracted to, and we may add up the contraction
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_reorder_info.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      logical, intent(out) ::
     &     possible
      logical, intent(in) ::
     &     modify_contr, set_reord_list
      integer, intent(in) ::
     &     idxop12
      type(contraction), intent(inout) ::
     &     contr
      type(reorder_info), intent(inout) ::
     &     reo_info
      integer, intent(inout) ::
     &     occ_vtx(ngastp,2,contr%nvtx)


      type(cntr_arc), pointer ::
     &     xarc(:), xarc_scr(:)
      type(cntr_vtx), pointer ::
     &     vertex(:)
      integer, pointer ::
     &     svertex(:) !, occ_shift(:,:)
      integer ::
     &     occ_shr(ngastp,2), occ_shl(ngastp,2),
     &     cnt_shr(ngastp,2), cnt_shl(ngastp,2)
      integer ::
     &     nxarc, ixarc, jxarc, ivtx1, ivtx2, 
     &     ixarc_shift, ixarc_sum, nxarc_new,
     &     maxreo, idx, idxsuper, ica, ica_vtx, hpvx

      integer, external ::
     &     imltlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'reorder_supvtx_x')
        write(luout,*) 'modify_contr: ',modify_contr
        write(luout,*) 'set_reord_list: ',set_reord_list
        write(luout,*) 'idxop12: ',idxop12
        write(luout,*) 'contr on input:'
        call prt_contr3(luout,contr,occ_vtx)
      end if

      possible = .true.

      if (set_reord_list) then
        maxreo = 2*(contr%nvtx - contr%nsupvtx)
c        reo_info%nreo = 0
c        allocate(reo_info%reo(maxreo))
        if (.not.associated(reo_info%reo))
     &       call quit(1,'reorder_supvtx_x',
     &       'call reorder_supvtx first')
      end if

      nxarc = contr%nxarc
      if (modify_contr) then
        xarc => contr%xarc
      else
        allocate(xarc_scr(contr%nxarc))
        xarc_scr = contr%xarc
        xarc => xarc_scr
      end if
      svertex => contr%svertex
      vertex =>  contr%vertex
      ! loop over xarcs
      do ixarc = 1, nxarc
        ! deleted xarc?
        if (xarc(ixarc)%link(1).le.0) cycle
        ! loop over xarcs and look for identical link(1) (if jxarc<ixarc)
        !                      or for identical link(2) (if jxarc>ixarc)
        do jxarc = 1, nxarc
          ! deleted xarc? or the xarc itself?
          if (xarc(ixarc)%link(1).le.0.or.ixarc.eq.jxarc) cycle
c dbg
c          print *,'ixarc, jxarc: ',ixarc,jxarc
c dbg
          
          ! same external ?
          if (xarc(ixarc)%link(2).ne.xarc(jxarc)%link(2)) cycle

          ! from which vertex?
          ivtx1 = xarc(ixarc)%link(1)
          ivtx2 = xarc(jxarc)%link(1)

c dbg
c          print *,'ivtx1,2: ',ivtx1,ivtx2
c dbg
          ! same super vertex?
          if (svertex(ivtx1).ne.svertex(ivtx2)) cycle
c dbg
c          print *,'survived b)'
c dbg

          occ_shr = 0
          occ_shl = 0
          cnt_shr = 0
          cnt_shl = 0
          do ica = 1, 2
            ! which C/A shall we look at in the associated vertices?
            ica_vtx = ica
            do hpvx = 1, ngastp
              if (xarc(ixarc)%occ_cnt(hpvx,ica).gt.0 .and.
     &            xarc(jxarc)%occ_cnt(hpvx,ica).gt.0) then
                if (xarc(ixarc)%occ_cnt(hpvx,ica).eq.
     &                        occ_vtx(hpvx,ica_vtx,ivtx1)) then
                  ! ixarc fully contracted to ivtx1
                  ! -> move jxarc here
c dbg
c                  print *,'left case (x)'
c dbg

                  cnt_shl(hpvx,ica) = xarc(jxarc)%occ_cnt(hpvx,ica)
                  occ_shl(hpvx,ica_vtx) = xarc(jxarc)%occ_cnt(hpvx,ica)

                else if (xarc(jxarc)%occ_cnt(hpvx,ica).eq.
     &                        occ_vtx(hpvx,ica_vtx,ivtx2)) then
                  ! jxarc fully contracted to ivtx2
                  ! -> move ixarc here
c dbg
c                  print *,'right case (x)'
c dbg

                  cnt_shr(hpvx,ica) = xarc(ixarc)%occ_cnt(hpvx,ica)
                  occ_shr(hpvx,ica_vtx) = xarc(ixarc)%occ_cnt(hpvx,ica)

c                else if (contr%nsupvtx.eq.2) then
cc dbg
c                  print *,'special case'
cc dbg
c                  
c                  cnt_shl(hpvx,ica) = xarc(jxarc)%occ_cnt(hpvx,ica)
c                  occ_shl(hpvx,ica_vtx) = xarc(jxarc)%occ_cnt(hpvx,ica)
c
c                  occ_shr(hpvx,ica_vtx) = occ_vtx(hpvx,ica_vtx,ivtx1)
c     &                 - xarc(ixarc)%occ_cnt(hpvx,ica)
c
                else
c dbg
c                  print *,'skipping difficult reo (x)'
c dbg
                  possible = .false.
                  return
c                  call quit(1,'reorder_supvtx','not yet')
                  
                end if
              end if

            end do
          end do
c dbg
c          print *,'CNT SHL:'
c          call wrt_occ(luout,cnt_shl)
c          print *,'CNT SHR:'
c          call wrt_occ(luout,cnt_shr)
c          print *,'OCC SHL:'
c          call wrt_occ(luout,occ_shl)
c          print *,'OCC SHR:'
c          call wrt_occ(luout,occ_shr)
c dbg

          ! update contractions
          xarc(ixarc)%occ_cnt =
     &         xarc(ixarc)%occ_cnt + cnt_shl - cnt_shr
          xarc(jxarc)%occ_cnt =
     &         xarc(jxarc)%occ_cnt - cnt_shl + cnt_shr

          ! update super-vertex occupation
          if (modify_contr) then
            occ_vtx(1:ngastp,1:2,ivtx1) =
     &           occ_vtx(1:ngastp,1:2,ivtx1) + occ_shl - occ_shr
            occ_vtx(1:ngastp,1:2,ivtx2) =
     &           occ_vtx(1:ngastp,1:2,ivtx2) - occ_shl + occ_shr
          end if

          if (set_reord_list) then
            if (iocc_nonzero(occ_shr)) then
              reo_info%nreo = reo_info%nreo+1
              idx = reo_info%nreo
              if (idx.gt.maxreo)
     &             call quit(1,'reorder_supvtx','unexpected event')
              idxsuper = svertex(ivtx1)
              reo_info%reo(idx)%idxsuper = idxsuper
              ! flag whether this is the result vertex of prev. binary contr.
              reo_info%reo(idx)%is_bc_result =
     &             idxop12.eq.vertex(ivtx1)%idx_op
              ! which components of super-vertex
              reo_info%reo(idx)%to =
     &             imltlist(idxsuper,svertex,ivtx2,1)
              reo_info%reo(idx)%from =
     &             imltlist(idxsuper,svertex,ivtx1,1)
              ! shifted occupation
              reo_info%reo(idx)%occ_shift = occ_shr
            end if
            if (iocc_nonzero(occ_shl)) then
              reo_info%nreo = reo_info%nreo+1
              idx = reo_info%nreo
              if (idx.gt.maxreo)
     &             call quit(1,'reorder_supvtx','unexpected event')
              idxsuper = svertex(ivtx1)
              reo_info%reo(idx)%idxsuper = idxsuper
              ! flag whether this is the result vertex of prev. binary contr.
              reo_info%reo(idx)%is_bc_result =
     &             idxop12.eq.vertex(ivtx1)%idx_op
              ! which components of super-vertex
              reo_info%reo(idx)%to =
     &             imltlist(idxsuper,svertex,ivtx1,1)
              reo_info%reo(idx)%from =
     &             imltlist(idxsuper,svertex,ivtx2,1)
              ! shifted occupation
              reo_info%reo(idx)%occ_shift = occ_shl
            end if
              
          end if

        end do

      end do

      ! remove deleted or zero xarcs:
      if (modify_contr) then
        nxarc_new = 0
        do ixarc = 1, nxarc
          if (xarc(ixarc)%link(1).le.0 .or.
     &        iocc_zero(xarc(ixarc)%occ_cnt) ) cycle
          nxarc_new = nxarc_new+1
          if (nxarc_new.lt.ixarc) xarc(nxarc_new) = xarc(ixarc)
        end do
        contr%nxarc = nxarc_new
      end if

      if (.not.modify_contr) deallocate(xarc_scr)

      if (ntest.ge.100) then
        write(luout,*) 'contr at the end of reorder_supvtx'
        if (.not.modify_contr) write(luout,*) 'should not have changed!'
        call prt_contr3(luout,contr,occ_vtx)
      end if

      end
