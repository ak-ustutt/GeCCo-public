*----------------------------------------------------------------------*
      subroutine reorder_supvtx(possible,
     &     modify_contr,set_reord_list,reo_before,reo_info,
     &     contr,occ_vtx,rstr_vtx,idxop12,orb_info)
*----------------------------------------------------------------------*
*     check whether a single vertex is contracted to more than one 
*     component of a super vertex
*     if this is the case, we may shift the involved occupations to
*     the super-vertex closest to the vertex they are commonly 
*     contracted to, and we may add up the contraction
*
*     basically, this process is of relevance only if both contractions
*     refer to indices of the same C/A H/P/V/X (corresponds to a partial
*     symmetrization of the intermediate represented by the super 
*     vertex); it will reduce the amount of storage necessary for the
*     intermediate
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_reorder_info.h'
      include 'def_orbinf.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      logical, intent(out) ::
     &     possible
      logical, intent(in) ::
     &     modify_contr, set_reord_list, reo_before
      integer, intent(in) ::
     &     idxop12
      type(contraction), intent(inout) ::
     &     contr
      type(reorder_info), intent(inout) ::
     &     reo_info
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(inout) ::
     &     occ_vtx(ngastp,2,contr%nvtx)
      integer, intent(inout) ::
     &     rstr_vtx(2,orb_info%ngas,2,2,contr%nvtx)

      type(cntr_arc), pointer ::
     &     arc(:), arc_scr(:)
      type(cntr_vtx), pointer ::
     &     vertex(:)
      integer, pointer ::
     &     svertex(:), hpvxgas(:,:)
      integer ::
     &     occ_shr(ngastp,2), occ_shl(ngastp,2),
     &     cnt_shr(ngastp,2), cnt_shl(ngastp,2),
     &     rstr_shr(2,orb_info%ngas,2,2), rstr_shl(2,orb_info%ngas,2,2)
      integer ::
     &     narc, iarc, jarc, ivtx0, ivtx1, ivtx2, iprim, isuper,
     &     iarc_shift, iarc_sum, narc_new,
     &     maxreo, idx, idxsuper, ica, ica_vtx, hpvx,
     &     iblk, ivtx, idxnew, ireo, ngas, nspin
      logical ::
     &     renamed(contr%nvtx), modified
      logical, pointer ::
     &     reo_generated(:)
      integer, pointer ::
     &     arc_involved(:)

      integer, external ::
     &     imltlist

c dbg
c      print *,'reorder_supvtx: on input nreo ',reo_info%nreo
c dbg
      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'reorder_supvtx')
        write(lulog,*) 'modify_contr: ',modify_contr
        write(lulog,*) 'set_reord_list: ',set_reord_list
        write(lulog,*) 'idxop12: ',idxop12
        write(lulog,*) 'contr on input:'
        call prt_contr3(lulog,contr,occ_vtx)
      end if

      if (orb_info%nspin.gt.1)
     &     call quit(1,'reorder_supvtx','adapt for nspin>1')
      nspin = orb_info%nspin
      ngas = orb_info%ngas
      hpvxgas => orb_info%ihpvgas

      possible = .true.

      maxreo = 2*contr%nvtx
      if (set_reord_list.and..not.associated(reo_info%reo)) then
        reo_info%nreo = 0
        allocate(reo_info%reo(maxreo),reo_info%nca_vtx(contr%nvtx))
        reo_info%nvtx_contr = contr%nvtx
        call set_nca_vtx(reo_info%nca_vtx,occ_vtx,contr%nvtx)
      end if
      allocate(reo_generated(maxreo),arc_involved(maxreo))
      reo_generated = .false.
      arc_involved = -1

      narc = contr%narc
      if (modify_contr) then
        arc => contr%arc
      else
        allocate(arc_scr(contr%narc))
        arc_scr = contr%arc
        arc => arc_scr
      end if
      svertex => contr%svertex
      vertex =>  contr%vertex
      ! loop over arcs
      do iarc = 1, narc
        ! deleted arc?
        if (arc(iarc)%link(1).le.0) cycle
        ! loop over arcs and look for identical link(1) (if jarc<iarc)
        !                      or for identical link(2) (if jarc>iarc)
        do jarc = 1, narc
          ! deleted arc? or the arc itself?
          if (arc(iarc)%link(1).le.0.or.iarc.eq.jarc) cycle

c dbg
c          print *,'iarc, jarc: ',iarc,jarc
c dbg
          
          ! primitive vertex
          if (iarc.gt.jarc) iprim = 1
          if (iarc.lt.jarc) iprim = 2
          ! primitive vertex is the same?
c dbg
c          print *,'iprim = ',iprim
c          print *,'arc(iarc)%link(iprim),arc(jarc)%link(iprim):',
c     &         arc(iarc)%link(iprim),arc(jarc)%link(iprim)
c dbg
          if (arc(iarc)%link(iprim).ne.arc(jarc)%link(iprim)) cycle

          ! assumed super vertex
          if (iarc.gt.jarc) isuper = 2
          if (iarc.lt.jarc) isuper = 1

          ivtx0 = arc(iarc)%link(iprim)
          ivtx1 = arc(iarc)%link(isuper)
          ivtx2 = arc(jarc)%link(isuper)

c fix -- do not consider this case:
c          if (vertex(ivtx1)%idx_op.ne.idxop12) cycle

c dbg
c          print *,'iprim,isuper:',iprim,isuper
c          print *,'ivtx0,1,2: ',ivtx0,ivtx1,ivtx2
c dbg
          ! ensure ivtx0 < ivtx1 and ivtx0 < ivtx2 
          !     or ivtx0 > ivtx1 and ivtx0 > ivtx2
          if (.not. ( (ivtx0.lt.ivtx1.and.ivtx0.lt.ivtx2) .or.
     &                (ivtx0.gt.ivtx1.and.ivtx0.gt.ivtx2) ) ) cycle
c dbg
c          print *,'survived a)'
c dbg

          ! same super vertex?
          if (svertex(ivtx1).ne.svertex(ivtx2)) cycle
c dbg
c          print *,'survived b)'
c dbg

          ! do contractions have non-vanishing overlap?
c          if (iocc_zero(iocc_overlap(arc(iarc)%occ_cnt,.false.,
c     &                               arc(jarc)%occ_cnt,.false.))) cycle

          occ_shr = 0
          occ_shl = 0
          cnt_shr = 0
          cnt_shl = 0
          do ica = 1, 2
            ! which C/A shall we look at in the associated vertices?
            ica_vtx = ica
            if (ivtx0.lt.ivtx1) ica_vtx = 3-ica
            do hpvx = 1, ngastp
              if (arc(iarc)%occ_cnt(hpvx,ica).gt.0 .and.
     &            arc(jarc)%occ_cnt(hpvx,ica).gt.0) then
                if (arc(iarc)%occ_cnt(hpvx,ica).eq.
     &                        occ_vtx(hpvx,ica_vtx,ivtx1)) then
                  ! iarc fully contracted to ivtx1
                  ! -> move jarc here
c dbg
c                  print *,'left case'
c dbg

                  cnt_shl(hpvx,ica) = arc(jarc)%occ_cnt(hpvx,ica)
                  occ_shl(hpvx,ica_vtx) = arc(jarc)%occ_cnt(hpvx,ica)

                else if (arc(jarc)%occ_cnt(hpvx,ica).eq.
     &                        occ_vtx(hpvx,ica_vtx,ivtx2)) then
                  ! jarc fully contracted to ivtx2
                  ! -> move iarc here
c dbg
c                  print *,'right case'
c dbg

                  cnt_shr(hpvx,ica) = arc(iarc)%occ_cnt(hpvx,ica)
                  occ_shr(hpvx,ica_vtx) = arc(iarc)%occ_cnt(hpvx,ica)

                else if (contr%nsupvtx.eq.2) then
                ! former treatment (outcommented by cmh) lead to problems
                ! in case where 2 njoined=2 intermediates were involved
                ! (leading to a njoined=2 result)
c dbg
c                  print *,'special case: just going ahead'
c dbg
                  
cmh                  cnt_shr(hpvx,ica) = arc(iarc)%occ_cnt(hpvx,ica)
cmh                  occ_shr(hpvx,ica_vtx) = arc(iarc)%occ_cnt(hpvx,ica)

cmh                  occ_shr(hpvx,ica_vtx) = occ_vtx(hpvx,ica_vtx,ivtx1)
cmh     &                 - arc(iarc)%occ_cnt(hpvx,ica)

                else
c dbg
c                  print *,'skipping difficult reo'
c dbg
cmh                  possible = .false.
cmh               go ahead for now, maybe it will become possible after
cmh               some other contractions. If it remains not allowed,
cmh               this will be checked by allowed_contr
cmh                  return
c                  call quit(1,'reorder_supvtx','not yet')
                  
                end if
              end if

            end do
          end do
c dbg
c          print *,'CNT SHL:'
c          call wrt_occ(lulog,cnt_shl)
c          print *,'CNT SHR:'
c          call wrt_occ(lulog,cnt_shr)
c          print *,'OCC SHL:'
c          call wrt_occ(lulog,occ_shl)
c          print *,'OCC SHR:'
c          call wrt_occ(lulog,occ_shr)
c dbg

          ! update contractions
          arc(iarc)%occ_cnt =
     &         arc(iarc)%occ_cnt + cnt_shl - cnt_shr
          arc(jarc)%occ_cnt =
     &         arc(jarc)%occ_cnt - cnt_shl + cnt_shr

          ! update super-vertex occupation
          if (modify_contr) then
            call fit_restr(rstr_shr,occ_shr,
     &           rstr_vtx(1,1,1,1,ivtx1),hpvxgas,ngas)
            call fit_restr(rstr_shl,occ_shl,
     &           rstr_vtx(1,1,1,1,ivtx2),hpvxgas,ngas)
c dbg
c            print *,'vtx # 1:'
c            call wrt_occ_rstr(6,ivtx1,occ_vtx(1,1,ivtx1),
c     &                        rstr_vtx(1,1,1,1,ivtx1),ngas,nspin)
c            print *,'vtx # 2:'
c            call wrt_occ_rstr(6,ivtx2,occ_vtx(1,1,ivtx2),
c     &                        rstr_vtx(1,1,1,1,ivtx2),ngas,nspin)
c            print *,'occ_shr'
c            call wrt_occ_rstr(6,0,occ_shr,
c     &                        rstr_shr,ngas,nspin)
c            print *,'occ_shl'
c            call wrt_occ_rstr(6,0,occ_shl,
c     &                        rstr_shl,ngas,nspin)
c dbg
c            rstr_vtx(1:2,1:ngas,1:2,1:2,1:nspin,ivtx1) =
c     &           rstr_vtx(1:2,1:ngas,1:2,1:2,1:nspin,ivtx1)
            rstr_vtx(1:2,1:ngas,1:2,1:2,ivtx1) =
     &           rstr_vtx(1:2,1:ngas,1:2,1:2,ivtx1)
     &           + rstr_shl - rstr_shr
c            rstr_vtx(1:2,1:ngas,1:2,1:2,1:nspin,ivtx2) =
c     &           rstr_vtx(1:2,1:ngas,1:2,1:2,1:nspin,ivtx2)
            rstr_vtx(1:2,1:ngas,1:2,1:2,ivtx2) =
     &           rstr_vtx(1:2,1:ngas,1:2,1:2,ivtx2)
     &           - rstr_shl + rstr_shr
            occ_vtx(1:ngastp,1:2,ivtx1) =
     &           occ_vtx(1:ngastp,1:2,ivtx1) + occ_shl - occ_shr
            occ_vtx(1:ngastp,1:2,ivtx2) =
     &           occ_vtx(1:ngastp,1:2,ivtx2) - occ_shl + occ_shr
c dbg
c            print *,'new vtx # 1:'
c            call wrt_occ_rstr(6,ivtx1,occ_vtx(1,1,ivtx1),
c     &                        rstr_vtx(1,1,1,1,ivtx1),ngas,nspin)
c            print *,'new vtx # 2:'
c            call wrt_occ_rstr(6,ivtx2,occ_vtx(1,1,ivtx2),
c     &                        rstr_vtx(1,1,1,1,ivtx2),ngas,nspin)
c dbg
          end if

          if (set_reord_list) then
            if (iocc_nonzero(occ_shr)) then
              reo_info%nreo = reo_info%nreo+1
              reo_generated(reo_info%nreo) = .true.
              arc_involved(reo_info%nreo) = iarc
              idx = reo_info%nreo
              if (idx.gt.maxreo) then
                write(lulog,*) 'idx,maxreo: ',idx,maxreo
                call quit(1,'reorder_supvtx','unexpected event')
              end if
              idxsuper = svertex(ivtx1)
              reo_info%reo(idx)%idxsuper = idxsuper
              reo_info%reo(idx)%idxop_ori =
     &             vertex(ivtx1)%idx_op
              reo_info%reo(idx)%iblkop_ori =
     &             vertex(ivtx1)%iblk_op
              reo_info%reo(idx)%dagger_ori =
     &             vertex(ivtx1)%dagger
              ! flag whether this is the result vertex of prev. binary contr.
              reo_info%reo(idx)%is_bc_result =
     &             idxop12.eq.vertex(ivtx1)%idx_op
              ! flag whether reo should occur before contraction
              reo_info%reo(idx)%reo_before   = reo_before
              ! which components of super-vertex
              reo_info%reo(idx)%to =
     &             imltlist(idxsuper,svertex,ivtx2,1)
              reo_info%reo(idx)%from =
     &             imltlist(idxsuper,svertex,ivtx1,1)
              reo_info%reo(idx)%to_vtx   = ivtx2
              reo_info%reo(idx)%from_vtx = ivtx1
              ! shifted occupation
              reo_info%reo(idx)%occ_shift = occ_shr
              reo_info%reo(idx)%shift_i0 = .false.
            end if
            if (iocc_nonzero(occ_shl)) then
              reo_info%nreo = reo_info%nreo+1
              reo_generated(reo_info%nreo) = .true.
              arc_involved(reo_info%nreo) = jarc
              idx = reo_info%nreo
              if (idx.gt.maxreo) then
                write(lulog,*) 'idx,maxreo: ',idx,maxreo
                call quit(1,'reorder_supvtx','unexpected event')
              end if
              idxsuper = svertex(ivtx1)
              reo_info%reo(idx)%idxsuper = idxsuper
              reo_info%reo(idx)%idxop_ori =
     &             vertex(ivtx1)%idx_op
              reo_info%reo(idx)%iblkop_ori =
     &             vertex(ivtx1)%iblk_op
              reo_info%reo(idx)%dagger_ori =
     &             vertex(ivtx1)%dagger
              ! flag whether this is the result vertex of prev. binary contr.
              reo_info%reo(idx)%is_bc_result =
     &             idxop12.eq.vertex(ivtx1)%idx_op
              ! flag whether reo should occur before contraction
              reo_info%reo(idx)%reo_before   = reo_before
              ! which components of super-vertex
              reo_info%reo(idx)%to =
     &             imltlist(idxsuper,svertex,ivtx1,1)
              reo_info%reo(idx)%from =
     &             imltlist(idxsuper,svertex,ivtx2,1)
              reo_info%reo(idx)%to_vtx   = ivtx1
              reo_info%reo(idx)%from_vtx = ivtx2
              ! shifted occupation
              reo_info%reo(idx)%occ_shift = occ_shl
              reo_info%reo(idx)%shift_i0 = .false.
            end if
              
          end if

        end do

      end do

      ! remove deleted or zero arcs:
      if (modify_contr) then
        modified = .false.
        narc_new = 0
        do iarc = 1, narc
          if (arc(iarc)%link(1).le.0 .or.
     &        iocc_zero(arc(iarc)%occ_cnt) ) cycle
          narc_new = narc_new+1
          if (narc_new.lt.iarc) arc(narc_new) = arc(iarc)
        end do
        contr%narc = narc_new
        modified = narc.ne.narc_new

        ! a quickie: new intermediate
        renamed = .false.
        idxnew = idxop12-1
        do ireo = 1, reo_info%nreo
          if (.not.reo_generated(ireo)) cycle

          modified = .true.
          
          reo_info%reo(ireo)%idxop_new = idxop12 ! default
          reo_info%reo(ireo)%dagger_new = .false. ! default
          idxsuper = reo_info%reo(idx)%idxsuper
          iblk = 0
          do ivtx = 1, contr%nvtx
            if (renamed(ivtx)) cycle ! do not rename twice
            if (svertex(ivtx).ne.idxsuper) cycle
            if (reo_before.or.vertex(ivtx)%idx_op.ne.idxop12) then
              renamed(ivtx) = .true.
              vertex(ivtx)%idx_op = idxnew
              iblk = iblk+1
              vertex(ivtx)%iblk_op = iblk
              vertex(ivtx)%dagger = .false.
            end if
          end do
          if (iblk.ne.0) reo_info%reo(ireo)%idxop_new = idxnew
          if (iblk.ne.0) idxnew = idxnew-1

          if (contr%index_info)
     &         call reorder_string_info(contr,
     &              reo_info%reo(ireo)%from_vtx,
     &              reo_info%reo(ireo)%to_vtx,
     &              arc_involved(ireo),
     &              reo_info%reo(ireo)%occ_shift,
     &              .false.)
          
        end do

        ! set 0-contraction, if necessary
        ! as they were removed above
        call check_disconnected(contr)

        if (contr%index_info.and.modified)
     &       call update_string_info(contr)
        
      end if

      if (.not.modify_contr) deallocate(arc_scr)

      deallocate(reo_generated,arc_involved)

      if (ntest.ge.100) then
        write(lulog,*) 'contr at the end of reorder_supvtx: nreo = ',
     &       reo_info%nreo
        if (.not.modify_contr) write(lulog,*) 'should not have changed!'
        call prt_contr3(lulog,contr,occ_vtx)
      end if

      end
