*----------------------------------------------------------------------*
      subroutine reorder_supvtx(possible,
     &     modify_contr,set_reord_list,reo_before,reo_info,
     &     contr,occ_vtx,idxop12)
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
      integer, intent(inout) ::
     &     occ_vtx(ngastp,2,contr%nvtx)


      type(cntr_arc), pointer ::
     &     arc(:), arc_scr(:)
      type(cntr_vtx), pointer ::
     &     vertex(:)
      integer, pointer ::
     &     svertex(:) !, occ_shift(:,:)
      integer ::
     &     occ_shr(ngastp,2), occ_shl(ngastp,2),
     &     cnt_shr(ngastp,2), cnt_shl(ngastp,2)
      integer ::
     &     narc, iarc, jarc, ivtx0, ivtx1, ivtx2, iprim, isuper,
     &     iarc_shift, iarc_sum, narc_new,
     &     maxreo, idx, idxsuper, ica, ica_vtx, hpvx,
     &     iblk, ivtx, idxnew, ireo
      logical, pointer ::
     &     reo_generated(:)

      integer, external ::
     &     imltlist

c dbg
c      print *,'reorder_supvtx: on input nreo ',reo_info%nreo
c dbg
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'reorder_supvtx')
        write(luout,*) 'modify_contr: ',modify_contr
        write(luout,*) 'set_reord_list: ',set_reord_list
        write(luout,*) 'idxop12: ',idxop12
        write(luout,*) 'contr on input:'
        call prt_contr3(luout,contr,occ_vtx)
      end if

      possible = .true.

      maxreo = 2*contr%nvtx
      if (set_reord_list.and..not.associated(reo_info%reo)) then
        reo_info%nreo = 0
        allocate(reo_info%reo(maxreo),reo_info%nca_vtx(contr%nvtx))
        reo_info%nvtx_contr = contr%nvtx
        call set_nca_vtx(reo_info%nca_vtx,occ_vtx,contr%nvtx)
      end if
      allocate(reo_generated(maxreo))
      reo_generated = .false.

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
c dbg
c                  print *,'special case'
c dbg
                  
                  cnt_shl(hpvx,ica) = arc(jarc)%occ_cnt(hpvx,ica)
                  occ_shl(hpvx,ica_vtx) = arc(jarc)%occ_cnt(hpvx,ica)

                  occ_shr(hpvx,ica_vtx) = occ_vtx(hpvx,ica_vtx,ivtx1)
     &                 - arc(iarc)%occ_cnt(hpvx,ica)

                else
c dbg
c                  print *,'skipping difficult reo'
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
          arc(iarc)%occ_cnt =
     &         arc(iarc)%occ_cnt + cnt_shl - cnt_shr
          arc(jarc)%occ_cnt =
     &         arc(jarc)%occ_cnt - cnt_shl + cnt_shr

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
              reo_generated(reo_info%nreo) = .true.
              idx = reo_info%nreo
              if (idx.gt.maxreo) then
                write(luout,*) 'idx,maxreo: ',idx,maxreo
                call quit(1,'reorder_supvtx','unexpected event')
              end if
              idxsuper = svertex(ivtx1)
              reo_info%reo(idx)%idxsuper = idxsuper
              reo_info%reo(idx)%idxop_ori =
     &             vertex(ivtx1)%idx_op
              reo_info%reo(idx)%iblkop_ori =
     &             vertex(ivtx1)%iblk_op
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
            end if
            if (iocc_nonzero(occ_shl)) then
              reo_info%nreo = reo_info%nreo+1
              reo_generated(reo_info%nreo) = .true.
              idx = reo_info%nreo
              if (idx.gt.maxreo) then
                write(luout,*) 'idx,maxreo: ',idx,maxreo
                call quit(1,'reorder_supvtx','unexpected event')
              end if
              idxsuper = svertex(ivtx1)
              reo_info%reo(idx)%idxsuper = idxsuper
              reo_info%reo(idx)%idxop_ori =
     &             vertex(ivtx1)%idx_op
              reo_info%reo(idx)%iblkop_ori =
     &             vertex(ivtx1)%iblk_op
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
            end if
              
          end if

        end do

      end do

      ! remove deleted or zero arcs:
      if (modify_contr) then
        narc_new = 0
        do iarc = 1, narc
          if (arc(iarc)%link(1).le.0 .or.
     &        iocc_zero(arc(iarc)%occ_cnt) ) cycle
          narc_new = narc_new+1
          if (narc_new.lt.iarc) arc(narc_new) = arc(iarc)
        end do
        contr%narc = narc_new

        ! a quickie: new intermediate
        idxnew = idxop12-1
        do ireo = 1, reo_info%nreo
          if (.not.reo_generated(ireo)) cycle
          reo_info%reo(ireo)%idxop_new = idxop12 ! default
          idxsuper = reo_info%reo(idx)%idxsuper
          iblk = 0
          do ivtx = 1, contr%nvtx
            if (svertex(ivtx).ne.idxsuper) cycle
            if (vertex(ivtx)%idx_op.ne.idxop12) then
              vertex(ivtx)%idx_op = idxnew
              iblk = iblk+1
              vertex(ivtx)%iblk_op = iblk
            end if
          end do
          if (iblk.ne.0) reo_info%reo(ireo)%idxop_new = idxnew
          if (iblk.ne.0) idxnew = idxnew-1
        end do

        ! set 0-contraction, if necessary
        ! as they were removed above
        call check_disconnected(contr)

      end if

      if (.not.modify_contr) deallocate(arc_scr)

      deallocate(reo_generated)

      if (ntest.ge.100) then
        write(luout,*) 'contr at the end of reorder_supvtx: nreo = ',
     &       reo_info%nreo
        if (.not.modify_contr) write(luout,*) 'should not have changed!'
        call prt_contr3(luout,contr,occ_vtx)
      end if

      end
