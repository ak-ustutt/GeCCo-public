*----------------------------------------------------------------------*
      subroutine reorder_supvtx(possible,
     &     modify_contr,set_reord_list,reo_info,
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
     &     maxreo, idx, idxsuper, ica, ica_vtx, hpvx

      integer, external ::
     &     imltlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'reorder_supvtx')
        write(luout,*) 'modify_contr: ',modify_contr
        write(luout,*) 'set_reord_list: ',set_reord_list
        write(luout,*) 'idxop12: ',idxop12
        write(luout,*) 'contr on input:'
        call prt_contr3(luout,contr,occ_vtx)
      end if

      possible = .true.

      if (set_reord_list) then
        maxreo = 2*(contr%nvtx - contr%nsupvtx)
        reo_info%nreo = 0
        allocate(reo_info%reo(maxreo))
      end if

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
                  print *,'left case'
c dbg

                  cnt_shl(hpvx,ica) = arc(jarc)%occ_cnt(hpvx,ica)
                  occ_shl(hpvx,ica_vtx) = arc(jarc)%occ_cnt(hpvx,ica)

                else if (arc(jarc)%occ_cnt(hpvx,ica).eq.
     &                        occ_vtx(hpvx,ica_vtx,ivtx2)) then
                  ! jarc fully contracted to ivtx2
                  ! -> move iarc here
c dbg
                  print *,'right case'
c dbg

                  cnt_shr(hpvx,ica) = arc(iarc)%occ_cnt(hpvx,ica)
                  occ_shr(hpvx,ica_vtx) = arc(iarc)%occ_cnt(hpvx,ica)

                else if (contr%nsupvtx.eq.2) then
c dbg
                  print *,'special case'
c dbg
                  
                  cnt_shl(hpvx,ica) = arc(jarc)%occ_cnt(hpvx,ica)
                  occ_shl(hpvx,ica_vtx) = arc(jarc)%occ_cnt(hpvx,ica)

                  occ_shr(hpvx,ica_vtx) = occ_vtx(hpvx,ica_vtx,ivtx1)
     &                 - arc(iarc)%occ_cnt(hpvx,ica)

                else
c dbg
                  print *,'skipping difficult reo'
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
      end if

      if (.not.modify_contr) deallocate(arc_scr)

      if (ntest.ge.100) then
        write(luout,*) 'contr at the end of reorder_supvtx'
        if (.not.modify_contr) write(luout,*) 'should not have changed!'
        call prt_contr3(luout,contr,occ_vtx)
      end if

      end
