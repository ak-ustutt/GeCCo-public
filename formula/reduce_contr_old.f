*----------------------------------------------------------------------*
      subroutine reduce_contr_old(contr,occ_vtx,
     &     iarc_red,idxop_new,ivtx_new,
     &     njoined_res,
     &     update_ori,ivtx_ori,iarc_ori,
     &     update_info,irestr_vtx,info_vtx,irestr_res,orb_info)
*----------------------------------------------------------------------*
*     successor of reduce_graph:
*     generate reduced multiple contraction after contracting vertices
*     connected by iarc_red
*     the new name for the vertex is given by idx_op_new
*     the new fused vertex number is given by ivtx_new
*     if (update_ori):
*       return the original (fused) vertex names on ivtx_ori
*       and the primary original arc number on iarc_ori
*     if (update_ori):
*       update the restriction and Ms/IRREP info arrays, as well
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_orbinf.h'
      include 'ifc_operators.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 100

      type(contraction), intent(inout) ::
     &     contr
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(inout) ::
     &     occ_vtx(ngastp,2,*)
      logical, intent(in) ::
     &     update_ori, update_info
      integer, intent(inout) ::
     &     ivtx_ori(*), iarc_ori(*),
     &     irestr_vtx(2,orb_info%ngas,2,2,contr%nvtx+njoined_res),
     &     info_vtx(2,contr%nvtx+njoined_res)
      integer, intent(in) ::
     &     iarc_red, idxop_new, ivtx_new, njoined_res, 
     &     irestr_res(2,orb_info%ngas,2,2,njoined_res)

      logical ::
     &     merge
      integer ::
     &     nvtx, narc, narc_new, nsupvtx, ngas,
     &     ivtx, jvtx, ivtx1, ivtx2, jvtx1, jvtx2, kvtx1, kvtx2,
     &     iarc, jarc, ilist, jlist, len_list, nmvleft, idx_merge,
     &     iarc_prm, isupvtx1, isupvtx2, idum, ms_new, gm_new
      integer ::
     &     arc_list(contr%narc), svmap(contr%nvtx),
     &     topomap(contr%nvtx,contr%nvtx), ireo(contr%nvtx),
     &     imvleft(contr%nvtx), icyc(contr%nvtx,2)
      integer, pointer ::
     &     svertex_ori(:), occ_vtx_ori(:,:,:),
     &     ivtx_ori_ori(:), info_vtx_ori(:,:),
     &     irestr_vtx_ori(:,:,:,:,:)
      type(cntr_vtx), pointer ::
     &     vertex_ori(:)
      integer, pointer ::
     &     occ_cnt(:,:), joined(:,:), svertex(:)
      type(cntr_arc), pointer ::
     &     arc(:)
      type(cntr_vtx), pointer ::
     &     vertex(:)

      integer, external ::
     &     idx_merge_vtx1vtx2
c      logical, external ::
c     &     merge_vtx1vtx2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &       'reduce_contr at your service')
        write(luout,*) 'contraction to be carried out: ',iarc_red
        write(luout,*) 'contr on entry:'
        call prt_contr3(luout,contr,occ_vtx(1,1,njoined_res+1))
        if (update_ori) then
          write(luout,*) 'ivtx_ori:'
          write(luout,'(3x,5i14)') ivtx_ori(1:contr%nvtx)
          write(luout,*) 'iarc_ori:'
          write(luout,'(3x,5i14)') iarc_ori(1:contr%narc)
        end if
      end if

      nsupvtx = contr%nsupvtx
      nvtx = contr%nvtx
      narc = contr%narc

      ngas = orb_info%ngas

      arc => contr%arc
      vertex => contr%vertex
      svertex => contr%svertex
      joined => contr%joined

      ! super vertex? find further primitive vertices ...
      ! and further contractions between the supervertices
      call get_associated_arcs(arc_list,len_list,iarc_red,contr)

      ! complain, if iarc_red was not the first of these
      if (arc_list(1).ne.iarc_red)
     &       call quit(1,'reduce_contr',
     &       'multi-arc contraction: expecting index of first arc '//
     &       'on iarc_red!')

      ! get super-vertex map for result
      if (njoined_res.eq.1) then
        svmap(1:nvtx) = 1
      else
        call svmap4contr(svmap,contr,occ_vtx,njoined_res)
      end if

      ! get topology info
      call topomap4contr(1,topomap,idum,idum,idum,contr,
     &                       occ_vtx(1,1,njoined_res+1))

      ! init reordering array
      do ivtx = 1, nvtx
        ireo(ivtx) = ivtx
      end do

      isupvtx1 = contr%svertex(contr%arc(iarc_red)%link(1))
      isupvtx2 = contr%svertex(contr%arc(iarc_red)%link(2))

      ! loop over primitive contractions
      do ilist = 1, len_list
        iarc_prm = arc_list(ilist)
c dbg
        print *,'current primitive arc: ',iarc_prm
        print *,'ilist, len: ',ilist,len_list
        print *,'current contr:'
        call prt_contr3(luout,contr,occ_vtx(1,1,njoined_res+1))
c dbg

        ! first round: generate new vertices
        ! contraction info
        ivtx1 = contr%arc(iarc_prm)%link(1)
        ivtx2 = contr%arc(iarc_prm)%link(2)
        occ_cnt =>  contr%arc(iarc_prm)%occ_cnt

c        merge = merge_vtx1vtx2(ivtx1,ivtx2,svertex,svmap,topomap,nvtx)
        idx_merge = idx_merge_vtx1vtx2(ivtx1,ivtx2,isupvtx1,isupvtx2,
     &       nmvleft,imvleft,svertex,svmap,topomap,nvtx)
        merge = idx_merge.gt.0
c dbg
        print *,'<<<< idx_merge: ',idx_merge
        print *,'imvleft:   ',imvleft(1:nmvleft)
c dbg

        if (update_info) then
          ms_new = info_vtx(1,ivtx1+njoined_res)+
     &             info_vtx(1,ivtx2+njoined_res)
          gm_new = multd2h(info_vtx(2,ivtx1+njoined_res),
     &                     info_vtx(2,ivtx2+njoined_res))
        end if

        if (merge) then
          ! generate occupation of new vertex 
          !  (+njoined as first entry on occ_vtx is result vertex)
          occ_vtx(1:ngastp,1:2,ivtx1+njoined_res) =
     &         occ_vtx(1:ngastp,1:2,ivtx1+njoined_res) +
     &         occ_vtx(1:ngastp,1:2,ivtx2+njoined_res)
     &         - occ_cnt - iocc_dagger(occ_cnt)
          vertex(ivtx1)%idx_op = idxop_new
          vertex(ivtx1)%iblk_op = 1
          svertex(ivtx1) = nsupvtx+1
          ! modify arcs
          do jarc = 1, narc
            jvtx1 = contr%arc(jarc)%link(1)
            jvtx2 = contr%arc(jarc)%link(2)
            if (jvtx1.eq.ivtx2)
     &           contr%arc(jarc)%link(1) = ivtx1
            if (jvtx2.eq.ivtx2)
     &           contr%arc(jarc)%link(2) = ivtx1
          end do
          ! mark for deletion
          vertex(ivtx2)%idx_op = 0
          if (update_ori)
     &         ivtx_ori(ivtx1) = ivtx_new
          if (update_info) then
            info_vtx(1,ivtx1+njoined_res) = ms_new
            info_vtx(2,ivtx1+njoined_res) = gm_new
            ! the preliminary solution:
            ! to be fixed for super-vertices
            !  find out to which result vertex the current vertex
            !  actually contributes
            if (njoined_res.eq.1) then
              call fit_restr(irestr_vtx(1,1,1,1,ivtx1+njoined_res),
     &                     occ_vtx(1,1,ivtx1+njoined_res),irestr_res,
     &                     orb_info%ihpvgas,ngas)
            else
              call dummy_restr(irestr_vtx(1,1,1,1,ivtx1+njoined_res),
     &           occ_vtx(1,1,ivtx1+njoined_res),1,orb_info%ihpvgas,ngas)
            end if
          end if
          ! update reodering array
          call update_reo(ireo,nvtx,ivtx1,ivtx2,
     &                    idx_merge,imvleft,nmvleft)

        else
          occ_vtx(1:ngastp,1:2,ivtx1+njoined_res) =
     &         occ_vtx(1:ngastp,1:2,ivtx1+njoined_res) +
     &         - occ_cnt
          occ_vtx(1:ngastp,1:2,ivtx2+njoined_res) =
     &         occ_vtx(1:ngastp,1:2,ivtx2+njoined_res)
     &         - iocc_dagger(occ_cnt)
          vertex(ivtx1)%idx_op = idxop_new
          vertex(ivtx1)%iblk_op = 1
          vertex(ivtx2)%idx_op = idxop_new
          vertex(ivtx2)%iblk_op = 1
          svertex(ivtx1) = nsupvtx+1
          svertex(ivtx2) = nsupvtx+1

          if (update_ori) then
            ivtx_ori(ivtx1) = ivtx_new
            ivtx_ori(ivtx2) = ivtx_new
          end if

          if (update_info) then
            info_vtx(1,ivtx1+njoined_res) = ms_new
            info_vtx(2,ivtx1+njoined_res) = gm_new
            info_vtx(1,ivtx2+njoined_res) = ms_new
            info_vtx(2,ivtx2+njoined_res) = gm_new
            ! the preliminary solution (2x):
            if (njoined_res.eq.1) then
              call fit_restr(irestr_vtx(1,1,1,1,ivtx1+njoined_res),
     &                     occ_vtx(1,1,ivtx1+njoined_res),irestr_res,
     &                     orb_info%ihpvgas,ngas)
              call fit_restr(irestr_vtx(1,1,1,1,ivtx2+njoined_res),
     &                     occ_vtx(1,1,ivtx2+njoined_res),irestr_res,
     &                     orb_info%ihpvgas,ngas)
            else
              call dummy_restr(irestr_vtx(1,1,1,1,ivtx1+njoined_res),
     &                     occ_vtx(1,1,ivtx1+njoined_res),1,
     &                     orb_info%ihpvgas,ngas)
              call dummy_restr(irestr_vtx(1,1,1,1,ivtx2+njoined_res),
     &                     occ_vtx(1,1,ivtx2+njoined_res),1,
     &                     orb_info%ihpvgas,ngas)
            end if
          end if
        end if
        
        ! mark arc for deletion
        contr%arc(iarc_prm)%link(1) = 0

        if (ivtx2.le.ivtx1) then
          write(luout,*) 'ivtx1, ivtx2: ',ivtx1,ivtx2
          call quit(1,'reduce_contr','did not expect ivtx2.le.ivtx1 !')
        end if

      end do ! ilist

      ! second round: delete old vertices and arcs
      ! 1) reoder vertices

      allocate(vertex_ori(nvtx), occ_vtx_ori(ngastp,2,nvtx),
     &         svertex_ori(nvtx))
      vertex_ori = vertex
      occ_vtx_ori(1:ngastp,1:2,1:nvtx) = occ_vtx(1:ngastp,1:2,
     &     njoined_res+1:njoined_res+nvtx)
      svertex_ori(1:nvtx) = svertex(1:nvtx)
      if (update_ori) then
        allocate(ivtx_ori_ori(nvtx))
        ivtx_ori_ori = ivtx_ori(1:nvtx)
      end if
      if (update_info) then
        allocate(info_vtx_ori(2,nvtx),irestr_vtx_ori(2,ngas,2,2,nvtx))
        info_vtx_ori(1:2,1:nvtx) = info_vtx(1:2,
     &                                  njoined_res+1:njoined_res+nvtx)
        irestr_vtx_ori(1:2,1:ngas,1:2,1:2,1:nvtx) =
     &       irestr_vtx(1:2,1:ngas,
     &                          1:2,1:2,njoined_res+1:njoined_res+nvtx)
      end if

c dbg
      print *,'ireo: ',ireo(1:nvtx)
c dbg
      jvtx = 0
c      ireo(1:nvtx) = 0
      do ivtx = 1, nvtx   ! loop over old vertices
        if (vertex(ivtx)%idx_op.eq.0) cycle
        jvtx = jvtx+1
c        ireo(ivtx) = jvtx
c        if (jvtx.lt.ivtx) then 
        vertex(ireo(ivtx)) = vertex_ori(ivtx)
          occ_vtx(1:ngastp,1:2,ireo(ivtx)+njoined_res) =
     &         occ_vtx_ori(1:ngastp,1:2,ivtx)
          svertex(ireo(ivtx)) = svertex_ori(ivtx)
          if (update_ori)
     &         ivtx_ori(ireo(ivtx)) = ivtx_ori_ori(ivtx)
          if (update_info)
     &         info_vtx(1:2,ireo(ivtx)+njoined_res) =
     &         info_vtx_ori(1:2,ivtx)
          if (update_info)
     &         irestr_vtx(1:2,1:ngas,1:2,1:2,ireo(ivtx)+njoined_res)=
     &         irestr_vtx_ori(1:2,1:ngas,1:2,1:2,ivtx)
c        end if
      end do
      contr%nvtx = jvtx

      deallocate(vertex_ori,svertex_ori,occ_vtx_ori)
      if (update_ori) deallocate(ivtx_ori_ori)
      if (update_info) deallocate(info_vtx_ori,irestr_vtx_ori)
c dbg
      print *,'after reo of vertices'
        call prt_contr3(luout,contr,occ_vtx(1,1,2))
c dbg

      ! update vertex numbers on arcs:
      do iarc = 1, narc         ! loop over old arcs
        if (contr%arc(iarc)%link(1).le.0) cycle
        jvtx1 = contr%arc(iarc)%link(1)
        jvtx2 = contr%arc(iarc)%link(2)
        contr%arc(iarc)%link(1) = ireo(jvtx1) ! convert to new
        contr%arc(iarc)%link(2) = ireo(jvtx2) !   numbering
      end do
c dbg
      print *,'after rename of vertices in arcs'
        call prt_contr3(luout,contr,occ_vtx(1,1,2))
c dbg
      
      ! remove and join arcs:
      narc_new = 0              ! counter for new arcs
      do iarc = 1, narc         ! loop over old arcs
        if (contr%arc(iarc)%link(1).le.0) cycle

        narc_new = narc_new+1
        jvtx1 = contr%arc(iarc)%link(1) ! get involved vertices
        jvtx2 = contr%arc(iarc)%link(2) !  (old numbers)
        
        if (jvtx1.gt.jvtx2) then
          call quit(1,'reduce_contr','this should not happen')
        end if

        do jarc = iarc+1, contr%narc ! add arcs with same vertices
          kvtx1 = contr%arc(jarc)%link(1)
          kvtx2 = contr%arc(jarc)%link(2)
          if (jvtx1.eq.kvtx1.and.jvtx2.eq.kvtx2) then
            contr%arc(iarc)%occ_cnt =
     &           contr%arc(iarc)%occ_cnt +
     &           contr%arc(jarc)%occ_cnt
            contr%arc(jarc)%link(1) = 0 ! mark for deletion
          else if (jvtx1.eq.kvtx2.and.jvtx2.eq.kvtx1) then
            call quit(1,'reduce_contr','this should not happen (2)')
          end if            
        end do

c        contr%arc(iarc)%link(1) = ireo(jvtx1) ! convert to new
c        contr%arc(iarc)%link(2) = ireo(jvtx2) !   numbering

        if (iarc.gt.narc_new) then ! copy if necessary
          contr%arc(narc_new) = contr%arc(iarc)
          if (update_ori)
     &         iarc_ori(narc_new) = iarc_ori(iarc)
        end if
      end do
      contr%narc = narc_new

      call update_svtx4contr(contr)

      if (ntest.ge.100) then
        write(luout,*) 'contr on exit:'
        call prt_contr3(luout,contr,occ_vtx(1,1,2))
        if (update_ori) then
          write(luout,*) 'ivtx_ori:'
          write(luout,'(3x,5i14)') ivtx_ori(1:contr%nvtx)
          write(luout,*) 'iarc_ori:'
          write(luout,'(3x,5i14)') iarc_ori(1:contr%narc)
        end if
      end if
      
      return
      end
