*----------------------------------------------------------------------*
      subroutine reduce_contr(contr,occ_vtx,
     &     iarc_red,idxop_new,ivtx_new,
     &     update_ori,ivtx_ori,iarc_ori)
*----------------------------------------------------------------------*
*     successor of reduce_graph:
*     generate reduced multiple contraction after contracting vertices
*     connected by iarc_red
*     the new name for the vertex is given by idx_op_new
*     the new fused vertex number is given by ivtx_new
*     if (update_ori):
*       return the original (fused) vertex names on ivtx_ori
*       and the primary original arc number on iarc_ori
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(inout) ::
     &     contr
      integer, intent(inout) ::
     &     occ_vtx(ngastp,2,*)
      logical, intent(in) ::
     &     update_ori
      integer, intent(inout) ::
     &     ivtx_ori(*), iarc_ori(*)
      integer, intent(in) ::
     &     iarc_red, idxop_new, ivtx_new

      integer ::
     &     nvtx, narc, narc_new,
     &     ivtx, jvtx, ivtx1, ivtx2, jvtx1, jvtx2, kvtx1, kvtx2,
     &     iarc, jarc
      integer, pointer ::
     &     occ_cnt(:,:)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &       'reduce_contr at your service')
        write(luout,*) 'contraction to be carried out: ',iarc_red
        write(luout,*) 'contr on entry:'
        call prt_contr3(luout,contr,occ_vtx)
        if (update_ori) then
          write(luout,*) 'ivtx_ori:'
          write(luout,'(3x,5i14)') ivtx_ori(1:contr%nvtx)
          write(luout,*) 'iarc_ori:'
          write(luout,'(3x,5i14)') iarc_ori(1:contr%narc)
        end if
      end if

      ! nsvtx = contr%nsupvtx
      nvtx = contr%nvtx
      narc = contr%narc

      ! super vertex? find further primitive vertices ...
      ! and further contractions between the supervertices

      ! complain, if iarc_red was not the first of these

      ! loop over primitive contractions

      ! first round: generate new vertices
      ! contraction info
      ivtx1 = contr%arc(iarc_red)%link(1)
      ivtx2 = contr%arc(iarc_red)%link(2)
      occ_cnt =>  contr%arc(iarc_red)%occ_cnt

      ! check whether vertices can be joined
      ! we may join vertices if neither of them connects to an
      ! operator vertex located between them
      ! unless the we have neighbouring vertices ...
c      if (ivtx2-ivtx1.eq.1) then
c        join = .true.
c      else
c        ! ... we must have a look at the topology
c        allocate(topomap(nvtx,nvtx))
c
c        call topomap4contr(1,topomap,idum,idum,idum,contr,occ_vtx(1,1,nj+1))
c
c        join = .true.
c        do ivtx = ivtx1+1, ivtx2-1
c          join = join.and.topomap(ivtx,ivtx2).eq.0
c        end do
c
c        deallocate(topomap)
c
c      end if

      ! generate occupation of new vertex 
      !  (+1 as first entry on occ_vtx is result vertex)
      occ_vtx(1:ngastp,1:2,ivtx1+1) =
     &       occ_vtx(1:ngastp,1:2,ivtx1+1) +
     &       occ_vtx(1:ngastp,1:2,ivtx2+1)
     &       - occ_cnt - iocc_dagger(occ_cnt)
      contr%vertex(ivtx1)%idx_op = idxop_new
      contr%vertex(ivtx1)%iblk_op = 1

      ! modify arcs
      do jarc = 1, narc
        jvtx1 = contr%arc(jarc)%link(1)
        jvtx2 = contr%arc(jarc)%link(2)
        if (jvtx1.eq.ivtx2)
     &       contr%arc(jarc)%link(1) = ivtx1
        if (jvtx2.eq.ivtx2)
     &       contr%arc(jarc)%link(2) = ivtx1
      end do

c      ! mark for deletion
c      contr%vertex(ivtx2)%idx_op = -1
c      contr%arc(iarc_red(idx))%link(1) = -1

      if (ivtx2.le.ivtx1) then
        write(luout,*) 'ivtx1, ivtx2: ',ivtx1,ivtx2
        call quit(1,'reduce_contr','did not expect ivtx2.le.ivtx1!')
      end if

      if (update_ori) then
        ! delete ivtx2 in ivtx_ori
        ivtx_ori(ivtx1) = ivtx_new
        do ivtx = ivtx2, nvtx-1
          ivtx_ori(ivtx) = ivtx_ori(ivtx+1)
        end do
      end if

      ! second round: delete old vertices and arcs
      do ivtx = ivtx2, nvtx-1 ! loop over old vertices
        contr%vertex(ivtx) = contr%vertex(ivtx+1)
        occ_vtx(1:ngastp,1:2,ivtx+1) =
     &       occ_vtx(1:ngastp,1:2,ivtx+2)
      end do
      contr%nvtx = nvtx-1

      contr%arc(iarc_red)%link(1) = -1
      narc_new = 0  ! counter for new arcs
      do iarc = 1, narc  ! loop over old arcs
        if (contr%arc(iarc)%link(1).lt.0) cycle

        narc_new = narc_new+1
        jvtx1 = contr%arc(iarc)%link(1) ! get involved vertices
        jvtx2 = contr%arc(iarc)%link(2) !  (old numbers)
        
        ! bring into correct sequence
        if (jvtx1.gt.jvtx2) then
          contr%arc(iarc)%link(1) = jvtx2
          contr%arc(iarc)%link(2) = jvtx1
          jvtx1 = contr%arc(iarc)%link(1)
          jvtx2 = contr%arc(iarc)%link(2)          
          contr%arc(iarc)%occ_cnt =
     &         iocc_dagger(contr%arc(iarc)%occ_cnt)
        end if

        do jarc = iarc+1, contr%narc ! add arcs with same vertices
          kvtx1 = contr%arc(jarc)%link(1)
          kvtx2 = contr%arc(jarc)%link(2)
          if (jvtx1.eq.kvtx1.and.jvtx2.eq.kvtx2) then
            contr%arc(iarc)%occ_cnt =
     &           contr%arc(iarc)%occ_cnt +
     &           contr%arc(jarc)%occ_cnt
            contr%arc(jarc)%link(1) = -1 ! mark for deletion
          else if (jvtx1.eq.kvtx2.and.jvtx2.eq.kvtx1) then
            contr%arc(iarc)%occ_cnt =
     &           contr%arc(iarc)%occ_cnt +
     &           iocc_dagger(contr%arc(jarc)%occ_cnt)
            contr%arc(jarc)%link(1) = -1 ! mark for deletion
          end if            
        end do
        if (jvtx1.gt.ivtx2)
     &       contr%arc(iarc)%link(1) = jvtx1-1 ! convert to new
        if (jvtx2.gt.ivtx2)
     &       contr%arc(iarc)%link(2) = jvtx2-1 !   numbering
        if (iarc.gt.narc_new) then     ! copy if necessary
          contr%arc(narc_new) = contr%arc(iarc)
          if (update_ori)
     &         iarc_ori(narc_new) = iarc_ori(iarc)
        end if
      end do
      contr%narc = narc_new

      if (ntest.ge.100) then
        write(luout,*) 'contr on exit:'
        call prt_contr3(luout,contr,occ_vtx)
        if (update_ori) then
          write(luout,*) 'ivtx_ori:'
          write(luout,'(3x,5i14)') ivtx_ori(1:contr%nvtx)
          write(luout,*) 'iarc_ori:'
          write(luout,'(3x,5i14)') iarc_ori(1:contr%narc)
        end if
      end if

      return
      end
