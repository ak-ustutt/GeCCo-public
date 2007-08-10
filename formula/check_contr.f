*----------------------------------------------------------------------*
      logical function check_contr(contr,proto)
*----------------------------------------------------------------------*
*     check whether the contraction contr contains all elements
*     required by the proto-contraction proto
*     NOTE: the vertices must be in identical order, so check before 
*        calling canon_contr which might reorder contr
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     mxfound = 10

      type(contraction), intent(in) ::
     &     contr, proto

      logical ::
     &     found, visited(proto%narc)
      integer ::
     &     nvtx, narc, nparc, ivtx, iarc, iparc, jparc,
     &     ivtx1, ivtx2, idx1, idx2, right_left,
     &     nproto, ncontr
      type(cntr_vtx), pointer ::
     &     vertex(:), pvertex(:)
      type(cntr_arc), pointer ::
     &     arc(:), parc(:)
      integer, pointer ::
     &     pocc(:,:)

      integer ::
     &     pocc_found(ngastp,2,mxfound),
     &     occ_found(ngastp,2,mxfound)

      logical, external ::
     &     check_occ_partition

c dbg
c      print *,'in check_contr'
c dbg
      check_contr = .false.
      
      ! all vertices should be identical
      ! make check *before* reordering
      nvtx = contr%nvtx
      if (proto%nvtx.ne.nvtx) return
c dbg
c      print *,'check 1 OK'
c dbg

      vertex => contr%vertex
      pvertex => proto%vertex
      do ivtx = 1, nvtx
        if (vertex(ivtx)%idx_op.ne.pvertex(ivtx)%idx_op) return
        if (vertex(ivtx)%iblk_op.ne.pvertex(ivtx)%iblk_op) return
      end do
c dbg
c      print *,'check 2 OK'
c dbg

      narc = contr%narc
      nparc = proto%narc

      if (nparc.eq.0) then
        check_contr = .true.
        return
      end if
c dbg
c      print *,'need more checks'
c dbg

      arc => contr%arc
      parc => proto%arc

      ! look at prototype arcs 
      ! round A: is everything connected that must be connected
      !          by (x y -1) proto-arcs
      !          is everything disconnected that must be disconnected
      !          by (x y [0]) proto-arcs
      do iparc = 1, nparc
        ivtx1 = parc(iparc)%link(1)
        ivtx2 = parc(iparc)%link(2)
        pocc => parc(iparc)%occ_cnt
        ! no connection allowed
        if (iocc_zero(pocc)) then
          found = .false.
          do iarc = 1, narc
            found = found.or.
     &           (arc(iarc)%link(1).eq.ivtx1 .and.
     &            arc(iarc)%link(2).eq.ivtx2 .and.
     &            iocc_nonzero(arc(iarc)%occ_cnt) )
          end do
c dbg
c          print *,ivtx1,ivtx2,' must not be connected -> ',.not.found
c dbg
          if (found) return

        ! a connection must be present
        else if (pocc(1,1).lt.0) then
          found = .false.
          do iarc = 1, narc
            found = found.or.
     &           (arc(iarc)%link(1).eq.ivtx1 .and.
     &            arc(iarc)%link(2).eq.ivtx2 .and.
     &            iocc_nonzero(arc(iarc)%occ_cnt) )
          end do
c dbg
c          print *,ivtx1,ivtx2,' must be connected -> ',found
c dbg
          if (.not.found) return
        end if
      end do
c dbg
c      print *,'final checks'
c dbg

      ! round B: look for right-side/left-side connections required by
      !          ( x 0 [C1] ), (x 0 [C2]) proto arcs:
      !          vertex x must be connected to vertices on the right/left
      !          with [C1] [C2] or partitions thereof
      do right_left = 1, 2
        if (right_left.eq.1) then
          idx1 = 1
          idx2 = 2
        else
          idx1 = 2
          idx2 = 1
        end if
        visited(1:nparc) = .false.
        do iparc = 1, nparc
          if (visited(iparc)) cycle
          visited(iparc) = .true.
          ivtx1 = parc(iparc)%link(idx1)
          ivtx2 = parc(iparc)%link(idx2)
          pocc_found(1:ngastp,1:2,1) = parc(iparc)%occ_cnt
          if (ivtx2.eq.0) then
            ! find further proto-arcs with same ivtx1
            nproto = 1
            do jparc = iparc+1, nparc
              if (visited(jparc)) cycle
              if (parc(jparc)%link(idx1).eq.ivtx1.and.
     &             iocc_nonzero(parc(jparc)%occ_cnt)) then
                nproto = nproto + 1
                if (nproto.gt.mxfound)
     &             call quit(1,'check_contr','increase mxfound')
                pocc_found(1:ngastp,1:2,nproto) = parc(jparc)%occ_cnt
              end if
            end do
            ! find all arcs in actual contraction with ivtx1
            ncontr = 0
            do iarc = 1, narc
              if (arc(iarc)%link(idx1).eq.ivtx1.and.
     &             iocc_nonzero(arc(iarc)%occ_cnt)) then
                ncontr = ncontr+1
                if (ncontr.gt.mxfound)
     &             call quit(1,'check_contr','increase mxfound')
                occ_found(1:ngastp,1:2,ncontr) = arc(iarc)%occ_cnt
              end if            
            end do
            if (.not.
     &           check_occ_partition(pocc_found,nproto,
     &                               occ_found,ncontr)) return
          end if
        end do
      end do
c dbg
c      print *,'at the end: T'
c dbg

      check_contr = .true.

      return
      end
