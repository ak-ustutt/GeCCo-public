*----------------------------------------------------------------------*
      logical function contr_in_contr(contra,contrb)
*----------------------------------------------------------------------*
*     check whether contraction A on contra is contained in contraction
*     B on contrb
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(in) ::
     &     contra, contrb

      logical ::
     &     check_commute, found, c_left, c_right
      integer ::
     &     nvtxa, nvtxb, narca, narcb,
     &     idx_a, idx_b, idx_first, idx_last, idx, ivtx1, ivtx2

      integer, pointer ::
     &     idx_op1, iblk_op1,
     &     iocc_a(:,:), ivtxa4ivtxb(:)
      logical, pointer ::
     &     assigned(:)

      logical, external ::
     &     iocc_equal

      if (ntest.ge.100) then
        write(luout,*) '--------------------------'
        write(luout,*) 'output from contr_in_contr'
        write(luout,*) '--------------------------'
        write(luout,*) 'vertices A/B: ',contra%nvtx,contrb%nvtx
        write(luout,*) 'arcs A/B:     ',contra%narc,contrb%narc
      end if
      
      ! pessimistic start
      contr_in_contr = .false.

      ! number of vertices and arcs of B must be greater than those of A
      ! else we can directly say good-bye
      if (contrb%nvtx.lt.contra%nvtx.and.
     &    contrb%narc.lt.contra%narc) return

      nvtxa = contra%nvtx
      nvtxb = contrb%nvtx
      narca = contra%narc
      narcb = contrb%narc

      ! operators of B which interrupt the sequence of A
      ! must be checked to commute with remaining sequence
      check_commute = .false.
c      allocate(check_commute(nvtxb))
      ! array: for number of vertex in B give number of vertex in A
      allocate(ivtxa4ivtxb(nvtxb))
      ivtxa4ivtxb(1:nvtxb) = 0
      ! 1) compare vertices
      idx_b = 0
      idx_first = 0
      idx_last = 0
      do idx_a = 1, nvtxa
        ! get operator idx and block
        idx_op1  => contra%vertex(idx_a)%idx_op
        iblk_op1 => contra%vertex(idx_a)%iblk_op
        found = .false.
        ! compare with B
        do while(idx_b.lt.nvtxb)
          idx_b = idx_b+1
          if (idx_op1 .eq.contrb%vertex(idx_b)%idx_op .and.
     &        iblk_op1.eq.contrb%vertex(idx_b)%iblk_op) then
            ivtxa4ivtxb(idx_b) = idx_a
            if (idx_first.eq.0) idx_first = idx_b
            idx_last = idx_b
            found = .true.
            exit
          end if
          ! looking for first member? continue
          if (idx_a.eq.1) cycle

          ! some other operator comes in between
          ! we must ensure, that this operator commutes with all
          ! members of A appearing to the right of this operator
          check_commute = .true.

        end do
        if (.not.found) exit
      end do

      if (ntest.ge.100) then
        write(luout,*) 'after vertex compare: found = ',found
        write(luout,*) ' ivtxa4ivtxb: ',ivtxa4ivtxb(1:nvtxb)
        write(luout,*) ' idx_first, idx_last = ',idx_first, idx_last
        write(luout,*) ' check_commute: ',check_commute
      end if

      ! not all vertices found: good-bye
      if (.not.found) return

      ! 1a) check for commutativity, if necessary
      if (check_commute) then
        b_loop: do idx_b = idx_first, idx_last
          if (ivtxa4ivtxb(idx_b).ne.0) cycle
          ! the current vertex in B does not appear in A
          ! but interrupts sequence of A
          ! check that it either is only connected to the left
          ! or to the right:
          c_left = .false.
          c_right = .false.
          ! loop over arcs of B
          do idx = 1, narcb
            ! does it appear as first link index?
            if (contrb%arc(idx)%link(1).eq.idx_b) then
              ! is second operator contained in A?
              c_right = ivtxa4ivtxb(contrb%arc(idx)%link(2)).ne.0
            ! does it appear as second link index?
            else if (contrb%arc(idx)%link(2).eq.idx_b) then
              ! is first operator contained in A?
              c_left = ivtxa4ivtxb(contrb%arc(idx)%link(1)).ne.0
            end if
            ! if contracted to left and right, we cannot commute
            if (c_left.and.c_right) exit b_loop
          end do

        end do b_loop

        if (ntest.ge.100) then
          write(luout,*) 'c_left, c_right: ',c_left, c_right
        end if

        ! test failed: good-bye
        if (c_left.and.c_right) return

      end if

      if (narca.eq.0) contr_in_contr = .true.

      ! 2) compare arcs
      ! allocate an array to remember which vertices/arcs of B 
      ! were already assigned to vertices/arcs of A
      allocate(assigned(narcb))
      assigned(1:narcb) = .false.
      ! loop over arcs of A
      a_arcs: do idx_a = 1, narca
        ivtx1 = contra%arc(idx_a)%link(1)
        ivtx2 = contra%arc(idx_a)%link(2)
        iocc_a => contra%arc(idx_a)%occ_cnt
        found = .false.
        ! loop over arcs of B
        b_arcs: do idx_b = 1, narcb
          ! only those which are unassigned yet
          if (assigned(idx_b)) cycle
          ! sequence of Op1, Op2 matters !
          if (ivtx1.ne.ivtxa4ivtxb(contrb%arc(idx_b)%link(1)))
     &         cycle b_arcs
          if (ivtx2.ne.ivtxa4ivtxb(contrb%arc(idx_b)%link(2)))
     &         cycle b_arcs
          if (.not.iocc_equal(iocc_a,.false.,
     &         contrb%arc(idx_b)%occ_cnt,.false.))
     &         cycle b_arcs
          ! if we survived, the arcs are equal
          assigned(idx_b) = .true.
          found = .true.
          exit b_arcs
        end do b_arcs
        contr_in_contr = found
        ! if .not.found, A is not contained in B, so exit
        if (.not.found) exit a_arcs
      end do a_arcs

      if (ntest.ge.100) then
        write(luout,*) 'after comparing arcs: ',contr_in_contr
      end if

      deallocate(assigned,ivtxa4ivtxb)

      return
      end
