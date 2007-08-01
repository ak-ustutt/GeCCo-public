*----------------------------------------------------------------------*
      logical function contr_in_contr(contra,contrb,op_info)
*----------------------------------------------------------------------*
*     check whether contraction A on contra is contained in contraction
*     B on contrb
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 0

      type(contraction), intent(in) ::
     &     contra, contrb
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     extended_check, found, c_left, c_right, shift
      integer ::
     &     nvtxa, nvtxb, narca, narcb,
     &     idx_a, idx_b, idx_first, idx_last, idx, jdx, ivtx1, ivtx2,
     &     ninter, njoined, ioff, iarc
      integer ::
     &     occ_test(ngastp,2)

      type(operator), pointer ::
     &     opres
      integer, pointer ::
     &     idx_op1, iblk_op1,
     &     iocc_a(:,:), ivtxa4ivtxb(:),
     &     conn_tab(:,:), cnt_b(:,:,:),
     &     occres(:,:,:)
      logical, pointer ::
     &     assigned(:)

c      logical, external ::
c     &     iocc_equal

      if (ntest.ge.100) then
        write(luout,*) '--------------------------'
        write(luout,*) 'output from contr_in_contr'
        write(luout,*) '--------------------------'
        write(luout,*) 'vertices A/B: ',contra%nvtx,contrb%nvtx
        write(luout,*) 'arcs A/B:     ',contra%narc,contrb%narc
      end if
      
      ! pessimistic start
      contr_in_contr = .false.

      ! number of vertices and arcs of B must be >= than those of A
      ! else we can directly say good-bye
      if (contrb%nvtx.lt.contra%nvtx.and.
     &    contrb%narc.lt.contra%narc) return

      nvtxa = contra%nvtx
      nvtxb = contrb%nvtx
      narca = contra%narc
      narcb = contrb%narc

      ! operators of B which interrupt the sequence of A
      ! must be checked to commute with remaining sequence
      extended_check = .false.
c      allocate(extended_check(nvtxb))
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
          ! we must ensure, that this operator correctly connects with
          ! members of A
          extended_check = .true.

        end do
        if (.not.found) exit
      end do

      if (ntest.ge.150) then
        write(luout,*) 'after vertex compare: found = ',found
        write(luout,*) ' ivtxa4ivtxb: ',ivtxa4ivtxb(1:nvtxb)
        write(luout,*) ' idx_first, idx_last = ',idx_first, idx_last
        write(luout,*) ' extended_check: ',extended_check
      end if

      ! not all vertices found: good-bye
      if (.not.found) then
        deallocate(ivtxa4ivtxb)
        return
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

      deallocate(assigned)

      ! exit if the arc-test failed
      if (.not.contr_in_contr) then
        deallocate(ivtxa4ivtxb)
        return
      end if

      ! finally: extended check of connectivity, if necessary
      if (extended_check) then

        if (ntest.ge.100)
     &       write(luout,*) 'running extended check ...'

        if (ntest.ge.150) then
          write(luout,*) 'A:'
          call prt_contr2(luout,contra,op_info)
          write(luout,*) 'B:'
          call prt_contr2(luout,contrb,op_info)
        end if

        ! set up an auxiliary array for connectivity
        allocate(conn_tab(nvtxb,nvtxb))
        conn_tab(1:nvtxb,1:nvtxb) = 0
        do idx = 1, narcb
          ivtx1 = contrb%arc(idx)%link(1)
          ivtx2 = contrb%arc(idx)%link(2)
          conn_tab(ivtx1,ivtx2) = idx
          conn_tab(ivtx2,ivtx1) = idx
        end do

        if (ntest.ge.150) then
          write(luout,*) 'the connectivity table'
          do idx = 1, nvtxb
            write(luout,'(3x,10i5)') conn_tab(1:nvtxb,idx)
          end do
        end if

        if (ntest.ge.150) then
          write(luout,*) 'before shifting (first,last=',
     &         idx_first,idx_last,')'
          write(luout,'(3x,10i5)') ivtxa4ivtxb(idx_first:idx_last)
        end if

        ! check whether B vertices can be moved to the left
        b_loop1: do idx_b = idx_first, idx_last
          if (ivtxa4ivtxb(idx_b).ne.0) cycle
          ! the current vertex in B does not appear in A
          ! but interrupts sequence of A

          ! does it connect to any vertex on the left?
          c_left = .false.
          do idx = idx_first, idx_b-1
            if (ivtxa4ivtxb(idx).eq.-1) cycle  ! ignore removed vertices
            if (conn_tab(idx,idx_b).ne.0) then
              c_left = .true.
              exit
            end if
          end do
          if (c_left) ivtxa4ivtxb(idx_b) = -2  ! mark as critical vertex
          if (.not.c_left) ivtxa4ivtxb(idx_b) = -1 ! remove from list

        end do b_loop1

        if (ntest.ge.150) then
          write(luout,*) 'after shifing left'
          write(luout,'(3x,10i5)') ivtxa4ivtxb(idx_first:idx_last)
        end if

        ! change critical vertices back to 0 for second round:
        do idx_b = idx_first, idx_last
          if (ivtxa4ivtxb(idx_b).eq.-2) ivtxa4ivtxb(idx_b) = 0
        end do

        ! for remaining B-vertices:
        ! check whether B vertices can be moved to the right
        b_loop2: do idx_b = idx_last, idx_first, -1
          if (ivtxa4ivtxb(idx_b).ne.0) cycle

          ! does it connect to any vertex on the right?
          c_right = .false.
          do idx = idx_last, idx_b+1, -1
            if (ivtxa4ivtxb(idx).eq.-1) cycle  ! ignore removed vertices
            if (conn_tab(idx,idx_b).ne.0) then
              c_right = .true.
              exit
            end if
          end do
          if (c_right) ivtxa4ivtxb(idx_b) = -2  ! mark as critical vertex
          if (.not.c_right) ivtxa4ivtxb(idx_b) = -1 ! remove from list

        end do b_loop2

        if (ntest.ge.150) then
          write(luout,*) 'after shifing right'
          write(luout,'(3x,10i5)') ivtxa4ivtxb(idx_first:idx_last)
        end if

        ! change critical vertices back to 0 and count
        ninter = 0
        do idx_b = idx_first, idx_last
          if (ivtxa4ivtxb(idx_b).eq.-2) then
            ivtxa4ivtxb(idx_b) = 0
            ninter = ninter+1
          end if
        end do

        if (ntest.ge.150) then
          write(luout,*) 'ninter: ', ninter
        end if

        ! if all B operators can be moved outside, we are done
        if (ninter.eq.0) then
          contr_in_contr = .true.
          deallocate(conn_tab,ivtxa4ivtxb)
          return
        end if

        ! else: check whether this is compatible with 
        ! the result type of A

        ! resort to pessimistic view
        contr_in_contr = .false.

        opres => op_info%op_arr(contra%idx_res)%op
        njoined = opres%njoined

        if (ntest.ge.150) then
          write(luout,*) 'njoined: ', njoined
        end if

        ! no way ...
        if (njoined.eq.1) then
          contr_in_contr = .false.
          deallocate(conn_tab,ivtxa4ivtxb)
          return
        end if

        ioff = (contra%iblk_res-1)*njoined
        allocate(occres(ngastp,2,njoined))
        occres(1:ngastp,1:2,1:njoined) =
     &       opres%ihpvca_occ(1:ngastp,1:2,ioff+1:ioff+njoined)

        if (ntest.ge.150) then
          write(luout,*) 'occres: '
          call wrt_occ_n(luout,occres,njoined)
        end if

        ! loop over critical B vertices
        ! add up connectivities, if operators can be joined
        allocate(cnt_b(ngastp,2,ninter)) ! maximum
        ninter = 0  ! reset counter for insertion places
        idx_a = 0
        do idx_b = idx_first, idx_last
          if (ivtxa4ivtxb(idx_b).ne.0) cycle
          ! can we shift this vertex to the last one?
          shift = .false.
          if (idx_a.gt.0) then
            shift = .true.
            do idx = idx_a+1, idx_b-1
              shift = shift.and.conn_tab(idx,idx_b).eq.0
            end do
          end if
          ! if not, we must require one further inserting place into 
          ! the A intermediate
          if (.not.shift) then
            ninter = ninter+1
            cnt_b(1:ngastp,1:2,ninter) = 0
          end if
          do idx = idx_first, idx_last
            if (ivtxa4ivtxb(idx).le.0) cycle
            iarc = conn_tab(idx,idx_b)
            if (iarc.eq.0) cycle
            if (contrb%arc(iarc)%link(1).eq.idx_b) then
              cnt_b(1:ngastp,1:2,ninter) =
     &           cnt_b(1:ngastp,1:2,ninter) +
     &           contrb%arc(iarc)%occ_cnt(1:ngastp,1:2)
            else
              cnt_b(1:ngastp,1:2,ninter) =
     &           cnt_b(1:ngastp,1:2,ninter) +
     &           iocc_dagger(contrb%arc(iarc)%occ_cnt(1:ngastp,1:2))
            end if
          end do
          idx_a = idx_b

        end do

        if (ntest.ge.150) then
          write(luout,*) 'cnt_b (ninter = ',ninter,')'
          call wrt_occ_n(luout,cnt_b,ninter)
        end if

        ! in fact this test works only for njoined.eq.2
        if (njoined.ne.2)
     &       call quit(1,'contr_in_contr','not yet completely general')
        if (ninter.lt.njoined) then
          contr_in_contr = .true.
          do idx = 1, ninter

            occ_test = 0
            ! max. connections from below (EX)
            do jdx = idx+1, njoined
              occ_test = occ_test +
     &             iocc_xdn(1,occres(1:ngastp,1:2,jdx))
            end do

            ! max. connections from above (DX)
            do jdx = 1, idx
              occ_test = occ_test +
     &             iocc_xdn(2,occres(1:ngastp,1:2,jdx))
            end do

            contr_in_contr = contr_in_contr .and.
     &           iocc_bound('<=',cnt_b(1:ngastp,1:2,idx),.false.,
     &                           occ_test,               .false.)

            if (ntest.ge.150) then
              write(luout,*) 'idx, contr_in_contr: ',idx, contr_in_contr
              write(luout,*) 'cnt_b: '
              call wrt_occ(luout,cnt_b(1,1,idx))
              write(luout,*) 'occ_test: '
              call wrt_occ(luout,occ_test)
            end if
            
          end do

        end if

        deallocate(cnt_b,conn_tab,occres)

      end if

      if (ntest.ge.100)
     &     write(luout,*) 'final: ',contr_in_contr

      deallocate(ivtxa4ivtxb)

      return
      end
