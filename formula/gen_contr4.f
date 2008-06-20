*----------------------------------------------------------------------*
      subroutine gen_contr4(strict,
     &                      form_list,proto_main,
     &                      occ_vtx,ol_map,op_info)
*----------------------------------------------------------------------*
*     generate all contraction originating from a proto-contraction
*     on input: proto_main - proto-contraction, i.e. a contraction with
*                         - resulting operator
*                         - all vertices in desired sequence
*                           (including "open line vertices" see 
*                            description of occ_vtx; they can be
*                            identified from their operator number which
*                            must be the same as the result operator)
*                         - all arcs with connections which are fixed
*                           (including fixed external lines, given as
*                            usual arcs, however. They will be converted
*                            to xarcs, as soon as the actual contraction
*                            is generated).
*                         - all arcs which must be non-zero (i.e.
*                              the corresp. vertices must be connected)
*               occ_vtx - for convenience: occupations of all vertices
*                         in the sequence conforming with proto_main
*                         the open line vertices have to be set to the
*                         adjoint of the corresponding open line occu-
*                         pation such that we can treat it as if it
*                         were an operator vertex
*               ol_map  - array which tells us whether the given
*                         vertex is no open line (0)
*                         or to which vertex of the result
*                         it contributes (1...njoined_result)
*               op_info - needed for debugging only
*     on output: all new contractions appended to form-list
*
*     extension to more general targets
*
*     andreas, march 2008
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
c      include 'def_operator.h'
      include 'def_formula_item.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      
      logical, intent(in) ::
     &     strict
      type(formula_item), intent(in), target ::
     &     form_list
      type(contraction), intent(in) ::
     &     proto_main
      integer ::
     &     occ_vtx(ngastp,2,proto_main%nvtx),
c dbg
     &     occ_test(ngastp,2,proto_main%nvtx),
c dbg
     &     ol_map(proto_main%nvtx)
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     form_pnt
      integer ::
     &     nvtx, ivtx, ij_tgt
      integer ::
     &     zero_occ(ngastp,2)
      logical ::
     &     fix_vtx(proto_main%nvtx)

      type(contraction) ::
     &     proto

      if (ntest.ge.100) then
        write(luout,*) '===================='
        write(luout,*) ' gen_contr speaking'
        write(luout,*) '===================='
        write(luout,*) ' version: 4'
c        call prt_contr2(luout,proto_main,op_info)
        call prt_contr3(luout,proto_main,occ_vtx)
      end if

      if (abs(proto_main%fac).lt.1d-20) then
        call quit(1,'gen_contr',
     &       'called with zero or near-zero prefactor')
      end if

      ! constant definitions of multi-dim. array do not work on
      ! all compilers:
      zero_occ = 0

      form_pnt => form_list

      nvtx = proto_main%nvtx

      ! get a copy of the proto-type for further processing
      call init_contr(proto)
      call copy_contr(proto_main,proto)
      call strip_contr(proto,1) ! remove the (x 0 [C]) type arcs
c dbg
c      print *,'the stripped proto:'
c      call prt_contr3(luout,proto,occ_vtx)
c dbg
      do ivtx = 1, nvtx
        fix_vtx(ivtx) = ol_map(ivtx).ne.0
      end do

      ! call recursive kernel
      ivtx = 1
      call gen_contr_rec(ivtx,proto)

      ! release the allocated space for proto
      call dealloc_contr(proto)

      return

      contains

      recursive subroutine gen_contr_rec(ivtx,proto_in)

      implicit none


      integer, intent(in) ::
     &     ivtx
      type(contraction), intent(in) ::
     &     proto_in

      logical ::
     &     init, ok, must_connect(nvtx), must_not_connect(nvtx), reo
      integer ::
     &     nvtx_prev, jvtx, kvtx, iarc, narc_old, ieqvfac, ierr
      integer ::
     &     ivtx_reo(nvtx), ivtx_reo2(nvtx), ol_map_reo(nvtx),
     &     occ_ol_prev(ngastp,2), occ_ol_rem(ngastp,2),
     &     occ_ol_vtx(ngastp,2,nvtx),
     &     occ_ol_scr(ngastp,2,nvtx),
     &     occ_ex(ngastp,2), occ_dx(ngastp,2),
     &     occ_cnt2prev(ngastp,2),
     &     occ_cnt2prev_min(ngastp,2),
     &     occ_cnt2prev_max(ngastp,2),
     &     occ_cnt2rem_max(ngastp,2),
     &     occ_sum(ngastp,2)
      type(contraction) ::
     &     proto, proto_new
      logical, external ::
     &     next_part_connection, next_dist2, check_contr

      integer, allocatable ::
     &     occ_conn(:,:,:)

      if (ntest.ge.100) then
        write(luout,*) '-------------------------------------------'
        write(luout,*) ' recursive part, ivtx = ',ivtx
        write(luout,*) '-------------------------------------------'
      end if

      ! get working copies of proto-contraction
      call init_contr(proto)
      call copy_contr(proto_in,proto)

      ! get currently unconnected vertex lines
      ! _vtx:  per vertex
      ! _prev: summed over 1:ivtx-1
      ! _rem : summed over ivtx+1:nvtx
      call gen_contr4_unconn(occ_ol_prev,occ_ol_rem,occ_ol_vtx,
     &     ivtx,proto,occ_vtx)
      if (ntest.ge.100) then
        write(luout,*) 'current proto-contraction:'
        call prt_contr3(luout,proto,occ_vtx)
        write(luout,*) 'occ_ol_vtx:'
        call wrt_occ_n(luout,occ_ol_vtx,nvtx)
        write(luout,*) 'occ_ol_prev,occ_ol_rem:'
        call wrt_occ2(luout,occ_ol_prev,occ_ol_rem)
      end if

      ! get unconnected excitation and deexcitation part of operator
      occ_ex = iocc_xdn(1,occ_ol_vtx(1:ngastp,1:2,ivtx))
      occ_dx = iocc_xdn(2,occ_ol_vtx(1:ngastp,1:2,ivtx))

      ! max. possible contractions with previous operators
      occ_cnt2prev_max = iocc_overlap(occ_ex,.true.,
     &     iocc_xdn(2,occ_ol_prev),.false.)
      ! max. possible contractions with remaining operators
      occ_cnt2rem_max = iocc_overlap(occ_dx,.false.,
     &     iocc_xdn(1,occ_ol_rem),.true.)

      if (ntest.ge.100) then
        write(luout,*) '[EX],[DX],[C2prev],[C2rem]'
        call wrt_occ4(luout,occ_ex,occ_dx,
     &                      occ_cnt2prev_max,occ_cnt2rem_max)
      end if

      ! validity checks
      ! (a) do not overshoot the number of excitations
      !   [EX_current] + [EX_previous] - [CNT2prev_max] 
      !         element-wise <= 0
      occ_sum = occ_ex + iocc_xdn(1,occ_ol_prev)
     &     - iocc_dagger(occ_cnt2prev_max)

      ok = iocc_bound('<=',occ_sum,.false.,zero_occ,.false.) 
c dbg
c      call wrt_occ2(luout,occ_sum,occ_tgt_ex)
c      print *,'ok 1: ',ok
c dbg
      if (ok) then
        ! (b) do not overshoot the number of deexcitations in target
        ! [DX_current] - [CNT2rem_max] element-wise <= 0
        occ_sum = occ_dx - occ_cnt2rem_max

        ok = iocc_bound('<=',occ_sum,.false.,zero_occ,.false.)
c dbg
c      call wrt_occ(luout,occ_sum)
c        print *,'ok 2: ',ok
c dbg
      end if

      if (ok) then
        ! set up array: which operators must be connected with
        ! current one?
        must_connect(1:ivtx-1) = .false.
        do iarc = 1, proto%narc
          if (proto%arc(iarc)%link(2).eq.ivtx .and.
     &         proto%arc(iarc)%occ_cnt(1,1).lt.0)
     &         must_connect(proto%arc(iarc)%link(1)) = .true.
        end do

        ! which operators must not be connect with current one?
        ! (any connections between operators for which connection 
        ! was fixed on input, as we here are dealing with the
        ! additional connections, only. If a 0-contraction was
        ! given on input, any connection between these operators
        ! is avoided effectively.
        must_not_connect(1:ivtx-1) = .false.
        do iarc = 1, proto%narc
          if (proto%arc(iarc)%link(2).eq.ivtx .and.
     &        (strict.and.proto%arc(iarc)%occ_cnt(1,1).ge.0 .or.
     &    .not.strict.and.iocc_zero(proto%arc(iarc)%occ_cnt)) )
     &         must_not_connect(proto%arc(iarc)%link(1)) = .true.
        end do

        ! a must_connect overrides must_not_connect
        do jvtx = 1, ivtx-1
          must_not_connect(jvtx) = must_not_connect(jvtx).and.
     &                            .not.must_connect(jvtx)
        end do

        ! also, we should avoid contractions between open lines:
        if (ol_map(ivtx).ne.0) then
          do jvtx = 1, ivtx-1
            must_not_connect(jvtx) =
     &           must_not_connect(jvtx).or.ol_map(jvtx).ne.0
          end do
        end if

        ! count number of previous operators to which connection is
        ! possible
        ok = .true.
        nvtx_prev = 0
        do jvtx = 1, ivtx-1
          if (iocc_nonzero(occ_ol_vtx(1:ngastp,1:2,jvtx))) then
            nvtx_prev = nvtx_prev+1
            ivtx_reo(nvtx_prev) = jvtx
          else if (must_connect(jvtx)) then
            ok = .false.
          end if
        end do
c dbg
c        print *,'ok 3: ',ok
c        print *,'must_connect:     ',must_connect(1:ivtx-1)
c        print *,'must_not_connect: ',must_not_connect(1:ivtx-1)
c dbg

        if (ok) then
          ! (we might set up occ_cnt2prev_min here)

          ! loop over all possible connections
          ! we start with the unconnected case
          occ_cnt2prev_min = 0
          allocate(occ_conn(ngastp,2,max(1,nvtx_prev)))
          occ_cnt2prev = occ_cnt2prev_min
          cnn_dist: do
            if (ntest.ge.100) then
              write(luout,*) 'occ_cnt2prev: (nvtx_prev = ',nvtx_prev,')'
              call wrt_occ(luout,occ_cnt2prev)
            end if

            init = .true.
            connect: do while(next_part_connection(occ_conn,
     &           init,occ_cnt2prev,max(1,nvtx_prev),
     &           occ_cnt2prev_min,occ_cnt2prev))
              init = .false.

              if (ntest.ge.100) then
                write(luout,*) 'current partition:'
                call wrt_occ_n(luout,occ_conn,max(1,nvtx_prev))
              end if

              ! check connections
              ! (a) do connections fit to vertices?
              ok = .true.
              do jvtx = 1, nvtx_prev
                ok = ok.and.
     &               iocc_bound('<=',
c     &               occ_conn(1:ngastp,1:2,ivtx_reo(jvtx)),.false.,
c     &               occ_ol_vtx(1:ngastp,1:2,jvtx),.false.)
     &               occ_conn(1:ngastp,1:2,jvtx),.false.,
     &               occ_ol_vtx(1:ngastp,1:2,ivtx_reo(jvtx)),.false.)
              end do
c dbg
c              print *,'check (a): ',ok
c dbg

              ! (b) are all operators connected that must be connected?
              if (ok) then
                do jvtx = 1, nvtx_prev
                  if (must_connect(ivtx_reo(jvtx)))
     &               ok = ok.and.
     &               iocc_nonzero(occ_conn(1:ngastp,1:2,jvtx))
                end do
              end if
c dbg
c              print *,'check (b): ',ok
c dbg

              ! (c) no connections between operators for which
              !  connection was fixed (exept we generated a 
              ! 0-contraction, which will be ignored below anyway)
              if (ok) then
                do jvtx = 1, nvtx_prev
                  if (must_not_connect(ivtx_reo(jvtx)))
     &                 ok = ok.and..not.
     &                     iocc_nonzero(occ_conn(1:ngastp,1:2,jvtx))
                end do
              end if
c dbg
c              print *,'check (c): ',ok
c dbg

              ! set up prototype
              if (ok) then
                ! create a new instance, as some of the information
                ! on proto must be retained for further "connect" loops 
                call init_contr(proto_new)
                call copy_contr(proto,proto_new)

                narc_old = proto_new%narc
                do jvtx = 1, nvtx_prev
                  kvtx = ivtx_reo(jvtx)
                  ! loop over existing arcs and replace proto-arc
                  ! by actual connection
                  ok = .false.
                  do iarc = 1, narc_old
                    if (proto_new%arc(iarc)%link(1).eq.kvtx .and.
     &                  proto_new%arc(iarc)%link(2).eq.ivtx) then
c dbg
c                      print *,'(a) add arc at ',iarc
c                      call wrt_occ(luout,occ_conn(1:ngastp,1:2,jvtx))
c dbg
                      ! only if this was a proto-connection (indicated
                      ! by a negative number on first entry)
                      ! else: do nothing, as we have above checked that
                      ! only an additional 0-connection remains
c                      if (proto_new%arc(iarc)%occ_cnt(1,1).lt.0)
c     &                     proto_new%arc(iarc)%occ_cnt =
c     &                     occ_conn(1:ngastp,1:2,jvtx)
                      ! do not completely understand that behaviour,
                      ! trying instead:
                      if (proto_new%arc(iarc)%occ_cnt(1,1).lt.0) then
                        proto_new%arc(iarc)%occ_cnt =
     &                     occ_conn(1:ngastp,1:2,jvtx)
                      else
                        proto_new%arc(iarc)%occ_cnt =
     &                       proto_new%arc(iarc)%occ_cnt +
     &                       occ_conn(1:ngastp,1:2,jvtx)
                      end if
                      ok = .true.
                    end if
                  end do
                  
                  ! in a few cases two or more contractions between
                  ! the same vertices end up on different arcs:
                  ! merge these and adapt narc accordingly
                  if (proto_new%narc.gt.0)
     &               call contr_clean_arcs(proto_new%arc,proto_new%narc)

                  if (.not.ok.and.
     &                 iocc_nonzero(occ_conn(1:ngastp,1:2,jvtx))) then
                    call resize_contr(proto_new,nvtx,
     &                   proto_new%narc + 1,0,0)
                    proto_new%narc = proto_new%narc + 1
                    iarc = proto_new%narc
c dbg
c                      print *,'(b) add arc at ',iarc
c                      call wrt_occ(luout,occ_conn(1:ngastp,1:2,jvtx))
c dbg
                    proto_new%arc(iarc)%link(1) = kvtx
                    proto_new%arc(iarc)%link(2) = ivtx
                    proto_new%arc(iarc)%occ_cnt =
     &                     occ_conn(1:ngastp,1:2,jvtx)
                  end if
                end do

                if (ivtx.lt.nvtx) then
                  ! next recursion level
                  call gen_contr_rec(ivtx+1,proto_new)!,
c     &               occ_ol_prev,occ_ol_rem,occ_ol_vtx)
                  if (ntest.ge.100) then
                    write(luout,*) '- back in level ',ivtx
                  end if
                else
                  ! final check for validity:
                  ! no (de)excitations do remain?
                  ! [EX_current] + [EX_prev] - [CNT2prev_max^+] == 0
                  occ_sum = occ_ex + iocc_xdn(1,occ_ol_prev) -
     &                 iocc_dagger(occ_cnt2prev)
                  ok = .not.iocc_nonzero(occ_sum)
c dbg
c                  call wrt_occ4(luout,occ_ex,iocc_xdn(1,occ_ol_prev),
c     &                 iocc_dagger(occ_cnt2prev),occ_tgt_ex)
c                  call wrt_occ(luout,occ_sum)
c                  print *,'EX : ',ok
c dbg                  
                  if (ok) then
                    ! [DX_current] + [DX_prev] - [CNT2prev_max] == 0
                    occ_sum = occ_dx + iocc_xdn(2,occ_ol_prev) -
     &                   occ_cnt2prev
                    ok = .not.iocc_nonzero(occ_sum)
c dbg
c                    call wrt_occ4(luout,occ_dx,iocc_xdn(2,occ_ol_prev),
c     &                   occ_tgt_dx,occ_sum)
c                    print *,'DX : ',ok
c dbg                  
                  end if

                  if (ok) then

                    ! make topological analysis
                    call topo_contr(ieqvfac,reo,ivtx_reo2,
     &                   proto_new,occ_vtx,fix_vtx)

                    ok = ieqvfac.gt.0
c dbg
c                    print *,'check ieqvfac: ',ok,ieqvfac
c dbg
                  end if

                  ! very last check:
                  ! all (0 x [C]) and (x 0 [C]) prototype arcs
                  ! must be present (as such or partitioned)
                  if (ok) then
c dbg
c                    print *,'+----+'
c                    call prt_contr2(luout,proto_new,op_info)
c                    call prt_contr2(luout,proto_main,op_info)
c                    print *,'+----+'
c dbg
                    ok = check_contr(proto_new,proto_main)

                  end if
c dbg
c                  print *,'final: ',ok
c dbg

                  if (ok) then

                    if (ntest.ge.100) then
                      write(luout,*) 'new contraction found!'
                    end if

                    ! ensure canonical order
                    call canon_contr(proto_new,reo,ivtx_reo2)
                    
                    ! reorder ol_map as well
                    call perm_mult(ol_map_reo,ol_map,ivtx_reo2,nvtx)

                    ! apply prefactor from equivalent permutation of vertices
                    proto_new%fac = proto_new%fac/dble(ieqvfac)

                    ! store contraction in formula list
                    ! and advance formula pointer
                    call new_formula_item(
     &                   form_pnt,command_add_contribution,
     &                   proto_new%idx_res)

                    ! transform to actual contraction form
                    ! i.e. remove external line vertices (if any)
c dbg
c                    do kvtx = 1, nvtx
c                      occ_test(1:ngastp,1:2,kvtx) =
c     &                     occ_vtx(1:ngastp,1:2,ivtx_reo2(kvtx))
c                    end do
c dbg
                    call proto2contr(form_pnt%contr,proto_new,
     &                   ol_map_reo,occ_test,op_info)

                    form_pnt => form_pnt%next

                  end if
                end if

                call dealloc_contr(proto_new)
              end if

            end do connect
            if (.not.next_dist2(occ_cnt2prev,2*ngastp,
     &           occ_cnt2prev_min,occ_cnt2prev_max,+1)) exit cnn_dist
          end do cnn_dist
          deallocate(occ_conn)

        end if
      end if
      call dealloc_contr(proto)
      
      if (ntest.ge.100) then
        write(luout,*) 'returning from level ',ivtx
      end if

      return
      end subroutine

      subroutine gen_contr4_unconn(occ_ol_prev,occ_ol_rem,occ_ol_vtx,
     &     ivtx,proto,occ_vtx)

      implicit none

      integer, intent(in) ::
     &     ivtx,
     &     occ_vtx(ngastp,2,*)
      integer, intent(out) ::
     &     occ_ol_prev(ngastp,2), occ_ol_rem(ngastp,2),
     &     occ_ol_vtx(ngastp,2,*)
      type(contraction), intent(in) ::
     &     proto

      type(cntr_arc), pointer ::
     &     arc(:)

      integer ::
     &     idx1, idx2, jvtx, iarc

      ! initialize with vertex occupations
      do jvtx = 1, nvtx
        occ_ol_vtx(1:ngastp,1:2,jvtx) = occ_vtx(1:ngastp,1:2,jvtx)
      end do

      ! loop over all arcs and remove contractions
      arc => proto%arc
      do iarc = 1, proto%narc
        ! ignore proto-type arcs
        if (arc(iarc)%occ_cnt(1,1).lt.0) cycle
        idx1 = arc(iarc)%link(1)
        idx2 = arc(iarc)%link(2)
        ! -- dto. --
        if (idx1.le.0.or.idx2.le.0) cycle
        occ_ol_vtx(1:ngastp,1:2,idx1) = occ_ol_vtx(1:ngastp,1:2,idx1)
     &       - arc(iarc)%occ_cnt
        occ_ol_vtx(1:ngastp,1:2,idx2) = occ_ol_vtx(1:ngastp,1:2,idx2)
     &       - iocc_dagger(arc(iarc)%occ_cnt)
      end do

      if (ivtx.le.0) return

      ! accumulate all vertices up to ivtx-1:
      occ_ol_prev = 0
      do jvtx = 1, ivtx-1
        occ_ol_prev = occ_ol_prev + occ_ol_vtx(1:ngastp,1:2,jvtx)
      end do

      ! accumulate all vertices from ivtx+1 to nvtx:
      occ_ol_rem = 0
      do jvtx = ivtx+1, nvtx
        occ_ol_rem = occ_ol_rem + occ_ol_vtx(1:ngastp,1:2,jvtx)
      end do

      return
      end subroutine

      end

