*----------------------------------------------------------------------*
      subroutine join_contr(contr_abc,contr_ac,contr_b,
     &     idxop_abc,iblk_abc,op_info)
*----------------------------------------------------------------------*
*     join contractions ac and b into one with result occupation abc
*     
*     will work only, if this leads to a unique result
*     should be replaced by a routine to handle contractions of
*     two (or more) operators
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 000

      type(contraction), intent(out) ::
     &     contr_abc

      type(contraction), intent(in) ::
     &     contr_ac, contr_b
      integer, intent(in) ::
     &     idxop_abc, iblk_abc
      type(operator_info), intent(in) ::
     &     op_info
 
      integer ::
     &     nvtx_abc, nvtx_ac, nvtx_a, nvtx_b, nvtx_c,
     &     narc_abc, narc_ac, narc_b, 
     &     idx, idx_abc, iarc,
     &     nconnect_a, nconnect_c,
     &     nconnect_b_dx, nconnect_b_ex
      integer ::
     &     iocc(ngastp,2), jocc(ngastp,2),
     &     iocc_a_dx(ngastp,2), iocc_b_dx(ngastp,2),
     &     iocc_b_ex(ngastp,2), iocc_c_ex(ngastp,2)
     
      integer, pointer ::
     &     ivtx_a(:), ivtx_b_dx(:), ivtx_b_ex(:), ivtx_c(:),
     &     ivtx_ac_reo(:), ivtx_b_reo(:)

      type(cntr_arc), pointer ::
     &     arc(:)

      if (ntest.ge.100) then
        write(luout,*) '====================' 
        write(luout,*) ' This is join_contr'
        write(luout,*) '====================' 
        write(luout,*) 'joining: AC, B'
        call prt_contr2(luout,contr_ac,op_info%op_arr)
        call prt_contr2(luout,contr_b,op_info%op_arr)
      end if

      nvtx_ac = contr_ac%nvtx
      narc_ac = contr_ac%narc
      ! count all vertices up to first "0" entry
      nvtx_a = 0
      do idx = 1, nvtx_ac
        if (contr_ac%vertex(idx)%idx_op.ne.0) then
          nvtx_a = nvtx_a+1
c          vtxlist(idx) = idx
        else
          exit
        end if
      end do

      nvtx_c = nvtx_ac-nvtx_a-1

      if (contr_b%nvtx.le.0)
     &     call quit(1,'join_contr',
     &     'inserted contraction fragment must be non-empty')

      nvtx_b = contr_b%nvtx
      narc_b = contr_b%narc

      nvtx_abc = nvtx_a+nvtx_b+nvtx_c
      narc_abc = narc_ac+narc_b+nvtx_a*(nvtx_b-1)+nvtx_c*(nvtx_b-1)

      if (ntest.ge.1000) then
        write(luout,*) 'nvtx_a, nvtx_b, nvtx_c: ',nvtx_a,nvtx_b,nvtx_c
        write(luout,*) 'narc_ac, narc_b, narc_abc: ',
     &       narc_ac, narc_b, narc_abc
      end if

      call resize_contr(contr_abc,nvtx_abc,narc_abc,0)

      if (nvtx_ac.gt.0) allocate(ivtx_ac_reo(nvtx_ac))
      if (nvtx_b.gt.0) allocate(ivtx_b_reo(nvtx_b))

      ! set prefactor
      contr_abc%fac = contr_ac%fac*contr_b%fac
      ! set result
      contr_abc%idx_res = idxop_abc
      contr_abc%iblk_res = iblk_abc

      ! set vertices and reordering arrays
      ! collect vertices from A
      do idx = 1, nvtx_a
        contr_abc%vertex(idx) = contr_ac%vertex(idx)
        ivtx_ac_reo(idx) = idx
      end do
      ! collect vertices from B
      idx_abc = nvtx_a
      do idx = 1, nvtx_b
        idx_abc = idx_abc+1
        contr_abc%vertex(idx_abc) = contr_b%vertex(idx)
        ivtx_b_reo(idx) = idx_abc        
      end do
      ! collect vertices from C
      do idx = nvtx_a+2,nvtx_ac
        idx_abc = idx_abc+1
        contr_abc%vertex(idx_abc) = contr_ac%vertex(idx)
        ivtx_ac_reo(idx) = idx_abc
      end do

      iarc = 0
      ! add all arcs from A and C, except the external ones
      do idx = 1, narc_ac
        if (contr_ac%arc(idx)%link(1).eq.0.or.
     &      contr_ac%arc(idx)%link(2).eq.0) cycle
        iarc = iarc + 1
        contr_abc%arc(idx)%link(1) =
     &       ivtx_ac_reo(contr_ac%arc(idx)%link(1))
        contr_abc%arc(idx)%link(2) =
     &       ivtx_ac_reo(contr_ac%arc(idx)%link(2))
        contr_abc%arc(idx)%occ_cnt =
     &       contr_ac%arc(idx)%occ_cnt
      end do
      ! add all arcs from B
      idx_abc = iarc !narc_ac
      do idx = 1, narc_b
        idx_abc = idx_abc+1
        iarc = iarc + 1
        contr_abc%arc(idx_abc)%link(1) =
     &       ivtx_b_reo(contr_b%arc(idx)%link(1))
        contr_abc%arc(idx_abc)%link(2) =
     &       ivtx_b_reo(contr_b%arc(idx)%link(2))
        contr_abc%arc(idx_abc)%occ_cnt =
     &       contr_b%arc(idx)%occ_cnt
      end do

      contr_abc%nvtx = nvtx_abc
      ! preliminary setting
      contr_abc%narc = iarc
      contr_abc%nfac = 0

      ! ===========================================================
      ! add new arcs:
      ! ===========================================================
      if (nvtx_a.gt.0) allocate(ivtx_a(nvtx_a))
      if (nvtx_b.gt.0) allocate(ivtx_b_ex(nvtx_b))
      if (nvtx_b.gt.0) allocate(ivtx_b_dx(nvtx_b))
      if (nvtx_c.gt.0) allocate(ivtx_c(nvtx_c))

      ! loop over all A vertices and
      ! get a list of vertices with unconnected deexcitations
      nconnect_a = 0
      iocc_a_dx = 0
c      do idx = 1, nvtx_a
c        ! get unconnected part ...
c        call get_unconnected4vertex(iocc,idx,contr_abc,op_info)
c        ! ... but deexcitation part only:
c        iocc = iocc_xdn(2,iocc)
c        if (iocc_nonzero(iocc)) then
c          nconnect_a = nconnect_a+1
c          ivtx_a(nconnect_a) = idx
c          iocc_a_dx = iocc_a_dx + iocc
c        end if
c      end do
      do idx = 1, narc_ac
        if (contr_ac%arc(idx)%link(2).eq.0) then
          nconnect_a = nconnect_a+1
          iocc_a_dx = iocc_a_dx + contr_ac%arc(idx)%occ_cnt
        end if
      end do

      ! loop over B vertices and
      ! get lists of vertices with unconnected excitations 
      ! and deexcitations (separately)
      nconnect_b_ex = 0
      nconnect_b_dx = 0
      iocc_b_ex = 0
      iocc_b_dx = 0
      do idx = nvtx_a+1, nvtx_a+nvtx_b
        call get_unconnected4vertex(iocc,idx,contr_abc,op_info)
        jocc = iocc_xdn(1,iocc)
        iocc = iocc_xdn(2,iocc)
        if (iocc_nonzero(jocc)) then
          nconnect_b_ex = nconnect_b_ex+1
          iocc_b_ex = iocc_b_ex + jocc
          ivtx_b_ex(nconnect_b_ex) = idx
        end if
        if (iocc_nonzero(iocc)) then
          nconnect_b_dx = nconnect_b_dx+1
          iocc_b_dx = iocc_b_dx + iocc
          ivtx_b_dx(nconnect_b_dx) = idx
        end if
      end do
      ! loop over C vertices and 
      ! get a list of vertices with unconnected excitations
      nconnect_c = 0
      iocc_c_ex = 0
c      do idx = nvtx_a+nvtx_b+1, nvtx_a+nvtx_b+nvtx_c
c        ! get unconnected part ...
c        call get_unconnected4vertex(iocc,idx,contr_abc,op_info)
c        ! ... but excitation part only:
c        iocc = iocc_xdn(1,iocc)
c        if (iocc_nonzero(iocc)) then
c          nconnect_c = nconnect_c+1
c          iocc_c_ex = iocc_c_ex + iocc
c          ivtx_c(nconnect_c) = idx
c        end if
c      end do
      do idx = 1, narc_ac
        if (contr_ac%arc(idx)%link(1).eq.0) then
          nconnect_c = nconnect_c+1
          iocc_c_ex = iocc_c_ex + iocc_dagger(contr_ac%arc(idx)%occ_cnt)
        end if
      end do
      
      ! -----------------------------------------------
      ! finally, make connections between A and B ....
      ! -----------------------------------------------
      arc => contr_abc%arc
      !iarc = narc_ac + narc_b 
      
      if (nconnect_a*nconnect_b_ex.gt.0) then

        ! one A vertex
        if (nconnect_a.eq.1) then
          ! get maximum connection
          iocc = iocc_overlap(iocc_a_dx,.false.,iocc_b_ex,.true.)
          if (nconnect_b_ex.eq.1) then
            iarc = iarc+1
            arc(iarc)%link(1) = ivtx_a(1)
            arc(iarc)%link(2) = ivtx_b_ex(1)
            arc(iarc)%occ_cnt = iocc
          else
            ! not equal to iocc_b_ex?
            if (.not.iocc_equal(iocc,.false.,iocc_b_ex,.true.))
     &           call quit(1,'join_contr',
     &       'cannot handle non-unique re-combinations (1)')
            ! make new arcs
            do idx = 1, nconnect_b_ex
              iarc = iarc+1
              arc(iarc)%link(1) = ivtx_a(1)
              arc(iarc)%link(2) = ivtx_b_ex(idx)
              call get_unconnected4vertex
     &             (iocc,ivtx_b_ex(idx),contr_abc,op_info)
              iocc = iocc_xdn(1,iocc)
              arc(iarc)%occ_cnt = iocc_dagger(iocc)
            end do
          end if
        else if (nconnect_b_ex.eq.1) then
          ! we simply reuse the external connections stored on contr_ac
          do idx = 1, narc_ac
            if (contr_ac%arc(idx)%link(2).eq.0) then
              iarc = iarc+1
              arc(iarc)%link(1) = ivtx_ac_reo(contr_ac%arc(idx)%link(1))
              arc(iarc)%link(2) = ivtx_b_ex(1)
              arc(iarc)%occ_cnt = contr_ac%arc(idx)%occ_cnt
            end if
          end do
c          ! get maximum connection
c          iocc = iocc_overlap(iocc_a_dx,.false.,iocc_b_ex,.true.)
c          ! A is more than one node, so we currently must insist
c          ! that all A lines are connected to get a unique result:
c          if (.not.iocc_equal(iocc_a_dx,.false.,iocc,.false.))
c     &       call quit(1,'join_contr',
c     &       'cannot handle non-unique re-combinations (2)')
c          ! make new arcs
c          do idx = 1, nconnect_a
c            iarc = iarc+1
c            arc(iarc)%link(1) = ivtx_a(idx)
c            arc(iarc)%link(2) = ivtx_b_ex(1)
c            call get_unconnected4vertex
c     &           (iocc,ivtx_a(idx),contr_abc,op_info)
c            iocc = iocc_xdn(2,iocc)
c            arc(iarc)%occ_cnt = iocc
c          end do
        
        else if (nconnect_a.gt.1.and.nconnect_b_ex.gt.1) then
          call quit(1,'join_contr',
     &       'cannot handle non-unique re-combinations (3)')
        end if

      end if
      ! -----------------------------------------------
      ! and between B and C ....
      ! -----------------------------------------------

      if (nconnect_b_dx*nconnect_c.gt.0) then

        ! one B vertex
        if (nconnect_b_dx.eq.1) then
          ! we may simply reuse the external connection on contr_ac
          do idx = 1, narc_ac
            if (contr_ac%arc(idx)%link(1).eq.0) then
              iarc = iarc+1
              arc(iarc)%link(1) = ivtx_b_dx(1)
              arc(iarc)%link(2) = ivtx_ac_reo(contr_ac%arc(idx)%link(2))
              arc(iarc)%occ_cnt = contr_ac%arc(idx)%occ_cnt
            end if
          end do
c          ! get maximum connection
c          iocc = iocc_overlap(iocc_b_dx,.false.,iocc_c_ex,.true.)
c          if (nconnect_c.eq.1) then
c            iarc = iarc+1
c            arc(iarc)%link(1) = ivtx_b_dx(1)
c            arc(iarc)%link(2) = ivtx_c(1)
c            arc(iarc)%occ_cnt = iocc
c          else
c          ! not equal to iocc_c?
c          if (.not.iocc_equal(iocc,.false.,iocc_c_ex,.true.))
c     &       call quit(1,'join_contr',
c     &       'cannot handle non-unique re-combinations (4)')
c          ! make new arcs
c          do idx = 1, nconnect_c
c            iarc = iarc+1
c            arc(iarc)%link(1) = ivtx_b_dx(1)
c            arc(iarc)%link(2) = ivtx_c(idx)
c            call get_unconnected4vertex
c     &           (iocc,ivtx_c(idx),contr_abc,op_info)
c            iocc = iocc_xdn(1,iocc)
c            arc(iarc)%occ_cnt = iocc_dagger(iocc)
c          end do
c        end if
        else if (nconnect_c.eq.1) then
          ! get maximum connection
          iocc = iocc_overlap(iocc_b_dx,.false.,iocc_c_ex,.true.)
          ! B is more than one node, so we currently must insist
          ! that all B lines are connected to get a unique result:
          if (.not.iocc_equal(iocc_b_dx,.false.,iocc,.false.))
     &       call quit(1,'join_contr',
     &       'cannot handle non-unique re-combinations (5)')
          ! make new arcs
          do idx = 1, nconnect_b_dx
            iarc = iarc+1
            arc(iarc)%link(1) = ivtx_b_dx(idx)
            arc(iarc)%link(2) = ivtx_c(1)
            call get_unconnected4vertex
     &           (iocc,ivtx_b_dx(idx),contr_abc,op_info)
            iocc = iocc_xdn(2,iocc)
            arc(iarc)%occ_cnt = iocc
          end do
        
        else if (nconnect_a.gt.1.and.nconnect_b_ex.gt.1) then
          call quit(1,'join_contr',
     &       'cannot handle non-unique re-combinations (6)')
        end if

      end if

      ! set number of arcs
      if (iarc.gt.narc_abc) call quit(1,'join_contr','iarc.gt.narc_abc')
      contr_abc%narc = iarc
      contr_abc%nfac = 0

      ! release reordering arrays
      if (nvtx_a) deallocate(ivtx_a)
      if (nvtx_b) deallocate(ivtx_b_ex)
      if (nvtx_b) deallocate(ivtx_b_dx)
      if (nvtx_c) deallocate(ivtx_c)
      if (nvtx_ac.gt.0) deallocate(ivtx_ac_reo)
      if (nvtx_b.gt.0) deallocate(ivtx_b_reo)

      if (ntest.ge.1000) then
        write(luout,*) 'generated contraction (bef. reorder):'
        call prt_contr2(luout,contr_abc,op_info%op_arr)
      end if

      ! enforce law and order
      call canon_contr(contr_abc)

      if (ntest.ge.100) then
        write(luout,*) 'generated contraction:'
        call prt_contr2(luout,contr_abc,op_info%op_arr)
      end if
      
      return
      end
