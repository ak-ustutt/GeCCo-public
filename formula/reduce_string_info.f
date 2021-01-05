*----------------------------------------------------------------------*
      subroutine reduce_string_info(contr_red,contr,
     &              nvtx, nvtx_new,arc_list,nlist,ireo_vtx)
*----------------------------------------------------------------------*
*
*     update index info on contr_red
*     - remove contracted indices
*     - reorder vertices
*
*     contr is the original contraction, contr_red has been reduced by
*     reduce_contr2()
*
*     this routine is setup for being called withing reduce_contr2()
*
*     andreas, dec. 2020
*     
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      character(len=18), parameter ::
     &     i_am = 'reduce_string_info'
      integer, parameter ::
     &     ntest = 100

      type(contraction), intent(in) ::
     &     contr
      type(contraction), intent(inout) ::
     &     contr_red
      integer, intent(in) ::
     &     nvtx, nvtx_new, nlist,
     &     arc_list(nlist), ireo_vtx(nvtx)

      type(string_element), pointer ::
     &     string_tmp(:)
      integer, pointer ::
     &     ireo_arc(:)

      logical ::
     &     on_list
      integer ::
     &     nidx, nidx_new, ivtx, idx, jdx, ilist, iarc, jarc,
     &     narc, narc_red, sign

      logical, external ::
     &     iocc_equal

!     make a copy
      nidx = contr%nidx
      allocate(string_tmp(nidx))

      string_tmp(1:nidx) = contr%contr_string(1:nidx)

      if (ntest>=100) then
        write(lulog,*) 'initial'
        call print_string_idx(string_tmp,nidx,nvtx)
      end if

!     we need a reodering array for the arcs
      narc = contr%narc
      narc_red = contr_red%narc

      allocate(ireo_arc(narc))
      ireo_arc = 0
      do iarc = 1, narc
        on_list = .false.
        do ilist = 1, nlist
          on_list = on_list.or.iarc.eq.arc_list(ilist)
        end do
        if (on_list) cycle      ! ignore arcs that are removed anyway
        do jarc = 1, narc_red
          if (ireo_vtx(contr%arc(iarc)%link(1))
     &         .eq.contr_red%arc(jarc)%link(1) .and.
     &        ireo_vtx(contr%arc(iarc)%link(2))
     &         .eq.contr_red%arc(jarc)%link(2)
c     it seems unnecessary to compare this (and can be wrong if arc were merged
c     &         .and.
c     &        iocc_equal(contr%arc(iarc)%occ_cnt,.false.,
c     &               contr_red%arc(jarc)%occ_cnt,.false.)
     &         ) then
            ireo_arc(iarc) = jarc
          end if
        end do
      end do
      if (ntest>=100) then
        write(lulog,'(1x,"arc reo: ",10i4)') ireo_arc(1:narc)
      end if
      
      nidx_new = nidx
!     mark all elements removed by contraction
      do idx = 1, nidx
        if (string_tmp(idx)%ext) cycle
        on_list = .false.
        do ilist = 1, nlist
          on_list = on_list.or.arc_list(ilist).eq.string_tmp(idx)%cnt
        end do
        if (on_list) then
          string_tmp(idx)%del = .true.
          nidx_new = nidx_new-1
        end if
      end do
      if (ntest>=100) then
        write(lulog,*) 'arcs removed, nidx_new = ',nidx_new
        call print_string_idx(string_tmp,nidx,nvtx)
        write(lulog,*) 'now reordering:'
        write(lulog,'(1x,"vtx reo: ",10i4)') ireo_vtx
      end if

!     assign new vertex and arc number
      do idx = 1, nidx
!        write(lulog,*) idx, string_tmp(idx)%vtx,
!     &       ireo_vtx(string_tmp(idx)%vtx)
        string_tmp(idx)%vtx = ireo_vtx(string_tmp(idx)%vtx)
        if (.not.string_tmp(idx)%ext)
     &       string_tmp(idx)%cnt = ireo_arc(string_tmp(idx)%cnt)
      end do
      if (ntest>=100) then
        write(lulog,'(1x,"new vtx: ",10i4)') string_tmp(1:nidx)%vtx
      end if
      
!     init new string_info on contr_red
      if (contr_red%index_info) then
        if (associated(contr_red%contr_string))
     &       deallocate(contr_red%contr_string)
        if (associated(contr_red%result_string))
     &       deallocate(contr_red%contr_string)
      end if
      allocate(contr_red%contr_string(nidx_new))
      jdx = 0
c dbg
c      write(lulog,*) 'nvtx_new,nidx,nidx_new: '
c     &     ,nvtx_new,nidx,nidx_new
c dbg
      do ivtx = 1, nvtx_new
        do idx = 1, nidx
          if (.not.string_tmp(idx)%del.and.
     &         string_tmp(idx)%vtx==ivtx) then
            jdx = jdx+1
c dbg
c            write(lulog,*) ivtx,idx,jdx
c dbg           
            contr_red%contr_string(jdx) = string_tmp(idx)
          end if
        end do
      end do
      contr_red%nidx = nidx_new
      
!     we sort all vertices, no implications for overall sign of diagram
      idx = 1
      sign = 0
      do ivtx = 1, nvtx_new
        if (idx>nidx_new) exit
        nidx = 0
        jdx = idx
        do while ( contr_red%contr_string(idx)%vtx==ivtx )
          nidx = nidx+1
          idx = idx+1
          if (idx>nidx_new) exit
        end do
c     dbg
c        write(lulog,*) 'ivtx,idx,jdx,nidx: ',ivtx,idx,jdx,nidx
c     dbg
        call sort_string(sign,contr_red%contr_string(jdx),nidx)
      end do
             
      deallocate(string_tmp,ireo_arc)

!     also copy final result info
      if (contr%nxidx>0) then
        allocate(contr_red%result_string(contr%nxidx))
        contr_red%result_string = contr%result_string
        contr_red%nxidx = contr%nxidx
      end if
      contr_red%total_sign = contr%total_sign

      contr_red%index_info = .true.
      
      end

     
