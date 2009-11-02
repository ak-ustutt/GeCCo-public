*----------------------------------------------------------------------*
      subroutine get_arc_list(arc_list,len_list,contr,orb_info)
*----------------------------------------------------------------------*
*     get list of (non-redundant) arcs, ordered according to
*     contraction strength (descending)
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     len_list, arc_list(*)
      type(contraction), intent(in), target ::
     &     contr
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     iarc, igas, idx, isvarc, ica, hpvx, idx_sv1, idx_sv2
      integer ::
     &     weight(ngastp)
      type(cntr_arc), pointer ::
     &     arc(:)
      integer, pointer ::
     &     narc, nsym, svertex(:), ihpvgas(:,:), igassh(:,:),
     &     sv_arc(:,:), cnt_strength(:)

      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'get_arc_list')

      weight(1:ngastp) = 0
      ihpvgas => orb_info%ihpvgas
      igassh  => orb_info%igassh
      nsym => orb_info%nsym
      do igas = 1, orb_info%ngas
        weight(ihpvgas(igas,1)) = weight(ihpvgas(igas,1))+
     &       sum(igassh(1:nsym,igas))
      end do

      if (ntest.ge.100) then
        write(luout,*) 'weights: ',weight(1:ngastp)
      end if

      narc => contr%narc

      allocate(sv_arc(2,narc),cnt_strength(narc))

      arc => contr%arc
      svertex => contr%svertex
      len_list = 0
      do iarc = 1, narc
        idx_sv1 = svertex(arc(iarc)%link(1))
        idx_sv2 = svertex(arc(iarc)%link(2))
        idx = 0
cmh  if there are still null arcs within a supervertix, try this:
cmh        if (idx_sv1.eq.idx_sv2.and.
cmh     &      sum(arc(iarc)%occ_cnt(1:ngastp,1:2)).eq.0) cycle
        do isvarc = 1, len_list
          if ((sv_arc(1,isvarc).eq.idx_sv1.and.
     &         sv_arc(2,isvarc).eq.idx_sv2) .or.
     &        (sv_arc(1,isvarc).eq.idx_sv2.and.
     &         sv_arc(2,isvarc).eq.idx_sv1)) then
            idx = isvarc
          end if
        end do
        if (idx.eq.0) then
          len_list = len_list+1
          idx = len_list
          sv_arc(1,idx) = idx_sv1
          sv_arc(2,idx) = idx_sv2
          cnt_strength(idx) = 0
          arc_list(idx) = iarc
        end if
        do ica = 1, 2
          do hpvx = 1, ngastp
            cnt_strength(idx) = cnt_strength(idx) +
     &           arc(iarc)%occ_cnt(hpvx,ica)*weight(hpvx)
          end do
        end do
      end do

      if (ntest.ge.100) then
        write(luout,*) 'len_list = ',len_list
        write(luout,*) 'raw:'
        do idx = 1, len_list
          write(luout,'(4x,3i5)') idx,cnt_strength(idx),arc_list(idx)
        end do
      end if

      call idxsort(cnt_strength,arc_list,len_list,-1)

      if (ntest.ge.100) then
        write(luout,*) 'sorted:'
        do idx = 1, len_list
          write(luout,'(4x,3i5)') idx,cnt_strength(idx),arc_list(idx)
        end do
      end if

      deallocate(sv_arc,cnt_strength)

      return
      end
