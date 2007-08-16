*----------------------------------------------------------------------*
      logical function merge_check(contr,ivtxsuper1,ivtxsuper2,
     &     ivtx1,ivtx2,iocc_ex1,iocc_ex2,topomap)
*----------------------------------------------------------------------*
*     check whether we REALLY may merge the given primitive vertices 
*     (i.e. antisymmetrize equal open lines)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     ivtxsuper1,ivtxsuper2,
     &     ivtx1, ivtx2, iocc_ex1(ngastp,2), iocc_ex2(ngastp,2),
c     &     iocc_res(ngastp,2,njoined),
     &     topomap(contr%nvtx,contr%nvtx)

      logical ::
     &     cnt_dag13, cnt_dag23
      integer ::
     &     nvtx, narc, iarc, ica, hpvx, nstr1, nstr2, ivtx3,
     &     ivtx_lookup13_1, ivtx_lookup13_2,
     &     ivtx_lookup23_1, ivtx_lookup23_2, ncnt1, ncnt2
      type(cntr_arc), pointer ::
     &     arc(:)

      nvtx = contr%nvtx
      narc = contr%narc
      arc => contr%arc

      merge_check = .true.

      outer_loop: do ica = 1, 2
        do hpvx = 1, ngastp
          ! first check: any non-zero entries for current hpvx?
          if (iocc_ex1(hpvx,ica).eq.0.or.
     &        iocc_ex2(hpvx,ica).eq.0) cycle

          nstr1 = iocc_ex1(hpvx,ica)
          nstr2 = iocc_ex2(hpvx,ica)
          ! if so: check whether these indices are going to be fully
          ! contracted to the same operator (or whether they are
          ! open lines)
          ! loop over remaining operators in contraction:
          vtx_loop: do ivtx3 = 1, nvtx
            if (ivtx3.eq.ivtx1.or.ivtx3.eq.ivtx2) cycle
            ! look at topomap
            if (topomap(ivtx3,ivtx1).eq.0.and.topomap(ivtx3,ivtx2).eq.0)
     &           cycle
            ! must not belong to current supervertices:
            if (contr%svertex(ivtx3).eq.ivtxsuper1.or.
     &          contr%svertex(ivtx3).eq.ivtxsuper2) cycle
            ! look for actual connection to this node
            ivtx_lookup13_1 = min(ivtx1,ivtx3)
            ivtx_lookup13_2 = max(ivtx1,ivtx3)
            cnt_dag13 = ivtx1.gt.ivtx3
            ivtx_lookup23_1 = min(ivtx2,ivtx3)
            ivtx_lookup23_2 = max(ivtx2,ivtx3)
            cnt_dag23 = ivtx2.gt.ivtx3
            ncnt1 = 0
            ncnt2 = 0
            do iarc = 1, narc
              if (arc(iarc)%link(1).eq.ivtx_lookup13_1.and.
     &            arc(iarc)%link(2).eq.ivtx_lookup13_2) then
                ncnt1 = arc(iarc)%occ_cnt(hpvx,ica)
                if (cnt_dag13) ncnt1 = arc(iarc)%occ_cnt(hpvx,3-ica)
c                found1 = .true.
              else if
     &           (arc(iarc)%link(1).eq.ivtx_lookup23_1.and.
     &            arc(iarc)%link(2).eq.ivtx_lookup23_2) then
                ncnt2 = arc(iarc)%occ_cnt(hpvx,ica)
                if (cnt_dag23) ncnt2 = arc(iarc)%occ_cnt(hpvx,3-ica)
c                found2 = .true.
              end if
            end do
            ! fully contracted? OK, go on with checks
c dbg
c            print *,'ivtx:',ivtx1,ivtx2,ivtx3
c            print *,'ncnt:',ncnt1,ncnt2
c            print *,'nstr:',nstr1,nstr2
c dbg            
            if (ncnt1.eq.nstr1.and.ncnt2.eq.nstr2)
     &           exit vtx_loop
            ! uncontracted? still OK, as well, try next operator
            if (ncnt1.eq.0.and.ncnt2.eq.0) cycle vtx_loop

            ! partially contracted: do not merge:
            merge_check = .false.
            exit outer_loop

          end do vtx_loop

        end do
      end do outer_loop

c dbg
c      print *,'merge_check = ',merge_check,ivtx1,ivtx2
c dbg

      return
      end
