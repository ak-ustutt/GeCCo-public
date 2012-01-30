*----------------------------------------------------------------------*
      subroutine get_vtx_in_contr(ivtx,idxop,adjop,njoined,vtxst,contr)
*----------------------------------------------------------------------*
*     return the first occurrence of operator idxop's vertices in contr
*     version of vtx_in_contr for super-vertices
*     you should still use vtx_in_contr() for testing whether the
*     operator is contained at all
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     idxop, njoined, vtxst
      integer, intent(out) ::
     &     ivtx(njoined)
      logical, intent(in) ::
     &     adjop

      logical ::
     &     first
      integer ::
     &     ijoin, jvtx, iblk_next

      first = .true.
      ijoin = 0
      do jvtx = vtxst, contr%nvtx
        if (first.and.contr%vertex(jvtx)%idx_op.eq.idxop.and.
     &               (contr%vertex(jvtx)%dagger.eqv.adjop)) then
          ivtx(1) = jvtx
          ijoin = 1
          if (ijoin.eq.njoined) exit
          iblk_next = contr%vertex(jvtx)%iblk_op+1
          first = .false.
        else if (contr%vertex(jvtx)%idx_op.eq.idxop .and.
     &          (contr%vertex(jvtx)%dagger.eqv.adjop) .and.
     &           contr%vertex(jvtx)%iblk_op.eq.iblk_next) then
          ijoin = ijoin+1
          ivtx(ijoin) = jvtx
          if (ijoin.eq.njoined) exit
          iblk_next = iblk_next+1
        end if
      end do

      if (ijoin.ne.njoined) then
        write(luout,*) 'idxop:        ',idxop
        write(luout,*) 'adjop:        ',adjop
        write(luout,*) 'ijoin,njoined:', ijoin,njoined
        do jvtx = 1, contr%nvtx
          write(luout,*) jvtx,contr%vertex(jvtx)%idx_op,
     &                        contr%vertex(jvtx)%iblk_op
        end do
        call quit(1,'get_vtx_in_contr',
     &     'did not find all vertices for operator')
      end if

      return
      end
