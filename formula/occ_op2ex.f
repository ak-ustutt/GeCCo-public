*----------------------------------------------------------------------*
      subroutine occ_op2ex(iocc_ex,iocc_cnt,merge_map,
     &                     set_cnt,set_map,ld_map,
     &                     iop1or2,iocc_op,njoined,isupvtx_op,
     &                     contr,arc_list,len_list)
*----------------------------------------------------------------------*
*     return the remaining ("external") occupation of super-vertex
*     isupvtx_op; for my convenience you have to provide the readily
*     extracted occupations of all primitive vertices on iocc_op
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     set_map, set_cnt
      integer, intent(in) ::
     &     ld_map, njoined, len_list,
     &     iop1or2,
     &     iocc_op(ngastp,2,njoined),isupvtx_op,arc_list(len_list)
      type(contraction), intent(in) ::
     &     contr
      integer, intent(inout) ::
     &     iocc_cnt(ngastp,2,len_list)
      integer, intent(out) ::
     &     iocc_ex(ngastp,2,njoined), merge_map(ld_map,2,njoined)

      integer ::
     &     iop_other, ijoin, ivtx, ilist, idx, jdx, iarc, ii
      integer, external ::
     &     imltlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'occ_op2ex')
        write(luout,*) 'set_cnt, set_map: ', set_cnt, set_map
        write(luout,*) 'isupvtx_op: ',isupvtx_op
        write(luout,*) 'iop1or2: ',iop1or2
        write(luout,*) 'OP:'
        call wrt_occ_n(luout,iocc_op,njoined)
      end if

      ! iop1or2 == 1: the current supervertex precedes the second, so
      if (iop1or2.eq.1) iop_other = 2
      ! else:
      if (iop1or2.eq.2) iop_other = 1

      ! we start out with the operator occupations
      iocc_ex(1:ngastp,1:2,1:njoined) =
     &     iocc_op(1:ngastp,1:2,1:njoined)
      
      if (set_map) then
        ! reset map
        merge_map = 0
      
        ! this is easy: each ex vertex is related with 
        ! the respective op vertex, so:
        do ijoin = 1, njoined
          merge_map(1,2,ijoin) = ijoin
        end do

      end if

      do ilist = 1, len_list
        iarc = arc_list(ilist)

        if (contr%svertex(contr%arc(iarc)%
     &        link(iop1or2)).eq.isupvtx_op) then
          ivtx = contr%arc(iarc)%link(iop1or2)
          if (set_cnt)
     &         iocc_cnt(1:ngastp,1:2,ilist) = contr%arc(iarc)%occ_cnt
        else if (contr%svertex(contr%arc(iarc)%
     &        link(iop_other)).eq.isupvtx_op) then
          ivtx = contr%arc(iarc)%link(iop_other)
          if (set_cnt)
     &         iocc_cnt(1:ngastp,1:2,ilist) =
     &         iocc_dagger(contr%arc(iarc)%occ_cnt)
        else 
          call quit(1,'occ_op2ex','arc does not contain supervertex?')
        end if

        ! get number of primitive supervertex (position in iocc_ex):
        ! count their number in svertex up to current position
        idx = imltlist(isupvtx_op,contr%svertex,ivtx,1)
c dbg
c        print *,'idx, ivtx: ',idx,ivtx
c dbg
        if (iop1or2.eq.1) then
          iocc_ex(1:ngastp,1:2,idx) = iocc_ex(1:ngastp,1:2,idx)
     &                              - iocc_cnt(1:ngastp,1:2,ilist)
        else
          iocc_ex(1:ngastp,1:2,idx) = iocc_ex(1:ngastp,1:2,idx)
     &                      - iocc_dagger(iocc_cnt(1:ngastp,1:2,ilist))
        end if

        if (set_map) then
          ! store merging info: the current contraction contributes
          ! to vertex idx of the original operator
          ! get a free entry:
          jdx = 1
          do while(merge_map(jdx,1,idx).gt.0.and.jdx.le.ld_map) 
            jdx = jdx+1
          end do
          if (jdx.gt.ld_map)
     &         call quit(1,'occ_op2ex','ld_map too small?')
          merge_map(jdx,1,idx) = ilist
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'EX:'
        call wrt_occ_n(luout,iocc_ex,njoined)
        write(luout,*) 'CNT:'
        call wrt_occ_n(luout,iocc_cnt,len_list)

        write(luout,*) 'merge map:'
        do ivtx = 1, njoined
          do ii = 1, 2
            write(luout,'(4x,i4," <-",i4," :",10i4)')
     &           ivtx,ii,merge_map(1:ld_map,ii,ivtx)
          end do
        end do
      end if

      return
      end
