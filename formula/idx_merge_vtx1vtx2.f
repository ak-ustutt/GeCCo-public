*----------------------------------------------------------------------*
      integer function idx_merge_vtx1vtx2(ivtx1,ivtx2,ivtxsup1,ivtxsup2,
     &     nmvleft,imvleft,
     &     svertex,svmap,topomap,nvtx)
*----------------------------------------------------------------------*
*     check whether vertices may be merged into one after a contraction
*     and return the position in the contraction where the merged 
*     vertex can be positions without violating the contraction 
*     sequence
*     imvleft(1:nmvleft) contains the numbers of vertices that must
*     be moved left of the resulting merged vertex
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ivtx1, ivtx2, ivtxsup1, ivtxsup2,
     &     nvtx, svmap(nvtx), topomap(nvtx,nvtx),
     &     svertex(nvtx)
      integer, intent(out) ::
     &     nmvleft, imvleft(*)

      logical ::
     &     merge, commute
      integer ::
     &     isupres1, ivtx, jvtx, ivtx1rm, ivtx2lm

      integer, external ::
     &     imltlist

      ! first of all: to which vertex of the result do they
      !  contribute open lines?
      merge = svmap(ivtx1).eq.0 .or.       ! either no contribution
     &        svmap(ivtx2).eq.0 .or.
     &        svmap(ivtx1).eq.svmap(ivtx2) ! or to the same vertex

      nmvleft = 0

      if (merge) then
c dbg
c        print *,'initial step OK'
c        print *,'ivtx1, ivtx2: ',ivtx1, ivtx2
c dbg        
        ! are we neighbours?
        if (ivtx2-ivtx1.eq.1) then
          ! for neighbours life is easy
          ivtx1rm = ivtx1
          ivtx2lm = ivtx2
        else
          ! if no:
          ! find leftmost position where the right vertex may commute to
          ivtx2lm = ivtx2
          do ivtx = ivtx2-1, ivtx1+1, -1
            if (topomap(ivtx,ivtx2).gt.0) exit
            ivtx2lm = ivtx
          end do
          
          ! find rightmost positition where left vertex may commute to
          ivtx1rm = ivtx1
          do ivtx = ivtx1+1, ivtx2lm-1
            if (topomap(ivtx1,ivtx).gt.0) exit
            ivtx1rm = ivtx
c            nmvleft = nmvleft+1
c            imvleft(nmvleft) = ivtx
          end do
        end if

c dbg
c        print *,'ivtx1rm, ivtx2lm: ',ivtx1rm, ivtx2lm
c dbg
        ! can we be neighbours?
        if (ivtx2lm-ivtx1rm.eq.1) then
          idx_merge_vtx1vtx2 = ivtx1rm
        else
          ! else we have to check whether the remaining vertices
          ! between the two nodes will be contracted with 
          ! vtx1 and vtx2 in the same step, i.e. if they belong to
          ! the same super-vertex node of the result:

          ! before, remove any remaining vertex between ivtx1rm and
          ! ivtx2lm which commutes to the left with all vertices up 
          ! to ivtx1rm
          do ivtx = ivtx1rm+1, ivtx2lm-1
            commute = .true.
            do jvtx = ivtx-1, ivtx1rm+1
              commute = commute.and.topomap(ivtx,jvtx).eq.0
            end do
            ! it must also commute with vtx1:
            commute = commute.and.topomap(ivtx,ivtx1).eq.0
            if (commute) then
              nmvleft = nmvleft+1
              imvleft(nmvleft) = ivtx
            end if
          end do

          ! super vertex to which ivtx1 and ivtx2 contribute
          ! (if decided yet, 0 if undecided)
          merge = .true.
          isupres1 = max(svmap(ivtx1),svmap(ivtx2))
          do ivtx = ivtx1rm+1, ivtx2lm-1
c dbg
c          print *,'ivtx = ',ivtx
c          print *,'isupres = ',isupres1
c          print *,'svertex(ivtx).eq.ivtxsup1:',
c     &         svertex(ivtx).eq.ivtxsup1
c          print *,'svertex(ivtx).eq.ivtxsup2:',
c     &         svertex(ivtx) .eq.ivtxsup2
c          print *,'svmap(ivtx).eq.isupres1:',
c     &         svmap(ivtx).eq.isupres1
c dbg
            merge = merge.and.(  ! either it commutes (see above):
     &           imltlist(ivtx,imvleft,nmvleft,1).gt.0 .or.
     &           ((svertex(ivtx).eq.ivtxsup1.or. ! or it is contracted
     &             svertex(ivtx).eq.ivtxsup2).and. 
     &              (svmap(ivtx).eq.isupres1.or. ! with same result vertex
     &               isupres1.eq.0)) )         ! (or still fully contracted)
            if (isupres1.eq.0) isupres1 = svmap(ivtx) ! set result vertex if
                                                    ! still undecided
          end do
          if (merge) then
            idx_merge_vtx1vtx2 = ivtx1rm
          else
            idx_merge_vtx1vtx2 = 0
          end if
        end if
      else  
        idx_merge_vtx1vtx2 = 0
      end if

c dbg
c      print *,'final: ',merge
c dbg
      
      return
      end
