*----------------------------------------------------------------------*
      logical function merge_vtx1vtx2(ivtx1,ivtx2,
     &     svertex,svmap,topomap,nvtx)
*----------------------------------------------------------------------*
*     check whether vertices may be merged into one after a contraction
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ivtx1, ivtx2, nvtx, svmap(nvtx), topomap(nvtx,nvtx),
     &     svertex(nvtx)

      logical ::
     &     merge
      integer ::
     &     isupvtx1, ivtx

      ! first of all: to which vertex of the result do they
      !  contribute open lines?
      merge = svmap(ivtx1).eq.0 .or.       ! either no contribution
     &        svmap(ivtx2).eq.0 .or.
     &        svmap(ivtx1).eq.svmap(ivtx2) ! or to the same vertex

      ! furthermore:
      ! we may merge vertices if neither of them connects to an
      ! operator vertex located between them
      ! either we have neighbouring vertices ...
      if (merge.and.ivtx2-ivtx1.ne.1) then
        ! ... or we must take a look at the current topology          
c dbg
c        print *,'testing ...'
c dbg
        merge = .true.
        isupvtx1 = max(svmap(ivtx1),svmap(ivtx2))
        do ivtx = ivtx1+1, ivtx2-1
c dbg
c          print *,'ivtx = ',ivtx
c          print *,'topomap(ivtx,ivtx2).eq.0:',topomap(ivtx,ivtx2).eq.0
c          print *,'svertex(ivtx1).eq.svertex(ivtx):',
c     &         svertex(ivtx1).eq.svertex(ivtx)
c          print *,'svertex(ivtx) .eq.svertex(ivtx2):',
c     &         svertex(ivtx) .eq.svertex(ivtx2)
c          print *,'svmap(ivtx).eq.isupvtx1:',
c     &         svmap(ivtx).eq.isupvtx1
c dbg
          merge = merge.and.
     &         (     topomap(ivtx,ivtx2).eq.0 ! either not connected
     &                      ! or vertex belongs to either supervertex ...
     &         .or.((svertex(ivtx1).eq.svertex(ivtx).or. 
     &               svertex(ivtx) .eq.svertex(ivtx2)).and. 
     &              (svmap(ivtx).eq.isupvtx1.or. ! with same result vertex
     &               isupvtx1.eq.0)) )         ! (or still fully contracted)
          if (isupvtx1.eq.0) isupvtx1 = svmap(ivtx) ! set result vertex if
                                                    ! still undecided
        end do
        
        ! if .false. we may test whether ivtx1 may be commuted to ivtx2:
        if (.not.merge) then
          merge = .true.
          do ivtx = ivtx1+1, ivtx2-1
            merge = merge.and.
     &           topomap(ivtx1,ivtx).eq.0 
          end do
        end if
            
      end if

c dbg
c      print *,'final: ',merge
c dbg
      merge_vtx1vtx2 = merge
      
      return
      end
