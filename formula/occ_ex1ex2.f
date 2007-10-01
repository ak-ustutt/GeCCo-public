*----------------------------------------------------------------------*
      subroutine occ_ex1ex2(iocc_ex1ex2,inum_ori,
     &                      iocc_ex1,iocc_ex2,
     &                      njoined1,njoined2,
     &                      ivtxsuper1,ivtxsuper2,
     &                      svertex,nvtx)
*----------------------------------------------------------------------*
*     assemble the two external occupations in the correct ordering
*     (according to the super-vertex sequence on svertex(nvtx))
*     report wether the entry on iocc_ex1ex2 was from either ex1
*     or ex2 (and the old number) on inum_ori
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, intent(in) ::
     &     njoined1, njoined2, nvtx,
     &     iocc_ex1(ngastp,2,njoined1),
     &     iocc_ex2(ngastp,2,njoined2),
     &     ivtxsuper1,ivtxsuper2,svertex(nvtx)
      integer, intent(out) ::
     &     iocc_ex1ex2(ngastp,2,njoined1+njoined2),
     &     inum_ori(2,njoined1+njoined2)

      integer ::
     &     idx1, idx2, idx12, ivtx

      idx1  = 0
      idx2  = 0
      idx12 = 0
      do ivtx = 1, nvtx
        if (svertex(ivtx).eq.ivtxsuper1) then
          idx1  = idx1 +1
          if (idx1.gt.njoined1)
     &         call quit(1,'occ_ex1ex2','njoined1 too small?')
          idx12 = idx12+1
          iocc_ex1ex2(1:ngastp,1:2,idx12) =
     &         iocc_ex1(1:ngastp,1:2,idx1)
          inum_ori(1,idx12) = 1
          inum_ori(2,idx12) = idx1
        else if (svertex(ivtx).eq.ivtxsuper2) then
          idx2  = idx2 +1
          if (idx2.gt.njoined2)
     &         call quit(1,'occ_ex1ex2','njoined2 too small?')
          idx12 = idx12+1
          iocc_ex1ex2(1:ngastp,1:2,idx12) =
     &         iocc_ex2(1:ngastp,1:2,idx2)
          inum_ori(1,idx12) = 2
          inum_ori(2,idx12) = idx2
        end if
      end do
c dbg
c      print *,'OCCEX1EX2: idx12 = ',idx12
c      print *,'svertex ',svertex(1:nvtx)
c      print *,'ivtxsuper1,ivtxsuper2:',ivtxsuper1,ivtxsuper2
c      print *,'inum_ori 1>',inum_ori(1,1:idx12)
c      print *,'inum_ori 2>',inum_ori(2,1:idx12)
c      call wrt_occ_n(6,iocc_ex1ex2,idx12)
c      if (idx12.ne.njoined1+njoined2) stop 'aha'
c dbg

      end
