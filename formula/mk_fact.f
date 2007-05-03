*----------------------------------------------------------------------*
      subroutine mk_fact(ifact,nfact,icntsq,iconn,narc,nvtx)
*----------------------------------------------------------------------*
*     a sequence of contractions is given on icntsq which is applied
*     to the graph described by the connectivity array iconn(3,narc)
*     for description of iconn see: form_fact
*     obtain the sequence of binary contractions on ifact(nfact)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'

      integer, intent(in) ::
     &     narc, nvtx, icntsq(narc), iconn(3,narc)
      integer, intent(out) ::
     &     nfact, ifact(4,*)
      
      integer ::
     &     idx, idxl, jdx, inum, narc_cur, narc_next, ivtxnew,
     &     icnt, ibase, ivtx1, ivtx2, nexp1, nexp2
      integer ::
     &     icoscr(3,narc), icoscr2(3,narc),
     &     iarceq(2,narc),
     &     idx_exp(nvtx)

      integer, external ::
     &     maxvtx, idxlist, int_pack, int_expand
      
      if (ntest.ge.100) then
        write(luout,*) '================='
        write(luout,*) ' mk_fact at work'
        write(luout,*) '================='
        write(luout,*) ' icntsq : ',icntsq(1:narc)
        write(luout,*) ' graph described by: '
        do idx = 1, narc
          write(luout,*)
     &         '  ',iconn(1,idx),':',iconn(2,idx),'-',iconn(3,idx)
        end do
      end if
      
      do idx = 1, narc
        iarceq(1,idx) = iconn(1,idx) 
        iarceq(2,idx) = iconn(1,idx) 
      end do
      idx = 1
      nfact = 0
      icoscr(1:3,1:narc) = iconn(1:3,1:narc)
      narc_cur = narc
      do while (narc_cur.gt.0)
        inum = icntsq(idx)  ! arc number form sequence list
        idxl = idxlist(inum,iarceq(1,1),narc,2)
        inum = iarceq(2,idxl) ! convert to equivalent arc
        ! find arc on current list
        idxl = idxlist(inum,icoscr(1,1),narc_cur,3)
        ! get equivalent arcs and join information into one number
        icnt = 0
        ibase = 1
        do jdx = 1, narc
          if (iarceq(2,jdx).eq.inum) then
            icnt = icnt+ibase*iarceq(1,jdx)
            ibase = ibase*(narc+1)
          end if
        end do
        if (idxl.gt.0) then
          ! code names of vertices into new vertex name:
          ivtx1 = icoscr(2,idxl)
          ivtx2 = icoscr(3,idxl)
          ! a) expand vertex names
          nexp1 = int_expand(ivtx1,nvtx+1,idx_exp)
          nexp2 = int_expand(ivtx2,nvtx+1,idx_exp(nexp1+1))
c dbg
          if (nexp1+nexp2.gt.nvtx) stop 'error trap -- nexp.gt.nvtx'
c dbg
          ! b) sort new vertex list
          call isort(idx_exp,nexp1+nexp2,+1)
          ! c) generate new name
          ivtxnew = int_pack(idx_exp,nexp1+nexp2,nvtx+1)
          call reduce_graph(icoscr2,narc_next,iarceq,
     &         inum,ivtxnew,icoscr,narc_cur,narc)
          nfact = nfact+1
          ifact(1,nfact) = icoscr(2,idxl)
          ifact(2,nfact) = icoscr(3,idxl)
          ifact(3,nfact) = ivtxnew
          ifact(4,nfact) = icnt
        end if
        
        idx = idx+1
        icoscr(1:3,1:narc_next) = icoscr2(1:3,1:narc_next)
        narc_cur = narc_next
      end do

      if (ntest.ge.100) then
        write(luout,*) 'resulting binary operation list:'
        do idx = 1, nfact
          icnt = ifact(4,idx)
          jdx = 1
          do
            iarceq(1,jdx) = mod(icnt,narc+1)
            icnt = icnt/(narc+1)
            if (icnt.eq.0) exit
            jdx = jdx+1
          end do
          write(luout,'(2x,i4,a,i4,a,i4,a,8i3)')
     &         ifact(1,idx),' * ',ifact(2,idx),
     &         ' -> ',ifact(3,idx),'  C:',iarceq(1,1:jdx)
        end do
      end if

      return
      end
