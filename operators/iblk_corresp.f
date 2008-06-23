*----------------------------------------------------------------------*
      integer function iblk_corresp(idx_next_poss,iblk1,op1,op2,dag2)
*----------------------------------------------------------------------*
*     find block of operator op2 that corresponds to block iblk1 of op1
*     if op2 has a different particle creation/annihilation rank, find
*     the occupations resulting from deleting C or A operators from
*     the occupation of op1, iblk1. Several possibilities may result
*     in cases where we have more than one HPVX occupied
*     
*     result: block of operator 2
*             if >=zero:     nothing found
*     idx_next_poss: if more than one possibility exists
*         on entry: choose the poss. to be considered in this call
*                   (1 on first call)
*         on exit:  -1 if no further poss. exists
*                   >1 (= number of next poss.) if several poss. ex.
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      type(operator), intent(in) ::
     &     op1, op2
      integer, intent(inout) ::      
     &     idx_next_poss
      integer, intent(in) ::      
     &     iblk1
      logical, intent(in) ::
     &     dag2

      integer ::
     &     nj1, nj2, idxblk1, idxblk2,
     &     pcr1, pcr2, idx, nposs, ij, igastp, ica

      integer, pointer ::
     &     occ_blk1(:,:,:), occ_blk2(:,:,:),
     &     occ_scr(:,:,:)
      
      integer, external ::
     &     rank_occ, iblk_occ, imltlist

      iblk_corresp = 0
      nj1 = op1%njoined
      nj2 = op2%njoined
      if (nj1.ne.nj2) return

      occ_blk1 => op1%ihpvca_occ
      occ_blk2 => op1%ihpvca_occ

      ! get number of particles created/annihilated
      pcr1 = rank_occ('C-A',occ_blk1,nj1)
      pcr2 = rank_occ('C-A',occ_blk1,nj2)
      
      idxblk1 = (iblk1-1)*nj1+1

      ! easy game
      if (pcr1.eq.pcr2) then
c dbg
        print *,'in easy section: ',idxblk1,dag2,trim(op2%name)
c dbg
        iblk_corresp = iblk_occ(occ_blk1(1,1,idxblk1),dag2,op2)
        idx_next_poss = -1
      else if (abs(pcr1-pcr2).eq.1) then
        if (pcr1-pcr2.eq.+1) ica = 1
        if (pcr1-pcr2.eq.-1) ica = 2

        nposs = ngastp
     &       - imltlist(0,occ_blk1(1:ngastp,ica,idxblk1:idxblk1-1+nj1),
     &                  ngastp*nj1,1)

        if (nposs.lt.idx_next_poss)
     &       call quit(1,'iblk_corresp','idx_next_poss out of range')

        allocate(occ_scr(ngastp,2,nj1))

        idx = 0
        this_loop: do ij = 1, nj1
          do igastp = 1, ngastp
            if (occ_blk1(igastp,ica,idxblk1-1+ij).eq.0) cycle
            idx = idx+1
            if (idx.ne.idx_next_poss) cycle
            occ_scr(1:ngastp,1:2,1:nj1) =
     &           occ_blk1(1:ngastp,1:2,idxblk1:idxblk1-1+nj1)
            occ_scr(igastp,ica,ij) =
     &           occ_scr(igastp,ica,ij) - 1
            iblk_corresp = iblk_occ(occ_scr,dag2,op2)
            exit this_loop
          end do
        end do this_loop

        deallocate(occ_scr)

        if (idx_next_poss.lt.nposs) then
          idx_next_poss = idx_next_poss+1
        else
          idx_next_poss = -1
        end if

      end if

      return
      end
