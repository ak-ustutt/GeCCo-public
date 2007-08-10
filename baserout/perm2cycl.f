*----------------------------------------------------------------------*
      subroutine perm2cycl(iperm,icycl,npgr)
*----------------------------------------------------------------------*
*     get cycle structure of permutation iperm, described by its
*     cycle structure (in icycl(1:npgr,1)) and its cycles (in 
*     icycl(1:npgr,2))
*----------------------------------------------------------------------*
      
      implicit none

      integer, intent(in) ::
     &     npgr, iperm(npgr)

      integer, intent(out) ::
     &     icycl(npgr,2)

      integer ::
     &     idx, jdx, istart, icnt, ipos, ncyc, imax, imin, idxmax
      integer ::
     &     icycraw(npgr), ib2cyc(npgr), lencyc(npgr)


      ! int cycle structure 
      icycl(1:npgr,1) = 0

      ! start with position of "1" in iperm
      do idx = 1, npgr
        if (iperm(idx).eq.1) then
          istart = idx
          exit
        end if
      end do

      idx = 1
      ncyc = 1
      ib2cyc(1:npgr) = 0 ! init "belongs to cycle" array

      ! loop through permutation
      loop1: do
        ! remember start value and follow cycle
        ipos = istart
        icnt = 1
        loop2: do
          icycraw(ipos) = icnt ! remember way through iperm
          ib2cyc(ipos) = ncyc  ! and to which cycle it belonged
          ipos = iperm(ipos)
          idx = idx+1
          if (ipos.eq.istart) exit loop2
          icnt = icnt+1
        end do loop2

        ! update cycle structure descriptor
        icycl(icnt,1) = icycl(icnt,1)+1 
        ! remember length of this cycle
        lencyc(ncyc) = icnt

        ! no more cycles?
        if (idx.gt.npgr) exit loop1

        ! increment cycle counter
        ncyc = ncyc+1

        ! get lowest next start value (next unassigned element of iperm)
        imin = npgr+1
        loop3: do jdx = 1, npgr
          if (ib2cyc(jdx).eq.0.and.iperm(jdx).lt.imin
     &       .and.iperm(jdx).gt.0) then             ! ib2cyc knows it ...
            imin = iperm(jdx)
            istart = jdx
          end if
        end do loop3
        if (imin.gt.npgr) exit loop1

      end do loop1
c dbg
      print *,'iperm  : ',iperm(1:npgr)
      print *,'icycraw: ',icycraw(1:npgr)
      print *,'ib2cyc : ',ib2cyc(1:npgr)
      print *,'lencyc : ',lencyc(1:ncyc)
c dbg
      
      ! resort to standard sequence (descending cycle order)
      jdx = 0
      do
        ! find maximum cycle length in lencyc
        imax = 0        
        do idx = 1, ncyc
          if (lencyc(idx).gt.imax) then
            imax = lencyc(idx)
            idxmax = idx
          end if
        end do
        if (imax.eq.0) exit

        lencyc(idxmax) = -imax
        
        do idx = 1, npgr
          if (ib2cyc(idx).eq.idxmax) then
            icycl(jdx+icycraw(idx),2) = iperm(idx)
          end if
        end do
        jdx = jdx+imax
        
      end do
c dbg
      print *,'->',icycl(1:npgr,1)
      print *,'  ',icycl(1:npgr,2)
c dbg

      end
