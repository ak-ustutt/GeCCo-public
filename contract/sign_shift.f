*----------------------------------------------------------------------*
      integer function sign_shift(occ_k,ivtx1,ivtx2,occ_i,njoined)
*----------------------------------------------------------------------*
*     sign for shifting string K characterized by occupation occ_k 
*     from vertex ivtx1 to vertex ivtx2 of super-vertex I, characterized
*     by occupation occ_i(,,njoined)
*     we always consider right-resolved restrictions, i.e.
*      I1;...;I2 -> I1',K;...;I2 -> I1';...;I2,K -> I1';...;I2'   or
*      I1;...;I2 -> I1;...;I2",K -> I1,K;...;I2" -> I1";...;I2"   
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     njoined, ivtx1, ivtx2,
     &     occ_k(ngastp,2), occ_i(ngastp,2,njoined)

      integer ::
     &     ivtx_l, ivtx_h, ivtx, nshift_c, nshift_a,
     &     nenclosed, n_1c, n_1a, n_2c, n_2a,
     &     n_lc, n_la, n_hc, n_ha, sign0

      ivtx_l = min(ivtx1,ivtx2)
      ivtx_h = max(ivtx1,ivtx2)
      
      ! elements in shift string
      nshift_c = sum(occ_k(1:ngastp,1))
      nshift_a = sum(occ_k(1:ngastp,2))

      ! elements in start vertex after removing shifted elements
      n_1c = sum(occ_i(1:ngastp,1,ivtx1)) - nshift_c
      n_1a = sum(occ_i(1:ngastp,2,ivtx1)) - nshift_a

      ! elements in target vertex
      n_2c = sum(occ_i(1:ngastp,1,ivtx2))
      n_2a = sum(occ_i(1:ngastp,2,ivtx2))

      ! assign 1, 2 to l(ow) and h(igh)
      ! the sign formula is valid for both left and right shifts then
      if (ivtx1.lt.ivtx2) then
        n_lc = n_1c
        n_la = n_1a
        n_hc = n_2c
        n_ha = n_2a
      else
        n_lc = n_2c
        n_la = n_2a
        n_hc = n_1c
        n_ha = n_1a
      end if

      ! count elements in enclosed strings:
      nenclosed = 0
      do ivtx = ivtx_l+1, ivtx_h-1
        nenclosed = nenclosed + sum(occ_i(1:ngastp,1:2,ivtx))
      end do

      !--------------------------------------------------------!
      ! count the elements passed by KC and KA to verify the   !
      ! sign formula. reverting the sequence of shifting, i.e. !
      ! either KA first or KA first, does not change the sign. !
      !                                                        !
      ! right: {IlC' KC IlA' KA} {enclosed} {IhC IhA }         !
      !              |       +----------------------|          !
      !              +--------------------------|              !
      !                                                        !
      ! left:  {IlC IlA } {enclosed} {IhC' KC IhA' KA}         !
      !            |-----------------------+       |           !
      !                |---------------------------+           !
      !--------------------------------------------------------!
      sign0 = mod(nshift_c*(n_la + nenclosed + n_hc) +
     &            nshift_a*(n_hc + n_ha + nenclosed)  , 2)

      sign_shift = 1

      if (sign0.eq.1) sign_shift = -1

      return
      end
