*----------------------------------------------------------------------*
      integer function sign_shift(occ_k,ivtx1,ivtx2,
     &     occ_i0,occ_k_tomove,occ_k_moved,nencl_pass,njoined)
*----------------------------------------------------------------------*
*     sign for shifting string K characterized by occupation occ_k 
*     from vertex ivtx1 to vertex ivtx2 of super-vertex I, characterized
*     by occupation occ_i(,,njoined) = occ_i0+occ_k_tomove,occ_k_moved
*     occ_k_tomove: simultaneous shifts not yet considered
*     occ_k_moved:  simultaneous shifts already considered
*     we always consider right-resolved restrictions, i.e.
*      I1;...;I2 -> I1',K;...;I2 -> I1';...;I2,K -> I1';...;I2'   or
*      I1;...;I2 -> I1;...;I2",K -> I1,K;...;I2" -> I1";...;I2"   
*
*     nencl_pass: number of CA on passive vertices that are passed by
*                 the presently considered operator
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     njoined, ivtx1, ivtx2, nencl_pass,
     &     occ_k(ngastp,2), occ_i0(ngastp,2,njoined),
     &     occ_k_tomove(ngastp,2,njoined), occ_k_moved(ngastp,2,njoined)

      integer ::
     &     ivtx_l, ivtx_h, ivtx, nshift_c, nshift_a,
     &     nenclosed, n_1c, n_1a, n_2c, n_2a,
     &     n_lc, n_la, n_hc, n_ha, nnym_la, nnym_lc, sign0, nnym_hc

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_func,'sign_shift')
        write(lulog,*) 'from, to ',ivtx1,ivtx2
        call wrt_occ(lulog,occ_k)
      end if
      
      ivtx_l = min(ivtx1,ivtx2)
      ivtx_h = max(ivtx1,ivtx2)
      
      ! elements in shift string
      nshift_c = sum(occ_k(1:ngastp,1))
      nshift_a = sum(occ_k(1:ngastp,2))

      ! elements in start vertex after removing shifted elements
      ! = the unchange part of the occupation +
      ! those elements which already moved here in
      ! a previous step (occ_k_moved)
      n_1c = sum(occ_i0(1:ngastp,1,ivtx1))
     &     + sum(occ_k_moved(1:ngastp,1,ivtx1))
      n_1a = sum(occ_i0(1:ngastp,2,ivtx1))
     &     + sum(occ_k_moved(1:ngastp,2,ivtx1))

      ! elements in target vertex
      ! see above
      n_2c = sum(occ_i0(1:ngastp,1,ivtx2))
     &     + sum(occ_k_moved(1:ngastp,1,ivtx2))
      n_2a = sum(occ_i0(1:ngastp,2,ivtx2))
     &     + sum(occ_k_moved(1:ngastp,2,ivtx2))

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
      ! includes those elements which either
      ! will move in succeeding step (occ_k_tomove)
      ! or those which already moved (occ_k_moved)
      nenclosed = nencl_pass
      do ivtx = ivtx_l+1, ivtx_h-1
        nenclosed = nenclosed +
     &       sum(occ_i0(1:ngastp,1:2,ivtx))+
     &       sum(occ_k_tomove(1:ngastp,1:2,ivtx))+
     &       sum(occ_k_moved(1:ngastp,1:2,ivtx))
      end do
c dbg
c      print *,'effective nenclosed: ',nenclosed
c dbg

      ! elements not yet moved at lower index:
      nnym_lc = sum(occ_k_tomove(1:ngastp,1,ivtx_l))
      nnym_la = sum(occ_k_tomove(1:ngastp,2,ivtx_l))
      ! elements not yet moved at higher index:
      nnym_hc = sum(occ_k_tomove(1:ngastp,1,ivtx_h))

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
      sign0 = mod(nshift_c*(n_la + nnym_lc + nnym_la + nenclosed + n_hc)
cmh     &         +  nshift_a*(n_hc + n_ha +    nnym_la + nenclosed)  , 2)
     &       + nshift_a*(n_hc + n_ha + nnym_hc + nnym_la + nenclosed),2)

      if (ntest.ge.100) then
        write(lulog,*) 'nshift_c, nshift_a: ',nshift_c,nshift_a
        write(lulog,*) 'n_la, nnym_lc, nnym_la, nenclosed, n_hc: ',
     &       n_la, nnym_lc, nnym_la, nenclosed, n_hc
        write(lulog,*) 'summed: ',n_la+nnym_lc+nnym_la+nenclosed+n_hc
        write(lulog,*) 'n_hc, n_ha, nnym_la, nenclosed:          ',
     &       n_hc, n_ha, nnym_la, nenclosed
        write(lulog,*) 'summed: ',n_hc+n_ha+nnym_la+nenclosed
      end if

      sign_shift = 1

      if (sign0.eq.1) sign_shift = -1

      return
      end
