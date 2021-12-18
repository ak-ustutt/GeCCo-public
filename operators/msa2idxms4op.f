*----------------------------------------------------------------------*
      pure integer function msa2idxms4op(msa,mst,occ_a,occ_c)
*----------------------------------------------------------------------*
*     convert Ms(A)-value (actually 2 times Ms) to index value of
*     operator block
*     e.g.  2, 0, -2  to 1, 2, 3 ....  in most cases
*     if (mst>0 or occ(A)-occ(C)<0) a number of ms(A) values is not
*     allowed (as the neccessary ms(C) cannot be realized) and thus 
*     skipped; we account for this offset by calculating (noting that
*     ms(A)max = occ(A); ms(C)max=occ(C))
*       ms(C)needed = ms(A)max + ms(total)
*       offset = (ms(C)needed - ms(C)max)/2
*     only positive offsets need to be considered as we enumerate
*     the blocks with ms(A)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     msa, mst, occ_a, occ_c

      integer ::
     &     offset

      msa2idxms4op = (occ_a - msa)/2 + 1
      offset = occ_a-occ_c+mst
      if ( offset .gt. 0 )
     &     msa2idxms4op = msa2idxms4op - offset/2

      return
      end
