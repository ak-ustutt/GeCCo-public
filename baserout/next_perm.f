*----------------------------------------------------------------------*
*     some tests to get an idea how to obtain permutations
*----------------------------------------------------------------------*
      logical function next_perm(iperm,nel)
*----------------------------------------------------------------------*
*     get next permutation
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     iperm(nel)
      integer, parameter ::
     &     way = 2
      logical, external ::
     &     next_perm_rec, next_perm_iter

      if (way.eq.1) next_perm = next_perm_rec(iperm,nel)
      if (way.eq.2) next_perm = next_perm_iter(iperm,nel)
      
      return
      end 

*----------------------------------------------------------------------*
      logical recursive function next_perm_rec(iperm,nel) result(lres)
*----------------------------------------------------------------------*
*     recursive kernel function
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     iperm(nel)

      logical ::
     &     succ
      integer ::
     &     ihlp, idx, jdx, imin, idxmin

      if (nel.eq.2) then
        if (iperm(1).lt.iperm(2)) then
          ihlp = iperm(1)
          iperm(1) = iperm(2)
          iperm(2) = ihlp
          succ = .true.
        else
          succ = .false.
        end if
      else
        ! try to find next permutation in substring
        succ = next_perm_rec(iperm(2),nel-1) 
        ! find next larger element in iperm and set as first
        if (.not.succ) then
          imin = huge(imin)
          idxmin = -1
          do idx = 1, nel  ! ??  2,nel ??
            if (iperm(idx).gt.iperm(1).and.iperm(idx).lt.imin) then
              imin = iperm(idx)
              idxmin = idx
            end if
          end do
          succ = idxmin.gt.0
          if (succ) then
            ! swap iperm(1) with next larger element:
            ihlp = iperm(1)
            iperm(1) = iperm(idxmin)
            iperm(idxmin) = ihlp
            ! sort remaining list (2...nel) 
            ! by insertion algorithm
            do idx = 3, nel
              ihlp = iperm(idx)
              jdx = idx-1
              do while (jdx.ge.2.and.iperm(jdx).gt.ihlp)
                iperm(jdx+1) = iperm(jdx)
                jdx = jdx-1
              end do
              iperm(jdx+1) = ihlp
            end do
            succ = .true.

          end if
        end if
      end if

      !next_perm_rec = succ
      lres = succ

      return
      end

*----------------------------------------------------------------------*
      logical function next_perm_iter(iperm,nel)
*----------------------------------------------------------------------*
*     iterative kernel function
*----------------------------------------------------------------------*      
      implicit none

      integer, intent(in) ::
     &     nel
      integer, intent(inout) ::
     &     iperm(nel)

      logical ::
     &     succ
      integer ::
     &     ihlp, idx, jdx, kdx, imin, idxmin

      succ = .false.
      main_loop : do idx = nel-1, 1, -1

        if (nel-idx.eq.1) then
          if (iperm(idx).lt.iperm(idx+1)) then
            ihlp = iperm(idx)
            iperm(idx) = iperm(idx+1)
            iperm(idx+1) = ihlp
            succ = .true.
            exit main_loop
          end if
        else
          imin = huge(imin)
          idxmin = -1
          do jdx = idx+1, nel
            if (iperm(jdx).gt.iperm(idx).and.iperm(jdx).lt.imin) then
              imin = iperm(jdx)
              idxmin = jdx
            end if
          end do
          succ = .false.
          if (idxmin.gt.0) then
            ihlp = iperm(idx)
            iperm(idx) = iperm(idxmin)
            iperm(idxmin) = ihlp
            ! sort remaining list
            ! by insertion algorithm
            do jdx = idx+2, nel
              ihlp = iperm(jdx)
              kdx = jdx-1
              do while (kdx.gt.idx.and.iperm(kdx).gt.ihlp)
                iperm(kdx+1) = iperm(kdx)
                kdx = kdx-1
              end do
              iperm(kdx+1) = ihlp
            end do
            succ = .true.

            exit main_loop
          end if
          
        end if

      end do main_loop
 
      next_perm_iter = succ

      return
      end
