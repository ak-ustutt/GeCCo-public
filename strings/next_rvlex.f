*----------------------------------------------------------------------*
      logical function next_rvlex(nel,idistr,ielmin,ielmax)
*----------------------------------------------------------------------*
*     next element in reverse lexical order
*
*     allow elements to vary within interval [ielmin(ipos):ielmax(ipos)]
*     (may be different for each position)
*     allow up to nel identical elements in string
*
*     purpose is generation of subspace and gamma distributions
*
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     nel, ielmin(nel), ielmax(nel)

      integer, intent(inout) ::
     &     idistr(nel)

      logical ::
     &     lnext
      integer ::
     &     ipos

      if (nel.eq.0) call quit(1,'next_rvlex','nel=0 error')
      ipos = 0
      lnext = .true.
      ! loop over elements
      ipos_loop: do
        ipos = ipos+1
        ! element lt next element, or last element
        if (ipos.eq.nel.or.idistr(ipos).lt.idistr(ipos+1)) then    
          ! check bounds
          if (idistr(ipos).lt.ielmax(ipos)) then
            idistr(ipos) = idistr(ipos)+1
            exit ipos_loop
          else if (ipos.eq.nel) then
            lnext = .false.
            exit ipos_loop
          end if
        end if
      end do ipos_loop

      if (lnext.and.ipos.gt.1) then
        idistr(1:ipos-1) = ielmin(1:ipos-1)
      end if

      next_rvlex = lnext

      return

      end
