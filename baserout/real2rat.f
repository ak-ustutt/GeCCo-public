*----------------------------------------------------------------------* 
      subroutine real2rat(num,den,rl)
*----------------------------------------------------------------------*
*     get best approximation of real number rl by rational number
*     rat = num/den
*
*     stupid brute force algorithm
*
*----------------------------------------------------------------------*
      implicit none

      integer, intent(out) ::
     &     num, den
      real(8), intent(in) ::
     &     rl

      integer, parameter ::
     &     maxtest = 40320  !   = 8!
      real(8), parameter ::
     &     thrsh = 1d-12

      integer ::
     &     numlw0, numup0, denlw0, denup0, denlw, denup, numlw, numup

      if (rl.gt.1d0) then
        numlw0 = int(rl)
        numup0 = numlw0+1
        den_loop: do den = 1, maxtest
          numlw = numlw0*den
          numup = numup0*den
          do num = numlw, numup
            if (abs(dble(num)/dble(den)-rl).lt.thrsh) exit den_loop
            if (abs(dble(num)/dble(den)-rl).gt.0d0) cycle
          end do
        end do den_loop
      else
        denlw0 = int(1d0/rl)
        denup0 = denlw0+1
        num_loop: do num = 1, maxtest
          denlw = max(1,denlw0*num)
          denup = denup0*num
          do den = denlw, denup
            if (abs(dble(num)/dble(den)-rl).lt.thrsh) exit num_loop
            if (abs(dble(num)/dble(den)-rl).lt.0d0) cycle
          end do
        end do num_loop
      end if

      if (.not.(abs(dble(num)/dble(den)-rl).lt.thrsh))
     &       call quit(1,'real2rat','too stupid for this')

      end

