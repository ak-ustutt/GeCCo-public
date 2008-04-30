      pure logical function zero_i8vec(ivec,nel,inc)

      implicit none

      integer, intent(in) ::
     &     nel, inc
      integer(8), intent(in) ::
     &     ivec(nel)
      
      integer ::
     &     iel, idx

      zero_i8vec = .true.
      idx = 1
      do iel = 1, nel
        zero_i8vec = zero_i8vec.and.ivec(idx).eq.0
        idx = idx+inc
      end do

      end
