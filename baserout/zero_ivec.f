      pure logical function zero_ivec(ivec,nel)

      implicit none

      integer, intent(in) ::
     &     nel, ivec(nel)
      
      integer ::
     &     iel

      zero_ivec = .true.
      do iel = 1, nel
        zero_ivec = zero_ivec.and.ivec(iel).eq.0
      end do

      end
