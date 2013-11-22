      subroutine occ2descr(descr,maxlen,occ,nj)
      !
      ! turn occupation into descriptor string
      ! if descr is too small, the output will be incomplete,
      ! but we do not quit or anything (assuming this is
      ! uncritical information output)
      !
      ! andreas, nov 2013
      !
      implicit none
      include 'opdim.h'

      integer, intent(in) :: 
     &    nj, occ(ngastp,2,nj), maxlen
     
      character(len=maxlen), intent(out) ::
     &    descr

      character, parameter ::
     &    symbol(1:4) = (/'H','P','V','X'/)

      integer ::
     &    ipos, ij, ica, itp, nn, kk

      ! init with spaces
      descr(1:maxlen) = ' '

      ! position counter
      ipos = 1

      ! loop over occ
      ij_loop: do ij = 1, nj
        do ica = 1, 2
          do itp = 1, ngastp
            nn = occ(itp,ica,ij)
            if (nn.eq.0) cycle
            do kk = 1, nn
              descr(ipos:ipos) = symbol(itp)
              ipos = ipos+1
              if (ipos.gt.maxlen) exit ij_loop
            end do
          end do
          if (ica.eq.1) then
            descr(ipos:ipos) = ','
            ipos = ipos+1
            if (ipos.gt.maxlen) exit ij_loop
          end if 
        end do
        if (ij.lt.nj) then
          descr(ipos:ipos) = ';'
          ipos = ipos+1
          if (ipos.gt.maxlen) exit ij_loop
        end if
      end do ij_loop

      end subroutine

