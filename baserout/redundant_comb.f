*----------------------------------------------------------------------*
      subroutine redundant_comb(comb,newcomb,redun,items,imax,
     &                                nreg)
*----------------------------------------------------------------------*
*     returns non-redundant combination on newcomb
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     imax(nreg), items, comb(items), redun(nreg*maxval(imax)),
     &     nreg

      integer, intent(out) ::
     &     newcomb(items)

      integer ::
     &     ii, counter(nreg*maxval(imax)), 
     &     temp(nreg*maxval(imax)), ireg, jj

      ! determine pattern of combination and register
      counter = 0
      ireg = 0
      do ii = 1,items
        counter(comb(ii)) = counter(comb(ii)) + 1
        if (ireg.eq.0) then
          ireg = (comb(ii)-1) / maxval(imax) + 1
        else if (ireg.ne.(comb(ii)-1)/maxval(imax)+1) then
          call quit(1,'redundant_comb','forbidden combination')
        end if
      end do

      temp = redun

      ! redefine redun pattern to allow for determination of newcomb
      do ii = 1,nreg
        do jj = (ii-1)*maxval(imax)+1,(ii-1)*maxval(imax)+imax(ii)
          temp(jj) = temp(temp(jj))
        end do
      end do

      ! determine non-redundant combination
      do ii = 1,items
        newcomb(ii) = temp(comb(ii))
      end do
      ! sort in ascending order
      call isort(newcomb,items,1)

      return
      end
