*----------------------------------------------------------------------*
      integer function sign_contr(cnt,i0,j0,n_enclosed,reversed)
*----------------------------------------------------------------------*
*     get sign of contraction:
*
*      {KC I0C I0A KA}{enclosed}{KA+ J0C J0A KC+}
*       |          |             |           |
*       |          +-------------+           |
*       +------------------------------------+
*
*     the sign for {IC IA} -> {KC I0C I0A KA} is calculated in
*     sign_hpvx (use mode 1)
*
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     cnt(ngastp,2), i0(ngastp,2), j0(ngastp,2),
     &     n_enclosed
      logical, intent(in) ::
     &     reversed

      integer ::
     &     ncntc, ncnta, ni0, nj0, ntrapo

      ncntc = sum(cnt(1:ngastp,1))
      ncnta = sum(cnt(1:ngastp,2))

      ni0   = sum(i0(1:ngastp,1:2))
      nj0   = sum(j0(1:ngastp,1:2))

      if (.not.reversed) then
        ! n(CNT(A)) * n(enclosed)
        ntrapo = mod(ncnta*n_enclosed,2)

        ! n(CNT(C)) * { n(enclosed) + n(I0) + n(J0)}
        ntrapo = mod(ntrapo + ncntc*(n_enclosed+ni0+nj0),2)
      else
        ! n(CNT(C)) * n(enclosed)
        ntrapo = mod(ncntc*n_enclosed,2)

        ! n(CNT(A)) * { n(enclosed) + n(I0) + n(J0)}
        ntrapo = mod(ntrapo + ncnta*(n_enclosed+ni0+nj0),2)
      end if

      sign_contr = 1
      if (ntrapo.eq.1) sign_contr = -1

      return
      end
