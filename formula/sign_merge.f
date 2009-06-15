*----------------------------------------------------------------------*
      integer function sign_merge(i0,j0,n_enclosed,reversed)
*----------------------------------------------------------------------*
*     get sign of merge:
*
*      {I0C I0A }{enclosed}{J0C J0A} -> {I0C J0C I0A J0A}{}{}
*          ^   ^            |   |
*          |   +----------------+
*          +----------------+
*
*     or if reversed
*
*      { I0C I0A }{enclosed}{J0C J0A} -> {J0C I0C J0A I0A}{}{}
*       ^   ^                |   |
*       |   +--------------------+
*       +--------------------+
*
*
*     the sign for {I0C H0C I0A J0A} -> {LC LA} is calculated in
*     sign_hpvx (use mode 2)
*
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'

      logical, intent(in) ::
     &     reversed
      integer, intent(in) ::
     &     i0(ngastp,2), j0(ngastp,2),
     &     n_enclosed

      integer ::
     &     ni0a, nj0a, ni0c, nj0c, ntrapo

      ni0c  = sum(i0(1:ngastp,1))
      ni0a  = sum(i0(1:ngastp,2))
      nj0c  = sum(j0(1:ngastp,1))
      nj0a  = sum(j0(1:ngastp,2))

      if (.not.reversed) then
        ! n(J0(A)) * n(enclosed)
        ntrapo = mod(nj0a*n_enclosed,2)

        ! n(J0(C)) * { n(enclosed) + n(I0(A))}
        ntrapo = mod(ntrapo + nj0c*(n_enclosed+ni0a),2)
      else
        ! n(J0(A)) * { n(enclosed) + n(I0(A))}
        ntrapo = mod(ntrapo + nj0a*(n_enclosed+ni0a),2)

        ! n(J0(C)) * { n(enclosed) + n(I0(A)) + n(I0(C))}
        ntrapo = mod(ntrapo + nj0c*(n_enclosed+ni0a+ni0c),2)
        
      end if

      sign_merge = 1
      if (ntrapo.eq.1) sign_merge = -1

      return
      end
