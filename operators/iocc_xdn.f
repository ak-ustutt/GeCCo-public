
*----------------------------------------------------------------------*
      function iocc_xdn(ixdn,iocc)
*----------------------------------------------------------------------*
*     return (dependent on ixdn) 1: e(x)citation, 2: (d)e-excitation, or
*     3: (n)on-excitation part of occupation
*     include interface file ifc_ioccfunc.inc in calling routines!
*     actually, treatment of active spaces/non-excitation is a bit
*     unclear, still .... should be rewritten anyway ....
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer ::
     &     iocc_xdn(ngastp,2)

      integer, intent(in) ::
     &     ixdn,
     &     iocc(ngastp,2)

      integer ::
     &     ica, ihpv

      iocc_xdn(1:ngastp,1:2) = 0
      do ica = 1, 2
        if (ixdn.eq.1.and.ica.eq.1) ihpv = ipart
        if (ixdn.eq.1.and.ica.eq.2) ihpv = ihole
        if (ixdn.eq.2.and.ica.eq.1) ihpv = ihole
        if (ixdn.eq.2.and.ica.eq.2) ihpv = ipart
        if (ixdn.eq.3) ihpv = 3
        iocc_xdn(ihpv,ica) = iocc(ihpv,ica)
      end do

      return
      end
