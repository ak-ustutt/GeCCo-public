
*----------------------------------------------------------------------*
      function iocc_xdn(ixdn,iocc)
*----------------------------------------------------------------------*
*     return (dependent on ixdn) 1: e(x)citation, 2: (d)e-excitation, or
*     3: (n)on-excitation part of occupation
*     include interface file ifc_ioccfunc.inc in calling routines!
*     actually, treatment of active spaces/non-excitation is a bit
*     unclear, still .... should be rewritten anyway ....
*     Modified to cope with external, virtual  orbitals. (GWR Apr.07)
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
      if(ixdn.eq.1)then
        iocc_xdn(ihole,2)=iocc(ihole,2)
        iocc_xdn(ipart,1)=iocc(ipart,1)
        iocc_xdn(iextr,1)=iocc(iextr,1)
      elseif(ixdn.eq.2)then
        iocc_xdn(ihole,1)=iocc(ihole,1)
        iocc_xdn(ipart,2)=iocc(ipart,2)
        iocc_xdn(iextr,2)=iocc(iextr,2)
      elseif(ixdn.eq.3)then
        iocc_xdn(ivale,1:2)=iocc(ivale,1:2)
      else
        call quit(1,'iocc_xdn','undefined part of matrix')
      endif
  
c      do ica = 1, 2
c        if (ixdn.eq.1.and.ica.eq.1) ihpv = ipart
c        if (ixdn.eq.1.and.ica.eq.2) ihpv = ihole
c        if (ixdn.eq.2.and.ica.eq.1) ihpv = ihole
c        if (ixdn.eq.2.and.ica.eq.2) ihpv = ipart
c        if (ixdn.eq.3) ihpv = 3
c        iocc_xdn(ihpv,ica) = iocc(ihpv,ica)
c      end do

      return
      end
