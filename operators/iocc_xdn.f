
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
     &     ica, ihpv, iscr(ngastp,2)

      iscr(1:ngastp,1:2) = 0

      if(ixdn.eq.1)then
        iscr(ihole,2)=iocc(ihole,2)
        iscr(ipart,1)=iocc(ipart,1)
        iscr(ivale,1)=iocc(ivale,1) ! treat val. like particle for now
        iscr(iextr,1)=iocc(iextr,1)
      elseif(ixdn.eq.2)then
        iscr(ihole,1)=iocc(ihole,1)
        iscr(ipart,2)=iocc(ipart,2)
        iscr(ivale,2)=iocc(ivale,2)
        iscr(iextr,2)=iocc(iextr,2)
      elseif(ixdn.eq.3)then
        iscr(ivale,1:2)=iocc(ivale,1:2)
      else
        call quit(1,'iocc_xdn','undefined value of ixdn')
      endif

      iocc_xdn = iscr

      return
      end
