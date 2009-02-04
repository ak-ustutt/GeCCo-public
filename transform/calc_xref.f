*----------------------------------------------------------------------*
      subroutine calc_xref(xref,xao,cmo,xhlf,nsym,nbas,nxbas,
     &                     ngas,hpvx_gas,idxcmo,mostnd)
*----------------------------------------------------------------------*
*
*     calculates reference expectation value for x-operator
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      real(8), intent(out) ::
     &     xref

      real(8), intent(in) ::
     &     cmo(*), xao(*)

      real(8), intent(inout) ::
     &     xhlf(*)

      integer, intent(in) ::
     &     nsym, nbas(nsym), nxbas(nsym), ngas, hpvx_gas(ngas),
     &     idxcmo(nsym,ngas,1), mostnd(2,nsym,ngas)

      integer ::
     &     isym, igas, ixao, icmoc, norb

      real(8), external ::
     &     ddot

      ixao = 1
      xref = 0d0

      do isym = 1, nsym

        if (ntest.ge.100)
     &       write(luout,*) 'isym: ',isym

        if (nbas(isym).eq.0) cycle

        do igas = 1, ngas

          if (ntest.ge.100)
     &       write(luout,*) 'igas,hpvx_gas(igas): ',igas,hpvx_gas(igas)

          if (hpvx_gas(igas).ne.1) cycle

          icmoc = idxcmo(isym,igas,1)
          norb = mostnd(2,isym,igas)-mostnd(1,isym,igas)+1

          if (ntest.ge.100)
     &         write(luout,*) 'norb,nbas(isym),ntao: ',
     &             norb,nbas(isym),nbas(isym)+nxbas(isym)

          if (norb.gt.0) then

            if (ntest.ge.100)
     &            write(luout,*) 'xao(ixao) = ',xao(ixao)
            if (ntest.ge.100)
     &            write(luout,*) 'cmo(icmoc) = ',cmo(icmoc)

            call dgemm('n','n',nbas(isym),norb,nbas(isym),
     &                1d0,xao(ixao),nbas(isym)+nxbas(isym),
     &                    cmo(icmoc),nbas(isym),
     &                0d0,xhlf,nbas(isym))

            if (ntest.ge.100)
     &            write(luout,*) 'xhlf(1) = ',xhlf(1)

            xref = xref + 2*ddot(nbas(isym)*norb,xhlf(1),1,cmo(icmoc),1)

            if (ntest.ge.100)
     &            write(luout,*) 'xref = ',xref

          end if
        end do
        ixao = ixao + (nbas(isym) + nxbas(isym))**2
      end do

      if (ntest.ge.100)
     &          write(luout,*) 'total xref = ',xref

      return
      end
