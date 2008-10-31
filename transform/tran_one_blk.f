*----------------------------------------------------------------------*
      subroutine tran_one_blk(xmo,xao,cmo,xhlf,
     &     isym,idxcmo,hpvx_c,hpvx_a,
     &     nbas,nxbas,mostnd,iad_gas,hpvx_gas,ngas,nsym)
*----------------------------------------------------------------------*
*     transform one-particle block (of an operator) to the
*     MO basis
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'multd2h.h'
      include 'opdim.h'
      include 'hpvxseq.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ngas, nsym, isym, hpvx_c, hpvx_a,
     &     idxcmo(nsym,ngas,2),
     &     iad_gas(ngas), hpvx_gas(ngas),
     &     mostnd(2,nsym,ngas), nbas(nsym), nxbas(nsym)
      real(8), intent(in) ::
     &     cmo(*), xmo(*)
      real(8), intent(inout) ::
     &     xao(*), xhlf(*)

      logical ::
     &     transp, xbas_c, xbas_a
      integer ::
     &     icmoc, icmoa, ixao, ixmo, igasc, igasa,
     &     isymc, isyma, nmoa, nmoc, naoa, naoc, ntaoa, ntaoc
      real(8) :: fac


      fac = 1d0

      ixao = 1
      ixmo = 1

      transp = hpvxseq(hpvx_c).gt.hpvxseq(hpvx_a)
      xbas_c = hpvx_c.eq.IEXTR
      xbas_a = hpvx_a.eq.IEXTR

      do isyma = 1, nsym
        isymc = multd2h(isyma,isym)

        if (ntest.ge.100)
     &       write(luout,*) 'isymc,isyma: ',isymc,isyma

        naoa = nbas(isyma)
        naoc = nbas(isymc)
        ntaoa = nbas(isyma)+nxbas(isyma)
        ntaoc = nbas(isymc)+nxbas(isymc)
 
        if (xbas_c) naoc = ntaoc
        if (xbas_a) naoa = ntaoa

        if (naoa*naoc.eq.0) cycle

        do igasa = 1, ngas
          if (iad_gas(igasa).ne.2) cycle
          if (hpvx_gas(igasa).ne.hpvx_a) cycle

          if (xbas_a) then
            icmoa  = idxcmo(isyma,igasa,2)
          else
            icmoa  = idxcmo(isyma,igasa,1)
          end if

          nmoa = mostnd(2,isyma,igasa)-mostnd(1,isyma,igasa)+1

          do igasc = 1, ngas
            if (iad_gas(igasc).ne.2) cycle
            if (hpvx_gas(igasc).ne.hpvx_c) cycle

            if (xbas_c) then
              icmoc = idxcmo(isymc,igasc,2)
            else
              icmoc = idxcmo(isymc,igasc,1)
            end if

            if (ntest.ge.100)
     &           write(luout,*)
     &           'igasc,igasa,icmoc,icmoa: ',igasc,igasa,icmoc,icmoa

            nmoc = mostnd(2,isymc,igasc)-mostnd(1,isymc,igasc)+1
          
            if (nmoa*nmoc.gt.0) then

              if (ntest.ge.100)
     &             write(luout,*) 'nmoa,nmoc,naoa,naoc,ntaoc: ',
     &                 nmoa,nmoc,naoa,naoc,ntaoc

              if (.not.transp) then

                if (ntest.ge.100)
     &              write(luout,*) 'no transp> xao(ixao) = ',xao(ixao)
                if (ntest.ge.100)
     &              write(luout,*) 'no transp> cmo(icmoc) = ',cmo(icmoc)

                call dgemm('t','n',nmoc,naoa,naoc,
     &                    1d0,cmo(icmoc),naoc,
     &                        xao(ixao),ntaoc,
     &                    0d0,xhlf,nmoc)

                if (ntest.ge.100)
     &              write(luout,*) 'no transp> xhlf(1) = ',xhlf(1)
                if (ntest.ge.100)
     &              write(luout,*) 'no transp> cmo(icmoa) = ',cmo(icmoa)

                call dgemm('n','n',nmoc,nmoa,naoa,
     &                   fac,xhlf,nmoc,
     &                       cmo(icmoa),naoa,
     &                   0d0,xmo(ixmo),nmoc)

              if (ntest.ge.100)
     &               write(luout,'(a,4i4)')
     &               'result block ',isymc,isyma,igasc,igasa
              if (ntest.ge.100)
     &             call wrtmat2(xmo(ixmo),nmoc,nmoa,nmoc,nmoa)

              else

                if (ntest.ge.100)
     &              write(luout,*) 'transp   > xao(ixao) = ',xao(ixao)
                if (ntest.ge.100)
     &              write(luout,*) 'transp   > cmo(icmoc) = ',cmo(icmoc)

                call dgemm('t','n',nmoc,naoa,naoc,
     &                    1d0,cmo(icmoc),naoc,
     &                        xao(ixao),ntaoc,
     &                    0d0,xhlf,nmoc)

                if (ntest.ge.100)
     &              write(luout,*) 'transp   > xhlf(1) = ',xhlf(1)
                if (ntest.ge.100)
     &              write(luout,*) 'transp   > cmo(icmoa) = ',cmo(icmoa)

                call dgemm('t','t',nmoa,nmoc,naoa,
     &                   fac,cmo(icmoa),naoa,
     &                       xhlf,nmoc,
     &                   0d0,xmo(ixmo),nmoa)

                if (ntest.ge.100)
     &               write(luout,'(a,4i4)')
     &               'result block ',isymc,isyma,igasc,igasa
                if (ntest.ge.100)
     &               call wrtmat2(xmo(ixmo),nmoa,nmoc,nmoa,nmoc)

              end if

            end if

            ixmo = ixmo + nmoc*nmoa
          end do
        end do

        ixao = ixao + ntaoc*ntaoa
      end do

      return
      end
