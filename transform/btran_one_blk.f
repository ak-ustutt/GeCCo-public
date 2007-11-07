*----------------------------------------------------------------------*
      subroutine btran_one_blk(xao,xmo,cmo,xhlf,
     &     isym,idxcmo,hpvx_c,hpvx_a,
     &     nbas,mostnd,iad_gas,hpvx_gas,ngas,nsym)
*----------------------------------------------------------------------*
*     back-transform one-particle block (of a density) to the
*     AO basis
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     ngas, nsym, isym, hpvx_c, hpvx_a,
     &     idxcmo(nsym,ngas),
     &     iad_gas(ngas), hpvx_gas(ngas),
     &     mostnd(2,nsym,ngas), nbas(nsym)
      real(8), intent(in) ::
     &     cmo(*), xmo(*)
      real(8), intent(inout) ::
     &     xao(*), xhlf(*)

      integer ::
     &     icmoc, icmoa, ixao, ixmo, igasc, igasa,
     &     isymc, isyma, nmoa, nmoc, naoa, naoc

c dbg
c      print *,'hpvx_ca: ',hpvx_c, hpvx_a
c dbg
      ixao = 1
      ixmo = 1

      do isyma = 1, nsym
        isymc = multd2h(isyma,isym)

        naoa = nbas(isyma)
        naoc = nbas(isymc)
c dbg
c        print *,'nbas: ',naoa,naoc
c dbg
        if (naoa*naoc.eq.0) cycle

        do igasa = 1, ngas
c dbg
c          print *,'igasa : ',igasa,iad_gas(igasa),hpvx_gas(igasa)
c dbg
          if (iad_gas(igasa).ne.2) cycle
          if (hpvx_gas(igasa).ne.hpvx_a) cycle

          icmoa = idxcmo(isyma,igasa)

          nmoa = mostnd(2,isyma,igasa)-mostnd(1,isyma,igasa)+1

          do igasc = 1, ngas
c dbg
c          print *,'igasc : ',igasc,iad_gas(igasc),hpvx_gas(igasc)
c dbg
            if (iad_gas(igasc).ne.2) cycle
            if (hpvx_gas(igasc).ne.hpvx_c) cycle

            icmoc = idxcmo(isymc,igasc)

            nmoc = mostnd(2,isymc,igasc)-mostnd(1,isymc,igasc)+1
c dbg
c            print *,'nmo: ',nmoa,nmoc
c dbg
          
            if (nmoa*nmoc.gt.0) then

c dbg
c            print *,'isym, igas: ',isym, igas
c              print *,'cmo c',icmoc
c              call wrtmat2(cmo(icmoc),naoc,nmoc,naoc,nmoc)
c              print *,'xmo ',ixmo
c              call wrtmat2(xmo(ixmo),nmoc,nmoa,nmoc,nmoa)
c dbg
              call dgemm('n','n',naoc,nmoa,nmoc,
     &                  1d0,cmo(icmoc),naoc,
     &                      xmo(ixmo),nmoc,
     &                  0d0,xhlf,naoc)
c dbg
c              print *,'xhlf '
c              call wrtmat2(xhlf,naoc,nmoa,naoc,nmoa)
c              print *,'cmo a',icmoa
c              call wrtmat2(cmo(icmoa),naoa,nmoa,naoa,nmoa)
c dbg

              call dgemm('n','t',naoc,naoa,nmoa,
     &                  1d0,xhlf,naoc,
     &                      cmo(icmoa),naoa,
     &                  1d0,xao(ixao),naoc)
c dbg
c              print *,'xao ',ixao
c              call wrtmat2(xao(ixao),naoc,naoa,naoc,naoa)
c dbg
            end if

            ixmo = ixmo + nmoc*nmoa
          end do
        end do

        ixao = ixao + naoc*naoa
      end do

      return
      end
