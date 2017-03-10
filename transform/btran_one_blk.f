*----------------------------------------------------------------------*
      subroutine btran_one_blk(xao,xmo,cmo,xhlf,
     &     isym,idxcmo,hpvx_c,hpvx_a,idxms,trplt,
     &     nbas,mostnd,iad_gas,hpvx_gas,ngas,nsym)
*----------------------------------------------------------------------*
*     back-transform one-particle block (of a density) to the
*     AO basis
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'
      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(in) ::
     &     ngas, nsym, isym, hpvx_c, hpvx_a,idxms,
     &     idxcmo(nsym,ngas),
     &     iad_gas(ngas), hpvx_gas(ngas),
     &     mostnd(2,nsym,ngas), nbas(nsym)
      logical, intent(in) :: 
     &     trplt
      real(8), intent(in) ::
     &     cmo(*), xmo(*)
      real(8), intent(inout) ::
     &     xao(*), xhlf(*)

      logical ::
     &     transp
      integer ::
     &     icmoc, icmoa, ixao, ixmo, igasc, igasa,
     &     isymc, isyma, nmoa, nmoc, naoa, naoc
      real(8) :: fac

      ixao = 1
      ixmo = 1

      fac = 1d0
      ! minus one for HH-block:
      !   C ... A
      ! propertiy integral would insert as  C A  and must be
      ! reorderes to A C for contraction with density --> fac=-1
      if (hpvx_c.eq.1.and.hpvx_a.eq.1) fac = -1d0
      if (trplt.and.idxms.eq.2) fac = -1d0*fac
      transp = hpvxseq(hpvx_c).gt.hpvxseq(hpvx_a)

      do isyma = 1, nsym
        isymc = multd2h(isyma,isym)

        naoa = nbas(isyma)
        naoc = nbas(isymc)
        if (naoa*naoc.eq.0) cycle

        do igasa = 1, ngas
          if (iad_gas(igasa).ne.2) cycle
          if (hpvx_gas(igasa).ne.hpvx_a) cycle

          icmoa = idxcmo(isyma,igasa)

          nmoa = mostnd(2,isyma,igasa)-mostnd(1,isyma,igasa)+1

          do igasc = 1, ngas
            if (iad_gas(igasc).ne.2) cycle
            if (hpvx_gas(igasc).ne.hpvx_c) cycle
            icmoc = idxcmo(isymc,igasc)

            nmoc = mostnd(2,isymc,igasc)-mostnd(1,isymc,igasc)+1
          
            if (nmoa*nmoc.gt.0) then

              if (.not.transp) then

                call dgemm('n','n',naoc,nmoa,nmoc,
     &                    1d0,cmo(icmoc),naoc,
     &                        xmo(ixmo),nmoc,
     &                    0d0,xhlf,naoc)

              else

                call dgemm('n','t',naoc,nmoa,nmoc,
     &                    1d0,cmo(icmoc),naoc,
     &                        xmo(ixmo),nmoa,
     &                    0d0,xhlf,naoc)

              end if

              call dgemm('n','t',naoc,naoa,nmoa,
     &                  fac,xhlf,naoc,
     &                      cmo(icmoa),naoa,
     &                  1d0,xao(ixao),naoc)

            end if

            ixmo = ixmo + nmoc*nmoa
          end do
        end do

        ixao = ixao + naoc*naoa
      end do

      return
      end
