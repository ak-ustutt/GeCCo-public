      subroutine set_csorbs_molpro(csorb,orb_info)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(inout), target ::
     &     orb_info
      logical, intent(out) ::
     &     csorb(orb_info%ntoob)

      integer, pointer ::
     &     hpvxgas(:,:), igasorb(:), ireost(:), ntoob

      integer ::
     &     imomo, imogc

      ntoob => orb_info%ntoob
      hpvxgas => orb_info%ihpvgas
      igasorb => orb_info%igasorb
      ireost => orb_info%ireost

      ! loop over orbitals in molpro order
      do imomo = 1, ntoob
        ! get orbital in gecco type ordering
        imogc = ireost(imomo)
        csorb(imomo) = hpvxgas(igasorb(imogc),1).eq.IHOLE
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'set_csorbs_molpro: csorb = '
        write(lulog,'(1x,5l2,1x,5l2,2x,5l2,1x,5l2)') csorb
      end if

      end
