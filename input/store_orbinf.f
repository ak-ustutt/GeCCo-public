*----------------------------------------------------------------------*
      subroutine store_orbinf(ffarch,orb_info)
*----------------------------------------------------------------------*
*     store essential info from orb_info on sequential file ffarch
*     (should be open and positioned)
*----------------------------------------------------------------------*
      implicit none

      include 'def_orbinf.h'
      include 'def_filinf.h'

      type(filinf), intent(inout) ::
     &     ffarch
      type(orbinf), intent(in) ::
     &     orb_info

      integer, parameter ::
     &     maxbuf = 512
      integer(4) ::
     &     buffer(maxbuf)
      integer(4) ::
     &     idx

      integer ::
     &     unit, nsym, ngas, nspin, igas, ispin

      unit = ffarch%unit
      if (unit.le.0) call quit(1,'store_orbinf',
     &     'archive file not open? '//trim(ffarch%name))

      write(unit) 'ORB_INFO'
      ! store dimensions (make sure that it is integer(4))
      buffer(1:6) =
     &            (/orb_info%nspin,orb_info%nsym,
     &              orb_info%ngas, orb_info%ntoob,
     &              orb_info%nbast,orb_info%caborb/) 
      write(unit) buffer(1:6)

      nsym = orb_info%nsym
      ngas = orb_info%ngas
      nspin = orb_info%nspin
      idx = 0
      buffer(idx+1:idx+nsym) = orb_info%nbas(1:nsym)
      idx = idx+nsym
      buffer(idx+1:idx+nsym) = orb_info%ntoobs(1:nsym)
      idx = idx+nsym
      do igas = 1, ngas
        buffer(idx+1:idx+nsym) = orb_info%igassh(1:nsym,igas)
        idx = idx+nsym
      end do
      buffer(idx+1:idx+ngas) = orb_info%iad_gas(1:ngas)
      idx = idx+ngas
      do ispin = 1, nspin
        buffer(idx+1:idx+ngas) = orb_info%ihpvgas(1:ngas,ispin)
        idx = idx+ngas
      end do
      buffer(idx+1:idx+nsym) = orb_info%cab_orb(1:nsym)
      idx = idx+nsym

      write(unit) buffer(1:idx)

      return
      end
