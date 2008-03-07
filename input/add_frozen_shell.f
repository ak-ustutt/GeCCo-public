*----------------------------------------------------------------------*
      subroutine add_frozen_shell(ishell,len,orb_info)
*----------------------------------------------------------------------*
*     split the first shell
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     len, ishell(len)
      type(orbinf), intent(inout) ::
     &     orb_info

      integer ::
     &     idx, ierr, nsym, ngas, nspin, igas, iprint
      integer, allocatable ::
     &     igassh_sv(:,:), iad_gas_sv(:), ihpvgas_sv(:,:)

      iprint = max(iprlvl,ntest)

      nsym  = orb_info%nsym
      ngas  = orb_info%ngas
      nspin = orb_info%nspin
      if (iprint.ge.50) then
        write(luout,*) '------------------'
        write(luout,*) ' add_frozen_shell'
        write(luout,*) '------------------'
        write(luout,*) ' shell to add: ',ishell(1:len)
        write(luout,*) ' defined shells: '
        do igas = 1, ngas
          write(luout,'(x,i3,2x,8i4)') igas,orb_info%igassh(1:nsym,igas)
        end do
      end if

      if (len.ne.nsym)
     &     call quit(0,'add_frozen_shell',
     &     'definition inconsistent with symmetry')

      if (ngas.lt.1)
     &     call quit(0,'add_frozen_shell',
     &     'no shells defined yet')

      ierr = 0
      do idx = 1, len
        if (ishell(idx).gt.orb_info%igassh(idx,1)) ierr = ierr+1
      end do
      if (ierr.gt.0) then
        write(luout,*) 'orbital to freeze: ',ishell(1:len)
        write(luout,*) ' current shell #1: ',orb_info%igassh(1:len,1)
        call quit(0,'add_frozen_shell','inconsistency')
      end if

      allocate(igassh_sv(nsym,ngas),iad_gas_sv(ngas),
     &     ihpvgas_sv(ngas,nspin))
      
      igassh_sv(1:nsym,1:ngas) = orb_info%igassh(1:nsym,1:ngas)
      iad_gas_sv(1:ngas) = orb_info%iad_gas(1:ngas)
      ihpvgas_sv(1:ngas,1:nspin) = orb_info%ihpvgas(1:ngas,1:nspin)

      ! re-allocate array
      deallocate(orb_info%igassh,orb_info%iad_gas,orb_info%ihpvgas)
      
      orb_info%ngas = orb_info%ngas+1
      ngas = ngas+1

      allocate(orb_info%igassh(nsym,ngas),orb_info%iad_gas(ngas),
     &     orb_info%ihpvgas(ngas,nspin))

      orb_info%igassh(1:nsym,1) = ishell(1:nsym)
      orb_info%igassh(1:nsym,2) = igassh_sv(1:nsym,1)-ishell(1:nsym)

      do igas = 3, ngas      
        orb_info%igassh(1:nsym,igas) = igassh_sv(1:nsym,igas-1)
      end do

      orb_info%iad_gas(1) = 1
      orb_info%iad_gas(2:ngas) = iad_gas_sv(1:ngas-1)

      orb_info%ihpvgas(1,1:nspin) = 1
      orb_info%ihpvgas(2:ngas,1:nspin) = ihpvgas_sv(1:ngas-1,1:nspin)

      if (iprint.ge.50) then
        write(luout,*) ' new shell definition: '
        do igas = 1, ngas
          write(luout,'(x,i3,2x,8i4)') igas,orb_info%igassh(1:nsym,igas)
        end do
      end if

      deallocate(igassh_sv,iad_gas_sv,ihpvgas_sv)

      end
