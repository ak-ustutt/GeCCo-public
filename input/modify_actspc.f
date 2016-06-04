*----------------------------------------------------------------------*
      subroutine modify_actspc(ishell,len,nactel,orb_info,mode)
*----------------------------------------------------------------------*
*     split the first shell
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest =  00

      integer, intent(in) ::
     &     len, ishell(len), mode, nactel
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
        write(lulog,*) '------------------'
        if (mode.eq.1) then
          write(lulog,*) ' modify_actspc: change inactive orbitals'
        else
          write(lulog,*) ' modify_actspc: change active orbitals'
        end if
        write(lulog,*) '------------------'
      end if

      if (len.ne.nsym)
     &     call quit(0,'modify_actspc',
     &     'definition inconsistent with symmetry')

      if (orb_info%ihpvgas(2,1).ne.3) then
        if (orb_info%ihpvgas(2,1).eq.1)
     &     call quit(0,'modify_actspc',
     &     '2nd GAS must be valence. define frozen shells afterwards.')

        ! no active shell yet: add an empty one!
        allocate(igassh_sv(nsym,ngas),iad_gas_sv(ngas),
     &       ihpvgas_sv(ngas,nspin))
        igassh_sv(1:nsym,1:ngas) = orb_info%igassh(1:nsym,1:ngas)
        iad_gas_sv(1:ngas) = orb_info%iad_gas(1:ngas)
        ihpvgas_sv(1:ngas,1:nspin) = orb_info%ihpvgas(1:ngas,1:nspin)
        deallocate(orb_info%igassh,orb_info%iad_gas,orb_info%ihpvgas)
        orb_info%ngas = orb_info%ngas+1
        ngas = ngas+1
        allocate(orb_info%igassh(nsym,ngas),orb_info%iad_gas(ngas),
     &       orb_info%ihpvgas(ngas,nspin))
        orb_info%igassh(1:nsym,1) = igassh_sv(1:nsym,1)
        orb_info%igassh(1:nsym,2) = 0
        orb_info%igassh(1:nsym,3:ngas) = igassh_sv(1:nsym,2:ngas-1)
        orb_info%iad_gas(1) = iad_gas_sv(1)
        orb_info%iad_gas(2) = 2
        orb_info%iad_gas(3:ngas) = iad_gas_sv(2:ngas-1)
        orb_info%ihpvgas(1,1:nspin) = ihpvgas_sv(1,1:nspin)
        orb_info%ihpvgas(2,1:nspin) = 3
        orb_info%ihpvgas(3:ngas,1:nspin) = ihpvgas_sv(2:ngas-1,1:nspin)
        deallocate(igassh_sv,iad_gas_sv,ihpvgas_sv)
      end if

      if (mode.ne.1.and.mode.ne.2)
     &     call quit(0,'modify_actspc',
     &     'mode must be 1 (inact) or 2 (act)')

      if (nactel.ge.0.and.nactel.ne.orb_info%nactel) then
        orb_info%nactel = nactel
        if (iprint.ge.50) write(lulog,'(x,a,i4)')
     &       'new number of active electrons: ',nactel
      end if

      do idx = 1, nsym
        if (ishell(idx).ge.0) then
          orb_info%igassh(idx,mode) = ishell(idx)
          orb_info%igassh(idx,3) = orb_info%ntoobs(idx)
     &         - orb_info%igassh(idx,1) - orb_info%igassh(idx,2)
        end if
      end do
      orb_info%nactorb = sum(ishell(1:nsym))

      if (iprint.ge.50) then
        write(lulog,'(x,a,8i4)') 'inactive occupied orbitals: ',
     &          orb_info%igassh(1:nsym,1)
        write(lulog,'(x,a,8i4)') 'CAS orbitals:               ',
     &          orb_info%igassh(1:nsym,2)
        write(lulog,'(x,a,8i4)') 'inactive virtual orbitals:  ',
     &          orb_info%igassh(1:nsym,3)
      end if

      end
