*----------------------------------------------------------------------*
      subroutine modify_actspc(ishell,len,nactel,orb_info,mode)
*----------------------------------------------------------------------*
*     split the first shell
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 100

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
        write(luout,*) '------------------'
        if (mode.eq.1) then
          write(luout,*) ' modify_actspc: change inactive orbitals'
        else
          write(luout,*) ' modify_actspc: change cas orbitals'
        end if
        write(luout,*) '------------------'
      end if

      if (len.ne.nsym)
     &     call quit(0,'modify_actspc',
     &     'definition inconsistent with symmetry')

      if (orb_info%ihpvgas(2,1).ne.3)
     &     call quit(0,'modify_actspc',
     &     '2nd GAS must be valence. define frozen shells afterwards.')

      if (mode.ne.1.and.mode.ne.2)
     &     call quit(0,'modify_actspc',
     &     'mode must be 1 (inact) or 2 (cas)')

      if (nactel.ge.0.and.nactel.ne.orb_info%nactel) then
        orb_info%nactel = nactel
        if (iprint.ge.50) write(luout,'(x,a,i4)')
     &       'new number of active electrons: ',nactel
      end if

      do idx = 1, nsym
        if (ishell(idx).ge.0) then
          orb_info%igassh(idx,mode) = ishell(idx)
          orb_info%igassh(idx,3) = orb_info%ntoobs(idx)
     &         - orb_info%igassh(idx,1) - orb_info%igassh(idx,2)
        end if
      end do

      if (iprint.ge.50) then
        write(luout,'(x,a,8i4)') 'inactive occupied orbitals: ',
     &          orb_info%igassh(1:nsym,1)
        write(luout,'(x,a,8i4)') 'CAS orbitals:               ',
     &          orb_info%igassh(1:nsym,2)
        write(luout,'(x,a,8i4)') 'inactive virtual orbitals:  ',
     &          orb_info%igassh(1:nsym,3)
      end if

      end
