*----------------------------------------------------------------------*
      subroutine set_orbinf(orb_info,hole_rv)
*----------------------------------------------------------------------*
*     set up orbital info arrays
*
*     symmetry ordering: 
*        orbitals are ordered according to symmetry only
*
*     type ordering:
*        orbitals are ordered first according to the shell they
*        belong to (and within shell according to their symmetry)
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(inout) ::
     &     orb_info
      logical, intent(in) ::
     &     hole_rv

      integer ::
     &     ngas, nsym, ntoob, nspin,
     &     iprint, idx, jdx, kdx, isym, igas, igastp, ispin,
     &     idxst, idxnd, ist, ind, inc, igasr, j, caborb, iloop, loop
      integer ::
     &     icount(ngastp)
      integer, allocatable ::
     &     iadscr(:),koffs(:)


      iprint = max(ntest,iprlvl)

      if (iprint.ge.100) then
        write(luout,*) '************'
        write(luout,*) ' set_orbinf'
        write(luout,*) '************'
      end if

      ! for convenience
      nspin = orb_info%nspin
      ngas = orb_info%ngas
      nsym = orb_info%nsym
      ntoob = orb_info%ntoob
      caborb=orb_info%caborb 

      if (hole_rv.and.nspin.gt.1) then
        call quit(1,'set_orbinf',
     &       'hole_rv is not compatible with nspin.gt.1')
      end if

      ! allocate some arrays on orb_info structure
      allocate(orb_info%ireots(ntoob+caborb),
     &     orb_info%ireost(ntoob+caborb),orb_info%igamorb(ntoob+caborb), 
     &     orb_info%igasorb(ntoob+caborb),orb_info%mostnd(2,nsym,ngas),
     &     orb_info%ngas_hpv(ngastp),orb_info%nactt_hpv(ngastp),
     &     orb_info%idx_gas(ngastp),orb_info%ioff_gas(ngastp),
     &     orb_info%gas_reo(ngas),orb_info%norb_hpv(ngastp,nspin))

      ! set the info arrays:

      ! set ngas_hpv (number of spaces per H/P/V space)
      ! set nactt_hpv (number of active orbitals per H/P/V)
      ! open shell with p/h splitting: we double-count!
      orb_info%ngas_hpv(1:ngastp) = 0
      orb_info%nactt_hpv(1:ngastp) = 0
      orb_info%norb_hpv(1:ngastp,1:nspin) = 0
      do ispin = 1, nspin
        icount(1:ngastp) = 0
        do igas = 1, ngas
          icount(orb_info%ihpvgas(igas,ispin)) =
     &         icount(orb_info%ihpvgas(igas,ispin))+1
        end do
        do igastp = 1, ngastp
          orb_info%ngas_hpv(igastp) =
     &         max(orb_info%ngas_hpv(igastp),icount(igastp))
        end do

        icount(1:ngastp) = 0
        do igas = 1, ngas
          if (orb_info%iad_gas(igas).ne.2) cycle
          do isym = 1, nsym
            icount(orb_info%ihpvgas(igas,ispin)) =
     &           icount(orb_info%ihpvgas(igas,ispin))+
     &           orb_info%igassh(isym,igas)
          end do
        end do
        do igastp = 1, ngastp
          orb_info%nactt_hpv(igastp) =
     &         max(orb_info%nactt_hpv(igastp),icount(igastp))
        end do

        do igas = 1, ngas
          do isym = 1, nsym
            orb_info%norb_hpv(orb_info%ihpvgas(igas,ispin),ispin) =
     &           orb_info%norb_hpv(orb_info%ihpvgas(igas,ispin),ispin)+
     &           orb_info%igassh(isym,igas)
          end do
        end do

      end do

      ! index and offset array from ngas_hpv:
c      idx = 1
cmh adapted for use of valence space
      orb_info%ioff_gas(1:ngas) = 0
      orb_info%idx_gas(1:ngas) = 0
      do igastp = 1, ngastp
        do idx = ngas,1,-1
          orb_info%ioff_gas(orb_info%ihpvgas(idx,1)) = idx-1
          orb_info%idx_gas(orb_info%ihpvgas(idx,1)) = idx
        end do
c        orb_info%ioff_gas(igastp) = idx-1
c        orb_info%idx_gas(igastp) = idx
c        idx = idx+orb_info%ngas_hpv(igastp)
      end do

      if (iprint.ge.100) then
        write(luout,*) 'ihpvgas:   ',orb_info%ihpvgas(1:ngas,1)
        if (orb_info%nspin.gt.1)
     &       write(luout,*) '           ',orb_info%ihpvgas(1:ngas,2)
        write(luout,*) 'ngas_hpv:  ',orb_info%ngas_hpv(1:ngastp)
        write(luout,*) 'nactt_hpv: ',orb_info%nactt_hpv(1:ngastp)
        write(luout,*) 'ioff_gas:  ',orb_info%ioff_gas(1:ngastp)
        write(luout,*) 'idx_gas:   ',orb_info%idx_gas(1:ngastp)
      end if

      do igas = 1, ngas
        orb_info%gas_reo(igas) = igas
      end do
      if (hole_rv) then
        igasr = 0
        do igas = orb_info%ngas_hpv(ihole), 1, -1          
          igasr = igasr + 1
          orb_info%gas_reo(igas) = igasr
        end do
      end if

      ! set up igamorb (IRREP per orbital in type ordering) and
      ! set up igasorb (shell per orbital in type ordering)      
      ! set up mostnd (start and end indices)      
      ! start with hole spaces
      if (hole_rv) then
        ! reverse counting
        ist = orb_info%ngas_hpv(ihole)
        ind = 1
        inc = -1
      else        
        ist = 1
        ind = orb_info%ngas_hpv(ihole)
        inc = +1
      end if
      idxst = 1
      igasr = 0
      do igas = ist, ind, inc
        igasr = igasr+1
        do isym = 1, nsym
          idxnd = idxst+orb_info%igassh(isym,igas)-1
          orb_info%mostnd(1:2,isym,igasr) = (/idxst,idxnd/)
          if (idxst.le.idxnd) then
            orb_info%igamorb(idxst:idxnd) = isym 
            orb_info%igasorb(idxst:idxnd) = igasr
          end if
          idxst = idxnd+1
        end do
      end do
      ! remaining spaces
      ist = max(ist,ind)+1
      ind = ngas
      do igas = ist, ind
        do isym = 1, nsym
          idxnd = idxst+orb_info%igassh(isym,igas)-1
          orb_info%mostnd(1:2,isym,igas) = (/idxst,idxnd/)
          if (idxst.le.idxnd) then
            orb_info%igamorb(idxst:idxnd) = isym 
            orb_info%igasorb(idxst:idxnd) = igas 
          end if
          idxst = idxnd+1
        end do
      end do

      if (iprint.ge.100) then
        write(luout,*) 'gas_reo:',orb_info%gas_reo(1:ngas)
        write(luout,*) 'igamorb:'
        call iwrtma(orb_info%igamorb,1,ntoob+caborb,1,ntoob+caborb)
        write(luout,*) 'igasorb:'
        call iwrtma(orb_info%igasorb,1,ntoob+caborb,1,ntoob+caborb)
        write(luout,*) 'mostnd:'
        do igas = 1, ngas
          write(luout,'(2x,i4,2x,8(x,2i4))')
     &         igas,orb_info%mostnd(1:2,1:nsym,igas)
        end do
      end if

      ! generate symmetry ordering -> type ordering mapping
        jdx = 0
        do isym = 1, nsym
          do igas = 1, ngas
            ! ignore cabs orbitals
            if(caborb.gt.0.and.igas.eq.ngas)cycle
            if (orb_info%ihpvgas(igas,1).eq.ihole.and.hole_rv) then
              igasr = orb_info%ngas_hpv(ihole)-igas+1
              ist = orb_info%mostnd(2,isym,igasr)
              ind = orb_info%mostnd(1,isym,igasr)
              inc = -1
            else
              ist = orb_info%mostnd(1,isym,igas)
              ind = orb_info%mostnd(2,isym,igas)
              inc = +1
            end if
            do idx = ist, ind, inc
              jdx = jdx+1
              orb_info%ireost(jdx) = idx
            end do
          end do
        end do

        ! append CABS orbitals (no reordering)
        do idx = ntoob+1, ntoob+caborb
          orb_info%ireost(idx) = idx
        end do

        ! generate reverse mapping
        do idx = 1, ntoob+caborb
          orb_info%ireots(orb_info%ireost(idx)) = idx 
        end do

      if (iprint.ge.100) then
        write(luout,*) 'ireost:'
        call iwrtma(orb_info%ireost,1,ntoob,1,ntoob+caborb)
        write(luout,*) 'ireots:'
        call iwrtma(orb_info%ireots,1,ntoob,1,ntoob+caborb)
      end if

      ! reverse iad_array
      if (hole_rv) then
        allocate(iadscr(ngas))
        iadscr(1:ngas) = orb_info%iad_gas(1:ngas)
        do igas = 1, orb_info%ngas_hpv(ihole)
          orb_info%iad_gas(igas) =
     &         iadscr(orb_info%ngas_hpv(ihole)-igas+1)
        end do
        deallocate(iadscr)
      end if

      if (iprint.ge.100) then
        write(luout,*) 'iad_gas:'
        call iwrtma(orb_info%iad_gas,ngas,1,ngas,1)
      end if

      return
      end
