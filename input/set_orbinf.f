*----------------------------------------------------------------------*
      subroutine set_orbinf(orb_info,env_type,hole_rv)
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
      character, intent(in) ::
     &     env_type*(*)
      logical, intent(in) ::
     &     hole_rv

      integer ::
     &     ngas, nsym, ntoob, nspin,
     &     iprint, idx, jdx, kdx, isym, igas, igastp, ispin, ipass,
     &     idxst, idxnd, ist, ind, inc, igasr, j, caborb, iloop, loop
      integer ::
     &     icount(ngastp)
      integer, allocatable ::
     &     iadscr(:),koffs(:)
      logical, pointer ::
     &     assigned(:)


      iprint = max(ntest,iprlvl)

      if (iprint.ge.100) then
        write(lulog,*) '************'
        write(lulog,*) ' set_orbinf'
        write(lulog,*) '************'
        write(lulog,*) ' hole_rv: ',hole_rv
        write(lulog,*) ' nspin  = ',orb_info%nspin
        write(lulog,*) ' ngas   = ',orb_info%ngas
        write(lulog,*) ' nsym   = ',orb_info%nsym
        write(lulog,*) ' ntoob  = ',orb_info%ntoob
        write(lulog,*) ' caborb = ',orb_info%caborb
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
        write(lulog,*) 'ihpvgas:   ',orb_info%ihpvgas(1:ngas,1)
        if (orb_info%nspin.gt.1)
     &       write(lulog,*) '           ',orb_info%ihpvgas(1:ngas,2)
        write(lulog,*) 'ngas_hpv:  ',orb_info%ngas_hpv(1:ngastp)
        write(lulog,*) 'nactt_hpv: ',orb_info%nactt_hpv(1:ngastp)
        write(lulog,*) 'ioff_gas:  ',orb_info%ioff_gas(1:ngastp)
        write(lulog,*) 'idx_gas:   ',orb_info%idx_gas(1:ngastp)
        write(lulog,*) 'igassh:'
        do igas = 1, ngas
          write(lulog,'(3x,8i4)') orb_info%igassh(1:nsym,igas)
        end do
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
        call quit(1,'set_orbinf','who needs hole_rv?? (1)')
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
      ist = ind+1
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
        write(lulog,*) 'gas_reo:',orb_info%gas_reo(1:ngas)
        write(lulog,*) 'igamorb:'
        call iwrtma(orb_info%igamorb,1,ntoob+caborb,1,ntoob+caborb)
        write(lulog,*) 'igasorb:'
        call iwrtma(orb_info%igasorb,1,ntoob+caborb,1,ntoob+caborb)
        write(lulog,*) 'mostnd:'
        do igas = 1, ngas
          write(lulog,'(2x,i4,2x,8(x,2i4))')
     &         igas,orb_info%mostnd(1:2,1:nsym,igas)
        end do
      end if

      ! generate symmetry ordering -> type ordering mapping
      select case(trim(env_type))
      case('dalton','DALTON','dalton_special','DALTON_SPECIAL',
     &     'dalton64','DALTON64','molpro_dump','MOLPRO_DUMP')
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
      case ('molpro_ifc','MOLPRO_IFC','virtual','VIRTUAL')
        ! the interface removes the core orbitals and shifts the numbering
        ! we just pretend that the symmetry ordering is done for core and other
        ! orbitals separately
        ! the "virtual molecule" is just handled analogous to that
        ! (although we will never do actual calculations with it)
        if (hole_rv) call quit(1,'set_orbinf','who needs hole_rv??')
        jdx = 0
        do ipass = 1,2 ! two passes: 1-deleted core, 2-others
          do isym = 1, nsym
            do igas = 1, ngas
              if (ipass.eq.1.and.orb_info%iad_gas(igas).ne.1) cycle
              if (ipass.eq.2.and.orb_info%iad_gas(igas).eq.1) cycle
              ! ignore cabs orbitals
              if(caborb.gt.0.and.igas.eq.ngas)cycle
              ist = orb_info%mostnd(1,isym,igas)
              ind = orb_info%mostnd(2,isym,igas)
              inc = +1
              do idx = ist, ind, inc
                jdx = jdx+1
                orb_info%ireost(jdx) = idx
              end do
            end do
          end do
        end do

        ! append CABS orbitals (no reordering)
        do idx = ntoob+1, ntoob+caborb
          orb_info%ireost(idx) = idx
        end do
      case ('cfour','CFOUR')
        ! actually: symmetry ordering means "the way the external program 
        !   has ordered the orbitals"; cfour orders kind-of-type-ordering like
        !   we try to get this info by interpreting the ext_gamorb array
        call make_ext2typ_reo(orb_info)
      case ('gamess','GAMESS')
        ! (We loop over all orbitals to avoid errors when the active
        !  space is modified such that orbitals must be reordered.
        !  The DALTON inferface above also needs to be changed
        !  if we want to use this nasty kind of tricks there as well)
        if (orb_info%n_bound_orbs.ne.ntoob.or.caborb.ne.0)
     &        call quit(1,'set_orbinf','need full orb-sym list!')
        allocate(assigned(ntoob))
        assigned = .false.
        do igas = 1, ngas
          do isym = 1, nsym
            if (orb_info%igassh(isym,igas).eq.0) cycle
            jdx = orb_info%mostnd(1,isym,igas)
            do idx = 1, ntoob
              if (assigned(idx)) cycle
              if (orb_info%isym_bound_orbs(idx).eq.isym) then
                orb_info%ireost(idx) = jdx
                jdx = jdx + 1
                assigned(idx) = .true.
                if (jdx.gt.orb_info%mostnd(2,isym,igas)) exit
              end if
            end do
          end do
        end do
        if (.not.all(assigned(1:ntoob)))
     &    call quit(1,'set_orbinf','unable to assign all orbitals')
        deallocate(assigned)
      case default
        call quit(1,'set_orbinf','unknown type '//trim(env_type))
      end select

        ! generate reverse mapping
        do idx = 1, ntoob+caborb
          orb_info%ireots(orb_info%ireost(idx)) = idx 
        end do

      if (iprint.ge.100) then
        write(lulog,*) 'ireost:'
        call iwrtma(orb_info%ireost,1,ntoob,1,ntoob+caborb)
        write(lulog,*) 'ireots:'
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
        write(lulog,*) 'iad_gas:'
        call iwrtma(orb_info%iad_gas,ngas,1,ngas,1)
      end if

      return
      end
