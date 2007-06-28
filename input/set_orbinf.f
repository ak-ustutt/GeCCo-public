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
      include 'explicit.h'

      integer, parameter ::
     &     ntest = 100

      type(orbinf), intent(inout) ::
     &     orb_info
      logical, intent(in) ::
     &     hole_rv

      integer ::
     &     ngas, nsym, ntoob, iprint, idx, jdx, kdx, isym, igas, igastp,
     &     idxst, idxnd, ist, ind, inc, igasr, j, caborb, iloop, loop
      integer, allocatable ::
     &     iadscr(:),koffs(:)


      iprint = max(ntest,iprlvl)

      if (iprint.ge.100) then
        write(luout,*) '************'
        write(luout,*) ' set_orbinf'
        write(luout,*) '************'
      end if

      ! for convenience
      ngas = orb_info%ngas
      nsym = orb_info%nsym
      ntoob = orb_info%ntoob
      caborb=orb_info%caborb 

      ! allocate some arrays on orb_info structure
      allocate(orb_info%ireots(ntoob-caborb),
     &     orb_info%ireost(ntoob-caborb),orb_info%igamorb(ntoob), 
     &     orb_info%igasorb(ntoob),orb_info%mostnd(2,nsym,ngas),
     &     orb_info%ngas_hpv(ngastp),orb_info%nactt_hpv(ngastp),
     &     orb_info%idx_gas(ngastp),orb_info%ioff_gas(ngastp))
      if(explicit)then
        allocate(orb_info%xreosym(ntoob),orb_info%xreotyp(ntoob),
     &       koffs(1:nsym))
      endif  

      ! set the info arrays:

      ! set ngas_hpv (number of spaces per H/P/V space)
      ! set nactt_hpv (number of active orbitals per H/P/V)
      orb_info%ngas_hpv(1:ngastp) = 0
      orb_info%nactt_hpv(1:ngastp) = 0
      do igas = 1, ngas
        orb_info%ngas_hpv(orb_info%ihpvgas(igas)) =
     &       orb_info%ngas_hpv(orb_info%ihpvgas(igas))+1
        if (orb_info%iad_gas(igas).ne.2) cycle
        do isym = 1, nsym
          orb_info%nactt_hpv(orb_info%ihpvgas(igas)) =
     &         orb_info%nactt_hpv(orb_info%ihpvgas(igas)) +
     &         orb_info%igassh(isym,igas)
        end do
      end do

      ! index and offset array from ngas_hpv:
      idx = 1
      do igastp = 1, ngastp
        orb_info%ioff_gas(igastp) = idx-1
        orb_info%idx_gas(igastp) = idx
        idx = idx+orb_info%ngas_hpv(igastp)
      end do

      if (iprint.ge.100) then
        write(luout,*) 'ngas_hpv:  ',orb_info%ngas_hpv(1:ngastp)
        write(luout,*) 'nactt_hpv: ',orb_info%nactt_hpv(1:ngastp)
        write(luout,*) 'ioff_gas:  ',orb_info%ioff_gas(1:ngastp)
        write(luout,*) 'idx_gas:   ',orb_info%idx_gas(1:ngastp)
      end if

c      ! If an R12 calculation is requested then must do an initial 
c      ! reordering, as the external orbitals are extracted separately
c      ! from the others. All other orbitals are written together in 
c      ! symmetry ordering (i.e. h/p/v in one list). It is necessary to
c      ! combine this list with that of the x-space (also symmetry 
c      ! ordered) before subsequent arrays can act on the total set.
c      if(explicit)then
c        idx=0
c        jdx=0
c        kdx=ntoob-orb_info%caborb
c        do isym=1,nsym
c          do j=1,orb_info%ntoobs(isym)-orb_info%cab_orb(isym)
c            idx=idx+1
c            orb_info%xreosym(idx)=jdx+j
c          enddo
c          jdx=jdx+orb_info%ntoobs(isym)-orb_info%cab_orb(isym)
c          do j=1,orb_info%cab_orb(isym)
c            idx=idx+1
c            orb_info%xreosym(idx)=kdx+j
c          enddo
c          kdx=kdx+orb_info%cab_orb(isym)
c        enddo  
c      ! Produce an array which will place all orbitals in type-order.
c        idx=0
c        koffs(1:nsym)=0
c        do igas=1,ngas-1
c          jdx=0
c          do isym=1,nsym
c            do j=1,orb_info%igassh(isym,igas)
c              idx=idx+1
c              orb_info%xreotyp(idx)=j+jdx+koffs(isym)
c            enddo
c            jdx=jdx+orb_info%ntoobs(isym)-orb_info%cab_orb(isym)
c            koffs(isym)=koffs(isym)+orb_info%igassh(isym,igas)
c          enddo
c        enddo
c        do j=1,orb_info%caborb
c          idx=idx+1
c          orb_info%xreotyp(idx)=j+ntoob-caborb
c        enddo          
c      endif  

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
        write(luout,*) 'igamorb:'
        call iwrtma(orb_info%igamorb,1,ntoob,1,ntoob)
        write(luout,*) 'igasorb:'
        call iwrtma(orb_info%igasorb,1,ntoob,1,ntoob)
        write(luout,*) 'mostnd:'
        do igas = 1, ngas
          write(luout,'(2x,i4,2x,8(x,2i4))')
     &         igas,orb_info%mostnd(1:2,1:nsym,igas)
        end do
      end if

      ! generate symmetry ordering -> type ordering mapping
      if(explicit)then
        loop=2
      else
        loop=1
      endif
      do iloop=1,loop
        jdx = 0
        do isym = 1, nsym
          do igas = 1, ngas
            if(iloop.eq.1.and.igas.eq.ngas)cycle
            if (orb_info%ihpvgas(igas).eq.ihole.and.hole_rv) then
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
              if(iloop.eq.1)then
                orb_info%ireost(jdx) = idx
              else
                orb_info%xreosym(jdx)=idx
              endif  
            end do
          end do
        end do

        ! generate reverse mapping
        if(iloop.eq.1)then
          do idx = 1, orb_info%ntoob-caborb
            orb_info%ireots(orb_info%ireost(idx)) = idx 
          end do
        else
          do idx = 1, orb_info%ntoob
            orb_info%xreotyp(orb_info%xreosym(idx)) = idx 
          end do
        endif  
      enddo  

      if (iprint.ge.100) then
        write(luout,*) 'ireost:'
        call iwrtma(orb_info%ireost,1,ntoob-caborb,1,ntoob-caborb)
        write(luout,*) 'ireots:'
        call iwrtma(orb_info%ireots,1,ntoob-caborb,1,ntoob-caborb)
        if(explicit)then
          write(luout,*)'xreosym:'
          call iwrtma(orb_info%xreosym,1,ntoob,1,ntoob)
          write(luout,*)'xreotyp:'
          call iwrtma(orb_info%xreotyp,1,ntoob,1,ntoob)
        endif  
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

c      write(luout,*)'End of set_orbinf'
c      stop

      return
      end
