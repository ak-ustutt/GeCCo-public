      subroutine read_env_dalton64(orb_info)
      ! the same as read_env_dalton, but for DALTON files with integer(8)
      implicit none
      
      include 'ioparam.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'ifc_baserout.h'
      include 'par_dalton.h'
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 100
      real(8), parameter ::
     &     thr_auto_freeze = -2.5d0

      type(orbinf) ::
     &     orb_info

      type(filinf) ::
     &     ffsir
      integer ::
     &     lusir, luerr, iprint, ngas, loop, i, j, caborb, idx, jdx,
     &     n_frozen, n_act, n_as1, n_as2, n_as3, nspin
      logical ::
     &     logaux

      integer, parameter ::
     &     mxsym = 8
      integer(8), parameter ::
     &     four8 = 4
      integer(8) ::
     &     istate,ispin,nactel,lsym,nsym,
     &     nisht,nasht,nocct,norbt,nbast,nxbast,nconf,nwopt,nwoph,
     &     ncdets, ncmot,nnashx,nnashy,nnorbt,n2orbt,
     &     nelmn1, nelmx1, nelmn3, nelmx3, mctype,
     &     nish(mxsym),nash(mxsym),norb(mxsym),nbas(mxsym),
     &     nfro(mxsym),nrhf(mxsym),
     &     muld2h(mxsym,mxsym), nas1(mxsym), nas2(mxsym), nas3(mxsym),
     &     auxbas(mxsym),linind(mxsym),totbas(mxsym)
      integer ::
     &     inactorb(mxsym),cas(mxsym),actel
      real(8) ::
     &     potnuc,emy,eactiv,emcscf

      real(8), pointer ::
     &     orb_en(:), orb_en_sort(:)
      integer(8), pointer ::
     &     orb_sym_i8(:)

      iprint = max(iprlvl,ntest)

      inquire(file='AUXBAS',exist=logaux)

      call file_init(ffsir,sirifc,ftyp_sq_unf,0)
      call file_open(ffsir)

      lusir = ffsir%unit
      luerr = luout
      rewind lusir
      call mollab('SIR IPH ',lusir,luerr)

      read (lusir) potnuc,emy,eactiv,emcscf,istate,ispin,nactel,lsym
      read (lusir) nisht,nasht,nocct,norbt,nbast,nconf,nwopt,nwoph,
     &     ncdets, ncmot,nnashx,nnashy,nnorbt,n2orbt,
     &     nsym,muld2h(1:mxsym,1:mxsym),nrhf(1:mxsym),nfro(1:mxsym),
     &     nish(1:mxsym),nash(1:mxsym),norb(1:mxsym),nbas(1:mxsym),
     &     nelmn1, nelmx1, nelmn3, nelmx3, mctype,
     &     nas1(1:mxsym), nas2(1:mxsym), nas3(1:mxsym)

      ! read orbital energies from SIRIFC and set up
      ! proposal for frozen core
      allocate(orb_en(norbt),orb_sym_i8(norbt))

      call mollab('TRCCINT ',lusir,luerr)
      read (lusir)
      read (lusir) orb_en(1:norbt), orb_sym_i8(1:norbt)
      
      ! scan through orbital energies and find bound orbital
      orb_info%n_bound_orbs = 0
      do idx = 1, norbt
        if (orb_en(idx).lt.0d0)
     &       orb_info%n_bound_orbs = orb_info%n_bound_orbs+1 
        
        if (orb_en(idx).lt.thr_auto_freeze)
     &       orb_info%n_freeze_rcmd = orb_info%n_freeze_rcmd+1

      end do
      allocate(orb_info%isym_bound_orbs(orb_info%n_bound_orbs),
     &     orb_en_sort(orb_info%n_bound_orbs))
      
      jdx = 0
      do idx = 1, norbt
        if (orb_en(idx).ge.0d0) cycle
        jdx = jdx+1
        orb_en_sort(jdx) = orb_en(idx)
        orb_info%isym_bound_orbs(jdx) =
     &       orb_sym_i8(idx)
      end do

      call idxsortx(orb_en_sort,orb_info%isym_bound_orbs,
     &                          orb_info%n_bound_orbs,+1)

      deallocate(orb_en,orb_en_sort,orb_sym_i8)

      ! If necessary, deal with auxiliary basis functions.
      caborb=0 
      linind(1:nsym)=0
      if(logaux)then
        open(file='AUXBAS',unit=999,status='old',form='formatted',
     &       access='sequential')
        rewind(999)
        read(999,*)
        do i=1,nsym
          read(999,'(5x,3i5)')nbas(i),linind(i),totbas(i)
          auxbas(i)=totbas(i)-nbas(i)
          loop=linind(i)*(2+totbas(i)/4)
          if(mod(totbas(i),four8).ne.0)loop=loop+linind(i)
          do j=1,loop
            read(999,*)
          enddo
          caborb=caborb+linind(i)
        enddo
        close(unit=999,status='keep')
      endif  

      if (iprint.ge.50) then
        write(luout,*) 'raw data from section "SIR IPH "'
        write(luout,*) 'potnuc = ',potnuc
        write(luout,*) 'emy    = ',emy
        write(luout,*) 'eactiv = ',eactiv
        write(luout,*) 'emcscf = ',emcscf
        write(luout,*) 'istate = ',istate
        write(luout,*) 'ispin  = ',ispin
        write(luout,*) 'nactel = ',nactel
        write(luout,*) 'lsym   = ',lsym
        write(luout,*) 'nsym   = ',nsym        
        write(luout,*) 'mctype = ',mctype
        
        write(luout,'(x,a,8i4)') 'nrhf   = ',nrhf(1:8)
        write(luout,'(x,a,8i4)') 'nfro   = ',nfro(1:8)
        write(luout,'(x,a,8i4)') 'nish   = ',nish(1:8)
        write(luout,'(x,a,8i4)') 'nash   = ',nash(1:8)
        write(luout,'(x,a,8i4)') 'norb   = ',norb(1:8)
        if(logaux)then
          write(luout,'(x,a,8i4)') 'cabs   = ',linind(1:8)
        endif           
        write(luout,'(x,a,8i4)') 'nbas   = ',nbas(1:8)
        if(logaux)then
          write(luout,'(x,a,8i4)') 'abas   = ',auxbas(1:8)
          write(luout,'(x,a,8i4)') 'tbas   = ',totbas(1:8)
        endif  
        write(luout,'(x,a,8i4)') 'nas1   = ',nas1(1:8)
        write(luout,'(x,a,8i4)') 'nas2   = ',nas2(1:8)
        write(luout,'(x,a,8i4)') 'nas3   = ',nas3(1:8)

        write(luout,'(x,a,2i8)') 'nnorbt,n2orbt: ',nnorbt,n2orbt

        write(luout,'(x,a)') 'sym_bound_orbs:'
        write(luout,'(x,10i4)') orb_info%isym_bound_orbs
        write(luout,'(x,a,i4)')  'n_freeze_rcmd: ',
     &       orb_info%n_freeze_rcmd
      end if

      call file_close_keep(ffsir)

c     We will just ignore if orbitals were frozen within Dalton run
c      n_frozen = i4elsum(nfro,nsym)
      n_act    = ielsum(nash,nsym)
      n_as1    = ielsum(nas1,nsym)
      n_as2    = ielsum(nas2,nsym)
      n_as3    = ielsum(nas3,nsym)

      nspin = 1
      ngas = 2
c      if (n_frozen.gt.0) ngas = ngas+1
      if(logaux) ngas=ngas+1
      orb_info%nactel = nactel
      orb_info%nactorb = n_act
      orb_info%lsym = lsym
      orb_info%imult = ispin
      orb_info%ims = 1 - mod(ispin,2) ! assume low-spin case as default
      if (n_act.gt.0) then
        ! test whether this can be treated as a simple
        ! high spin open shell case:
        if (.false..and.nactel.eq.n_act.and.ispin.eq.nactel+1) then
          ! we should check the symmetry here ...
          write(luout,*) 'high-spin valence shell detected'
          ngas = ngas+1
          nspin = 2
        else
          write(luout,*) 'valence shell is not high-spin ...'
c let pass to enable CAS calculations
          ngas = ngas+1
c          call quit(1,'read_env_dalton',
c     &     'not adapted to non-high-spin open shell cases')
        end if
      end if
c      if (n_as1.gt.0.or.n_as2.gt.0.or.n_as3.gt.0)
c     &     call quit(1,'read_env_dalton',
c     &     'not adapted to RAS orbital spaces')

      nbast  = sum(nbas(1:nsym))
      nxbast = sum(auxbas(1:nsym))

      orb_info%nspin = nspin
      orb_info%nsym = nsym
      orb_info%ngas = ngas
      orb_info%ntoob = norbt
      orb_info%nbast = nbast
      orb_info%nxbast = nxbast
      orb_info%caborb = caborb

      ! use the data to initialize orb_info structure            
      allocate(orb_info%nbas(nsym),
     &     orb_info%ntoobs(nsym),orb_info%igassh(nsym,ngas),
     &     orb_info%iad_gas(ngas),orb_info%ihpvgas(ngas,nspin))
      allocate(orb_info%cab_orb(nsym),
     &              orb_info%nxbas(nsym) )
      

      orb_info%nbas(1:nsym) = nbas(1:nsym)
      orb_info%ntoobs(1:nsym) = norb(1:nsym)
      if (n_frozen.gt.0)
     &   write(luout,*) 'We will just ignore DALTON''s frozen orbitals'
c      if (n_frozen.gt.0) then
c        ! not sure whether DALTON's frozen orbitals are what we want
c        call quit(1,'read_env_dalton',
c     &       'check first what DALTON means by freezing orbitals!')
c        if(logaux)then
c          orb_info%iad_gas(1:ngas) = (/1,2,2,2/)
c          orb_info%ihpvgas(1:ngas,1) = (/1,1,2,4/)
c        else
c          orb_info%iad_gas(1:ngas) = (/1,2,2/)
c          orb_info%ihpvgas(1:ngas,1) = (/1,1,2/)
c        endif  
c        orb_info%igassh(1:nsym,1) = nfro(1:nsym)
c        orb_info%igassh(1:nsym,2) = nish(1:nsym)-nfro(1:nsym)
c        orb_info%igassh(1:nsym,3) = norb(1:nsym)-nish(1:nsym)
c        if(logaux)then
c          orb_info%igassh(1:nsym,4) = linind(1:nsym)
c        endif  
c      else
        if (nspin.eq.1.and.n_act.eq.0) then
          if(logaux)then
            orb_info%iad_gas(1:ngas) = (/2,2,2/)
            orb_info%ihpvgas(1:ngas,1) = (/1,2,4/)
          else  
            orb_info%iad_gas(1:ngas) = (/2,2/)
            orb_info%ihpvgas(1:ngas,1) = (/1,2/)
          endif
          orb_info%igassh(1:nsym,1) = nish(1:nsym)
          orb_info%igassh(1:nsym,2) = norb(1:nsym)-nish(1:nsym)
          if(logaux)then
            orb_info%igassh(1:nsym,3) = linind(1:nsym)
          endif
        else if (nspin.eq.1) then
          if(logaux)then
            orb_info%iad_gas(1:ngas) = (/2,2,2,2/)
            orb_info%ihpvgas(1:ngas,1) = (/1,3,2,4/)
          else
            orb_info%iad_gas(1:ngas) = (/2,2,2/)
            orb_info%ihpvgas(1:ngas,1) = (/1,3,2/)
          endif
          orb_info%igassh(1:nsym,1) = nish(1:nsym)
          orb_info%igassh(1:nsym,2) = nash(1:nsym)
          orb_info%igassh(1:nsym,3) 
     &                  = norb(1:nsym)-nish(1:nsym)-nash(1:nsym)
          if(logaux)then
            orb_info%igassh(1:nsym,4) = linind(1:nsym)
          endif
        else
          if(logaux)then
            orb_info%iad_gas(1:ngas) = (/2,2,2/)
            orb_info%ihpvgas(1:ngas,1) = (/1,1,2,4/)
            orb_info%ihpvgas(1:ngas,2) = (/1,2,2,4/)
          else  
            orb_info%iad_gas(1:ngas) = (/2,2/)
            orb_info%ihpvgas(1:ngas,1) = (/1,1,2/)
            orb_info%ihpvgas(1:ngas,2) = (/1,2,2/)
          endif
          orb_info%igassh(1:nsym,1) = nish(1:nsym)
          orb_info%igassh(1:nsym,2) = nash(1:nsym)
          orb_info%igassh(1:nsym,3) =
     &         norb(1:nsym)-nish(1:nsym)-nash(1:nsym)
          if(logaux)then
            orb_info%igassh(1:nsym,4) = linind(1:nsym)
          endif
        end if
c      end if
      if(logaux)then
        orb_info%cab_orb(1:nsym)=linind(1:nsym)
        orb_info%nxbas(1:nsym)=auxbas(1:nsym)
      else
        orb_info%cab_orb(1:nsym)=0
        orb_info%nxbas(1:nsym)=0
      endif 

      return
      end
