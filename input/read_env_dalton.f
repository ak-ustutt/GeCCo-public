      subroutine read_env_dalton(orb_info)

      implicit none
      
      include 'ioparam.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'ifc_baserout.h'
      include 'par_dalton.h'

      integer, parameter ::
     &     ntest = 100

      type(orbinf) ::
     &     orb_info

      type(filinf) ::
     &     ffsir
      integer ::
     &     lusir, iprint, ngas, loop, i, j, caborb
      logical ::
     &     have_frozen, have_act, have_as1, have_as2, have_as3, logaux

      integer, parameter ::
     &     mxsym = 8
      ! DALTON writes integer*4, so we must take care of that
      integer*4 ::
     &     istate,ispin,nactel,lsym,nsym,
     &     nisht,nasht,nocct,norbt,nbast,nconf,nwopt,nwoph,
     &     ncdets, ncmot,nnashx,nnashy,nnorbt,n2orbt,
     &     nelmn1, nelmx1, nelmn3, nelmx3, mctype,
     &     nish(mxsym),nash(mxsym),norb(mxsym),nbas(mxsym),
     &     nfro(mxsym),nrhf(mxsym),
     &     muld2h(mxsym,mxsym), nas1(mxsym), nas2(mxsym), nas3(mxsym),
     &     auxbas(mxsym),linind(mxsym),totbas(mxsym)
      real(8) ::
     &     potnuc,emy,eactiv,emcscf

      iprint = max(iprlvl,ntest)

      inquire(file='AUXBAS',exist=logaux)

      call file_init(ffsir,sirifc,ftyp_sq_unf,0)
      call file_open(ffsir)

      lusir = ffsir%unit
      rewind lusir
      call mollab('SIR IPH ',lusir,luout)

      read (lusir) potnuc,emy,eactiv,emcscf,istate,ispin,nactel,lsym
      read (lusir) nisht,nasht,nocct,norbt,nbast,nconf,nwopt,nwoph,
     &     ncdets, ncmot,nnashx,nnashy,nnorbt,n2orbt,
     &     nsym,muld2h(1:mxsym,1:mxsym),nrhf(1:mxsym),nfro(1:mxsym),
     &     nish(1:mxsym),nash(1:mxsym),norb(1:mxsym),nbas(1:mxsym),
     &     nelmn1, nelmx1, nelmn3, nelmx3, mctype,
     &     nas1(1:mxsym), nas2(1:mxsym), nas3(1:mxsym)

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
          if(mod(totbas(i),4).ne.0)loop=loop+linind(i)
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

        write(luout,'(x,a,2i4)') 'nnorbt,n2orbt: ',nnorbt,n2orbt
      end if

      call file_close_keep(ffsir)

      have_frozen = i4elsqsum(nfro,nsym).gt.0
      have_act    = i4elsqsum(nash,nsym).gt.0
      have_as1    = i4elsqsum(nas1,nsym).gt.0
      have_as2    = i4elsqsum(nas1,nsym).gt.0
      have_as3    = i4elsqsum(nas1,nsym).gt.0

      ngas = 2
      if (have_frozen) ngas = ngas+1
      if(logaux) ngas=ngas+1
      if (have_act.or.have_as1.or.have_as2.or.have_as3)
     &     call quit(1,'read_env_dalton',
     &     'not adapted to active orbitals')

      orb_info%nsym = nsym
      orb_info%ngas = ngas
      orb_info%ntoob = norbt
      orb_info%nbast = nbast
      orb_info%caborb = caborb

      ! use the data to initialize orb_info structure            
      allocate(orb_info%nbas(nsym),
     &     orb_info%ntoobs(nsym),orb_info%igassh(nsym,ngas),
     &     orb_info%iad_gas(ngas),orb_info%ihpvgas(ngas))
      if(logaux)allocate(orb_info%cab_orb(nsym))

      orb_info%nbas(1:nsym) = nbas(1:nsym)
      orb_info%ntoobs(1:nsym) = norb(1:nsym)
      if (have_frozen) then
        if(logaux)then
          orb_info%iad_gas(1:ngas) = (/1,2,2,2/)
          orb_info%ihpvgas(1:ngas) = (/1,1,2,4/)
        else
          orb_info%iad_gas(1:ngas) = (/1,2,2/)
          orb_info%ihpvgas(1:ngas) = (/1,1,2/)
        endif  
        orb_info%igassh(1:nsym,1) = nfro(1:nsym)
        orb_info%igassh(1:nsym,2) = nrhf(1:nsym)-nfro(1:nsym)
        orb_info%igassh(1:nsym,3) = norb(1:nsym)-nrhf(1:nsym)
        if(logaux)then
          orb_info%igassh(1:nsym,4) = linind(1:nsym)
        endif  
      else
        if(logaux)then
          orb_info%iad_gas(1:ngas) = (/2,2,2/)
          orb_info%ihpvgas(1:ngas) = (/1,2,4/)
        else  
          orb_info%iad_gas(1:ngas) = (/2,2/)
          orb_info%ihpvgas(1:ngas) = (/1,2/)
        endif
        orb_info%igassh(1:nsym,1) = nrhf(1:nsym)
        orb_info%igassh(1:nsym,2) = norb(1:nsym)-nrhf(1:nsym)
        if(logaux)then
          orb_info%igassh(1:nsym,3) = linind(1:nsym)
        endif  
      end if
      if(logaux)then
        orb_info%cab_orb(1:nsym)=linind(1:nsym)
      endif 

      return
      end
