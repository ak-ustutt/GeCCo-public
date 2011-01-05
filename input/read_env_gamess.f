      subroutine read_env_gamess(orb_info)

      implicit none
      
      include 'ioparam.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'ifc_baserout.h'
      include 'par_gamess.h'
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 100
      real(8), parameter ::
     &     thr_auto_freeze = -2.5d0

      type(orbinf) ::
     &     orb_info

      type(filinf) ::
     &     ffdict

      integer, parameter ::
     &     mxsym = 8
      integer ::
     &     ludict, iprint, ngas, idx, jdx, nspin, nsym, isym,
     &     mxirrep, lsym, nactel, ispin, nstate, igroup,
     &     nfro(mxsym), nrhf(mxsym), nash(mxsym), norb(mxsym)
      logical ::
     &     isgroup(8)
      character(3) ::
     &     groups(8), irreps(8)

      ! 64 bit version of GAMESS writes integer*8
      integer(8) ::
     &     irecst,ioda(len_da),ifilen(len_da),is,ipk,
     &     nfroz,nocct,n_act,norbs,norbt,num1
      real(8) ::
     &     potnuc,emcscf,sz,s2,etot2,eerd,e1,e2,ven,vee,epot,
     &     eelct,ekin,estate(mxrt),statn,ecore

      real(8), pointer ::
     &     orb_en(:)
      integer(8), pointer ::
     &     orb_sym(:)
      character(len=8), pointer ::
     &     orb_lab(:)

      groups(1:8) = (/'C1','Ci','Cs','C2','C2v','C2h','D2','D2h'/)

      iprint = max(iprlvl,ntest)

      ! open direct access file DICTNRY
      call file_init(ffdict,dictnry,ftyp_da_unf,irecln)
      call file_open(ffdict)

      ludict = ffdict%unit

      ! read first record: information about where to find which record
      read (ludict,rec=1) irecst,ioda,ifilen,is,ipk
c dbg
c      print *,'first record from DICTNRY. irecst,is,ipk:',irecst,is,ipk
c      print *,'          ii        ioda        nrec      ifilen'
c      print *,'------------------------------------------------'
c      do idx=1,len_da
c       if (ioda(idx).ge.0)
c     &     write(luout,'(x,4i12)') idx, ioda(idx),
c     &                             ifilen(idx)/irecln+1, ifilen(idx)
c      end do
c dbgend

      ! read energy quantities
      read (ludict,rec=ioda(6)) potnuc,eelct,etot2,sz,s2,ecore,emcscf,
     &                          eerd,e1,e2,ven,vee,epot,ekin,
     &                          estate,statn
      ispin = 2*nint(s2)+1
      if (nint(sz).ne.0) then
        print *,'sz = ',sz
        call warn('read_env_gamess',
     &           'can we deal with Sz.ne.0?')
      end if
      nstate = nint(statn)

      ! read number of orbital parameters
      read (ludict,rec=ioda(5)) nfroz,nocct,n_act,norbs,norbt,num1
      ! unpack number of electrons and irrep
      lsym = mod(num1,1000)
      nactel = num1/1000 - 2*nocct !what if nfroz>0 ???

      ! read orbital energies and set up proposal for frozen core
      allocate(orb_en(norbt),orb_sym(norbt),orb_lab(norbt))
      read (ludict,rec=ioda(17)) orb_en(1:norbt)
      read (ludict,rec=ioda(262)) orb_sym(1:norbt)
      read (ludict,rec=ioda(324)) orb_lab(1:norbt)
c      ! now get nsym (lowest possible 2^n number):
c      mxirrep = maxval(orb_sym(1:norbt))
c      nsym = mxsym
c      do while (nsym/2.ge.mxirrep)
c        nsym = nsym/2
c      end do

      ! For a strange (?) reason, GAMESS may mess up the irrep numbers
      ! for the orbitals. So let's restore them by help of the labels.
      ! First try to recognize point group
      isgroup(1:8) = .true. ! stands for: C1,Ci,Cs,C2,C2v,C2h,D2,D2h
      do idx = 1, norbt
        select case(trim(orb_lab(idx)))
        case('A') ! not Ci,Cs,C2v,C2h,D2h
          isgroup(2:3) = .false.
          isgroup(5:6) = .false.
          isgroup(8) = .false.
        case('AG','AU') ! not C1,Cs,C2,C2v,D2
          isgroup(1) = .false.
          isgroup(3:5) = .false.
          isgroup(7) = .false.
        case('A''','A''''') ! only Cs
          isgroup(1:2) = .false.
          isgroup(4:8) = .false.
        case('B') ! only C2
          isgroup(1:3) = .false.
          isgroup(5:8) = .false.
        case('A1','A2') ! only C2v
          isgroup(1:4) = .false.
          isgroup(6:8) = .false.
        case('B1','B2') ! only C2v or D2
          isgroup(1:4) = .false.
          isgroup(6) = .false.
          isgroup(8) = .false.
        case('BU','BG') ! only C2h
          isgroup(1:5) = .false.
          isgroup(7:8) = .false.
        case('B3') ! only D2
          isgroup(1:6) = .false.
          isgroup(8) = .false.
        case('B1G','B2G','B3G','B1U','B2U','B3U') ! only D2h
          isgroup(1:7) = .false.
        case default
          call quit(1,'read_env_gamess',
     &              'unknown irrep label: '//trim(orb_lab(idx)))
        end select
      end do
      isym = 0
      do idx = 1, 8
        if (isgroup(idx)) then
          isym = isym + 1
          igroup = idx
        end if
      end do
      if (isym.ne.1) call quit(1,'read_env_gamess',
     &       'Could not unambiguously identify point group!')

      ! irrep numbers as in GAMESS manual (for $DET):
      irreps(1:8) = '   '
      select case(igroup)
      case(1)
        nsym = 1
        irreps = (/'A  ','   ','   ','   ','   ','   ','   ','   '/)
      case(2)
        nsym = 2
        irreps = (/'AG ','AU ','   ','   ','   ','   ','   ','   '/)
      case(3)
        nsym = 2
        irreps = (/'A''','A''''','   ','   ','   ','   ','   ','   '/)
      case(4)
        nsym = 2
        irreps = (/'A  ','B  ','   ','   ','   ','   ','   ','   '/)
      case(5)
        nsym = 4
        irreps = (/'A1 ','A2 ','B1 ','B2 ','   ','   ','   ','   '/)
      case(6)
        nsym = 4
        irreps = (/'AG ','BU ','BG ','AU ','   ','   ','   ','   '/)
      case(7)
        nsym = 4
        irreps = (/'A  ','B1 ','B2 ','B3 ','   ','   ','   ','   '/)
      case(8)
        nsym = 8
        irreps = (/'AG ','B1G','B2G','B3G','AU ','B1U','B2U','B3U'/)
      case default
        call quit(1,'read_env_gamess','impossible.')
      end select

      do idx = 1, norbt
        do isym = 1, nsym
          if (trim(orb_lab(idx)).eq.trim(irreps(isym))) then
            orb_sym(idx) = isym
            exit
            if (isym.eq.nsym) call quit(1,'read_env_gamess',
     &               'There must be an error in this routine')
          end if
        end do
      end do
      deallocate(orb_lab)
c dbg
c      print *,'   orbital  symmetry    energy'
c      print *,'------------------------------'
c      do idx = 1, norbt
c        write(luout,'(x,2i10,f10.6)') idx,orb_sym(idx),orb_en(idx)
c      end do
c dbgend

      ! though irrep indices may differ from other codes
      ! (e.g. for C2v irrep order is A1 A2 B1 B2)
      ! the multiplication table should be the same as in multd2h.h.

      nfro = 0
      nrhf = 0
      nash = 0
      norb = 0
      do idx = 1, norbt
        isym = orb_sym(idx)
        norb(isym) = norb(isym) + 1
        if (idx.le.nfroz) nfro(isym) = nfro(isym) + 1
        if (idx.le.nocct) nrhf(isym) = nrhf(isym) + 1
        if (idx.gt.nocct.and.idx.le.nocct+n_act)
     &             nash(isym) = nash(isym) + 1
      end do

      ! GAMESS sorts orbitals according to energy so it's easy
      ! to set up irrep array for orbitals.
      ! we set it up for all orbitals since we need this information
      ! in set_orbinf
      ! find bound orbitals:
      orb_info%n_bound_orbs = norbt
      do idx = 1, norbt
        if (orb_en(idx).ge.thr_auto_freeze) exit
        orb_info%n_freeze_rcmd = orb_info%n_freeze_rcmd+1
      end do

      ! just copy from orb_sym:
      allocate(orb_info%isym_bound_orbs(orb_info%n_bound_orbs))
      orb_info%isym_bound_orbs(1:orb_info%n_bound_orbs)
     &      = orb_sym(1:orb_info%n_bound_orbs)

      deallocate(orb_en,orb_sym)

      if (iprint.ge.50) then
        write(luout,*) 'data from DICTNRY:'
        write(luout,*) 'potnuc= ',potnuc
        write(luout,*) 'eelct = ',eelct
        write(luout,*) 'etot2 = ',etot2
        write(luout,*) 'emcscf= ',emcscf
        write(luout,*) 'ecore = ',ecore
        write(luout,*) 'estate= ',estate(1:nstate)
        write(luout,*) 'nstate= ',nstate
        write(luout,*) 'ispin = ',ispin
        write(luout,*) 'nactel= ',nactel
        write(luout,*) 'lsym  = ',lsym
        write(luout,*) 'symm. = ',groups(igroup)
        write(luout,*) 'nsym  = ',nsym
        write(luout,*) 'nfroz = ',nfroz
        write(luout,*) 'nocct = ',nocct
        write(luout,*) 'n_act = ',n_act
        write(luout,*) 'norbs = ',norbs
        write(luout,*) 'norbt = ',norbt
        write(luout,'(x,a,8i4)') 'nfro   = ',nfro(1:8)
        write(luout,'(x,a,8i4)') 'nrhf   = ',nrhf(1:8)
        write(luout,'(x,a,8i4)') 'nash   = ',nash(1:8)
        write(luout,'(x,a,8i4)') 'norb   = ',norb(1:8)
        write(luout,'(x,a)') 'sym_bound_orbs:'
        write(luout,'(x,10i4)') orb_info%isym_bound_orbs
        write(luout,'(x,a,i4)')  'n_freeze_rcmd: ',
     &       orb_info%n_freeze_rcmd
      end if

      call file_close_keep(ffdict)

      nspin = 1
      ngas = 2
      if (nfroz.gt.0) ngas = ngas+1
      orb_info%nactel = nactel
      orb_info%lsym = lsym
      orb_info%imult = ispin
      if (n_act.gt.0) then
        ! test whether this can be treated as a simple
        ! high spin open shell case:
        if (nactel.eq.n_act.and.ispin.eq.nactel+1) then
          ! we should check the symmetry here ...
          write(luout,*) 'high-spin valence shell detected'
          ngas = ngas+1
          nspin = 2
        else
          write(luout,*) 'valence shell is not high-spin ...'
c let pass to enable CAS calculations
          ngas = ngas+1
        end if
      end if

      orb_info%nspin = nspin
      orb_info%nsym = nsym
      orb_info%ngas = ngas
      orb_info%ntoob = norbt
      orb_info%nbast = norbt ! number of BF = number of orbs (correct?)
      orb_info%nxbast = 0
      orb_info%caborb = 0

      ! use the data to initialize orb_info structure            
      allocate(orb_info%nbas(nsym),
     &     orb_info%ntoobs(nsym),orb_info%igassh(nsym,ngas),
     &     orb_info%iad_gas(ngas),orb_info%ihpvgas(ngas,nspin))
      allocate(orb_info%cab_orb(nsym),
     &              orb_info%nxbas(nsym) )
      

      orb_info%nbas(1:nsym) = norb(1:nsym) ! sym of BF = sym of orbs (?)
      orb_info%ntoobs(1:nsym) = norb(1:nsym)
      if (nfroz.gt.0) then
        ! not sure whether GAMESS' frozen orbitals are what we want
        call quit(1,'read_env_gamess',
     &       'check first what GAMESS means by freezing orbitals!')
        orb_info%iad_gas(1:ngas) = (/1,2,2/)
        orb_info%ihpvgas(1:ngas,1) = (/1,1,2/)
        orb_info%igassh(1:nsym,1) = nfro(1:nsym)
        orb_info%igassh(1:nsym,2) = nrhf(1:nsym)-nfro(1:nsym)
        orb_info%igassh(1:nsym,3) = norb(1:nsym)-nrhf(1:nsym)
      else
        if (nspin.eq.1.and.n_act.eq.0) then
          orb_info%iad_gas(1:ngas) = (/2,2/)
          orb_info%ihpvgas(1:ngas,1) = (/1,2/)
          orb_info%igassh(1:nsym,1) = nrhf(1:nsym)
          orb_info%igassh(1:nsym,2) = norb(1:nsym)-nrhf(1:nsym)
        else if (nspin.eq.1) then
          orb_info%iad_gas(1:ngas) = (/2,2,2/)
          orb_info%ihpvgas(1:ngas,1) = (/1,3,2/)
          orb_info%igassh(1:nsym,1) = nrhf(1:nsym)
          orb_info%igassh(1:nsym,2) = nash(1:nsym)
          orb_info%igassh(1:nsym,3) 
     &                  = norb(1:nsym)-nrhf(1:nsym)-nash(1:nsym)
        else
          orb_info%iad_gas(1:ngas) = (/2,2/)
          orb_info%ihpvgas(1:ngas,1) = (/1,1,2/)
          orb_info%ihpvgas(1:ngas,2) = (/1,2,2/)
          orb_info%igassh(1:nsym,1) = nrhf(1:nsym)
          orb_info%igassh(1:nsym,2) = nash(1:nsym)
          orb_info%igassh(1:nsym,3) =
     &         norb(1:nsym)-nrhf(1:nsym)-nash(1:nsym)
        end if
      end if
      orb_info%cab_orb(1:nsym)=0
      orb_info%nxbas(1:nsym)=0

      return
      end
