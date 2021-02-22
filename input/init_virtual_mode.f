      subroutine init_virtual_mode(orb_info)

      implicit none

      include 'ioparam.h'
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'ifc_baserout.h'
      include 'par_dalton.h'
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 00
      character(len=18), parameter ::
     &     i_am = 'init_virtual_mode'

      type(orbinf) ::
     &     orb_info

      type(filinf) ::
     &     ffintf
      integer ::
     &     luintf, iprint, ngas, idxgas, idelim

      integer, parameter ::
     &     mxsym = 8
      character(len=256) :: line, intfile

      real(8) ::
     &     ecore
      integer ::
     &     nirr, nel, isym, mult, mem_ext, istate, refstate,
     &     norbs(8), nocc(8), ncore(8),
     &     nclosed(8)
      integer ::
     &     ninact(8), nact(8), nvirt(8)
      logical ::
     &     error, rd_intfile, rd_nirr, rd_nel, rd_sym, rd_mult,
     &     rd_norbs, rd_nocc, rd_ncore, rd_closed, !rd_irrtyp
     &     rd_ecore, rd_istate, rd_refstate


      if (ntest.ge.10) write(lulog,*) 'entered '//i_am

      iprint = max(iprlvl,ntest)

      call file_init(ffintf,'virtual_molecule.dat',ftyp_sq_frm,0)
      call file_open(ffintf)
      luintf = ffintf%unit

      ! initialize to zero
      ! (saveguard for ill-defined nirr)
      norbs(1:8)=0; nocc(1:8)=0; ncore(1:8)=0; nclosed(1:8)=0;
      ninact(1:8)=0; nact(1:8)=0; nvirt(1:8)=0

      ! set defaults
      mem_ext = -1 ! use defined memory of GeCCo
      rd_intfile = .false.
      rd_nirr = .false.
      rd_nel = .false.
      rd_sym = .false.
      rd_istate = .false.
      rd_refstate = .false.
      rd_mult = .false.
      rd_norbs = .false.
      rd_nocc = .false.
      rd_ncore = .false.
      rd_closed = .false.
      rd_ecore = .false.

      ! loop over file and read info

      do
        read(luintf,'(a)',end=99) line

        idelim = min(len_trim(line)+1,index(line,' '))-1
        select case(line(1:idelim))
        case('nirrep')
          read (line(idelim+1:),*) nirr
          rd_nirr = .true.
        case('nelec')
          read (line(idelim+1:),*) nel
          rd_nel = .true.
        case('refsym')
          read (line(idelim+1:),*) isym
          rd_sym = .true.
        case('refmult')
          read (line(idelim+1:),*) mult
          rd_mult = .true.
        case('norbs')
          if (rd_nirr) then
            read (line(idelim+1:),*) norbs(1:nirr)
            rd_norbs = .true.
          end if
        case('nocc')
          if (rd_nirr) then
            read (line(idelim+1:),*) nocc(1:nirr)
            rd_nocc = .true.
          end if
        case('ncore')
          if (rd_nirr) then
            read (line(idelim+1:),*) ncore(1:nirr)
            rd_ncore = .true.
          end if
        case('nclosed')
          if (rd_nirr) then
            read (line(idelim+1:),*) nclosed(1:nirr)
            rd_closed = .true.
          end if
        end select

      end do

 99   continue

      error = .false.
      if (.not.rd_nirr) then
        write(lulog,*) 'Could not read number of IRREPs'
        error = .true.
      end if
      if (.not.rd_nel) then
        write(lulog,*) 'Could not read number of electrons'
        error = .true.
      end if
      if (.not.rd_sym) then
        write(lulog,*)
     &       'Could not read reference state IRREP, assuming 1'
        isym = 1
      end if
      if (.not.rd_mult) then
        write(lulog,*) 'Could not read reference state multiplicity'
        error = .true.
      end if
      if (.not.rd_norbs) then
        write(lulog,*) 'Could not read number of orbitals/IRREP'
        error = .true.
      end if
      if (.not.rd_ncore) then
        write(lulog,*)
     &       'Could not read number of core orbitals/IRREP, assuming 0'
        if (rd_nirr) ncore(1:nirr) = 0
      end if
      if (.not.rd_closed) then
        write(lulog,*) 'Could not read number of closed orbitals/IRREP'
        error = .true.
      end if
      if (.not.rd_nocc) then
        write(lulog,*)'Could not read number of occupied orbitals/IRREP'
        error = .true.
      end if

      if (error) call quit(1,i_am,'Error(s) reading interface file!')

      ! get inactive, active and secondary shells
      ninact(1:nirr) = nclosed(1:nirr)-ncore(1:nirr)
      nact(1:nirr)   = nocc(1:nirr)-nclosed(1:nirr)
      nvirt(1:nirr)  = norbs(1:nirr)-nocc(1:nirr)

      if (nirr.gt.2.and.nirr.ne.4.and.nirr.ne.8) then
        call warn(i_am,'Defined nirrep is strange, trying to '
     &                  //'recover')
        if (nirr.eq.3) nirr=4
        if (nirr.gt.4) nirr=8
      end if

      !if (iprint.ge.1) then
        write(luout,'(1x,"GeCCo will use these orbital spaces:")')
        write(luout,'(1x," core      ",8i4)') ncore(1:nirr)
        write(luout,'(1x," inactive  ",8i4)') ninact(1:nirr)
        write(luout,'(1x," active    ",8i4)') nact(1:nirr)
        write(luout,'(1x," virtual   ",8i4)') nvirt(1:nirr)
      !end if

      ngas=1
      if (sum(ncore(1:nirr)).gt.0)  ngas = ngas+1
      if (sum(ninact(1:nirr)).gt.0) ngas = ngas+1
      if (sum(nact(1:nirr)).gt.0)   ngas = ngas+1


      ! set orb_info
      orb_info%nsym = nirr
      orb_info%ngas = ngas
      orb_info%nspin = 1    ! no UHF

      orb_info%nbast = sum(norbs(1:nirr))
      orb_info%ntoob = sum(norbs(1:nirr))
      orb_info%lsym  = isym
      orb_info%imult = mult
      orb_info%ims   = 1 - mod(mult,2) ! MS = 0 or 1/2

      orb_info%nxbast = 0 ! no CABS stuff so far
      orb_info%caborb = 0

      orb_info%nactel = nel - 2*sum(nclosed(1:nirr))
      orb_info%nactorb = sum(nact(1:nirr))

      orb_info%ncore_ext = sum(ncore(1:nirr))
      orb_info%name_intfile_ext = intfile

      allocate(orb_info%igassh(nirr,ngas),orb_info%iad_gas(ngas),
     &         orb_info%ihpvgas(ngas,1),orb_info%nbas(nirr),
     &         orb_info%ntoobs(nirr))
      allocate(orb_info%cab_orb(nirr),
     &              orb_info%nxbas(nirr) )

      idxgas = 0
      if (sum(ncore(1:nirr)).gt.0) then
        idxgas = idxgas+1
        orb_info%iad_gas(idxgas)=1
        orb_info%ihpvgas(idxgas,1)=IHOLE
        orb_info%igassh(1:nirr,idxgas) = ncore(1:nirr)
      end if
      if (sum(ninact(1:nirr)).gt.0) then
        idxgas = idxgas+1
        orb_info%iad_gas(idxgas)=2
        orb_info%ihpvgas(idxgas,1)=IHOLE
        orb_info%igassh(1:nirr,idxgas) = ninact(1:nirr)
      end if
      if (sum(nact(1:nirr)).gt.0) then
        idxgas = idxgas+1
        orb_info%iad_gas(idxgas)=2
        orb_info%ihpvgas(idxgas,1)=IVALE
        orb_info%igassh(1:nirr,idxgas) = nact(1:nirr)
      end if
      idxgas = idxgas+1
      orb_info%iad_gas(idxgas)=2
      orb_info%ihpvgas(idxgas,1)=IPART
      orb_info%igassh(1:nirr,idxgas) = nvirt(1:nirr)

      orb_info%nbas(1:nirr) = norbs(1:nirr)
      orb_info%ntoobs(1:nirr) = norbs(1:nirr)
      orb_info%cab_orb(1:nirr)=0
      orb_info%nxbas(1:nirr)=0

      orb_info%mem_ext = mem_ext

      if (ntest.ge.100) then
        write(lulog,*) 'Have set:'
        write(lulog,*) 'igassh:'
        do idxgas = 1, ngas
          write(lulog,'(1x,8i6)') orb_info%igassh(1:nirr,idxgas)
        end do
      end if

      call file_close_keep(ffintf)

      return
      end
