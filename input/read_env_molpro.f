      subroutine read_env_molpro(orb_info)

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
      character(len=20), parameter ::
     &     i_am = 'read_env_molpro'

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
     &     nirr, nel, isym, mult, mem_mpro,
     &     norbs(8), nocc(8), ncore(8),
     &     nclosed(8)
      integer ::
     &     ninact(8), nact(8), nvirt(8)
      logical ::
     &     error, rd_intfile, rd_nirr, rd_nel, rd_sym, rd_mult,
     &     rd_norbs, rd_nocc, rd_ncore, rd_closed, !rd_irrtyp
     &     rd_ecore


      iprint = max(iprlvl,ntest)

      call file_init(ffintf,'mpro_gecco_ifc.dat',ftyp_sq_frm,0)
      call file_open(ffintf)
      luintf = ffintf%unit

      read(luintf,'(a)') line

      if (line(16:19).ne.'v1.0') then
        write(lulog,*) 'Is not v1.0 type file:'
        write(lulog,*) trim(line)
        call quit(0,i_am,'wrong format of interface file?')
      end if

      ! set defaults
      mem_mpro = -1 ! use defined memory of GeCCo
      rd_intfile = .false.
      rd_nirr = .false.
      rd_nel = .false.
      rd_sym = .false.
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
        case('dumpfile') 
          read (line(idelim+1:),'(a)') intfile
          rd_intfile = .true.
        case('memory')
          read (line(idelim+1:),*) mem_mpro
        case('nirrep')
          read (line(idelim+1:),*) nirr
          rd_nirr = .true.
        case('nelec')
          read (line(idelim+1:),*) nel
          rd_nel = .true.
        case('ecore')
          read (line(idelim+1:),*) ecore
          rd_ecore = .true.
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
!        case('itype-array')
!          if (rd_norbs.and.rd_ncore) then
!            nn = sum(norbs(1:nirr))-sum(ncore(1:nirr))
!            read(luintf,end=99,*) irrtyp(1:nn)
!            rd_irrtyp = .true.
        end select

      end do

 99   continue

      error = .false.
      if (.not.rd_intfile) then
        call warn(i_am,
     &     'could not read dumpfile name, trying "moints"')
        intfile = 'moints'
      end if
      if (.not.rd_nirr) then
        write(lulog,*) 'Could not read number of IRREPs'
        error = .true.
      end if
      if (.not.rd_nel) then
        write(lulog,*) 'Could not read number of electrons'
        error = .true.
      end if
      if (.not.rd_ecore) then
        write(lulog,*) 'Could not read core energy'
        error = .true.
      end if
      if (.not.rd_sym) then
        write(lulog,*) 'Could not read reference state IRREP'
        error = .true.
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
        write(lulog,*) 'Could not read number of core orbitals/IRREP'
        error = .true.
      end if
      if (.not.rd_closed) then
        write(lulog,*) 'Could not read number of closed orbitals/IRREP'
        error = .true.
      end if
      if (.not.rd_nocc) then
        write(lulog,*)'Could not read number of occupied orbitals/IRREP'
        error = .true.
      end if
!      if (.not.rd_irrtyp) then
!        write(lulog,*) 'Could not read array with IRREPs of orbitals'
!        error = .true.
!      end if

      if (error) call quit(1,i_am,'Error(s) reading interface file!')

      ! get inactive, active and secondary shells
      ninact(1:nirr) = nclosed(1:nirr)-ncore(1:nirr)
      nact(1:nirr)   = nocc(1:nirr)-nclosed(1:nirr)
      nvirt(1:nirr)  = norbs(1:nirr)-nocc(1:nirr)

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

      orb_info%ncore_mpro = sum(ncore(1:nirr))
      orb_info%name_intfile_mpro = intfile

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

      orb_info%mem_mpro = mem_mpro

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
