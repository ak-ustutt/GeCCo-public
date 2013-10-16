      subroutine read_env_molpro_dump(orb_info)

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
     &     ntest = 100
      character(len=20), parameter ::
     &     i_am = 'read_env_molpro_dump'

      type(orbinf) ::
     &     orb_info

      type(filinf) ::
     &     ffdump
      integer ::
     &     ludump, iprint, ngas, ii, nirr

      integer, parameter ::
     &     mxsym = 8
      character(len=256) :: line

      ! the namelist to read
      integer, parameter ::
     &     maxorb=1000
      integer :: 
     &     norb, nelec, ms2, orbsym(maxorb), isym
      namelist /FCI/ norb, nelec, ms2, orbsym, isym


      iprint = max(iprlvl,ntest)

      call file_init(ffdump,'FCIDUMP',ftyp_sq_frm,0)
      call file_open(ffdump)
      ludump = ffdump%unit

      read(ludump,'(a)') line

      if (line(1:5).ne.' &FCI') then
        write(luout,*) 'File FCIDUMP does not start with "&FCI":'
        write(luout,*) trim(line)
        call quit(0,i_am,'wrong format of FCIDUMP?')
      end if
      rewind ludump
      ! read as namelist
      read(ludump,nml=fci)

      if (ntest.ge.100) then
        write(luout,*) 'Namelist FCI:'
        write(luout,'(1x,a,i5)') 'NORB =  ',norb
        write(luout,'(1x,a,i5)') 'NELEC = ',nelec
        write(luout,'(1x,a,i5)') 'MS2 =   ',ms2
        write(luout,'(1x,a,i5)') 'ISYM =  ',isym
        write(luout,'(1x,a)')   'ORBSYM:'
        write(luout,'(1x,20i3)') ORBSYM(1:NORB)  
      end if

C takes too long for large cases, we read enuc when importing
C the integrals ...
C      ! fast-forward to end-of-file
C      do
C        read(ludump,'(a)',end=1234) line
C      end do
C 1234 read(line,*) enuc, i1, i2, i3, i4
C
C      if (ntest.ge.100) write(luout,*) 'last rec: ',enuc, i1, i2, i3, i4

      ! analyze orbsym to get number of irreps
      nirr = 0
      do ii = 1, norb
        nirr = max(nirr,orbsym(ii))
      end do

      ! well, only 1, 2, 4, and 8 are allowed answers
      if (nirr.eq.3) nirr = 4
      if (nirr.gt.4) nirr = 8

      ! set orb_info
      orb_info%nsym = nirr
      orb_info%ngas = 1 ! no way to extract info on other spaces yet
      orb_info%nspin = 1
 
      orb_info%nbast = norb
      orb_info%ntoob = norb
      orb_info%lsym  = isym
      orb_info%ims   = ms2
      orb_info%nxbast = 0
      orb_info%caborb = 0      

      ! save this info here
      orb_info%nactel = nelec

      allocate(orb_info%igassh(nirr,1),orb_info%iad_gas(1),
     &         orb_info%ihpvgas(1,1),orb_info%nbas(nirr),
     &         orb_info%ntoobs(nirr))
      allocate(orb_info%cab_orb(nirr),
     &              orb_info%nxbas(nirr) )

      orb_info%iad_gas(1)=0
      orb_info%ihpvgas(1,1)=IPART
      orb_info%igassh(1:nirr,1) = 0
      do ii = 1, norb
        orb_info%igassh(orbsym(ii),1) = orb_info%igassh(orbsym(ii),1)+1
      end do
      orb_info%nbas(1:nirr) = orb_info%igassh(1:nirr,1)
      orb_info%ntoobs(1:nirr) = orb_info%igassh(1:nirr,1)
      orb_info%cab_orb(1:nirr)=0
      orb_info%nxbas(1:nirr)=0

      call file_close_keep(ffdump)

      return
      end
