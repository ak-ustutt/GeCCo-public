*----------------------------------------------------------------------*
      subroutine import_hamint_molpro_dump(hlist,str_info,orb_info)
*----------------------------------------------------------------------*
*     import one- and two-electron integrals from MOLPRO's FCIDUMP file
*     we need:
*      FCIDUMP for fock matrix
*      FCIDUMP for two-electron integrals
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'

      character(len=40), parameter ::
     &     i_am = 'import_hamint_molpro_dump'
      
      type(me_list), intent(inout) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(filinf) ::
     &     ffdump
      integer ::
     &     ludump, ifree, nfock
      character(len=256) ::
     &     line

      real(8) ::
     &     ecore
      real(8), pointer ::
     &     fock(:)

      ! the namelist for dummy read
      integer, parameter ::
     &     maxorb=1000
      integer ::
     &     norb, nelec, ms2, orbsym(maxorb), isym
      namelist /FCI/ norb, nelec, ms2, orbsym, isym

      ! open then DUMP file
      call file_init(ffdump,'FCIDUMP',ftyp_sq_frm,0)
      call file_open(ffdump)
      ludump = ffdump%unit

      read(ludump,'(a)') line

      if (line(1:5).ne.' &FCI') then
        write(lulog,*) 'File FCIDUMP does not start with "&FCI":'
        write(lulog,*) trim(line)
        call quit(0,i_am,'wrong format of FCIDUMP?')
      end if
      rewind ludump
      ! read as namelist
      read(ludump,nml=fci)

      ! allocate space for Fock matrix
      ifree = mem_setmark('import_hamint')
      nfock = (orb_info%ntoob+orb_info%caborb)*
     &        (orb_info%ntoob+orb_info%caborb+1)/2

      ifree = mem_alloc_real(fock,nfock,'raw_fock')

      fock(1:nfock) = 0d0

      ! now the dumpfile is positioned on the first record
      ! that contains a two-electron integral
      ! transfer control to H2 reader      
      call import_h2_molpro_dump(ludump,hlist,fock,nfock,
     &     str_info,orb_info)

      ! the previous routine has exited with ludump
      ! positioned on the first record that contains F
      call import_fock_molpro_dump(ludump,hlist,ecore,fock,nfock,
     &     str_info,orb_info)

      ifree = mem_flushmark('import_hamint')

      return
      end
