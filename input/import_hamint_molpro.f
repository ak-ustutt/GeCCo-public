*----------------------------------------------------------------------*
      subroutine import_hamint_molpro(hlist,str_info,orb_info)
*----------------------------------------------------------------------*
*     import one- and two-electron integrals from the file generated
*     by the molpro interface ...
*     initially, we fix the name as "moints"
*
*     andreas, jan 16
*
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

      character(len=35), parameter ::
     &     i_am = 'import_hamint_molpro'
      
      type(me_list), intent(inout) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(filinf) ::
     &     ffmoint
      integer ::
     &     lumoint, ifree, nfock, ifpos
      character(len=256) ::
     &     line

      real(8) ::
     &     ecore
      real(8), pointer ::
     &     fock(:)

      ! open then DUMP file
      call file_init(ffmoint,'moints',ftyp_st_unf,0)
      call file_open(ffmoint)
      lumoint = ffmoint%unit

      ! allocate space for Fock matrix
      ifree = mem_setmark('import_hamint')
      nfock = (orb_info%ntoob+orb_info%caborb)*
     &        (orb_info%ntoob+orb_info%caborb+1)/2

      ifree = mem_alloc_real(fock,nfock,'raw_fock')

      fock(1:nfock) = 0d0

      ! initially the moint is positioned on the first record
      ifpos = 1
      ! that contains a two-electron integral
      ! transfer control to H2 reader      
      call import_h2_molpro_ifc(lumoint,ifpos,hlist,fock,nfock,
     &     str_info,orb_info)

      ! the previous routine has exited with the file 
      ! positioned on the first record that contains F
      call import_fock_molpro_ifc(lumoint,ifpos,hlist,ecore,fock,nfock,
     &     str_info,orb_info)

      ifree = mem_flushmark('import_hamint')

      call file_close_keep(ffmoint)

      return
      end
