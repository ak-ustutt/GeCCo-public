*----------------------------------------------------------------------*
      subroutine import_cmo_molpro(ffcmo,orb_info)
*----------------------------------------------------------------------*
*     read MO-coefficients from SIRIFC (sirius interface file)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'
      include 'multd2h.h'
      include 'par_molpro.h'

      integer, parameter ::
     &     ntest = 100

      type(orbinf), intent(in), target ::
     &     orb_info
      type(filinf), intent(inout) ::
     &     ffcmo

      type(filinf) ::
     &     ffsir
      logical ::
     &     closeit, ok
      real(8) ::
     &     cpu0,sys0,wall0,cpu,sys,wall
      integer ::
     &     lusir, isym, ifree, igas, igasr,
     &     ncmo, imo, imo_reo, ioff, ioff_reo, len,
     &     ioff_sym(8), ioff_mo(8),
     &     nmo, nao, i, j, ij, ji, ncoef

      real(8), pointer ::
     &     cmo(:), cmo_tmp(:), cmo_reo(:)

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('import_cmo_molpro')

      ! get buffers (SAO basis -> symmetry blocked)
      ncmo = 0
      imo = 0
      do isym = 1, orb_info%nsym
        ioff_sym(isym) = ncmo
        ioff_mo(isym) =  imo
        ncmo = ncmo +
     &       orb_info%ntoobs(isym)*orb_info%nbas(isym)
        imo = imo + orb_info%ntoobs(isym)
      end do
      
      ! buffer for CMO-matrix as read from SIRIFC
      ifree = mem_alloc_real(cmo_tmp,ncmo,'cmo_ori_tmp')
      ifree = mem_alloc_real(cmo,ncmo,'cmo_ori')
      ! buffer for reordered CMO-matrix
      ifree = mem_alloc_real(cmo_reo,ncmo,'cmo_reo')

      inquire(file=f_coeff,exist=ok)
      if (.not.ok) call quit(0,'import_cmo_molpro',
     &       'did not find the CMOMOL file for the coeffs')

      ! open files
      call file_init(ffsir,f_coeff,ftyp_sq_frm,0)
      call file_open(ffsir)

      lusir = ffsir%unit
      rewind lusir
      
      if (ffcmo%unit.le.0) then
        call file_open(ffcmo)
        closeit = .true.
      else
        closeit = .false.
      end if

      ! skip 'BEGIN DATA' 
      read(lusir,*)
      ! skip orbital type
      read(lusir,*)

      ! and finally: the coefficient matrix
      ! in square matrix form
      read (lusir,*,end=3,err=6) cmo_tmp(1:ncmo)
      !read (lusir,*) cmo_tmp(1:ncmo)


      ! important for molpro coefficient file:
      ! molpro store the data in a adjoint format than that of dalton
      ! here the cmo matrix is made adjoint after reading it from the file

      ncoef = 0
      do isym = 1, orb_info%nsym
        nao = orb_info%nbas(isym)
        nmo = orb_info%ntoobs(isym)
        do i = 1, nmo
          do j = 1, nao
            ij = (i-1)* nmo + j
            ji = (j-1)* nao + i
            cmo(ji + ncoef) = cmo_tmp(ij + ncoef)
          end do
        end do
        ncoef = ncoef +
     &       orb_info%ntoobs(isym)*orb_info%nbas(isym)
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'CMO (original):'
        call wr_blkmat(cmo,orb_info%nbas,orb_info%ntoobs,
     &                     orb_info%nsym,0)
      end if

      call file_close_keep(ffsir)

      ! and reorder
      ioff_reo = 0
      do imo_reo = 1, orb_info%ntoob
        isym = orb_info%igamorb(imo_reo)
        imo = orb_info%ireots(imo_reo) - ioff_mo(isym)
        len  = orb_info%nbas(isym)
        ioff = ioff_sym(isym) + (imo-1)*len
        cmo_reo(ioff_reo+1:ioff_reo+len) = cmo(ioff+1:ioff+len)
        ioff_reo = ioff_reo+len
      end do

      if (ntest.ge.100) then
        ioff_reo = 1
        do igas = 1, orb_info%ngas
          igasr = orb_info%gas_reo(igas)
          write(lulog,*) 'orbital shell # ',
     &         igas,orb_info%igassh(1:orb_info%nsym,igasr)
          call wr_blkmat(cmo_reo(ioff_reo),orb_info%nbas,
     &                                     orb_info%igassh(1,igasr),
     &                                     orb_info%nsym,0)          
          do isym = 1, orb_info%nsym
            ioff_reo = ioff_reo + orb_info%nbas(isym)*
     &                            orb_info%igassh(isym,igasr)
          end do
        end do
      end if

! 114 call quit(0,'import_cmo_molpro','reading CMOMOL ended early')

! 116 call quit(0,'import_cmo_molpro','error in reading from CMOMOL')
      
      ! write to file
      call put_vec(ffcmo,cmo_reo,1,ncmo)

      if (closeit)
     &     call file_close_keep(ffcmo)

      ifree = mem_flushmark('import_cmo_molpro')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.10) 
     &     call prtim(lulog,'time in cmo import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return

    3 call quit(0,'import_cmo_molpro','reading CMOMOL ended early')

    6 call quit(0,'import_cmo_molpro','error in reading from CMOMOL')
      end
