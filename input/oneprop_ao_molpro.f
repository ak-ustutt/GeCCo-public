*----------------------------------------------------------------------*
      subroutine oneprop_ao_molpro(ffdao,dens,orb_info)
*----------------------------------------------------------------------*
*     contract 1-particle density matrix (on ffdao) with available
*     integrals from the DALTON environment
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'par_molpro.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 100

      type(filinf), intent(inout) ::
     &     ffdao
      type(me_list), intent(in) ::
     &     dens
      type(orbinf), intent(in) ::
     &     orb_info

      character ::
     &     label*8
      logical ::
     &     closeit, ok
      integer ::
     &     ifree, nfull, nblkd, isym, jsym, ierr,
     &     naoint, i, j, ij, ji, nao_i, nao_j
      type(filinf) ::
     &     ffprop
      real(8) ::
     &     xnorm, xtrace
      character(len=256) ::
     &     line
      real(8), pointer ::
     &     dao(:), xint_raw(:), xint_blk(:)

      real(8), external ::
     &     dnrm2, ddot
      
      ifree = mem_setmark('1prop_molpro')

      nblkd = 0
      do isym = 1, orb_info%nsym
        jsym = multd2h(isym,dens%gamt)
        nblkd = nblkd+orb_info%nbas(isym)*orb_info%nbas(jsym)
      end do

c dbg
c      print *,'integrals:'
c dbg      
      ifree = mem_alloc_real(dao,nblkd,'dao')
      ifree = mem_alloc_real(xint_raw,nblkd,'int_raw')
      ifree = mem_alloc_real(xint_blk,nblkd,'int_blk')

      closeit = .false.
      if (ffdao%unit.le.0) then
        call file_open(ffdao)
        closeit = .true.
      end if

      call get_vec(ffdao,dao,1,nblkd)

      if (closeit) call file_close_keep(ffdao)

      inquire(file=aoproper,exist=ok)
      if (.not.ok) call quit(0,'oneprop_ao_molpro',
     &       'did not find any integral file')

      ! open file with property integrals
      call file_init(ffprop,aoproper,ftyp_sq_frm,0)
      call file_open(ffprop)
      ! loop over list of property integrals
      rewind ffprop%unit

c dbg
c        print *,'LABEL: ',trim(label)
c dbg

      read(ffprop%unit,*)
      read(ffprop%unit,'(a)') line
      label = line(10:17)

      read(ffprop%unit,*,end=3,err=6) xint_raw(1:nblkd)

      naoint = 0

      do isym = 1, orb_info%nsym
        jsym = multd2h(isym,dens%gamt)
        nao_i = orb_info%nbas(isym)
        nao_j = orb_info%nbas(jsym)
        do i = 1, nao_j
          do j = 1, nao_i
            ij = (i-1)* nao_j + j
            ji = (j-1)* nao_i + i
            xint_blk(ji + naoint) = xint_raw(ij + naoint)
          end do
        end do
        naoint = naoint +
     &       orb_info%nbas(jsym)*orb_info%nbas(isym)
      end do

      if (ntest.ge.100) then 
        write(lulog,*) 'AO integrals (original):'
        call wr_blkmat(xint_blk,orb_info%nbas,orb_info%ntoobs,
     &                     orb_info%nsym,0)
      end if

      ! check, whether integral block is nonzero:
      xnorm = dnrm2(nblkd,xint_blk,1)
!     if (xnorm.lt.1d-12) 

      xtrace = ddot(nblkd,xint_blk,1,dao,1)
      xtrace = -1.0*xtrace ! multiplied with -1, unlike dalton, as the property integrals from
                           ! molpro comes with the opposite sign 

      write(lulog,'(2x,">>> ",a," : ",g20.10)') trim(label),xtrace


      call file_close_keep(ffprop)

      ifree = mem_flushmark('1prop_molpro')

      return

    3 call quit(0,'import_cmo_molpro','reading property integrals 
     &             ended early')

    6 call quit(0,'import_cmo_molpro','error in reading property
     &  integrals')
      end

