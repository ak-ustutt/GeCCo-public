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
     &     ntest = 00

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
     &     ifree, nfull, nblkd, isym, jsym, ierr, sym,
     &     naoint, i, j, ij, ji, nao_i, nao_j
      type(filinf) ::
     &     ffprop
      real(8) ::
     &     xnorm, xtrace
      real(8), pointer ::
     &     dao(:), xint_raw(:), xint_blk(:)

      real(8), external ::
     &     dnrm2, ddot
      logical, external ::
     &     next_proper
      
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

      if (ntest.ge.100) write(lulog,*) 'scanning ',aoproper
c dbg
c        print *,'LABEL: ',trim(label)
c dbg

      do while(next_proper(ffprop%unit,sym,label))

         if (ntest.ge.100) write(lulog,*)'label: ',trim(label),
     &                                   ' sym: ',sym

         if (sym.ne.dens%gamt) cycle

         read(ffprop%unit,*,end=3,err=6) xint_raw(1:nblkd)

         ! not clear, why ....
c         naoint = 0
c
c         do isym = 1, orb_info%nsym
c           jsym = multd2h(isym,dens%gamt)
c           nao_i = orb_info%nbas(isym)
c           nao_j = orb_info%nbas(jsym)
c           write(lulog,*) 'isym, jsym, nao_i, nao_j, naoint:',
c     &          isym, jsym, nao_i, nao_j, naoint
c           do i = 1, nao_i
c             do j = 1, nao_j
c               ij = (i-1)* nao_j + j
c               ji = (j-1)* nao_i + i
c               xint_blk(ji + naoint) = xint_raw(ij + naoint)
c               write(lulog,*) 'resorting: ',ij,' to ',ji
c             end do
c           end do
c           naoint = naoint +
c     &          orb_info%nbas(jsym)*orb_info%nbas(isym)
c         end do
         ! seems to be in correct order ...
         xint_blk = xint_raw

         if (ntest.ge.100) then 
           write(lulog,*) 'AO integrals (original):'
           call wr_blkmat2(xint_blk,orb_info%nbas,orb_info%nbas,
     &                        orb_info%nsym,dens%gamt,0)
         end if

         ! check, whether integral block is nonzero:
         xnorm = dnrm2(nblkd,xint_blk,1)
!        if (xnorm.lt.1d-12) 

         xtrace = ddot(nblkd,xint_blk,1,dao,1)
         xtrace = -1.0*xtrace ! multiplied with -1, unlike dalton, as the property integrals from
                           ! molpro comes with the opposite sign 

         write(lulog,'(2x,">>> ",a," : ",g20.10)') trim(label),xtrace
         if (lulog.ne.luout)
     &       write(luout,'(2x,">>> ",a," : ",g20.10)') 
     &                                             trim(label),xtrace

      end do

      call file_close_keep(ffprop)

      ifree = mem_flushmark('1prop_molpro')

      return

    3 call quit(0,'import_cmo_molpro','reading property integrals 
     &             ended early')

    6 call quit(0,'import_cmo_molpro','error in reading property
     &  integrals')
      end


      logical function next_proper(lu,sym,label)

      implicit none

      include "stdunit.h"

      integer, intent(in) :: lu
      integer, intent(out) :: sym
      character(len=8), intent(out) :: label
     
      integer ipos
      character(len=256) line

      do
         read(lu,*,end=12) line
         if (line(1:10)=='BEGIN_DATA') then
            read(lu,'(a)',end=12) line
            label(1:8) = ' '
            label = line(10:17)
            ipos = index(line,'SYMMETRY=')
            if (ipos<=0) goto 12
            read(line(ipos+9:),*) sym
            next_proper = .true.
            return
         end if
      end do

 12   next_proper = .false.
      return

      end
