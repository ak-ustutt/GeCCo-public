*----------------------------------------------------------------------*
      subroutine oneprop_ao_dalton(ffdao,dens,orb_info)
*----------------------------------------------------------------------*
*     contract 1-particle density matrix (on ffdao) with available
*     integrals from the DALTON environment
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'par_dalton.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'multd2h.h'

      type(filinf), intent(inout) ::
     &     ffdao
      type(me_list), intent(in) ::
     &     dens
      type(orbinf), intent(in) ::
     &     orb_info

      character ::
     &     label*8
      logical ::
     &     closeit
      integer ::
     &     ifree, nfull, nblkd, isym, jsym, ierr
      type(filinf) ::
     &     ffprop
      real(8) ::
     &     xnorm, xtrace
      real(8), pointer ::
     &     dao(:), xint_raw(:), xint_blk(:)

      real(8), external ::
     &     dnrm2, ddot
      
      ifree = mem_setmark('1prop_dalton')

      nfull = (orb_info%nbast+1)*orb_info%nbast/2
      nblkd = 0
      do isym = 1, orb_info%nsym
        jsym = multd2h(isym,dens%gamt)
        nblkd = nblkd+orb_info%nbas(isym)*orb_info%nbas(jsym)
      end do

c dbg
c      print *,'integrals:'
c dbg      

      ifree = mem_alloc_real(dao,nblkd,'dao')
      ifree = mem_alloc_real(xint_raw,nfull,'int_raw')
      ifree = mem_alloc_real(xint_blk,nblkd,'int_blk')

      closeit = .false.
      if (ffdao%unit.le.0) then
        call file_open(ffdao)
        closeit = .true.
      end if

      call get_vec(ffdao,dao,1,nblkd)

      if (closeit) call file_close_keep(ffdao)

      ! open file with property integrals
      call file_init(ffprop,aoproper,ftyp_sq_unf,0)
      call file_open(ffprop)
      ! loop over list of property integrals
      rewind ffprop%unit
      do 
        call next_mollab(label,ffprop%unit,ierr)
        if (label.eq.eof_lab.or.ierr.lt.0) exit
c dbg
c        print *,'LABEL: ',trim(label)
c dbg
        read(ffprop%unit) xint_raw(1:nfull)
        call get_symblk_rank1(xint_blk,xint_raw,dens%gamt,1d0,
     &       orb_info%nbast,orb_info%nbas,orb_info%nsym)

        ! check, whether integral block is nonzero:
        xnorm = dnrm2(nblkd,xint_blk,1)
        if (xnorm.lt.1d-12) cycle

        xtrace = ddot(nblkd,xint_blk,1,dao,1)

        write(lulog,'(2x,">>> ",a," : ",g20.12)') trim(label),xtrace

      end do

      call file_close_keep(ffprop)

      ifree = mem_flushmark('1prop_dalton')

      return
      end

