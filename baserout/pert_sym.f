*----------------------------------------------------------------------*
      integer function pert_sym(int_name,orb_info)
*----------------------------------------------------------------------*
*     determines IRREP of a given perturbation direction
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'multd2h.h'
      include 'def_filinf.h'
      include 'par_dalton.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      character(len=8), intent(in) ::
     &     int_name

      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     isym, gamma, jsym, nao_blk, nao_full, ifree, luaoprop,
     &     luerror, len_blk(8)

      real(8) ::
     &     norm(orb_info%nsym)

      real(8), pointer ::
     &     ao_full(:), ao_blk(:)

      character(len=1) ::
     &     irr1, irr2

      type(filinf) ::
     &     ffaoprop

      logical ::
     &     irrep_found

      ! buffer for AO-matrix as read from AOPROPER
      nao_full = (orb_info%nbast+orb_info%nxbast+1)*
     &           (orb_info%nbast+orb_info%nxbast)/2
      ifree = mem_alloc_real(ao_full,nao_full,'ao_full')

      luerror = luout

      ! open files
      call file_init(ffaoprop,aoproper,ftyp_sq_unf,0)
      call file_open(ffaoprop)

      luaoprop = ffaoprop%unit
      rewind luaoprop

      call mollab(int_name,luaoprop,luerror)

      ! read matrix in upper triangular form
      read (luaoprop) ao_full(1:nao_full)

      call file_close_keep(ffaoprop)

      ! get norm for all IRREPS
      norm = 0d0
      do gamma = 1, orb_info%nsym

        ! get buffers (SAO basis -> symmetry blocked)
        nao_blk = 0
        do isym = 1, orb_info%nsym
          jsym = multd2h(isym,gamma)
          nao_blk = nao_blk
     &            + (orb_info%nbas(isym)+orb_info%nxbas(isym))*
     &              (orb_info%nbas(jsym)+orb_info%nxbas(jsym))
        end do

        ! buffer for extracted matrix in symmetry blocked form
        ifree = mem_alloc_real(ao_blk,nao_blk,'ao_blk')

        ! reorder
        len_blk(1:orb_info%nsym) =
     &      orb_info%nbas(1:orb_info%nsym)+
     &      orb_info%nxbas(1:orb_info%nsym)
        call reo_full2sym(ao_blk,ao_full,orb_info%nbast+orb_info%nxbast,
     &       len_blk,orb_info%nsym,gamma,dble(1))

        ! calculate norm
        do isym = 1,nao_blk
          norm(gamma) = norm(gamma) + ao_blk(isym)**2
        end do

        ! clear buffer for extracted matrix in symmetry blocked form
        ifree = mem_dealloc('ao_blk')
      end do

      ! clear buffer for AO matrix as read from AOPROPER
      ifree = mem_dealloc('ao_full')

      ! pert_sym is only irrep with a non-zero norm
      irrep_found = .false.
      pert_sym = 0
      do gamma = 1, orb_info%nsym
        if (norm(gamma).gt.1d-12.and..not.irrep_found) then
          pert_sym = gamma
          irrep_found = .true.
        else if (norm(gamma).gt.1d-12) then
          write(irr1,'(i1)') pert_sym
          write(irr2,'(i1)') gamma
          call quit(1,'pert_sym','cannot tell if '//int_name//
     &       ' belongs to irrep '//irr1//' or '//irr2)
        end if
      end do
      if (.not.irrep_found) call quit(1,'pert_sym',
     &       'no irrep found for '//int_name)

      return
      end
