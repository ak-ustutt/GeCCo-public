*----------------------------------------------------------------------*
      integer function pert_sym_dalton(int_name,orb_info)
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
     &     luerror, len_blk(8), psign

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

      select case(trim(int_name))
      case ('XDIPLEN','YDIPLEN','ZDIPLEN')
        psign = 1
      case ('XDIPVEL','YDIPVEL','ZDIPVEL',
     &      'XANGMOM','YANGMOM','ZANGMOM')
        psign = -1
      case default
        call quit(1,'pert_sym_dalton','DALTON: cannot handle list_type"'
     &       //trim(int_name)//'"')
      end select

      ! buffer for AO-matrix as read from AOPROPER
      nao_full = (orb_info%nbast+orb_info%nxbast+1)*
     &           (orb_info%nbast+orb_info%nxbast)/2
      ifree = mem_alloc_real(ao_full,nao_full,'ao_full')

      luerror = lulog

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
     &       len_blk,orb_info%nsym,gamma,dble(psign))

        ! calculate norm
        do isym = 1,nao_blk
          norm(gamma) = norm(gamma) + ao_blk(isym)**2
        end do

        ! clear buffer for extracted matrix in symmetry blocked form
        ifree = mem_dealloc('ao_blk')
      end do

      ! clear buffer for AO matrix as read from AOPROPER
      ifree = mem_dealloc('ao_full')

      ! pert_sym_dalton is only irrep with a non-zero norm
      irrep_found = .false.
      pert_sym_dalton = 0
      do gamma = 1, orb_info%nsym
        if (norm(gamma).gt.1d-12.and..not.irrep_found) then
          pert_sym_dalton = gamma
          irrep_found = .true.
        else if (norm(gamma).gt.1d-12) then
          write(irr1,'(i1)') pert_sym_dalton
          write(irr2,'(i1)') gamma
          call quit(1,'pert_sym_dalton','cannot tell if '//int_name//
     &       ' belongs to irrep '//irr1//' or '//irr2)
        end if
      end do
      if (.not.irrep_found) call quit(1,'pert_sym_dalton',
     &       'no irrep found for '//int_name)

      return

      end
*----------------------------------------------------------------------*
      integer function pert_sym_molpro(int_name,orb_info)
*----------------------------------------------------------------------*
*     Symmetry of the perturbation has been obtained in a different way
*     for integrals from MOLPRO interfaces than that from DALTON
*     interface. In MOLPRO the property integrals are written in blocks
*     rather than in triangular form like the integrals in DALTON.
*     Pradipta, 2017
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'par_molpro.h'
      include 'def_orbinf.h'

      character(len=8), intent(in) ::
     &     int_name

      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     luaoprop, luerr

      character(len=48) :: b
      character(len=1) :: c

      logical :: ok
      type(filinf) ::
     &     ffaoprop

      luerr = lulog

      inquire(file=aoproper,exist=ok)
      if (.not.ok) call quit(0,'pert_symbol_molpro',
     &       'did not find any property file')

      ! open files
      call file_init(ffaoprop,aoproper,ftyp_sq_frm,0)
      call file_open(ffaoprop)

      luaoprop = ffaoprop%unit
      rewind luaoprop

      ! here we will read each line to match the prop label
      ! and then will simply get the symmetry of the perturbation
      ! written in the same line

    1 read (luaoprop,9,end=3,err=6) b
      if (b(10:16).ne.int_name(1:7)) go to 1
      if (luerr.lt.0) luerr = 0
      c=b(46:46)
      read(c,'(I1)') pert_sym_molpro
      call file_close_keep(ffaoprop)
      return
    9 format(a48)
 
    3 if (luerr.lt.0) then
         luerr = -1
         return
      else
         write(luerr,4)int_name,luaoprop
         call quit(0,'pert_sym_molpro','molecule label not found 
     &             on file')
      end if
      call file_close_keep(ffaoprop)
    4 format(/' MOLECULE label ',A8,
     *        ' not found on unit',I4)
 
    6 if (luerr.lt.0) then
         luerr = -2
         return
      else
         write (luerr,7) luaoprop,int_name
         call quit(0,'pert_sym_molpro', 'error reading file')
      end if
    7 format(/' error reading unit',I4,
     *       /T22,'when searching for label ',A8)
      call file_close_keep(ffaoprop)
      

      end

