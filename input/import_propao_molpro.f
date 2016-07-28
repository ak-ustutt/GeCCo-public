*----------------------------------------------------------------------*
      subroutine import_propao_molpro(ffao,label,gamma,psym,orb_info)
*----------------------------------------------------------------------*
*     read property integral "label" from AOPROPER file
*     we have to provide the symmetry block (gamma) to be extracted
*     and the symmetry with respect to ij<->ji permutation (psym,
*     can be 1 or -1)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'
      include 'multd2h.h'
      include 'par_dalton.h'

      integer, parameter ::
     &     ntest = 100

      type(orbinf), intent(in), target ::
     &     orb_info
      type(filinf), intent(inout) ::
     &     ffao
      character(len=*), intent(in) ::
     &     label
      integer, intent(in) ::
     &     gamma, psym

      type(filinf) ::
     &     ffaoprop
      logical ::
     &     closeit,ok
      real(8) ::
     &     cpu0,sys0,wall0,cpu,sys,wall
      integer ::
     &     luaoprop, isym, jsym, ifree, luerror,
     &     nao_blk, nao_full,nao_i,nao_j,i,j,ij,ji,
     &     naoint,len_blk(8)
      character(8) ::
     &     propfile

      real(8), pointer ::
     &     ao_full(:), ao_blk(:)

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('import_ao_molpro')

      ! get buffers (SAO basis -> symmetry blocked)
      nao_blk = 0
      do isym = 1, orb_info%nsym
        jsym = multd2h(isym,gamma)
        nao_blk = nao_blk
     &          + orb_info%nbas(isym)*
     &            orb_info%nbas(jsym)
      end do
      
      ! buffer for AO-matrix as read from AOPROPER
      ifree = mem_alloc_real(ao_full,nao_blk,'ao_full')
      ! buffer for extracted matrix in symmetry blocked form
      ifree = mem_alloc_real(ao_blk,nao_blk,'ao_blk')

      inquire(file='DIPZ',exist=ok)
      if (.not.ok) call quit(0,'import_propao_molpro',
     &       'did not find any DIPZ file')

      propfile = 'DIPZ'

      ! open files
      call file_init(ffaoprop,propfile,ftyp_sq_frm,0)
      call file_open(ffaoprop)

      luaoprop = ffaoprop%unit
      rewind luaoprop
      
      if (ffao%unit.le.0) then
        call file_open(ffao)
        closeit = .true.
      else
        closeit = .false.
      end if

      read (luaoprop,*)
      read (luaoprop,*)
      ! read matrix in upper triangular form
      read (luaoprop,*) ao_full(1:nao_blk)

      call file_close_keep(ffaoprop)

      naoint = 0

      do isym = 1, orb_info%nsym
        jsym = multd2h(isym,gamma)
        nao_i = orb_info%nbas(isym)
        nao_j = orb_info%nbas(jsym)
        do i = 1, nao_j
          do j = 1, nao_i
            ij = (i-1)* nao_j + j
            ji = (j-1)* nao_i + i
            ao_blk(ji + naoint) = ao_full(ij + naoint)
          end do
        end do
        naoint = naoint +
     &       orb_info%nbas(jsym)*orb_info%nbas(isym)
      end do

      len_blk(1:orb_info%nsym) =
     &    orb_info%nbas(1:orb_info%nsym)

      if (ntest.ge.100) then
        write(lulog,*) 'Property integral: AO (blocked):'
        call wr_blkmat2(ao_blk,len_blk,len_blk,
     &                     orb_info%nsym,gamma,0)
      end if
      
      ! write to file
      call put_vec(ffao,ao_blk,1,nao_blk)

      if (closeit)
     &     call file_close_keep(ffao)

      ifree = mem_flushmark('import_ao_molpro')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.10) 
     &     call prtim(lulog,'time in property(ao) import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
