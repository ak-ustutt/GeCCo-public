*----------------------------------------------------------------------*
      subroutine import_propao_molpro(ffao,label,gamma,orb_info)
*----------------------------------------------------------------------*
*     read one-electron integrals from AOPROPER file in MOLPRO format
*
*     label: describe the operator using the following sintax:
*            <labelM>:<labelD>
*            where these are the conventions for Molpro and Dalton
*            see get_oper_info.f and one_el_op_convention.txt
*
*     gamma: symmetry block to be extracted
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'
      include 'multd2h.h'
      include 'par_molpro.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in), target ::
     &     orb_info
      type(filinf), intent(inout) ::
     &     ffao
      character(len=*), intent(in) ::
     &     label
      integer, intent(in) ::
     &     gamma

      type(filinf) ::
     &     ffaoprop
      logical ::
     &     closeit,ok
      real(8) ::
     &     cpu0,sys0,wall0,cpu,sys,wall
      integer ::
     &     luaoprop, isym, jsym, ifree,
     &     nao,
     &     len_blk(8),
     &     gammaN, i_gammaN, iMD
      character(len=100)
     &     :: line
      real(8), pointer ::
     &     ao_full(:)


      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('import_ao_molpro')

      ! get buffers (SAO basis -> symmetry blocked)
      ! integrals corresponding to the F12 calculations are not
      ! considered here
      nao = 0
      do isym = 1, orb_info%nsym
        jsym = multd2h(isym,gamma)
        nao = nao
     &       + orb_info%nbas(isym)*orb_info%nbas(jsym)
      end do
      
      ! buffer for AO-matrix as read from AOPROPER
      ifree = mem_alloc_real(ao_full,nao,'ao_full')

      inquire(file=aoproper,exist=ok)
      if (.not.ok) call quit(0,'import_propao_molpro',
     &       'did not find any '//aoproper//' file')

      ! Get integrals
      call file_init(ffaoprop,aoproper,ftyp_sq_frm,0)
      call file_open(ffaoprop)
c     rewind ffaoprop%unit
      
      gammaN = -1
!     Find the position of label in the file and get gammaN
      iMD = index( label, ':')
      do while (gammaN .EQ. -1)
        read (ffaoprop%unit, '(a100)', end=3, err=6) line
        if (line(1:1).EQ.'#') then
          if  (index(line,'MATRIX '//label(1:iMD-1)) .NE. 0 .OR.
     &         index(line,'MATRIX '//trim(label(iMD+1:))) .NE. 0) then
            i_gammaN = index(line,'SYMMETRY=') + 9
            read(line(i_gammaN:i_gammaN), '(I1)') gammaN
          end if
        end if
      end do

      if (gammaN.NE.gamma) then
         write(lulog, '("Inconsistent irrep for operator ",A,".")')
     &        trim(label)
         call quit(0,'import_propao_molpro',
     &        'Inconsistent irrep for operator '//trim(label)//'.')
      end if

      ! read matrix in blocks form
      read (ffaoprop%unit,*) ao_full(1:nao)

      call file_close_keep(ffaoprop)

      len_blk(1:orb_info%nsym) =
     &    orb_info%nbas(1:orb_info%nsym)

      if (ntest.ge.100) then
        write(lulog,*) '----------------------------------'
        write(lulog,*) 'One electron integrals of '//trim(label)//
     &        ' in AO (blocked):'
        call wr_blkmat2(ao_full,len_blk,len_blk,
     &                     orb_info%nsym,gamma,0)
      end if
      
!     Write to file
      if (ffao%unit.le.0) then
        call file_open(ffao)
        closeit = .true.
      else
        closeit = .false.
      end if
      call put_vec(ffao,ao_full,1,nao)
      if (closeit)
     &     call file_close_keep(ffao)

      !ifree = mem_flushmark('import_ao_molpro')
      ifree = mem_flushmark()

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.10) 
     &     call prtim(lulog,'time in property(ao) import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return

 3    write(lulog, 7) trim(label)
      call quit(0,'import_propao_molpro',
     &     'One electron operator label not found: end of file.')
c
 6    write(lulog, 7) trim(label)
      call quit(0,'import_propao_molpro',
     &     'One electron operator label not found: error.')
c
 7    format(/' Operator label ', A ,
     &     ' not found in the one electron integrals file'//
     &     ' Did you request it in the molpro input?')

      end
