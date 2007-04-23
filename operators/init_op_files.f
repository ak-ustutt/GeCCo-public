*----------------------------------------------------------------------*
      subroutine init_op_files(ffops,ops,nops)
*----------------------------------------------------------------------*
*     basic initialization of files containing operator matrix elements
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'par_globalmarks.h'
      include 'ifc_memman.h'
      include 'ioparam.h'
      include 'def_operator.h'
      include 'def_filinf.h'

      integer, intent(in) ::
     &     nops

      type(operator), intent(in) ::
     &     ops(nops)
      type(filinf), intent(inout) ::
     &     ffops(nops)

      integer ::
     &     iop, ifree, nbuff

      do iop = 1, nops
        ! basic initialization
        call file_init(ffops(iop),
     &       'op_'//trim(ops(iop)%name)//'_elements.da',ftyp_da_unf,
     &       lblk_da)
        ! incore options
        ! preliminary version:
        if (trim(ops(iop)%name).eq.'H') then
          ! we directly use the knowledge that blocks 1-5 should be
          ! incore
          ffops(iop)%buffered = .true.
          ! make sure, we are in the appropriate memory section
          call mem_pushmark()
          ifree = mem_gotomark(op_files)

          ifree = mem_alloc_int(ffops(iop)%incore,
     &         ops(iop)%n_occ_cls,trim(ops(iop)%name)//'_incore')
          ffops(iop)%incore(1:5) = +1 ! dirty
          ffops(iop)%incore(6:ops(iop)%n_occ_cls) = -1
          nbuff = ops(iop)%off_op_occ(6) ! dirty
          ifree = mem_alloc_real(ffops(iop)%buffer,
     &         nbuff,trim(ops(iop)%name)//'_buffer')

          call mem_popmark()

        end if
      end do

      return
      end
