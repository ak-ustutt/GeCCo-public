*----------------------------------------------------------------------*
      subroutine init_op_files(ffops,ops,nops)
*----------------------------------------------------------------------*
*     basic initialization of files containing operator matrix elements
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
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
     &     iop, ifree, nbuff, iocc

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
          ifree = mem_alloc_int(ffops(iop)%idxrec,
     &         ops(iop)%n_occ_cls,trim(ops(iop)%name)//'_idxrec')
          nbuff = 0
          do iocc = 1, ops(iop)%n_occ_cls
            ! allocate buffer for 1-hamiltonian only
            if (max(ops(iop)%ica_occ(1,iocc),
     &              ops(iop)%ica_occ(2,iocc)).gt.1 .or.
     &          iextr.gt.0.and. ! R12: ignore currently
     &          ops(iop)%ihpvca_occ(iextr,1,iocc)+
     &          ops(iop)%ihpvca_occ(iextr,2,iocc).gt.0 ) then
              ffops(iop)%incore(iocc) = -1
              ffops(iop)%idxrec(iocc) = -1
            else
              ffops(iop)%idxrec(iocc) = nbuff
              ffops(iop)%incore(iocc) = ops(iop)%len_op_occ(iocc)
              nbuff = nbuff + ops(iop)%len_op_occ(iocc)
            end if
          end do
          ffops(iop)%nbuffer = nbuff
          ifree = mem_alloc_real(ffops(iop)%buffer,
     &         nbuff,trim(ops(iop)%name)//'_buffer')

          call mem_popmark()

        end if
      end do

      return
      end
