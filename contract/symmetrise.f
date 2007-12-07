      subroutine symmetrise(fac,ffin,op_in,ffout,op_out,nocc_cls,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to symmetrise an operator matrix:
*
*     O^{symm}_{ij} = fac*(O_{ij}+O_{ji}) = O^{symm}_{ji}
*
*     ffin and op_in are the input, asymmetric matrices
*     ffout and op_out represent the output, symmetric matrices.
*     It is assumed that the two operators have the same shape.
*     GWR November 2007
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'def_operator_array.h'
      include 'def_file_list.h'
      include 'def_file_array.h'
      include 'def_operator_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 100

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(operator), intent(in) ::
     &     op_in, op_out
      type(filinf), intent(in) ::
     &     ffin
      type(filinf), intent(inout) ::
     &     ffout
      integer, intent(in) ::
     &     nocc_cls
      real(8), intent(in) ::
     &     fac

      logical ::
     &     bufin, bufout
      integer ::
     &     ifree, nblk, nbuff, iocc_cls, idxmsa,
     &     msmax, msa, igama, idx, jdx, ngam, len_gam_ms, ioff

      logical ::
     &     loop(nocc_cls)
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      integer, external ::
     &     idx_oplist2

      call file_open(ffout)

      ! Check whether files are buffered.
      bufin = .false.
      if(ffin%buffered) bufin = .true. 
      bufout = .false.
      if(ffout%buffered) bufout = .true.

      ifree = mem_setmark('symmetrise')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Loop array.
      loop(1:nocc_cls) = .false.

      ! Allocations made to maximum block length to save time.
      if(.not.bufin)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls 
          if(op_in%formal_blk(iocc_cls))
     &         loop(iocc_cls) = .true.

          nbuff = nbuff + op_in%len_op_occ(iocc_cls)
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
        call get_vec(ffin,buffer_in,1,nbuff)

      else
        if(ntest.ge.100)
     &       write(luout,*)'Symmetrise: input not incore'
        buffer_in => ffin%buffer(1:)
      endif

      if(.not.bufout)then
        nbuff=0
        do iocc_cls = 1, nocc_cls 
          nbuff = nbuff + op_out%len_op_occ(iocc_cls)
        enddo
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
        buffer_out(1:nbuff) = 0d0
      else
        if(ntest.ge.100)
     &       write(luout,*)'Symmetrise: output not incore'
        buffer_out => ffout%buffer(1:)
      endif

      ! Loop over occupation class.
      iocc_loop: do iocc_cls = 1, nocc_cls 

        ! Loop over Ms of annihilator string.
        idxmsa = 0
        msmax = 2
        msa_loop : do msa = msmax, -msmax, -2

          idxmsa = idxmsa+1
      
          ! Loop over Irrep of annihilator string.
          igama_loop: do igama =1, ngam

            len_gam_ms = int(sqrt(dble(op_out%
     &           len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa))))

            ioff = op_out%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)

            idx_loop: do idx = 1,len_gam_ms
              jdx_loop: do jdx = 1,len_gam_ms
            
                buffer_out((idx-1)*len_gam_ms+jdx+ioff) = 
     &             fac * (buffer_in((idx-1)*len_gam_ms+jdx+ioff) +
     &             buffer_in((jdx-1)*len_gam_ms+idx+ioff))
            
              enddo jdx_loop  
            enddo idx_loop
          
          enddo igama_loop
          
        enddo msa_loop
      enddo iocc_loop

      if(.not.bufout)then
        call put_vec(ffout,buffer_out,1,nbuff)
      endif  

      call file_close_keep(ffout)

      ifree = mem_flushmark('symmetrise')

      return
      end
