      subroutine invert(ffinp,op_inp,ffinv,op_inv,nocc_cls,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to invert a 2e-operator matrix. 
*     Currently only works for operators with a single occupancy block.
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
     &     op_inp, op_inv
      type(filinf), intent(in) ::
     &     ffinp
      type(filinf), intent(inout) ::
     &     ffinv
      integer, intent(in) ::
     &     nocc_cls

      logical ::
     &     bufin, bufout
      logical ::
     &     loop(nocc_cls)
      integer ::
     &     len_str, ifree, nbuff, idxmsa, idx_out, iocc_cls,
     &     msmax, msa, igama, idx, jdx, ngam,
     &     len_gam_ms, ioff, nrot
      real(8) ::
     &     min_eig, shift
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:), scratch(:,:), eigen_val(:),
     &     eigen_vec(:,:)

      integer, external ::
     &     idx_oplist2

      ! Check whether files are buffered.
      bufin = .false.
      if(ffinp%buffered) bufin = .true.
      bufout = .false.
      if(ffinv%buffered) bufout = .true.

      ifree = mem_setmark('invert')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Initialise loop array.
      loop(1:nocc_cls) = .false.

      ! Allocations made to maximum block length to save time.
      if(.not.bufin)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls
          if(op_inp%formal_blk(iocc_cls))
     &         loop(iocc_cls) = .true.

          nbuff = nbuff+op_inp%len_op_occ(iocc_cls)
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
        call get_vec(ffinp,buffer_in,1,nbuff)
      else
        if(ntest.ge.100)
     &       write(luout,*)'Invert: input not incore'
        buffer_in => ffinp%buffer(1:)
      endif

      if(.not.bufout)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls
          nbuff = nbuff + op_inv%len_op_occ(iocc_cls)
        enddo
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
        buffer_out(1:nbuff) = 0d0
      else
        if(ntest.ge.100)
     &       write(luout,*)'Symmetrise: output not incore'
        buffer_out => ffinv%buffer(1:)
      endif

      ! Do 2 loops.
      ! First finds the eigenvalues of the matrix and checks whether the
      ! matrix is positive definite. NB The matrix must be symmetric.
      ! Second does the inversion, but before it does so it adds the 
      ! absolute of the most negative eigenvalue to all diagonal elements.
      min_eig = 1d20
      out_loop: do idx_out = 1, 2

        if(idx_out.eq.2.and.min_eig.le.0d0)then
          write(luout,'(a,f12.6)')'Negative eigenvalue: ',min_eig
          write(luout,*) 'Matrix diagonal will be shifted.'
          shift = abs(min_eig)+0.5d0
        endif
        

        ! Loop over occupation class.
        iocc_loop: do iocc_cls = 1, nocc_cls
          if(loop(iocc_cls)) cycle iocc_loop

          ! Loop over Ms of annihilator string.
          idxmsa = 0
          msmax = 2
          msa_loop : do msa = msmax, -msmax, -2

            idxmsa = idxmsa+1
      
            ! Loop over Irrep of annihilator string.
            igama_loop: do igama =1, ngam

              len_gam_ms = int(sqrt(dble(op_inv%
     &           len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa))))

              ioff = op_inv%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)

              allocate(scratch(len_gam_ms,len_gam_ms))
              do idx = 1,len_gam_ms
                do jdx = 1,len_gam_ms
                  scratch(idx,jdx) =
     &                 buffer_in(ioff+(idx-1)*len_gam_ms+jdx)
                enddo
              enddo

              ! Call the inversion or eigenvalue routines.
              if(idx_out.eq.1)then
                if(len_gam_ms.gt.0)then
                  allocate(eigen_vec(len_gam_ms,len_gam_ms),
     &                 eigen_val(len_gam_ms))
                  call jacobi(scratch,len_gam_ms,len_gam_ms,eigen_val,
     &                 eigen_vec,nrot)

                  do idx = 1, len_gam_ms
                    if(eigen_val(idx).lt.min_eig)
     &                   min_eig = eigen_val(idx)
                  enddo
                  write(luout,*)'eigen_val: ',eigen_val
                    
                  deallocate(eigen_vec,eigen_val)
                endif
              else
                if(min_eig.eq.1d20)then
                  call quit(1,'invert','huge eigenvalues in matrix')
                elseif(min_eig.le.0d0)then
                  if(abs(min_eig).gt.1d20)then
                    call quit(1,'invert','huge negative eigenvalues')
                  else
                    do idx = 1,len_gam_ms
                      scratch(idx,idx) = scratch(idx,idx) + shift
                    enddo
                  endif
                endif

                call gaussj(scratch,len_gam_ms,len_gam_ms)

                do idx = 1,len_gam_ms
                  do jdx = 1,len_gam_ms
                    buffer_out((idx-1)*len_gam_ms+jdx+ioff) =
     &                   scratch(idx,jdx)                
                  enddo
                enddo
              endif

              deallocate(scratch)
                      
            enddo igama_loop
          
          enddo msa_loop

        enddo iocc_loop

      enddo out_loop

      if(.not.bufout)then
        call put_vec(ffinv,buffer_out,1,nbuff)
      endif  

      ifree = mem_flushmark('invert')

      return
      end
