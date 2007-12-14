      subroutine op_vec_mult(sgn,ffin,opin,vec,len_vec,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to multiply the matrix of an operator, with the vector,
*     vec. Vec is rearranged as an N*N matrix to fit the operator.
*     i.e. A.vec.... = vec(new)
*     GWR December 2007
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
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(filinf), intent(in) ::
     &     ffin
      type(operator), intent(in) ::
     &     opin
      integer, intent(in) ::
     &     len_vec
      real(8), intent(in) ::
     &     sgn 
      real(8), intent(inout) ::
     &     vec(len_vec)

      logical ::
     &     bufin, closeit, loop(opin%n_occ_cls)
      integer ::
     &     iops, nbuffered, idxbuf, iocc_cls ,nocc_cls, ibufloop,
     &     ifree, nbuff, idxmsa, msmax, msa, igama, idx, jdx, kdx, ngam,
     &     len_gam_ms, ioff
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:), scratch1(:,:), scratch2(:,:),
     &     scratch3(:,:)

      if(ntest.ge.100)then
        write(luout,*)'======================'
        write(luout,*)' Entered op_vec_mult. '
        write(luout,*)'======================'
      endif

      ifree = mem_setmark('op_vec_mult')

      ! Open the input file if necessary.
      closeit = .false.
      if(ffin%unit.le.0)then
        call file_open(ffin)
        closeit = .true.
      endif

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym
      nocc_cls = opin%n_occ_cls

      ! Check whether file is buffered.
      bufin = .false.
      if(ffin%buffered) bufin = .true.

      ! Allocations made to maximum block length to save time.
      loop(1:nocc_cls) = .false.
      if(.not.bufin)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls
          if(opin%formal_blk(iocc_cls))
     &         loop(iocc_cls) = .true.
          nbuff = nbuff + opin%len_op_occ(iocc_cls)
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
      endif

      iocc_loop: do iocc_cls = 1, nocc_cls
        if(loop(iocc_cls)) cycle iocc_loop

        ! Loop over Ms of annihilator string.
        idxmsa = 0
        msmax = 2
        msa_loop : do msa = msmax, -msmax, -2

          idxmsa = idxmsa+1
      
          ! Loop over Irrep of annihilator string.
          igama_loop: do igama = 1, ngam

            ! Set block length and offset.
            len_gam_ms = int(sqrt(dble(opin%
     &         len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa))))
            if(len_gam_ms.eq.0)cycle igama_loop

            ioff = opin%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)

            ! Allocate scratch space.
            allocate(scratch1(len_gam_ms,len_gam_ms),
     &           scratch2(len_gam_ms,len_gam_ms),
     &           scratch3(len_gam_ms,len_gam_ms))

            ! Get the necessary numbers from the buffers or files.
            if(.not.bufin)then
              call get_vec(ffin,buffer_in,1,nbuff)
            else
              buffer_in =>
     &             ffin%buffer(1:)
            endif

            ! Copy the operator blocks to the scratch arrays.
c            scratch1 = 0d0
c            scratch2 = 0d0
            do idx = 1,len_gam_ms
              do jdx = 1,len_gam_ms
                scratch1(idx,jdx) =
     &               buffer_in(ioff+(jdx-1)*len_gam_ms+idx)
                scratch2(idx,jdx) =
     &                 vec(ioff+(jdx-1)*len_gam_ms+idx)
              enddo
            enddo

            scratch3 = 0d0
            ! Do the multiplication. 
            do idx = 1,len_gam_ms
              do jdx = 1,len_gam_ms
                do kdx = 1,len_gam_ms
                  scratch3(idx,jdx) = scratch3(idx,jdx)+
     &                 scratch1(idx,kdx)*scratch2(kdx,jdx)
                enddo
              enddo
            enddo
            
            ! Copy the result back to the block of the right-operator.
            do idx = 1,len_gam_ms
              do jdx = 1,len_gam_ms

                vec((jdx-1)*len_gam_ms+idx+ioff) =
     &             sgn * scratch3(idx,jdx)

              enddo
            enddo

            deallocate(scratch1,scratch2,scratch3)
                      
          enddo igama_loop
          
        enddo msa_loop

      enddo iocc_loop

      ! Close the input file.
      if(closeit)
     &     call file_close_keep(ffin)

      ifree = mem_flushmark('op_vec_mult')

      return
      end
