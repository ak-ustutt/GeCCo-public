      subroutine op_mult2(sgn,ffin,ops_in,nops_in,ffout,op_out,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to multiply the matrices of nops_in operators, with the 
*     product operator, X:
*     i.e. A.B.C.... = X
*     Should be called by mult_op to ensure that operator files are 
*     opened and closed properly, and that they have the same shapes etc.
*     GWR November/December 2007
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
     &     ntest = 100

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(operator_array), intent(in) ::
     &     ops_in(nops_in)
      type(operator), intent(in) ::
     &     op_out
      type(file_array), intent(in) ::
     &     ffin(nops_in)
      type(filinf), intent(inout) ::
     &     ffout
      integer, intent(in) ::
     &     nops_in
      real(8), intent(in) ::
     &     sgn

      logical ::
     &     bufout
      logical ::
     &     bufin(nops_in), loop(op_out%n_occ_cls)
      integer ::
     &     iops, nbuffered, idxbuf, iocc_cls ,nocc_cls, ibufloop,
     &     ifree, nbuff, idxmsa, msmax, msa, igama, idx, jdx, kdx, ngam,
     &     len_gam_ms, ioff
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:), scratch1(:,:), scratch2(:,:),
     &     scratch3(:,:)

      if(ntest.ge.100)then
        write(luout,*)'=================='
        write(luout,*)' Entered op_mult. '
        write(luout,*)'=================='
        write(luout,*)'Output: ',trim(op_out%name)
      endif

      ifree = mem_setmark('op_mult')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym
      nocc_cls = op_out%n_occ_cls

      ! Check whether files are buffered.
      bufout = .false.
      if(ffout%buffered) bufout = .true.

      bufin(1:nops_in) = .false.
      nbuffered = 0
      idxbuf = 0
      do iops = 1, nops_in
        if(ffin(iops)%fhand%buffered) then
          bufin(iops) = .true.
          nbuffered = nbuffered + 1
        else
          if(idxbuf.eq.0) idxbuf = iops
        endif
      enddo

      ! Allocations made to maximum block length to save time.
      loop(1:nocc_cls) = .false.
      if(nbuffered.eq.0)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls
          if(ops_in(idxbuf)%op%formal_blk(iocc_cls))
     &         loop(iocc_cls) = .true.
          nbuff = nbuff + ops_in(idxbuf)%op%len_op_occ(iocc_cls)
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
      endif

      if(.not.bufout)then
        nbuff=0
        do iocc_cls = 1, nocc_cls
          ! Repeat this test in case input is non-buffered.
          if(op_out%formal_blk(iocc_cls))
     &         loop(iocc_cls) = .true.
          nbuff = nbuff + op_out%len_op_occ(iocc_cls)
        enddo
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
        buffer_out(1:nbuff) = 0d0
      else
        if(ntest.ge.100)
     &       write(luout,*)'Op_mult: output not incore'
        buffer_out => ffout%buffer(1:)
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
            len_gam_ms = int(sqrt(dble(op_out%
     &         len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa))))
            if(len_gam_ms.eq.0)cycle igama_loop

            ioff = op_out%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)

            ! Allocate scratch space.
            allocate(scratch1(len_gam_ms,len_gam_ms),
     &           scratch2(len_gam_ms,len_gam_ms),
     &           scratch3(len_gam_ms,len_gam_ms))

            iops_loop: do iops = 1, nops_in-1

              do ibufloop = 0,1

                ! Get the necessary numbers from the buffers or files.
                if(ibufloop.eq.1.or.iops.eq.1)then
                  if(.not.bufin(iops+ibufloop))then
                    call get_vec(ffin(iops+ibufloop)%fhand,
     &                   buffer_in,ioff+1,nbuff)
                  else
                    buffer_in =>
     &                   ffin(iops+ibufloop)%fhand%buffer(ioff+1:)
                  endif
                endif

                ! Copy the operator blocks to the scratch arrays.
c                scratch1 = 0d0
c                scratch2 = 0d0
                do idx = 1,len_gam_ms
                  do jdx = 1,len_gam_ms
                    if(ibufloop.eq.0)then
                      if(iops.eq.1)then
                        scratch1(idx,jdx) =
     &                       buffer_in((jdx-1)*len_gam_ms+idx)
                      else
                        scratch1(idx,jdx) = scratch3(idx,jdx)
                      endif
                    else
                      scratch2(idx,jdx) =
     &                     buffer_in((jdx-1)*len_gam_ms+idx)
                    endif
                  enddo
                enddo

              enddo

              scratch3 = 0d0
              ! Do the multiplication. 
              do idx = 1,len_gam_ms
                do jdx = 1,len_gam_ms
                  do kdx = 1,len_gam_ms
                    scratch3(idx,jdx) = scratch3(idx,jdx)+
     &                   sgn*scratch1(idx,kdx)*scratch2(kdx,jdx)
                  enddo
                enddo
              enddo

            enddo iops_loop
          
            ! Copy the result back to the block of the right-operator.
            do idx = 1,len_gam_ms
              do jdx = 1,len_gam_ms

                buffer_out((jdx-1)*len_gam_ms+idx+ioff) =
     &             scratch3(idx,jdx)

              enddo
            enddo

            deallocate(scratch1,scratch2,scratch3)
                      
          enddo igama_loop
          
        enddo msa_loop

      enddo iocc_loop

      ! Save the results if necessary.
      if(.not.bufout)then
        call put_vec(ffout,buffer_out,1,nbuff)
      endif  

      ifree = mem_flushmark('op_mult')

      return
      end
