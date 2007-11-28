      subroutine op_mult(sgn,ffleft,op_left,iblkleft,
     &     ffrig,op_rig,iblkrig,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to multiply the matrices of two operators, with the 
*     product operator replacing the right-hand operator:
*     i.e. A.B = B^new
*     Currently works for a single block of an operator.
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
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 100

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator), intent(in) ::
     &     op_left, op_rig
      type(filinf), intent(in) ::
     &     ffleft
      type(filinf), intent(inout) ::
     &     ffrig
      integer, intent(in) ::
     &     iblkleft, iblkrig
      real(8), intent(in) ::
     &     sgn

      logical ::
     &     bufleft, bufrig
      integer ::
     &     ifree, lblk, nblk, nbuff, ioff_left, ioff_rig, idxmsa,
     &     msmax, msa, igama, idx, jdx, kdx, ngam,
     &     len_gam_ms, ioff, len_blk_left, ioff_blk_left, len_blk_rig,
     &     ioff_blk_rig
      
      real(8), pointer ::
     &     buffer_left(:), buffer_rig(:), scratch1(:,:), scratch2(:,:),
     &     scratch3(:,:)

      if(ntest.ge.100)then
        write(luout,*)'=================='
        write(luout,*)' Entered op_mult. '
        write(luout,*)'=================='
      endif

      ! If the files are not already open, open them.
      if(ffrig%unit.le.0)then
        call file_open(ffrig)
      endif
      if(ffleft%unit.le.0)then
        call file_open(ffleft)
      endif

      ! Check whether files are buffered.
      bufleft = .false.
      if(ffleft%buffered) bufleft = ffleft%incore(iblkleft).gt.0
      bufrig = .false.
      if(ffrig%buffered) bufrig = ffrig%incore(iblkrig).gt.0

      ifree = mem_setmark('op_mult')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Find total block length.
      ioff_blk_left = op_left%off_op_occ(iblkleft)
      len_blk_left  = op_left%len_op_occ(iblkleft)
      ioff_blk_rig = op_rig%off_op_occ(iblkrig)
      len_blk_rig  = op_rig%len_op_occ(iblkrig)

      ! Allocations made to maximum block length to save time.
      if(.not.bufleft)then
        nbuff = len_blk_left
        ifree = mem_alloc_real(buffer_left,nbuff,'buffer_left')
        call get_vec(ffleft,buffer_left,ioff_blk_left+1,
     &       ioff_blk_left+len_blk_left)
      else
        if(ntest.ge.100)
     &       write(luout,*)'Op_mult: input not incore'
        buffer_left => ffleft%buffer(ioff_blk_left+1:)
      endif

      ! Extract the numerical values on the block for both operators.
      if(.not.bufrig)then
        nbuff=len_blk_rig
        ifree= mem_alloc_real(buffer_rig,nbuff,'buffer_rig')
        call get_vec(ffrig,buffer_rig,ioff_blk_rig+1,
     &       ioff_blk_rig+len_blk_rig)
      else
        if(ntest.ge.100)
     &       write(luout,*)'Op_mult: output not incore'
        buffer_rig => ffrig%buffer(ioff_blk_rig+1:)
      endif

      ! Loop over Ms of annihilator string.
      idxmsa = 0
      msmax = 2
      msa_loop : do msa = msmax, -msmax, -2

        idxmsa = idxmsa+1
        ! Usually have mst=0 operators => Ms(c)=Ms(a)
      
        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam

          len_gam_ms = int(sqrt(dble(op_rig%
     &         len_op_gmo(iblkrig)%gam_ms(igama,idxmsa))))

          if(len_gam_ms.ne.int(sqrt(dble(op_left%
     &         len_op_gmo(iblkleft)%gam_ms(igama,idxmsa)))))
     &         call quit(1,'op_mult','odd block lengths')

          ioff_rig = op_rig%off_op_gmo(iblkrig)%gam_ms(igama,idxmsa)-
     &         ioff_blk_rig
          ioff_left = op_left%off_op_gmo(iblkleft)%gam_ms(igama,idxmsa)-
     &         ioff_blk_left

          allocate(scratch1(len_gam_ms,len_gam_ms),
     &         scratch2(len_gam_ms,len_gam_ms),
     &         scratch3(len_gam_ms,len_gam_ms))

          ! Copy the operator blocks to the scratch arrays.
          do idx = 1,len_gam_ms
           do jdx = 1,len_gam_ms

              scratch1(idx,jdx) =
     &            buffer_left(ioff_left+(jdx-1)*len_gam_ms+idx)
              scratch2(idx,jdx) =
     &            buffer_rig(ioff_rig+(jdx-1)*len_gam_ms+idx)
                        
            enddo
          enddo

          ! Do the multiplication. 
          scratch3 = 0d0
          do idx = 1,len_gam_ms
            do jdx = 1,len_gam_ms
              do kdx = 1,len_gam_ms
                scratch3(idx,jdx) = scratch3(idx,jdx)+
     &               sgn*scratch1(idx,kdx)*scratch2(kdx,jdx)
              enddo
            enddo
          enddo

          ! Cpy the result back to the block of the right-operator.
          idx_loop: do idx = 1,len_gam_ms
            jdx_loop: do jdx = 1,len_gam_ms

             buffer_rig((jdx-1)*len_gam_ms+idx+ioff_rig) =
     &           scratch3(idx,jdx)

             enddo jdx_loop
           enddo idx_loop

           deallocate(scratch1,scratch2,scratch3)
                      
        enddo igama_loop
          
      enddo msa_loop

      ! Save the results if necessary.
      if(.not.bufrig)then
        call put_vec(ffrig,buffer_rig,ioff_blk_rig+1,
     &       ioff_blk_rig+len_blk_rig)
      endif  

      ! Close the files.
      call file_close_keep(ffleft)
      call file_close_keep(ffrig)

      ifree = mem_flushmark('op_mult')

      return
      end
