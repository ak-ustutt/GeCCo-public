      subroutine invert(fac,ffin,op_in,iblkin,ffout,op_out,iblkout,
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
      include 'multd2h.h'

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
     &     iblkin, iblkout
      real(8), intent(in) ::
     &     fac

      logical ::
     &     bufin, bufout
      integer ::
     &     len_str, idum, ifree, lblk, nblkmax,
     &     nblk, nbuff, ioffin, ioffout, idxst, idxnd, njoined,
     &     idoffin, idoffout, idxmsa,
     &     msmax, msa, msc, igama, igamc, idx, jdx, ngam,
     &     len_gam_ms, ioff, len_blk_in, ioff_blk_in, len_blk_out,
     &     ioff_blk_out
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:), scratch(:,:)

      integer, external ::
     &     idx_oplist2

      if(ffout%unit.le.0)then
        call file_open(ffout)
      endif

      ! Check whether files are buffered.
      bufin = .false.
      if(ffin%buffered) bufin = ffin%incore(iblkin).gt.0
      bufout = .false.
      if(ffout%buffered) bufout = ffout%incore(iblkout).gt.0

      ifree = mem_setmark('symmetrise')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Find total block length.
      ioff_blk_in = op_in%off_op_occ(iblkin)
      len_blk_in  = op_in%len_op_occ(iblkin)
      ioff_blk_out = op_out%off_op_occ(iblkout)
      len_blk_out  = op_out%len_op_occ(iblkout)

      ! Allocations made to maximum block length to save time.
      if(.not.bufin)then
        nbuff = len_blk_in
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
        call get_vec(ffin,buffer_in,ioff_blk_in+1,
     &       ioff_blk_in+len_blk_in)
      else
        if(ntest.ge.100)
     &       write(luout,*)'Symmetrise: input not incore'
        buffer_in => ffin%buffer(ioff_blk_in+1:)
      endif

      if(.not.bufout)then
        nbuff=len_blk_out
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
        buffer_out(1:nbuff) = 0d0
      else
        if(ntest.ge.100)
     &       write(luout,*)'Symmetrise: output not incore'
        buffer_out => ffout%buffer(ioff_blk_out+1:)
      endif

      ! Loop over Ms of annihilator string.
      idxmsa = 0
      msmax = 2
      msa_loop : do msa = msmax, -msmax, -2

        idxmsa = idxmsa+1
        msc = msa + op_out%mst
        ! Usually have mst=0 operators => Ms(c)=Ms(a)
      
        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam
          igamc = multd2h(igama,op_out%gamt)

          len_gam_ms = int(sqrt(dble(op_out%
     &         len_op_gmo(iblkout)%gam_ms(igama,idxmsa))))

          ioff = op_out%off_op_gmo(iblkout)%gam_ms(igama,idxmsa)-
     &         ioff_blk_out

          allocate(scratch(len_gam_ms,len_gam_ms))
          do idx = 1,len_gam_ms
           do jdx = 1,len_gam_ms

              scratch(idx,jdx) =
     &            buffer_in(ioff+(idx-1)*len_gam_ms+jdx)
            
            enddo
          enddo

          call gaussj(scratch,len_gam_ms,len_gam_ms)

          idx_loop: do idx = 1,len_gam_ms
            jdx_loop: do jdx = 1,len_gam_ms

             buffer_out((idx-1)*len_gam_ms+jdx+ioff) =
     &           scratch(idx,jdx)

             enddo jdx_loop
           enddo idx_loop

           deallocate(scratch)
                      
        enddo igama_loop
          
      enddo msa_loop

      if(.not.bufout)then
        call put_vec(ffout,buffer_out,ioff_blk_out+1,
     &       ioff_blk_out+len_blk_out)
      endif  

      call file_close_keep(ffout)

      ifree = mem_flushmark('symmetrise')

      return
      end
