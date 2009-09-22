      subroutine diagonal(fac,invert,mel_in,iblkin,mel_out,iblkout,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to extract the diagonal of an operator matrix with a 
*     switch to form the inverse of this if required. 
*     GWR November 2007
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_memman.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(me_list), intent(inout) ::
     &     mel_in, mel_out
      integer, intent(in) ::
     &     iblkin, iblkout
      real(8), intent(in) ::
     &     fac
      logical, intent(in) ::
     &     invert

      logical ::
     &     bufin, bufout
      integer ::
     &     len_str, idum, ifree, lblk, nblkmax,
     &     nblk, nbuff, ioffin, ioffout, idxst, idxnd, njoined,
     &     idoffin, idoffout, idxmsa,
     &     msmax, msa, msc, igama, igamc, idx, ngam,
     &     len_gam_ms, ioff, len_blk_in, ioff_blk_in, len_blk_out,
     &     ioff_blk_out
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      type(filinf), pointer ::
     &     ffin, ffout
      type(operator), pointer ::
     &     op_in, op_out

      integer, external ::
     &     idx_oplist2

      ffin  => mel_in%fhand
      ffout => mel_out%fhand
      op_in  => mel_in%op
      op_out => mel_out%op

      call file_open(ffin)
      call file_open(ffout)

      ! Check whether files are buffered.
      bufin = .false.
      if(ffin%buffered) bufin = ffin%incore(iblkin).gt.0
      bufout = .false.
      if(ffout%buffered) bufout = ffout%incore(iblkout).gt.0

      ifree = mem_setmark('symmetrise')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Find total block length.
      ioff_blk_in = mel_in%off_op_occ(iblkin)
      len_blk_in  = mel_in%len_op_occ(iblkin)
      ioff_blk_out = mel_out%off_op_occ(iblkout)
      len_blk_out  = mel_out%len_op_occ(iblkout)

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
        msc = msa + mel_out%mst
        ! Usually have mst=0 operators => Ms(c)=Ms(a)
      
        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam
          igamc = multd2h(igama,mel_out%gamt)

          len_gam_ms=int(sqrt(dble(mel_out%
     &         len_op_gmo(iblkout)%gam_ms(igama,idxmsa))))

          ioff=mel_out%off_op_gmo(iblkout)%gam_ms(igama,idxmsa)-
     &         ioff_blk_out

          idx_loop: do idx = 1,len_gam_ms
            
            if(.not.invert)then
              buffer_out((idx-1)*len_gam_ms+idx+ioff) = 
     &             fac * buffer_in((idx-1)*len_gam_ms+idx+ioff)
            else
              if(buffer_in((idx-1)*len_gam_ms+idx+ioff).ge.1d-10)then
                buffer_out((idx-1)*len_gam_ms+idx+ioff) = 
     &               fac / buffer_in((idx-1)*len_gam_ms+idx+ioff)
              endif  
            endif            
          enddo idx_loop
          
        enddo igama_loop
          
      enddo msa_loop

      if(.not.bufout)then
        call put_vec(ffout,buffer_out,ioff_blk_out+1,
     &       ioff_blk_out+len_blk_out)
      endif  

      call file_close_keep(ffin)
      call file_close_keep(ffout)

      ifree = mem_flushmark('symmetrise')

      return
      end
