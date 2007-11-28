      subroutine add_unity(fac,ffout,opout,iblkout,orb_info)
*----------------------------------------------------------------------*
*
*     Add antisymmetrised delta.
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 100

      real(8), intent(in)::
     &     fac
      type(filinf), intent(inout) ::
     &     ffout
      type(operator), intent(in) ::
     &     opout
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     iblkout

      logical ::
     &     bufout
      integer ::
     &     len_str, idum, ifree, lblk, nblkmax, nblk, nbuff,
     &     ioffin, ioffout, idxst, idxnd, njoined, join_off,
     &     idoffin, idoffout, idxmsa, msmax, msa, msc, igama, igamc,
     &     idx, ngam, len_gam_ms, ioff, len_blk, ioff_blk
      integer ::
     &     iocc(ngastp,2), opout_temp(ngastp,2)
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      logical, external ::
     &     iocc_equal, irestr_equal

      if (ntest.ge.100) then
        write(luout,*) '=========================='
        write(luout,*) ' add_unity messing around'
        write(luout,*) '=========================='
        write(luout,*) ' fac = ',fac
        write(luout,*) ' ffout: ',trim(ffout%name),
     &                   ' rec: ',ffout%current_record
        write(luout,*) ' opout: ',opout%name(1:len_trim(opout%name)),
     &       '  block: ',iblkout
      end if

      ! Set array with the necessary shape for the delta operator, then 
      ! compare it to the out operator.
      iocc(1:ngastp,1:2)=0
      iocc(ihole,1:2)=2

      njoined = opout%njoined
      join_off=(iblkout-1)*njoined
      opout_temp(1:ngastp,1:2)=0
      do idx=1,njoined
        opout_temp(1:ngastp,1:2) = opout_temp(1:ngastp,1:2)+
     &       opout%ihpvca_occ(1:ngastp,1:2,join_off+idx)
      enddo

c dbg  
c      if (.not.iocc_equal(iocc,.false.,opout_temp,opout%dagger)) then
c        call quit(1,'add_unity','output incompatible')
c      endif  
c dbg

      bufout = .false.
      if(ffout%buffered) bufout = ffout%incore(iblkout).gt.0

      ifree = mem_setmark('add_unity')

      ngam = orb_info%nsym

      ! Find total block length.
      ioff_blk = opout%off_op_occ(iblkout)
      len_blk  = opout%len_op_occ(iblkout)
!      len_str = int(sqrt(dble(opout%len_op_occ(iblkout))))
!     str_gmo(iblkout)%gam_ms(igama,idxmsa))))
      
      if(len_str.gt.ifree)
     &     call quit(1,'add_unity','insufficient space')
      
      ! Allocations made to maximum block length to save time.
      if(.not.bufout)then
        nbuff=len_blk
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
      else
        if(ntest.ge.100)
     &       write(luout,*)'Add_unity: Not incore'
        buffer_out => ffout%buffer(ioff_blk+1:)
      endif

      if(.not.bufout)then
        call get_vec(ffout,buffer_out,ioff_blk+1,ioff_blk+len_blk)
      endif  

      ! Loop over Ms of annihilator string.
      idxmsa = 0
      msmax = 2
      msa_loop : do msa = msmax, -msmax, -2

        idxmsa = idxmsa+1
        msc = msa + opout%mst
        ! Usually have mst=0 operators => Ms(c)=Ms(a)
      
        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam
          igamc = multd2h(igama,opout%gamt)

          len_gam_ms=int(sqrt(dble(opout%
     &         len_op_gmo(iblkout)%gam_ms(igama,idxmsa))))

          ioff=opout%off_op_gmo(iblkout)%gam_ms(igama,idxmsa)-
     &         ioff_blk
          idx_loop: do idx =1,len_gam_ms

            buffer_out((idx-1)*len_gam_ms+idx+ioff) = 1d0*fac+
     &         buffer_out((idx-1)*len_gam_ms+idx+ioff)

          enddo idx_loop
          
        enddo igama_loop
          
      enddo msa_loop

      if(.not.bufout)then
        call put_vec(ffout,buffer_out,ioff_blk+1,ioff_blk+len_blk)
      endif  

      ifree = mem_flushmark('add_unity')

      return
      end
