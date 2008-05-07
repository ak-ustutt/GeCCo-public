      subroutine add_unity(fac,mel_out,iblkout,orb_info)
*----------------------------------------------------------------------*
*
*     Add (antisymmetrised) unit operator
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      real(8), intent(in)::
     &     fac
      type(me_list), intent(inout) ::
     &     mel_out
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     iblkout

      logical ::
     &     bufout,closeit
      integer ::
     &     len_str, idum, ifree, lblk, nblkmax, nblk, nbuff,
     &     ioffin, ioffout, idxst, idxnd, njoined, join_off,
     &     idoffin, idoffout, idxmsa, msmax, msa, msc, igama, igamc,
     &     idx, ngam, len_gam_ms, ioff, len_blk, ioff_blk
      integer ::
     &     iocc(ngastp,2), opout_temp(ngastp,2)
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      type(filinf), pointer ::
     &     ffout
      type(operator), pointer ::
     &     opout

      logical, external ::
     &     iocc_equal, irestr_equal, occ_is_diag_blk

      ffout => mel_out%fhand
      opout => mel_out%op

      if (ntest.ge.100) then
        write(luout,*) '=========================='
        write(luout,*) ' add_unity messing around'
        write(luout,*) '=========================='
        write(luout,*) ' fac = ',fac
        write(luout,*) ' list = ',trim(mel_out%label)
        write(luout,*) ' ffout: ',trim(ffout%name),
     &                   ' rec: ',ffout%current_record
        write(luout,*) ' opout: ',opout%name(1:len_trim(opout%name)),
     &       '  block: ',iblkout
      end if

      closeit = .false.
      njoined = opout%njoined

      join_off=(iblkout-1)*njoined
      ! check whether the out operator is a diagonal block:
      if (.not.
     &     occ_is_diag_blk(opout%ihpvca_occ(1,1,join_off+1),njoined))
     &     return ! else we just tacitly return

      bufout = .false.
      if(ffout%buffered) bufout = ffout%incore(iblkout).gt.0

      ifree = mem_setmark('add_unity')

      ngam = orb_info%nsym

      ! Find total block length.
      ioff_blk = mel_out%off_op_occ(iblkout)
      len_blk  = mel_out%len_op_occ(iblkout)

      if(len_blk.gt.ifree)
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
        if(ffout%unit.le.0)then
          call file_open(ffout)
          closeit = .true.
        endif
        call get_vec(ffout,buffer_out,ioff_blk+1,ioff_blk+len_blk)
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
     &         ioff_blk
          idx_loop: do idx =1,len_gam_ms

            buffer_out((idx-1)*len_gam_ms+idx+ioff) = 1d0*fac+
     &         buffer_out((idx-1)*len_gam_ms+idx+ioff)

          enddo idx_loop
          
        enddo igama_loop
          
      enddo msa_loop

c dbg
      buffer_out(6) = 0.75d0
      buffer_out(7) = -0.25d0
      buffer_out(8) = -0.25d0
      buffer_out(9) = 0.75d0

c      buffer_out(37) = 0.75d0
c      buffer_out(45) = 0.75d0
c      buffer_out(78) = 0.75d0
c      buffer_out(85) = 0.75d0
c      buffer_out(92) = 0.75d0
c      buffer_out(99) = 0.75d0
c      buffer_out(106) = 0.75d0
c      buffer_out(113) = 0.75d0
c      buffer_out(114) = 0.75d0
c      buffer_out(121) = 0.75d0
c      buffer_out(128) = 0.75d0
c      buffer_out(135) = 0.75d0
c      buffer_out(142) = 0.75d0
c      buffer_out(149) = 0.75d0
c      buffer_out(150) = 0.75d0
c      buffer_out(157) = 0.75d0
c      buffer_out(164) = 0.75d0
c      buffer_out(171) = 0.75d0
c      buffer_out(178) = 0.75d0
c      buffer_out(185) = 0.75d0

c      buffer_out(38) = -0.25d0
c      buffer_out(44) = -0.25d0
c      buffer_out(80) = -0.25d0
c      buffer_out(87) = -0.25d0
c      buffer_out(90) = -0.25d0
c      buffer_out(97) = -0.25d0
c      buffer_out(107) = -0.25d0
c      buffer_out(112) = -0.25d0
c      buffer_out(116) = -0.25d0
c      buffer_out(123) = -0.25d0
c      buffer_out(126) = -0.25d0
c      buffer_out(133) = -0.25d0
c      buffer_out(143) = -0.25d0
c      buffer_out(148) = -0.25d0
c      buffer_out(151) = -0.25d0
c      buffer_out(156) = -0.25d0
c      buffer_out(166) = -0.25d0
c      buffer_out(173) = -0.25d0
c      buffer_out(176) = -0.25d0
c      buffer_out(183) = -0.25d0
c dbg

      if(.not.bufout)then
        call put_vec(ffout,buffer_out,ioff_blk+1,ioff_blk+len_blk)
      endif  
      if(closeit)
     &     call file_close_keep(ffout)

      ifree = mem_flushmark('add_unity')

      return
      end
