      subroutine symmetrise(fac,me_in,me_out,
     &     xnorm2_blk,
     &     op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to symmetrise an operator matrix:
*
*     O^{symm}_{ij} = fac*(O_{ij}+O_{ji}) = O^{symm}_{ji}
*
*     me_in are the input, asymmetric matrices
*     me_out represent the output, symmetric matrices.
*     It is assumed that the two operators have the same shape.
*     GWR November 2007
*
*     generalised and extended May 2008, andreas
*     
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 000

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(operator_info), intent(in) ::
     &     op_info
      type(me_list), intent(inout) ::
     &     me_in, me_out
      real(8), intent(in) ::
     &     fac
      real(8), intent(out) ::
     &     xnorm2_blk(*)

      logical ::
     &     bufin, bufout, open_close_in, open_close_out, same
      integer ::
     &     nocc_cls, njoined,
     &     ifree, nblk, nbuff, nbuff2, iblk, jblk, ioff_blk,
     &     idx, jdx, ngam, ioff, joff, lenblk
      real(8) ::
     &     fac_off, fac_dia, value

      real(8), pointer ::
     &     buffer_in(:), buffer_out(:), buffer_in2(:), buffer_out2(:)

      type(operator), pointer ::
     &     op_in, op_out
      type(filinf), pointer ::
     &     ffin, ffout

      real(8), external ::
     &     ddot
      logical, external ::
     &     occ_is_diag_blk, iocc_equal_n
      integer, external ::
     &     iblk_occ


      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'symmetrise')
        write(lulog,*) 'IN:  ',trim(me_in%label)
        write(lulog,*) 'OUT: ',trim(me_out%label)
        write(lulog,*) 'elements on input:'
      end if
      if (ntest.ge.1000) then
        call wrt_mel_file(lulog,5,
     &       me_in,
     &       1,me_in%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      ffin  => me_in%fhand
      ffout => me_out%fhand

      ! in and out on same file?
      same = associated(ffin,ffout)
      
      if (.not.same)
     &     call quit(1,'symmetrise','needed?')

      op_in  => me_in%op
      op_out => me_out%op

      nocc_cls = op_in%n_occ_cls
      njoined  = op_in%njoined

      if (nocc_cls.ne.op_out%n_occ_cls.or.
     &    njoined.ne.op_out%njoined)
     &     call quit(1,'symmetrise',
     &     'Input and output list do not have identical operators')
 
      open_close_in  = ffin%unit.le.0
      open_close_out = ffout%unit.le.0

      if (open_close_in ) call file_open(ffin )
      if (open_close_out.and..not.same) call file_open(ffout)

      ! Check whether files are buffered.
      bufin = .false.
      if(ffin%buffered) bufin = .true. 
      bufout = .false.
      if(ffout%buffered) bufout = .true.

      ifree = mem_setmark('symmetrise')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Allocations made to maximum block length
      if(.not.bufin)then
        nbuff = 0
        nbuff2 = 0
        do iblk = 1, nocc_cls 
          ioff_blk = (iblk-1)*njoined
          if(op_in%formal_blk(iblk))
     &         cycle
          if (occ_is_diag_blk(op_in%ihpvca_occ(1,1,ioff_blk+1),njoined))
     &         then
            nbuff = max(nbuff,me_in%len_op_occ(iblk))
          else
            ! need two buffers in case of off-diagonal blocks:
            nbuff = max(nbuff,me_in%len_op_occ(iblk)) 
            nbuff2 = max(nbuff2,me_in%len_op_occ(iblk)) 
          end if
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
        if (nbuff2.gt.0)
     &       ifree = mem_alloc_real(buffer_in2,nbuff,'buffer_in2')
      else
        if (nocc_cls.gt.1)
     &     call quit(1,'symmetrise','buffered file: not yet')
        buffer_in => ffin%buffer(1:)
      endif

      if(.not.bufout.and..not.same)then
        nbuff = 0
        nbuff2 = 0
        do iblk = 1, nocc_cls 
          ioff_blk = (iblk-1)*njoined
          if(op_out%formal_blk(iblk))
     &         cycle
          if (occ_is_diag_blk(op_out%ihpvca_occ(1,1,ioff_blk+1),
     &         njoined))
     &    then
            nbuff = max(nbuff,me_out%len_op_occ(iblk))
          else
            ! need two buffers in case of off-diagonal blocks:
            nbuff = max(nbuff,me_out%len_op_occ(iblk)) 
            nbuff2 = max(nbuff2,me_out%len_op_occ(iblk)) 
          end if
        enddo
        ifree = mem_alloc_real(buffer_out,nbuff,'buffer_out')
        if (nbuff2.gt.0)
     &       ifree = mem_alloc_real(buffer_out2,nbuff,'buffer_out2')
      else if (same) then
        buffer_out => buffer_in
        if (nbuff2.gt.0) buffer_out2 => buffer_in2
      else
        call quit(1,'symmetrise','buffered file 2: needed?')
        buffer_out => ffout%buffer(1:)
      endif

      ! factor for off-diagonal and diagonal
      fac_off = 0.5d0*fac
      fac_dia = fac

      if (ntest.ge.100) then
        write(lulog,*) 'size of buffer: ',nbuff,nbuff2
      end if

      ! quick exit for (rare) special case
      if (bufin.and.bufout.and.nocc_cls.eq.1
     &         .and.me_out%len_op_occ(1).eq.1) then
         buffer_out(1) = fac*buffer_in(1)
         ifree = mem_flushmark('symmetrise')
         return
      end if

      ! Loop over occupation classes.
      iocc_loop: do iblk = 1, nocc_cls 

        ioff_blk = (iblk-1)*njoined

        if (ntest.ge.100)
     &       write(lulog,*) 'iblk: ',iblk
        if (ntest.ge.100)
     &       call wrt_occ_n(lulog,op_in%ihpvca_occ(1,1,ioff_blk+1),
     &           njoined)

        if (.not.
     &       occ_is_diag_blk(op_in%ihpvca_occ(1,1,ioff_blk+1),njoined))
     &  then
          jblk = iblk_occ(op_in%ihpvca_occ(1,1,ioff_blk+1),
     &                     .true.,op_in,op_in%blk_version(iblk))
          if (jblk.le.0) then
            call wrt_occ_n(lulog,op_in%ihpvca_occ(1,1,ioff_blk+1),
     &           njoined)
            call quit(1,'symmetrise',
     &           'no adjoint found for this occupation')
          end if
            
          if (jblk.lt.iblk) cycle

          ioff = me_in%off_op_occ(iblk)
          joff = me_in%off_op_occ(jblk)
          lenblk = me_in%len_op_occ(iblk)

          if (lenblk.eq.0) cycle

          if (ntest.ge.100) then
            write(lulog,*) 'iblk, jblk: ',iblk, jblk
            write(lulog,*) 'ioff, joff, lenblk: ',ioff, joff, lenblk
            write(lulog,*) 'off-diagonal case'
          end if

          call get_vec(ffin,buffer_in,ioff+1,ioff+lenblk)
          call get_vec(ffin,buffer_in2,joff+1,joff+lenblk)
          call symmetrise_blk1blk2(buffer_out,buffer_out2,
     &                        buffer_in, buffer_in2,
     &                         fac_off,fac_dia,
     &                        me_in,iblk,jblk,
     &                        str_info,orb_info)

          call put_vec(ffout,buffer_out,ioff+1,ioff+lenblk)
          call put_vec(ffout,buffer_out2,joff+1,joff+lenblk)

          ! update norm^2
          xnorm2_blk(iblk) = ddot(lenblk,
     &         buffer_out,1,buffer_out,1)
          xnorm2_blk(jblk) = xnorm2_blk(iblk)

        else

          ioff = me_in%off_op_occ(iblk)
          lenblk = me_in%len_op_occ(iblk)

          if (ntest.ge.100) then
            write(lulog,*) 'iblk: ',iblk
            write(lulog,*) 'ioff, lenblk: ',ioff, lenblk
            write(lulog,*) 'diagonal case'
          end if

          call get_vec(ffin,buffer_in,ioff+1,ioff+lenblk)

          call symmetrise_blk1blk2(buffer_out,buffer_out,
     &                        buffer_in, buffer_in,fac_off,fac_dia,
     &                        me_in,iblk,iblk,
     &                        str_info,orb_info)

          call put_vec(ffout,buffer_out,ioff+1,ioff+lenblk)

          ! update norm^2
          xnorm2_blk(iblk) = ddot(lenblk,
     &         buffer_out,1,buffer_out,1)

        end if

      end do iocc_loop

      if (open_close_in ) call file_close_keep(ffin)
      if (open_close_out.and..not.same) call file_close_keep(ffout)

      if (ntest.ge.1000) then
        write(lulog,*) 'elements on output:'
        call wrt_mel_file(lulog,5,
     &       me_in,
     &       1,me_in%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      ifree = mem_flushmark('symmetrise')

      return
      end
