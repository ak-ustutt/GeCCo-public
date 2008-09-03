*----------------------------------------------------------------------*
      subroutine sym_ab_list(fac,me_in,me_out,
     &     xnorm2_blk,update_norm,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to symmetrize ME list wrt. MS-combinations
*     the symmetrization type (+1/-1) is given by me_out%absym
*     use fac=0.5d0 if both non-flipped and flipped parts have
*     been set before or if one of both has been pre-multiplied
*     by the appropriate symmetry factor
*     one might use fac=1d0, if flipped parts were omitted
*     previously
*
*     july 2008, andreas
*     
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'mdef_operator_info.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(operator_info), intent(in) ::
     &     op_info
      type(me_list), intent(inout) ::
     &     me_in, me_out
      real(8), intent(in) ::
     &     fac
      real(8), intent(out) ::
     &     xnorm2_blk(*)
      logical, intent(in) ::
     &     update_norm

      logical ::
     &     bufin, bufout, open_close_in, open_close_out, same
      integer ::
     &     nocc_cls, njoined, ioff_blk,
     &     ifree, nblk, nbuff, iblk, idisc_off_in, idisc_off_out,
     &     idx, ngam, ioff, lenblk
      real(8) ::
     &     fac_rel, fac_pre

      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      type(operator), pointer ::
     &     op_in, op_out
      type(filinf), pointer ::
     &     ffin, ffout

      real(8), external ::
     &     ddot

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'sym_ab_list')
        write(luout,*) 'IN:  ',trim(me_in%label)
        write(luout,*) 'OUT: ',trim(me_out%label)
      end if
      if (ntest.ge.1000) then
        write(luout,*) 'elements on input:'
        call wrt_mel_file(luout,5,
     &       me_in,
     &       1,me_in%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      if (me_out%absym.ne.+1.and.me_out%absym.ne.-1)
     &     call quit(1,'sym_ab_list',
     &     'absym.ne.+1.and.absym.ne.-1, list='//trim(me_out%label))
      ! factor for off-diagonal and diagonal
      fac_rel = dble(me_out%absym)
      fac_pre = fac

      ffin  => me_in%fhand
      ffout => me_out%fhand

      ! in and out on same file?
      same = associated(ffin,ffout)
      
      if (.not.same)
     &     call quit(1,'sym_ab_list','needed?')

      op_in  => me_in%op
      op_out => me_out%op

      nocc_cls = op_in%n_occ_cls
      njoined  = op_in%njoined

      if (nocc_cls.ne.op_out%n_occ_cls.or.
     &    njoined.ne.op_out%njoined)
     &     call quit(1,'sym_ab_list',
     &     'Input and output list do not have identical operators')
 
      open_close_in  = ffin%unit.le.0
      open_close_out = ffout%unit.le.0

      if (open_close_in ) call file_open(ffin )
      if (open_close_out.and..not.same) call file_open(ffout)

      ! Check whether files are buffered.
      bufin = ffin%buffered
      bufout = ffout%buffered

      ifree = mem_setmark('sym_ab_list')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Allocations made to maximum block length
      if (.not.bufin) then
        nbuff = 0
        do iblk = 1, nocc_cls 
          ioff_blk = (iblk-1)*njoined
          if(op_in%formal_blk(iblk))
     &         cycle
          nbuff = max(nbuff,me_in%len_op_occ(iblk))
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
      else
        call quit(1,'sym_ab_list','buffered file: not yet')
        buffer_in => ffin%buffer(1:)
      endif

      if (.not.bufout.and..not.same) then
        nbuff = 0
        do iblk = 1, nocc_cls 
          ioff_blk = (iblk-1)*njoined
          if(op_out%formal_blk(iblk))
     &         cycle
          nbuff = max(nbuff,me_out%len_op_occ(iblk))
        enddo
        ifree = mem_alloc_real(buffer_out,nbuff,'buffer_out')
      else if (same) then
        buffer_out => buffer_in
      else
        call quit(1,'sym_ab_list','buffered file: not yet')
        buffer_out => ffout%buffer(1:)
      endif

      if (ntest.ge.100) then
        write(luout,*) 'size of buffer: ',nbuff
      end if

      idisc_off_in  = ffin%length_of_record*(ffin%current_record-1)
      idisc_off_out = ffout%length_of_record*(ffout%current_record-1)

      ! Loop over occupation classes.
      iocc_loop: do iblk = 1, nocc_cls 

        if(op_out%formal_blk(iblk)) cycle

        ioff_blk = (iblk-1)*njoined

        ioff = me_in%off_op_occ(iblk)+idisc_off_in
        lenblk = me_in%len_op_occ(iblk)

        if (ntest.ge.100)
     &       write(luout,*) 'iblk, ioff, lenblk: ',iblk, ioff, lenblk

        call get_vec(ffin,buffer_in,ioff+1,ioff+lenblk)

        call sym_ab_blk(buffer_out,
     &                  buffer_in, fac_pre, fac_rel,
     &                  me_in,iblk,
     &                  str_info,strmap_info,orb_info)

        ioff = me_out%off_op_occ(iblk)+idisc_off_out
        call put_vec(ffout,buffer_out,ioff+1,ioff+lenblk)

        ! update norm^2
        if (update_norm) then
          xnorm2_blk(iblk) = ddot(lenblk,
     &         buffer_out,1,buffer_out,1)
        end if

      end do iocc_loop

      if (open_close_in ) call file_close_keep(ffin)
      if (open_close_out.and..not.same) call file_close_keep(ffout)

      if (ntest.ge.1000) then
        write(luout,*) 'elements on output:'
        call wrt_mel_file(luout,5,
     &       me_out,
     &       1,me_in%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      ifree = mem_flushmark('sym_ab_list')

      return
      end
