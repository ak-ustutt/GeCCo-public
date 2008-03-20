      subroutine symmetrise(fac,me_in,me_out,
     &     xnorm2_blk,
     &     op_info,orb_info)
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
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 100

      type(orbinf), intent(in) ::
     &     orb_info
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
     &     ifree, nblk, nbuff, iocc_cls, idxmsa,
     &     msmax, msa, igama, idx, jdx, ngam, len_gam_ms, ioff
      real(8) ::
     &     fac_off, fac_dia, value

      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      type(operator), pointer ::
     &     op_in, op_out
      type(filinf), pointer ::
     &     ffin, ffout

      real(8), external ::
     &     ddot
      logical, external ::
     &     occ_is_diag_blk, iocc_equal_n

      ffin  => me_in%fhand
      ffout => me_out%fhand

      ! in and out on same file?
      same = associated(ffin,ffout)
      
      op_in  => me_in%op
      op_out => me_out%op

      nocc_cls = op_in%n_occ_cls
      njoined  = op_in%njoined

      if (nocc_cls.ne.op_out%n_occ_cls.or.
     &    njoined.ne.op_out%njoined.or.
     &    .not.iocc_equal_n(op_in%ihpvca_occ,op_in%dagger,
     &                      op_out%ihpvca_occ,op_out%dagger,njoined))
     &     call quit(1,'symmetrise',
     &     'Input and output list do not have identical operators')

      ! well, currently for operators with one block only
      if (nocc_cls.gt.1)
     &     call quit(1,'symmetrise',
     &     'currently restricted to one-block operators')

      ! of course, the block must be a "diagonal" block (str(C)==str(A)):
      if (.not.occ_is_diag_blk(op_in%ihpvca_occ,njoined))
     &     call quit(1,'symmetrise',
     &     'the block cannot be symmetrised (not on diagonal)')

      if (me_out%off_op_gmox(1)%maxd.gt.1)
     &     call quit(1,'symmetrise',
     &         'more than 1 distribution - I am lost :-(')


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
        do iocc_cls = 1, nocc_cls 
          if(op_in%formal_blk(iocc_cls))
     &         cycle
          nbuff = nbuff + me_in%len_op_occ(iocc_cls)
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
        call get_vec(ffin,buffer_in,1,nbuff)

      else
        buffer_in => ffin%buffer(1:)
      endif

      if(.not.bufout.and..not.same)then
        nbuff=0
        do iocc_cls = 1, nocc_cls 
          nbuff = nbuff + me_out%len_op_occ(iocc_cls)
        enddo
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
c        buffer_out(1:nbuff) = 0d0
      else if (same) then
        buffer_out => buffer_in
      else
        buffer_out => ffout%buffer(1:)
      endif

      ! factor for off-diagonal and diagonal
      fac_off = 0.5d0*fac
      fac_dia = fac

      ! Loop over occupation classes.
      iocc_loop: do iocc_cls = 1, nocc_cls 

        ! Loop over Ms of annihilator string.
        idxmsa = 0
        msmax = 2
        msa_loop : do msa = msmax, -msmax, -2

          idxmsa = idxmsa+1
      
          ! Loop over Irrep of annihilator string.
          igama_loop: do igama =1, ngam

            ! a dirty quick fix to get the string length:
            len_gam_ms = int(sqrt(dble(me_out%
     &           len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa))))
            if (len_gam_ms.ne.
     &           me_out%ld_op_gmox(iocc_cls)%d_gam_ms(1,igama,idxmsa))
     &           call quit(1,'symmetrise','not true?')

            ioff = me_out%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)

            idx_loop: do idx = 1,len_gam_ms
              jdx_loop: do jdx = 1,idx-1
            
                value = 
     &             fac_off * (buffer_in((idx-1)*len_gam_ms+jdx+ioff) +
     &                        buffer_in((jdx-1)*len_gam_ms+idx+ioff))
                buffer_out((idx-1)*len_gam_ms+jdx+ioff) = value
                buffer_out((jdx-1)*len_gam_ms+idx+ioff) = value

              enddo jdx_loop
              buffer_out(idx*len_gam_ms+ioff) =
     &             fac_dia * buffer_in(idx*len_gam_ms+ioff)
            enddo idx_loop
          
          enddo igama_loop
          
        enddo msa_loop

        ! update norm^2
        xnorm2_blk(iocc_cls) = ddot(me_out%len_op_occ(iocc_cls),
     &       buffer_out(me_out%off_op_occ(iocc_cls)+1),1,
     &       buffer_out(me_out%off_op_occ(iocc_cls)+1),1)
      enddo iocc_loop

      if(.not.bufout)then
        call put_vec(ffout,buffer_out,1,nbuff)
      endif  

      if (open_close_in ) call file_close_keep(ffin)
      if (open_close_out.and..not.same) call file_close_keep(ffout)

      ifree = mem_flushmark('symmetrise')

      return
      end
