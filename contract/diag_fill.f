      subroutine diag_fill(me_dia,xdia,str_info,orb_info)
*-----------------------------------------------------------------------
*     Routine to compile the diagonal elements of a 2-electron operator.
*     This is done by simply copying the (kl,kl) element to the
*     (kl,mn) elements for all mn in the block. 
*
*     This is valid for approx. A, but will need to be more involved
*     later on. Currently works well for operators with only one block
*     but may need to check the routine if operators with more than one
*     are needed.

*     GWR November 2007
*-----------------------------------------------------------------------

      implicit none

      integer, parameter ::
     &     ntest = 00
      
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'
      include 'ifc_baserout.h'
      include 'hpvxseq.h'
      include 'multd2h.h'

      real(8), intent(in) ::
     &     xdia(*)
      type(me_list), intent(in) ::
     &     me_dia
      type(strinf), intent(in), target ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     ifree, ngam, nocc_cls, njoined, iocc_cls, ms, idxms, igam,
     &     nbuff, len_gam_ms, ioff_dia, ioff_out, idx, jdx
      logical ::
     &     buffered
      logical ::
     &     loop(me_dia%op%n_occ_cls)
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall
      real(8), pointer ::
     &     buffer(:)

      type(filinf), pointer ::
     &     ffdia
      type(operator), pointer ::
     &     op

      call atim_csw(cpu0,sys0,wall0)

      if(ntest.ge.100)then
        write(luout,*)'======================================='
        write(luout,*)' Formation of diagonal of 2-e operator '
        write(luout,*)'======================================='
      endif

      ifree = mem_setmark('diag_fill')

      op => me_dia%op
      ffdia => me_dia%fhand

      ! Define some parameters.
      ngam = orb_info%nsym
      nocc_cls = op%n_occ_cls

      ! Initialisations.
      loop(1:nocc_cls) = .false.
      njoined = op%njoined

      ! Open the output file
      call file_open(ffdia)
      buffered = ffdia%buffered
      
      ! Allocate memory for the output.
      nbuff = 0
      if(.not.buffered) then
        do iocc_cls = 1, nocc_cls
          if(op%formal_blk(iocc_cls))then
            loop(iocc_cls)=.true.
          endif
          nbuff = nbuff + me_dia%len_op_occ(iocc_cls)
        enddo
        if(nbuff.gt.ifree)
     &       call quit(1,'diag_fill','Inadequate memory')
        ifree = mem_alloc_real(buffer,nbuff,'buffer')
        buffer(1:nbuff) = 0d0
      else
        buffer => ffdia%buffer(1:)
      endif
      
      ioff_dia = 0
      ! Loop over the occupation classes of the operator.
      occ_loop: do iocc_cls = 1, nocc_cls
c        if(loop(iocc_cls)) cycle occ_loop

        ! Loop over Ms values.
        idxms = 0
        do ms = 2, -2, -2
          idxms = idxms+1

          ! Loop over symmetry blocks.
          do igam = 1, ngam

            len_gam_ms = int(sqrt(dble(me_dia%
     &           len_op_gmo(iocc_cls)%gam_ms(igam,idxms))))
            
            ioff_out = me_dia%off_op_gmo(iocc_cls)%gam_ms(igam,idxms)

            do idx = 1, len_gam_ms

              do jdx = 1, len_gam_ms
                if(.not.loop(iocc_cls))then

c                  if(.not.mp2)then
                    ! Copy diagonal elements to all elements of the row.
                    buffer(ioff_out+(jdx-1)*len_gam_ms+idx) =
     &                   xdia(ioff_dia+idx)
c                  else
c dbg
c                    ! Setting preconditioner to 1 for MP2-R12.
c                    buffer(ioff_out+(idx-1)*len_gam_ms+jdx) = 1d0
c                  endif
c dbg
                else
                  ! Set formal blocks to 0. 
                  ! ????? formal blocks should never appear in the
                  ! numerical part .... ????????
                  buffer(ioff_out+(idx-1)*len_gam_ms+jdx) = 0d0
                endif
              enddo
            enddo

            ioff_dia = ioff_dia + len_gam_ms

          enddo
        enddo
      enddo occ_loop

      if(.not.buffered)then
        call put_vec(ffdia,buffer,1,nbuff)
      endif

      call file_close_keep(ffdia)

      if(ntest.ge.1000)then
        call wrt_mel_file(luout,5,me_dia,1,
     &       nocc_cls,str_info,orb_info)
      endif

      call atim_csw(cpu,sys,wall)
      
      ifree = mem_flushmark()

      if (iprlvl.ge.5)
     &     call prtim(luout,'time in diag_fill ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
