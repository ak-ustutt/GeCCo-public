      subroutine diag_fill(ffdia,xdia,op,str_info,orb_info)
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
     &     ntest = 1000
      
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'
      include 'ifc_baserout.h'
      include 'hpvxseq.h'
      include 'multd2h.h'

      real(8), intent(in) ::
     &     xdia(*)
      type(operator), intent(in) ::
     &     op
      type(strinf), intent(in), target ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      type(filinf), intent(inout) ::
     &     ffdia

      integer ::
     &     ifree, ngam, nocc_cls, njoined, iocc_cls, ms, idxms, igam,
     &     nbuff, len_gam_ms, ioff_dia, ioff_out, idx, jdx
      logical ::
     &     buffered
      logical ::
     &     loop(op%n_occ_cls)
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall
      real(8), pointer ::
     &     buffer(:)

      call atim_csw(cpu0,sys0,wall0)

      if(ntest.ge.100)then
        write(luout,*)'======================================='
        write(luout,*)' Formation of diagonal of 2-e operator '
        write(luout,*)'======================================='
      endif

      ifree = mem_setmark('diag_fill')

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
            cycle
          else
            nbuff = nbuff + op%len_op_occ(iocc_cls)
          endif
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
        if(loop(iocc_cls)) cycle occ_loop

        ! Loop over Ms values.
        idxms = 0
        do ms = 2, -2, -2
          idxms = idxms+1

          ! Loop over symmetry blocks.
          do igam = 1, ngam

            len_gam_ms = int(sqrt(dble(op%
     &           len_op_gmo(iocc_cls)%gam_ms(igam,idxms))))
            
            ioff_out = op%off_op_gmo(iocc_cls)%gam_ms(igam,idxms)

            do idx = 1, len_gam_ms

              do jdx = 1, len_gam_ms
                ! Copy diagonal elements to all elements of the row.
c                buffer(ioff_out+(jdx-1)*len_gam_ms+idx) =
c     &               xdia(ioff_dia+idx)
c dbg
c                ! Column?
c                buffer(ioff_out+(idx-1)*len_gam_ms+jdx) =
c     &               xdia(ioff_dia+idx)
                buffer(ioff_out+(idx-1)*len_gam_ms+jdx) = 1d0
c dbg
              enddo
c dbg
c              ! Unit operator?
c              buffer(ioff_out+(idx-1)*len_gam_ms+idx) = 1d0
c dbg
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
        call wrt_op_file(luout,5,ffdia,op,1,
     &       nocc_cls,str_info,orb_info)
      endif

      call atim_csw(cpu,sys,wall)
      
      ifree = mem_flushmark()

      if (iprlvl.ge.5)
     &     call prtim(luout,'time in diag_fill ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
