*----------------------------------------------------------------------*
      subroutine set_r12intm(op,name,dagger,
     &     min_rank,max_rank,ncadiff,iformal,orb_info)
*----------------------------------------------------------------------*
*     set up intermediate for R12 with open lines within
*     settings analoguous to R12coeff
*     hpvx_mnmx,irestr are chosen appropriately
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_operator_array.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(operator), intent(inout) ::
     &     op
      character, intent(in) ::
     &     name*(*)
      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     min_rank, max_rank, ncadiff, iformal

      type(orbinf), intent(inout) ::
     &     orb_info

      integer ::
     &     idx
      integer ::
     &     occ0(ngastp,2)
      type(operator_array), pointer ::
     &     opscr(:)

      allocate(opscr(2))
      allocate(opscr(1)%op)
      allocate(opscr(2)%op)

      ! no external lines, so:
      occ0 = 0
      call set_uop(opscr(1)%op,'scr1',.false.,
     &     occ0,1,orb_info)

      ! the internal lines are defined by a coefficient-like 
      ! operator, so we use ...
      call set_r12c(opscr(2)%op,'scr2',dagger,
     &     min_rank,max_rank,ncadiff,iformal,orb_info)

      ! ... and define the intermediate:
      call set_gen_intermediate(op,name,
     &     opscr,2,orb_info)

      op%formal = .true.
      do idx = 1, op%n_occ_cls
        op%formal = op%formal.and.op%formal_blk(idx)
      end do

      call dealloc_operator(opscr(1)%op)
      call dealloc_operator(opscr(2)%op)
      deallocate(opscr(1)%op)
      deallocate(opscr(2)%op)
      deallocate(opscr)
      
      return
      
      end
