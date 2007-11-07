*----------------------------------------------------------------------*
      subroutine set_dens(op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,iformal,orb_info)
*----------------------------------------------------------------------*
*     set up density-like operator (minrank to maxrank)
*     settings analoguous to Hamiltonian
*     hpvx_mnmx,irestr are chosen appropriately
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_operator_array.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 100

      type(operator), intent(inout) ::
     &     op
      character, intent(in) ::
     &     name*(*)
      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     absym, casym, gamma, s2, ms,
     &     min_rank, max_rank, iformal

      type(orbinf), intent(inout) ::
     &     orb_info

      integer ::
     &     occ0(ngastp,2)
      type(operator_array), pointer ::
     &     opscr(:)

      allocate(opscr(2))
      allocate(opscr(1)%op)
      allocate(opscr(2)%op)

      ! density has no external lines, so:
      occ0 = 0
      call set_uop(opscr(1)%op,'scr1',.false.,0,0,1,1,0,
     &     occ0,1,orb_info)

      ! the internal lines are defined by a Hamiltonian-like
      ! operator, so we use ...
      call set_hop(opscr(2)%op,'scr2',dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,iformal,orb_info)

      ! ... and define the density:
      call set_gen_intermediate(op,name,
     &     opscr,2,orb_info)

      call dealloc_operator(opscr(1)%op)
      call dealloc_operator(opscr(2)%op)
      deallocate(opscr(1)%op)
      deallocate(opscr(2)%op)
      deallocate(opscr)
      
      return
      
      end
