*----------------------------------------------------------------------*
      subroutine set_unit_operator(op_info,orb_info)
*----------------------------------------------------------------------*
*     Set-up unit operator (an entirely formal object).
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'par_opnames_gen.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(operator), pointer ::
     &     op_pnt
      logical ::
     &     dagger
      character ::
     &     name*(len_opname)
      integer ::
     &     absym, casym, s2, ms, min_rank, max_rank, ncadiff,
     &     gamma, iarr(1), iformal, isim, idx
      integer ::
     &     ihpv_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2)

      integer, external ::
     &     idx_oplist2

      ! new entry: the unit operator
      call add_operator(op_unity,op_info)
      idx = idx_oplist2(op_unity,op_info)
      op_pnt => op_info%op_arr(idx)%op

      ! new entry
      name = op_unity
      dagger = .false.
      absym = 0
      casym = 0
      gamma = 1
      s2 = 0
      ms = 0
      min_rank = 0
      max_rank = 0
      ncadiff = 0
      iformal=0

      ihpv_mnmx = 0
      irestr    = 0

      call set_genop(op_pnt,name,optyp_operator,
     $     dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info)
 
      return

      end
