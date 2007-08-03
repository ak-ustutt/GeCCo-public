*----------------------------------------------------------------------*
      subroutine set_h_operators(op_info,orb_info,explicit)
*----------------------------------------------------------------------*
*     Hard-wired set-up of the Hamiltonian operator.
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'cc_routes.h'
      include 'ifc_input.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'par_opnames_gen.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      logical, intent(in) ::
     &     explicit
      type(orbinf) ::
     &     orb_info
      type(operator), pointer ::
     &     op_pnt, ham_pnt
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

      ! new entry: the Hamiltonian
      call add_operator(op_ham,op_info)
      idx = idx_oplist2(op_ham,op_info)
      op_pnt => op_info%op_arr(idx)%op

      ! new entry: the Hamiltonian
      name = op_ham
      dagger = .false.
      absym = 0
      casym = 0
      gamma = 1
      s2 = 0
      ms = 0
      min_rank = 0
      max_rank = 2
      ncadiff = 0
      iformal=2

      call set_hpvx_and_restr_for_h()

      call set_genop(op_pnt,name,optyp_operator,
     $     dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info)

      call get_argument_value('calculate.routes','simtraf',ival=isim)
      if (isim.gt.0) then
        call add_operator(op_hhat,op_info)
        idx = idx_oplist2(op_hhat,op_info)
        op_pnt => op_info%op_arr(idx)%op
        idx = idx_oplist2(op_ham,op_info)
        ham_pnt => op_info%op_arr(idx)%op
        call clone_operator(op_pnt,ham_pnt,orb_info)

      end if
 
      return
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      contains
*----------------------------------------------------------------------*
*     some loops extracted for better overview:
*----------------------------------------------------------------------*
      subroutine set_hpvx_and_restr_for_h()
*----------------------------------------------------------------------*

      implicit none

      integer ::
     &     ica, igastp, igas

      do ica = 1, 2
        do igastp = 1, ngastp-1
          if (orb_info%nactt_hpv(igastp).gt.0) then
            ihpv_mnmx(1,igastp,ica) = min_rank
            ihpv_mnmx(2,igastp,ica) = max_rank
          else
            ihpv_mnmx(1,igastp,ica) = 0
            ihpv_mnmx(2,igastp,ica) = 0
          end if
        end do
      end do
      if(explicit)then
        ihpv_mnmx(1,ngastp,1:2)=min_rank
        ihpv_mnmx(2,ngastp,1:2)=max_rank
      else
        ihpv_mnmx(1:2,ngastp,1:2)=0
      endif  

      irestr(1:2,1:orb_info%ngas,1:2,1:2) = 0
      do ica = 1, 2
        do igas = 1, orb_info%ngas
          irestr(1,igas,ica,1) = min_rank
          irestr(2,igas,ica,1) = max_rank
        end do
      end do

      return
      end subroutine set_hpvx_and_restr_for_h

      end
