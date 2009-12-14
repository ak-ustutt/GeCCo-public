*----------------------------------------------------------------------*
      subroutine interface_reo_info(reo_info,reord,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     interface routine to convert info from "reord" to "reo_info",
*     including some post-processing of "reo_info"
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_reorder_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      type(reorder_info), intent(inout) ::
     &     reo_info
      type(reorder), intent(in) ::
     &     reord
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     nj_in, nj_out, idum

      nj_out = reord%nj_out
      nj_in  = reord%nj_in

      reo_info%nreo = reord%nreo

      reo_info%n_op_reo = 1
      reo_info%sign_reo = reord%sign

      allocate(reo_info%iocc_reo(ngastp,2,reord%nreo),
     &         reo_info%iocc_opreo0(ngastp,2,nj_in))

      reo_info%iocc_reo = reord%occ_shift
      reo_info%iocc_opreo0 = reord%occ_op0

      if (nj_out.ne.nj_in)
     &     call quit(1,'interface_reo_info',
     &                 'nj_out.ne.nj_in: not yet')

      idum = 0
      call get_reo_info2(2,0,
     &     reord%occ_opout,reord%occ_opin,
     &     reord%rst_opout,reord%rst_opin,
     &     nj_out,idum,idum,
     &     reord%merge_stp1,reord%merge_stp1inv,
     &     reord%merge_stp2,reord%merge_stp2inv,
     &     idum,idum,idum,idum,idum,idum,
     &     reo_info,str_info,orb_info)

      return
      end
