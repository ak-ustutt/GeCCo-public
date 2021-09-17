*-----------------------------------------------------------------------*
      subroutine contr_get_eqv_line_factor(contr,op_info)
*-----------------------------------------------------------------------*
*
*     Count equivalent lines in the contraction and determine associated
*     factor. Not required for GeCCo itself as it contracts only over
*     non-redundant index tuples, but for ITF code it is needed
*
*     input/output: contr (set contr%eqvl_fact)
*     auxiliary input: op_info
*
*     andreas, march 2021
*
*-----------------------------------------------------------------------*
      Implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     ntest = 00
      character(len=25), parameter ::
     &     i_am = 'contr_get_eqv_line_factor'
      
      type(contraction), intent(inout) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     nvtx, narc, nxarc,
     &     iarc, ii, ica, igastp
      real(8) ::
     &     fact
      
      integer, pointer ::
     &     occ_vtx(:,:,:), occ_cnt(:,:), vtxseq(:,:)
      type(cntr_arc), pointer ::
     &     arc
      
      if (ntest.gt.0) call write_title(lulog,wst_dbg_subr,i_am)

      if (ntest.ge.100) then
        write(lulog,*) 'contraction: '
        call prt_contr2(lulog,contr,op_info)
      end if

      nvtx = contr%nvtx
      narc = contr%narc
      nxarc = contr%nxarc

      fact = 1.0
!     loop over contraction arcs
      do iarc = 1, narc
        arc => contr%arc(iarc)
        occ_cnt => arc%occ_cnt
        do igastp = 1, ngastp
          do ica = 1, 2
            if (occ_cnt(igastp,ica).gt.1) then
              ! update fact with 1/(n!)
              do ii = 2, occ_cnt(igastp,ica)
                fact = fact/dble(ii)
              end do
            end if
          end do
        end do
      end do

      contr%eqvl_fact = fact

      if (ntest.ge.100) then
        write(lulog,'(1x,"determined equiv. line factor = ",f20.12)')
     &       fact
      end if
        
      return
      
      end subroutine

      
      
