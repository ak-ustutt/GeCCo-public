*----------------------------------------------------------------------*
      subroutine check_itf_permute(contr,fl_fact)
*----------------------------------------------------------------------*
*     Check complete contraction for permutation factors
*     This information is used in the ITF translator
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      type(contraction), intent(inout) ::
     &     contr
      type(formula_item), target, intent(inout) ::
     &     fl_fact

      integer ::
     &     isum,           ! Sum of external indices
     &     icount,         ! Number of operators with odd number of external indices
     &     i, j, l         ! Loop index
      logical ::
     &     permute(ngastp,2)      ! True if there is a permutation factor
      type(formula_item), pointer ::
     &     fl_pnt

      permute = .false.

      ! Find the permutation factor
      ! If there is an odd number of external indices over multiple
      ! operators, then there is a factor
      do i = 1, ngastp
        do j = 1, 2
          icount = 0
          do l = 1, contr%nxarc
            isum = contr%xarc(l)%occ_cnt(i,j)
            if (mod(isum,2)>0) then
              icount = icount + 1
            end if
          end do
c     dbg
c          write(lulog,*) 'gas, icount: ',i, icount
c     dbg
          if (icount>1) then
            permute(i,j) = .true.
          end if
        end do
      end do

      ! Loop through all the binary contractions and set the perm flags
      fl_pnt => fl_fact
      do while (associated(fl_pnt%next))
         if (associated(fl_pnt%bcontr)) then
            fl_pnt%bcontr%perm = permute
         end if
         fl_pnt => fl_pnt%next
      end do


      return
      end
