*----------------------------------------------------------------------*
      logical function next_msgamdist_diag(first,
     &     msdis,gamdis,
     &     nn,occ,ms,gam,nsym)
*----------------------------------------------------------------------*
*     increment an entire set of Ms and IRREP distributions
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'multd2h.h'
      
      logical, intent(in) ::
     &     first
      integer, intent(in) ::
     &     nn, ms, gam, nsym,
     &     occ(nn)
      integer, intent(inout) ::
     &     msdis(nn),
     &     gamdis(nn)

      logical ::
     &     success

      logical, external ::
     &     first_msdistn, first_gamdistn,
     &     next_msdistn, next_gamdistn

      if (first) then
        ! initialize
        success = first_msdistn(msdis,ms,occ,nn)
        success = success.and.        
     &      first_gamdistn(gamdis,gam,nsym,nn)
      else          
        ! innermost index: increment IRREP distribution
        if (next_gamdistn(gamdis,gam,nsym,nn)) then
          success = .true.

        ! increment Ms distribution
        else if (next_msdistn(msdis,ms,occ,nn)) then
          ! it's the same call to first_gamdist as above, so
          ! actually the result *must* be true (not checked)
          success =
     &         first_gamdistn(gamdis,gam,nsym,nn)
        else
          success = .false.
        end if

      end if

      next_msgamdist_diag = success

      return
      end
