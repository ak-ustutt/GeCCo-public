*----------------------------------------------------------------------*
      logical function next_msgamdist2(first,
     &     msdis_c,msdis_a,gamdis_c,gamdis_a,
     &     nc,na,occ_c,occ_a,ms_c,ms_a,gam_c,gam_a,nsym,ms_fix)
*----------------------------------------------------------------------*
*     increment an entire set of Ms and IRREP distributions
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'multd2h.h'
      
      logical, intent(in) ::
     &     first, ms_fix
      integer, intent(in) ::
     &     nc, na, ms_a, ms_c, gam_a, gam_c, nsym,
     &     occ_c(nc), occ_a(na)
      integer, intent(inout) ::
     &     msdis_c(nc), msdis_a(na),
     &     gamdis_c(nc), gamdis_a(na)

      integer ::
     &     msa_tot, msc_tot, idx
      logical ::
     &     success

      logical, external ::
     &     first_msdistn, first_gamdistn,
     &     next_msdistn, next_gamdistn

      if (first) then
        ! initialize
        success = first_msdistn(msdis_c,ms_c,occ_c,nc)
        success = success.and.
     &       first_msdistn(msdis_a,ms_a,occ_a,na)
        success = success.and.        
     &      first_gamdistn(gamdis_c,gam_c,nsym,nc)
        success = success.and.        
     &      first_gamdistn(gamdis_a,gam_a,nsym,na)
      else          
        ! innermost index: increment IRREP distribution of C
        if (next_gamdistn(gamdis_c,gam_c,nsym,nc)) then
          success = .true.

        ! increment Ms distribution of C
        else if (next_msdistn(msdis_c,ms_c,occ_c,nc)) then
          ! it's the same call to first_gamdist as above, so
          ! actually the result *must* be true (not checked)
          success =
     &         first_gamdistn(gamdis_c,gam_c,nsym,nc)
        ! increment IRREP distribution of A
        else if (next_gamdistn(gamdis_a,gam_a,nsym,na))
     &         then
          success = first_msdistn(msdis_c,ms_c,occ_c,nc)
          success = success.and.
     &         first_gamdistn(gamdis_c,gam_c,nsym,nc)

        ! increment Ms distribution of A
        else if (next_msdistn(msdis_a,ms_a,occ_a,na)) then
          success =
     &         first_gamdistn(gamdis_a,gam_a,nsym,na)
          success = success.and.
     &         first_msdistn(msdis_c,ms_c,occ_c,nc)
          success = success.and.
     &         first_gamdistn(gamdis_c,gam_c,nsym,nc)

        else
          success = .false.
        end if

      end if

      ! Extra condition if we want to fix the c/a terms to have the 
      ! same spin distribution.      
      if(ms_fix)then
c dbg
c        print *,'na,nc',na,nc
c        print *,'msdis_a',msdis_a
c        print *,'msdis_c',msdis_c
c dbg
c        if(na.ne.nc) call quit(1,'next_msgamdist2',
c     &                         'na not equal to nc')
c        msa_tot = 0
c        msc_tot = 0
        do idx = 1, min(na,nc)
          success = success.and.(msdis_a(idx).eq.msdis_c(idx))
c          msa_tot = msa_tot + msdis_a(idx)
c          msc_tot = msc_tot + msdis_c(idx)
        enddo
c        success = success.and.(msa_tot.eq.msc_tot)
      endif          
      
      next_msgamdist2 = success

      return
      end
*----------------------------------------------------------------------*
*     a few aux-routines follow:
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      logical function next_msdistn(msdis,mst,msmxdis,n)
*----------------------------------------------------------------------*
*     function to generate MS distributions 
*     first index runs fastest, last index is fixed by mst
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     n, msmxdis(n), mst
      integer, intent(inout) ::
     &     msdis(n)

      integer ::
     &     i

      if (n.gt.1) then
        do
          next_msdistn = .false.
          do i = 1, n-1
            if (msdis(i).gt.-msmxdis(i)) then
              ! decrement if possible
              msdis(i) = msdis(i) - 2
              next_msdistn = .true.
              exit
            else
              ! else, set to max. value and try next index
              msdis(i) = msmxdis(i)
            end if
          end do
          ! if false: even (n-1) block was at min. Ms -> end
          if (.not.next_msdistn) exit
          ! last Ms is fixed by Ms(total)
          msdis(n) = mst - sum(msdis(1:n-1))
          ! is this a possible value? exit, else: go on ...
          if (abs(msdis(n)).le.msmxdis(n)) exit
        end do
      else
        next_msdistn = .false.
      end if

      return
      end
*----------------------------------------------------------------------*
      logical function first_msdistn(msdis,mst,msmxdis,n)
*----------------------------------------------------------------------*
*     function to generate first MS distribution 
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     n, msmxdis(n), mst
      integer, intent(inout) ::
     &     msdis(n)

      logical, external ::
     &     next_msdistn

      if (n.eq.0) then
        ! empty string can be Ms=0 only
        first_msdistn = mst.eq.0
        return
      end if

      first_msdistn = .true.
      ! start with max value for all entries
      if (n.gt.1) then
        msdis(1:n-1) = msmxdis(1:n-1)
        msdis(n) = mst-sum(msdis(1:n-1))
      else
        msdis(1) = mst
      end if
      ! if this is an acceptable value for msdis --> exit
      if (abs(msdis(n)).le.msmxdis(n)) return
      
      ! else scan for a better choice
      first_msdistn = next_msdistn(msdis,mst,msmxdis,n)

      return
      end

*----------------------------------------------------------------------*
      logical function next_gamdistn(gamdis,gamt,nsym,n)
*----------------------------------------------------------------------*
*     function to generate IRREP distributions
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     gamt, n, nsym
      integer, intent(inout) ::
     &     gamdis(n)

      integer ::
     &     gam_acc, i, i_last

      if (n.gt.1) then
        next_gamdistn = .false.
        do i = 1, n-1
          if (gamdis(i).lt.nsym) then
            ! increment if possible
            gamdis(i) = gamdis(i) + 1
            i_last = i
            next_gamdistn = .true.
            exit
          else
            ! reset to min. value and try next index
            gamdis(i) = 1
          end if
        end do
        ! if false: even (n-1) block was at max. IRREP -> end
        if (next_gamdistn) then
          ! fix last IRREP 
          ! (1 to i_last-1 were set to 1, so skip them)
          gam_acc = gamt
          do i = i_last, n-1
            gam_acc = multd2h(gam_acc,gamdis(i))
          end do
          gamdis(n) = gam_acc
        end if
      else
        next_gamdistn = .false.
      end if

      return
      end

*----------------------------------------------------------------------*
      logical function first_gamdistn(gamdis,gamt,nsym,n)
*----------------------------------------------------------------------*
*     function to generate first IRREP distribution
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     gamt, nsym, n
      integer, intent(inout) ::
     &     gamdis(n)

      if (n.eq.0) then
        ! empty string can be totally symmetric only
        first_gamdistn = gamt.eq.1
        return
      end if

      first_gamdistn = .true.
      ! easy: IRREP 1 for 1..n-1
      ! so IRREP gamt for n
      if (n.gt.1) gamdis(1:n-1) = 1 
      gamdis(n) = gamt

      return
      end
