*----------------------------------------------------------------------*
      subroutine get_genealogical_coeff(s2,ms2,nc,na,ncoup,
     &               coeff,cperm,aperm)
*----------------------------------------------------------------------*
*     retrieve the genealogical coupling coefficients for an
*     arbitrary operator consisting of a string of creators
*     followed by a string of annihilators for a given alpha/beta perm.
*
*     s2 = 2*S, ms2 = 2*Ms, nc = # of C operators, na = # of A operators
*     ncoup = # of possible coupling pathways = # of components
*     cperm, aperm: alpha/beta permutation for C and A string
*
*     matthias, apr 2013
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     s2, ms2, nc, na, ncoup, cperm(nc), aperm(na)
      real(8), intent(out) ::
     &     coeff(ncoup)

      integer ::
     &     icoup, nup, branchperm(nc+na), iel, is2, ims2, ic, ia
      real(8) ::
     &     sigma
      logical ::
     &     branchnext

      logical, external ::
     &     next_perm
    
      ! now come the loops through the branching diagram
      icoup = 0
      if (mod(nc+na+s2,2).ne.0)
     &   call quit(1,'get_genealogical_coeff','impossible S')
      nup = (nc+na+s2)/2 !total number of up-couplings
      branchperm(1:nup) = 0 ! initial up-couplings
      branchperm(nup+1:nc+na) = 1 !initial down-couplings
      branchnext = .true.
      branch_loop: do while (branchnext)
        ! we may not couple below S=0
        do iel = 1, nc+na
          if (sum(branchperm(1:iel)).gt.iel/2) then
            ! next branch
            branchnext = next_perm(branchperm,nc+na)
            cycle branch_loop
          end if
        end do
        icoup = icoup + 1
c dbg
c        write(lulog,'(a,20i2)') 'current branch:',branchperm
c dbgend

        ! walk along the path and couple everything together
        is2 = 0
        ims2 = 0
        iel = 0
        coeff(icoup) = 1d0
c dbg
c        print *,'icoup: ',icoup
c dbgend
        ! first C string
        do ic = 1, nc
          iel = iel + 1
          is2 = is2-2*branchperm(iel)+1 ! couple total spin up/down
          ims2 = ims2+cperm(ic) ! couple Ms up/down
          ! get coefficient
          sigma = 0.5d0*dble(cperm(ic)) ! alpha or beta
c dbg
c          print *,'get coeff. for is2,ims2,sigma=',is2,ims2,sigma
c dbgend
          if (branchperm(iel).eq.0) then ! up case
            coeff(icoup) = coeff(icoup) *
     &         sqrt((dble(is2)/2d0+sigma*dble(ims2))/dble(is2))
          else ! down case
            coeff(icoup) = -2d0*sigma*coeff(icoup) *
     &         sqrt((dble(is2)/2d0+1d0-sigma*dble(ims2))/
     &              (dble(is2)+2d0))
          end if
        end do
        ! now A string in reverse order
        do ia = na, 1, -1
          iel = iel + 1
          is2 = is2-2*branchperm(iel)+1 ! couple total spin up/down
          ims2 = ims2-aperm(ia) ! couple Ms down/up
          ! get coefficient
          sigma = -0.5d0*dble(aperm(ia)) ! beta or alpha
          if (branchperm(iel).eq.0) then ! up case
            coeff(icoup) = coeff(icoup) *
     &         sqrt((dble(is2)/2d0+sigma*dble(ims2))/dble(is2))
          else ! down case
            coeff(icoup) = -2d0*sigma*coeff(icoup) *
     &         sqrt((dble(is2)/2d0+1d0-sigma*dble(ims2))/
     &              (dble(is2)+2d0))
          end if
          ! beta? -> change sign!
          if (aperm(ia).eq.-1) coeff(icoup)
     &                        = -coeff(icoup)
        end do

        ! next branch
        branchnext = next_perm(branchperm,nc+na)
      end do branch_loop

      return
      end
