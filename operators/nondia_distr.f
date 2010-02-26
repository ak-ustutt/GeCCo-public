*----------------------------------------------------------------------*
      logical function nondia_distr(mapca,diag_idx,diag_ca,
     &                              ms_c,ms_a,gam_c,gam_a,
     &                              nc,msdiag,gamdiag)
*----------------------------------------------------------------------*
*     determines whether a ms/gamma distribution is diagonal and
*     whether the index tuple defining the diagonal has the given
*     symmetry and ms value.
*
*     matthias, feb 2010
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'multd2h.h'

      integer, intent(in) ::
     &     mapca(nc), diag_idx(nc), diag_ca(nc), ms_c(nc), ms_a(nc),
     &     gam_c(nc), gam_a(nc), nc, msdiag, gamdiag

      integer ::
     &     ms_diag, gam_diag, idx

      nondia_distr = .false.
      ms_diag = 0
      gam_diag = 1
      do idx = 1, nc
        nondia_distr = nondia_distr.or.
     &                 ms_a(mapca(idx)).ne.ms_c(idx).or.
     &                 gam_a(mapca(idx)).ne.gam_c(idx)
        if (diag_ca(idx).eq.1) then
          ms_diag = ms_diag + ms_c(diag_idx(idx))
          gam_diag = multd2h(gam_diag,gam_c(diag_idx(idx)))
        else
          ms_diag = ms_diag - ms_a(diag_idx(idx))
          gam_diag = multd2h(gam_diag,gam_a(diag_idx(idx)))
        end if
      end do
      nondia_distr = nondia_distr.or.
     &               (msdiag.ne.999.and.ms_diag.ne.msdiag).or.
     &               (gamdiag.ne.0.and.gam_diag.ne.gamdiag)

      return
      end
*----------------------------------------------------------------------*
