*----------------------------------------------------------------------*
      logical function next_part_connection(occ_conn,
     &     init,occ_dist,nvtx,!occ_ol_vtx,
     &     occ_min,occ_max)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'

      logical, intent(in) ::
     &     init
      integer, intent(out) ::
     &     occ_conn(ngastp,2,nvtx)
      integer, intent(in) ::
     &     nvtx, !occ_ol_vtx(ngastp,2,nvtx),
     &     occ_dist(ngastp,2), occ_min(ngastp,2), occ_max(ngastp,2)

      logical ::
     &     init_flag, ok
      integer ::
     &     ica, igastp
      integer ::
     &     iscr(nvtx)

      logical, external ::
     &     next_part_number

      init_flag = init
      ca_loop: do ica = 1, 2
        do igastp = 1, ngastp
          if (.not.init_flag) iscr(1:nvtx) = occ_conn(igastp,ica,1:nvtx)
          ok = next_part_number(init_flag,.false.,iscr,
     &         occ_dist(igastp,ica),nvtx,
     &         occ_min(igastp,ica),occ_max(igastp,ica))
          occ_conn(igastp,ica,1:nvtx) = iscr(1:nvtx) 
          if (.not.ok.and..not.init_flag) then
            if (ica.eq.2.and.igastp.eq.ngastp) exit ca_loop
            ok = next_part_number(.true.,.false.,iscr,
     &         occ_dist(igastp,ica),nvtx,
     &         occ_min(igastp,ica),occ_max(igastp,ica))
            occ_conn(igastp,ica,1:nvtx) = iscr(1:nvtx) 
          else if (.not.init_flag) then
            exit ca_loop
          end if

        end do
      end do ca_loop

      next_part_connection = ok

      return
      end
