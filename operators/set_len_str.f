*----------------------------------------------------------------------*
      subroutine set_len_str(len_str,nc,na,
     &                       graphs,
     &                       graph_c,idxms_c,gam_c,hpvx_c,
     &                       graph_a,idxms_a,gam_a,hpvx_a,
     &                       hpvxseq)
*----------------------------------------------------------------------*
*     set string lengths
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_graph.h'

      type(graph), intent(in), target ::
     &     graphs(*)
      integer, intent(in) ::
     &     nc, na,
     &     graph_c(nc),idxms_c(nc),gam_c(nc),hpvx_c(nc),
     &     graph_a(na),idxms_a(na),gam_a(na),hpvx_a(na),
     &     hpvxseq(ngastp)
      integer, intent(out) ::
     &     len_str(nc+na)

      integer ::
     &     idx_ca, hpvx, idx_hpvx, idx_c, idx_a
      
      if (na+nc.eq.0) return

      idx_ca = 0
      do idx_hpvx = 1, ngastp
        hpvx = hpvxseq(idx_hpvx)
        do idx_c = 1, nc
          if (hpvx_c(idx_c).eq.hpvx) then
            idx_ca = idx_ca + 1
            len_str(idx_ca) =
     &           graphs(graph_c(idx_c))%
     &           lenstr_gm(gam_c(idx_c),idxms_c(idx_c))
          end if
        end do
        do idx_a = 1, na
          if (hpvx_a(idx_a).eq.hpvx) then
            idx_ca = idx_ca + 1
            len_str(idx_ca) =
     &           graphs(graph_a(idx_a))%
     &           lenstr_gm(gam_a(idx_a),idxms_a(idx_a))
          end if
        end do
      end do

      end
