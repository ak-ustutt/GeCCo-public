*----------------------------------------------------------------------*
      subroutine set_len_str(len_str,nc,na,
     &                       graphs,
     &                       graph_c,idxms_c,gam_c,hpvx_c,
     &                       graph_a,idxms_a,gam_a,hpvx_a,
     &                       hpvxseq,resort)
*----------------------------------------------------------------------*
*     set string lengths
*     if (resort): resort to HPVX major sequence (default: C A major)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_graph.h'

      type(graph), intent(in), target ::
     &     graphs(*)
      logical, intent(in) ::
     &     resort
      integer, intent(in) ::
     &     nc, na,
     &     graph_c(nc),idxms_c(nc),gam_c(nc),hpvx_c(nc),
     &     graph_a(na),idxms_a(na),gam_a(na),hpvx_a(na),
     &     hpvxseq(ngastp)
      integer, intent(out) ::
     &     len_str(nc+na)

      integer ::
     &     idx_ca, hpvx, idx_hpvx, idx_c, idx_a, idx_high
      
      if (na+nc.eq.0) return
c dbg
c      print *,'gam_a: ',gam_a
c dbg

      idx_ca = 0
      idx_high = 1
      if (resort) idx_high = ngastp
      do idx_hpvx = 1, idx_high
        hpvx = hpvxseq(idx_hpvx)
        do idx_c = 1, nc
          if (hpvx_c(idx_c).eq.hpvx.or..not.resort) then
            idx_ca = idx_ca + 1
            len_str(idx_ca) =
     &           graphs(graph_c(idx_c))%
     &           lenstr_gm(gam_c(idx_c),idxms_c(idx_c))
          end if
        end do
        do idx_a = 1, na
          if (hpvx_a(idx_a).eq.hpvx.or..not.resort) then
            idx_ca = idx_ca + 1
            len_str(idx_ca) =
     &           graphs(graph_a(idx_a))%
     &           lenstr_gm(gam_a(idx_a),idxms_a(idx_a))
          end if
        end do
      end do

      end
