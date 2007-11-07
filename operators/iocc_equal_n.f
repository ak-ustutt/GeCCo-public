
*----------------------------------------------------------------------*
      logical function iocc_equal_n(iocc,dagi,jocc,dagj,njoined)
*----------------------------------------------------------------------*
*     compare two occupations and return .true. if equal
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
*----------------------------------------------------------------------*

      logical, intent(in) ::
     &     dagi, dagj
      integer, intent(in) ::
     &     njoined,
     &     iocc(ngastp,2,njoined), jocc(ngastp,2,njoined)

      integer ::
     &     ica, ica_i, ica_j, ihpv,
     &     ijoin, idx, jdx, inc, jnc

      iocc_equal_n = .true.

      if (.not.dagi) then
        idx = 1
        inc = 1
      else
        idx = njoined
        inc = -1
      end if

      if (.not.dagj) then
        jdx = 1
        jnc = 1
      else
        jdx = njoined
        jnc = -1
      end if

      outer_loop: do ijoin = 1, njoined

        do ica = 1,2
          ica_i = ica
          ica_j = ica
          if (dagi) ica_i = 3-ica
          if (dagj) ica_j = 3-ica
        
          do ihpv = 1, ngastp
            iocc_equal_n = iocc_equal_n.and.
     &           iocc(ihpv,ica_i,idx).eq.jocc(ihpv,ica_j,jdx)
          end do
          if (iocc_equal_n.eq..false.) exit outer_loop
        end do
        
        idx = idx+inc
        jdx = jdx+jnc
      end do outer_loop

      return
      end
