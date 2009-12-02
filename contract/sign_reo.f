*----------------------------------------------------------------------*
      integer function sign_reo(occ_opori,occ_op0,njoined,
     &     occ_reo,from_to,nreo,
     &     from_to_vtx,nca_vtx,is_op,nvtx)
*----------------------------------------------------------------------*
*     get sign for nreo (simultaneous!) reorderings given by occ_reo
*     occ_opori is the original occupation
*     occ_op0   is the original occupation less the reo occupations
*     extended to considering sign changes if intervening vertices with
*     odd occupation exist
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, intent(in) ::
     &     njoined, nreo, nvtx,
     &     occ_opori(ngastp,2,njoined), occ_op0(ngastp,2,njoined),
     &     occ_reo(ngastp,2,nreo), from_to(2,nreo),
     &     from_to_vtx(2,nreo), is_op(nvtx), nca_vtx(nvtx)

      integer ::
     &     occ_k_from(ngastp,2,njoined), occ_k_to(ngastp,2,njoined),
     &     occ_from_tmp(ngastp,2), occ_to_tmp(ngastp,2)

      logical ::
     &     error
      integer ::
     &     ireo, jreo, ivtx, nencl, ihpv, ica

      integer, external ::
     &     sign_shift, sign_hpvx


      ! savety check:
      error = .false.
      ireo_loop: do ireo = 1, nreo
        do jreo = ireo+1,nreo
          if (from_to(1,ireo).eq.from_to(1,jreo)) then
            do ihpv = 1, ngastp
              do ica = 1, 2
                if (occ_op0(ihpv,ica,from_to(1,ireo)).ne.0.and.
     &              occ_reo(ihpv,ica,ireo).ne.0.and.
     &              occ_reo(ihpv,ica,jreo).ne.0) then
                  error = .true.
                  exit ireo_loop
                end if
              end do
            end do
          end if
          if (from_to(2,ireo).eq.from_to(2,jreo)) then
            do ihpv = 1, ngastp
              do ica = 1, 2
                if (occ_op0(ihpv,ica,from_to(2,ireo)).ne.0.and.
     &              occ_reo(ihpv,ica,ireo).ne.0.and.
     &              occ_reo(ihpv,ica,jreo).ne.0) then
                  error = .true.
                  exit ireo_loop
                end if
              end do
            end do
          end if
        end do
      end do ireo_loop

      if (error) then
        write(luout,*) 'nreo = ',nreo
        call wrt_occ_n(luout,occ_reo,nreo)
        write(luout,*) 'from: ',from_to(1,1:nreo)
        write(luout,*) 'to:   ',from_to(2,1:nreo)
        write(luout,'(x,a,4i4)') 'conflict for ireo, jreo, ihpv, ica:',
     &       ireo, jreo, ihpv, ica
        call quit(1,'sign_reo','reo requires more than pairwise merge!')
      end if

      occ_k_from = 0
      occ_k_to = 0
      do ireo = 1, nreo
        occ_k_from(1:ngastp,1:2,from_to(1,ireo))
     &       = occ_reo(1:ngastp,1:2,ireo)
      end do
c dbg
c      print *,'initial occ_k_from'
c      call wrt_occ_n(luout,occ_k_from,njoined)
c dbg

      sign_reo = 1
      do ireo = 1, nreo

        occ_k_from(1:ngastp,1:2,from_to(1,ireo)) =
     &       occ_k_from(1:ngastp,1:2,from_to(1,ireo)) -
     &       occ_reo(1:ngastp,1:2,ireo)
c dbg
c        print *,'updated occ_k_from'
c        call wrt_occ_n(luout,occ_k_from,njoined)
c dbg

        ! unchanged part of "from" vertex is fix part + moved part:
        occ_from_tmp(1:ngastp,1:2) =
     &         occ_op0(1:ngastp,1:2,from_to(1,ireo))
     &       + occ_k_to(1:ngastp,1:2,from_to(1,ireo))
        ! same for "to" vertex:
        occ_to_tmp(1:ngastp,1:2) =
     &         occ_op0(1:ngastp,1:2,from_to(2,ireo))
     &       + occ_k_to(1:ngastp,1:2,from_to(2,ireo))

        ! count number of CA-op's on passive vertices between reo-vertices:
        nencl = 0
        do ivtx = min(from_to_vtx(1,ireo),from_to_vtx(2,ireo))+1,
     &            max(from_to_vtx(1,ireo),from_to_vtx(2,ireo))-1
          if (is_op(ivtx).gt.0) cycle
          nencl = nencl + nca_vtx(ivtx)
        end do
c dbg
c        if (mod(nencl,2).ne.0) write(luout,*) 'ODD nencl appeared!'
c dbg

        ! correct order before shifting is ...

        ! ... for "to" vertex: unchanged part / to-move part:
        sign_reo = sign_reo*sign_hpvx(2,
     &       occ_to_tmp(1,1),.false.,
     &       occ_k_from(1,1,from_to(2,ireo)),.false.)
c dbg
c        print *,'sign_reo after preparing "to" vertex   = ',sign_reo
c dbg

        ! ... for "from" vertex: unchanged part / reo part / to-move part:
        sign_reo = sign_reo*sign_hpvx(2,
     &       occ_from_tmp(1,1),.false.,
     &       occ_reo(1,1,ireo),.false.)
        ! add reo part to temporary "from" vertex
        occ_from_tmp(1:ngastp,1:2) = occ_from_tmp(1:ngastp,1:2)
     &                             + occ_reo(1:ngastp,1:2,ireo)
        sign_reo = sign_reo*sign_hpvx(2,
     &       occ_from_tmp(1,1),.false.,
     &       occ_k_from(1,1,from_to(1,ireo)),.false.)
c dbg
c        print *,'sign_reo after preparing "from" vertex = ',sign_reo
c dbg

        ! now calculate the sign for shifting
        sign_reo = sign_reo*sign_shift(
     &       occ_reo(1,1,ireo),from_to(1,ireo),from_to(2,ireo),
     &       occ_op0,occ_k_from,occ_k_to,nencl,njoined)
c dbg
c        print *,'sign_reo after shifting                = ',sign_reo
c dbg

        ! correct order after shifting is ...

        ! ... for "from" vertex: unchanged part / to-move part:
        ! first substract reo part again:
        occ_from_tmp(1:ngastp,1:2) = occ_from_tmp(1:ngastp,1:2)
     &                             + occ_reo(1:ngastp,1:2,ireo)
        sign_reo = sign_reo*sign_hpvx(2,
     &       occ_from_tmp(1,1),.false.,
     &       occ_k_from(1,1,from_to(1,ireo)),.false.)
c dbg
c        print *,'sign_reo after remerging "from" vertex = ',sign_reo
c dbg

        ! ... for "to" vertex: unchanged part / reo part / to-move part:
        sign_reo = sign_reo*sign_hpvx(2,
     &       occ_to_tmp(1,1),.false.,
     &       occ_reo(1,1,ireo),.false.)
        ! add reo part to temporary "to" vertex
        occ_to_tmp(1:ngastp,1:2) = occ_to_tmp(1:ngastp,1:2)
     &                           + occ_reo(1:ngastp,1:2,ireo)
        sign_reo = sign_reo*sign_hpvx(2,
     &       occ_to_tmp(1,1),.false.,
     &       occ_k_from(1,1,from_to(2,ireo)),.false.)
c dbg
c        print *,'sign_reo after remerging "to" vertex   = ',sign_reo
c dbg

        occ_k_to(1:ngastp,1:2,from_to(2,ireo)) =
     &       occ_k_to(1:ngastp,1:2,from_to(2,ireo)) +
     &       occ_reo(1:ngastp,1:2,ireo)
c dbg
c        print *,'updated occ_k_to'
c        call wrt_occ_n(luout,occ_k_to,njoined)
c dbg

      end do
c dbg
c      print *,'sign_reo(final) = ',sign_reo
c dbg

      return
      end
