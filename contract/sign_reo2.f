*----------------------------------------------------------------------*
      integer function sign_reo2(occ_opori,occ_op0,njoined,
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
     &     occ_k_from(ngastp,2,njoined), occ_k_to(ngastp,2,njoined)

      logical ::
     &     error
      integer ::
     &     ireo, jreo, ivtx, nencl

      integer, external ::
     &     sign_shift, sign_hpvx
      call quit(1,'sign_reo2','call to obsolete routine')


      ! savety check:
      error = .false.
      do ireo = 1, nreo
        do jreo = ireo+1,nreo
          if (from_to(1,ireo).eq.from_to(1,jreo) .or.
     &        from_to(2,ireo).eq.from_to(2,jreo)) then
            error = .true.
          end if
        end do
      end do

      if (error) then
        write(luout,*) 'nreo = ',nreo
        call wrt_occ_n(luout,occ_reo,nreo)
        write(luout,*) 'from: ',from_to(1,1:nreo)
        write(luout,*) 'to:   ',from_to(2,1:nreo)
c        call quit(1,'sign_reo','reo from or to same vertex occurred!')
        sign_reo2 = 0
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

      sign_reo2 = 1
      do ireo = 1, nreo

        occ_k_from(1:ngastp,1:2,from_to(1,ireo)) =
     &       occ_k_from(1:ngastp,1:2,from_to(1,ireo)) -
     &       occ_reo(1:ngastp,1:2,ireo)
c dbg
c        print *,'updated occ_k_from'
c        call wrt_occ_n(luout,occ_k_from,njoined)
c dbg
        ! count number of CA-op's on passive vertices between reo-vertices:
        nencl = 0
        do ivtx = min(from_to_vtx(1,ireo),from_to_vtx(2,ireo))+1,
     &            max(from_to_vtx(1,ireo),from_to_vtx(2,ireo))-1
          if (is_op(ivtx).gt.0) cycle
          nencl = nencl + nca_vtx(ivtx)
        end do
c dbg
        if (mod(nencl,2).ne.0) write(luout,*) 'ODD nencl appeared!'
c dbg

        sign_reo2 = sign_reo2*sign_shift(
     &       occ_reo(1,1,ireo),from_to(1,ireo),from_to(2,ireo),
     &       occ_op0,occ_k_from,occ_k_to,nencl,njoined)
c dbg
c        print *,'sign_reo2(1) = ',sign_reo2
c dbg
        sign_reo2 = sign_reo2*sign_hpvx(2,
     &       occ_op0(1,1,from_to(1,ireo)),.false.,
     &       occ_reo(1,1,ireo),.false.)
c dbg
c        print *,'sign_reo2(2a) = ',sign_reo2
c dbg
        sign_reo2 = sign_reo2*sign_hpvx(2,
     &       occ_op0(1,1,from_to(2,ireo)),.false.,
     &       occ_reo(1,1,ireo),.false.)
c dbg
c        print *,'sign_reo2(2b) = ',sign_reo2
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
c      print *,'sign_reo2(final) = ',sign_reo2
c dbg

      return
      end
