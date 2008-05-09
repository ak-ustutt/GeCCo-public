      subroutine reduce_fact_info(contr_red,contr,idxst,ireo)

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      type(contraction), intent(inout) ::
     &     contr_red
      integer, intent(in) ::
     &     ireo(contr_red%nvtx), idxst

      integer ::
     &     nfac, idx, jdx, ivtx1, ivtx2, iarc, jarc

      nfac = contr%nfac
      contr_red%nfac = nfac

      call resize_contr(contr_red,contr_red%nvtx,contr_red%narc,
     &                            contr_red%nxarc,contr%nfac)
      contr_red%inffac = contr%inffac
      
      do idx = idxst, nfac
        iarc = contr%inffac(4,idx)
        ! old vertices
        ivtx1 = contr%arc(iarc)%link(1)
        ivtx2 = contr%arc(iarc)%link(2)
        ! new vertices
        ivtx1 = ireo(ivtx1)
        ivtx2 = ireo(ivtx2)
        ! look for arc on reduced contraction
        jarc = 0
        do jdx = 1, contr_red%narc
          if (contr_red%arc(jdx)%link(1).eq.ivtx1 .and.
     &        contr_red%arc(jdx)%link(2).eq.ivtx2) then
            jarc = jdx
            exit
          end if
        end do
        if (jarc.eq.0) then
          write(luout,*) 'idx = ',idx
          write(luout,*) 'no appropriate arc on contr_red for ',
     &         ivtx1,ivtx2
          write(luout,*) 'ireo = ',ireo(1:contr%nvtx)
          call prt_contr3(luout,contr,-1)
          call prt_contr3(luout,contr_red,-1)
          call quit(1,'reduce_fact_info','buggy?')
        end if
        contr_red%inffac(4:5,idx) = jarc
      end do

      end
