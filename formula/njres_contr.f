      integer function njres_contr(contr)

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr

      type(cntr_arc), pointer ::
     &     xarc(:)
      integer ::
     &     iarc, nxarc

      xarc => contr%xarc
      nxarc = contr%nxarc
      
      njres_contr = 1
      do iarc = 1, nxarc
        njres_contr = max(njres_contr,xarc(iarc)%link(2))
      end do

      return
      end
