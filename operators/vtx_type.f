*----------------------------------------------------------------------*
      integer function vtx_type(op)
*----------------------------------------------------------------------*
*     analyze blocks of operator and return type as defined in 
*     def_operator.h
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_operator.h'
      include 'ifc_operators.h'

      type(operator), intent(in) ::
     &     op

      logical ::
     &     ex, dx, vl
      integer ::
     &     iblk, nblk
      integer, pointer ::
     &     hpvx_occ(:,:,:)

      hpvx_occ => op%ihpvca_occ
      nblk = op%n_occ_cls
      
      ex = .false.
      dx = .false.
      vl = .false.
      do iblk = 1, nblk
        ex = ex.or.iocc_nonzero(iocc_xdn(1,hpvx_occ(1:ngastp,1:2,iblk)))
        dx = dx.or.iocc_nonzero(iocc_xdn(2,hpvx_occ(1:ngastp,1:2,iblk)))
        vl = vl.or.iocc_nonzero(iocc_xdn(3,hpvx_occ(1:ngastp,1:2,iblk)))
      end do
      
      if (.not.ex.and..not.dx.and..not.vl) vtx_type = vtxtyp_scalar
      if (ex.and..not.dx.and..not.vl) vtx_type = vtxtyp_ph_ex
      if (dx.and..not.ex.and..not.vl) vtx_type = vtxtyp_ph_dx
      if (dx.and.ex.and..not.vl) vtx_type = vtxtyp_ph
      if (vl) vtx_type = vtxtyp_val

      return
      end
