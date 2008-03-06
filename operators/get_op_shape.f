      subroutine get_op_shape(op_int,op_shape,shape_name,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     Routine to examine the outermost vertices of an intermediate-type 
*     operator, op, and produce a single-vertex operator, op_shape, with
*     the same shape (upper and lower arcs).
*     GWR February 2008
*----------------------------------------------------------------------*

      implicit none 

      integer, parameter ::
     &     ntest = 100

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'

      type(operator), intent(in) ::
     &     op_int
      type(operator), intent(inout) ::
     &     op_shape
      character, intent(in) ::
     &     shape_name*(*)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      character(3), parameter ::
     &     opdum_temp  = '_T_'

      integer ::
     &     idx, jdx, iblk, iblk_off, idx_temp
      integer ::
     &     occ_def(ngastp,2)
      type(operator), pointer ::
     &     op_temp

      integer, external ::
     &     idx_oplist2


      if(ntest.ge.100)then
        write(luout,*) '--------------------------'
        write(luout,*) ' Shape of an intermediate '
        write(luout,*) '--------------------------'
      endif

      ! Define an operator based on the first and last vertices of the 
      ! intermediate's blocks, excitation and de-excitation parts 
      ! respectively.
      do idx = 1, op_int%n_occ_cls
        occ_def = 0
        iblk_off = (idx-1)*op_int%njoined
        do jdx = 1, 2
          if(jdx.eq.1)then
            iblk = iblk_off + 1
          else
            iblk = iblk_off + op_int%njoined
          endif

          occ_def = occ_def + iocc_xdn(jdx,
     &         op_int%ihpvca_occ(1:ngastp,1:2,iblk))
        enddo
        if(idx.eq.1)then
          ! Set the first, often only, block.
          call set_uop(op_shape,trim(shape_name),.false.,
     &         occ_def,1,orb_info)
        else
          ! Set the subsequent blocks then join to the first.
          call add_operator(opdum_temp,op_info)
          idx_temp = idx_oplist2(opdum_temp,op_info)
          op_temp => op_info%op_arr(idx_temp)%op
          call set_uop(op_temp,opdum_temp,.false.,
     &         occ_def,1,orb_info)
          call join_operator(op_shape,op_temp,orb_info)

          call del_operator(opdum_temp,op_info)
        endif
      enddo

      if(ntest.ge.100)then
        write(luout,*)'Intermediate shape: '
        call print_op_occ(luout,op_shape)
      endif

      return
      end
