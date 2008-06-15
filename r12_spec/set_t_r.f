      subroutine set_t_r(flist_t_r,bar,idxsop,idxtop,
     &     idxr12,idxc12,idxcpp12,
     &     r12op,r12fix,op_info)

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      type(formula_item), intent(inout), target ::
     &     flist_t_r

      logical, intent(in) ::
     &     bar
      integer, intent(in) ::
     &     idxsop, idxtop, idxr12, idxc12, idxcpp12,
     &     r12op,r12fix

      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     fl_t_r_pnt

      integer ::
     &     iblk_pxhp, iblk_xxhp, iblk_pxpp, iblk_xxpp,
     &     occ(ngastp,2)

      integer, external ::
     &     iblk_occ
      
      fl_t_r_pnt => flist_t_r
      call new_formula_item(fl_t_r_pnt,command_set_target_init,idxsop)
      fl_t_r_pnt => fl_t_r_pnt%next
      call expand_op_product(fl_t_r_pnt,idxsop,
     &       1d0,1,idxtop,-1,-1,
     &       0,0,.false.,op_info)
      do while(associated(fl_t_r_pnt%next))
        fl_t_r_pnt => fl_t_r_pnt%next
      end do
      
      if (.not.bar) then
        if(r12op.eq.0.and..not.r12fix)then
          call expand_op_product(fl_t_r_pnt,idxsop,
     &         1d0,2,(/idxc12,idxr12/),-1,-1,
     &         (/1,2/),1,.false.,op_info)
        else
          call expand_op_product(fl_t_r_pnt,idxsop,
     &         1d0,1,idxr12,-1,-1,
     &         0,0,.false.,op_info)
        endif
      else
        if(r12op.eq.0.and..not.r12fix)then
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,-idxr12,idxc12,idxsop/),
     &         (/1,2,3,1/),
     &         -1,-1,
     &         (/2,3/),1,
     &         0,0,
     &         0,0,op_info)
        else
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,3,2,
     &         (/idxsop,-idxr12,idxsop/),
     &         (/1,2,1/),
     &         -1,-1,
     &         0,0,
     &         0,0,
     &         0,0,op_info)
        endif
      end if

      if(r12op.eq.1.or.r12op.ge.3)then
        do while(associated(fl_t_r_pnt%next))
          fl_t_r_pnt => fl_t_r_pnt%next
        enddo
        
        ! find xp|hp and xx|hp blocks of R12
        occ = 0
        occ(IPART,1) = 1
        occ(IEXTR,1) = 1
        occ(IHOLE,2) = 1
        occ(IPART,2) = 1
        iblk_pxhp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op)
        occ = 0
        occ(IEXTR,1) = 2
        occ(IHOLE,2) = 1
        occ(IPART,2) = 1
        iblk_xxhp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op)
      

        if (.not.bar) then
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxr12,idxc12,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
     &         (/1,iblk_pxhp,1,1/),(/0,iblk_pxhp,0,0/),
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        else
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxc12,-idxr12,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
     &         (/1,1,iblk_pxhp,1/),(/0,0,iblk_pxhp,0/),
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        end if
        do while(associated(fl_t_r_pnt%next))
          fl_t_r_pnt => fl_t_r_pnt%next
        enddo
        if (.not.bar) then
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxr12,idxc12,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
     &         (/1,iblk_xxhp,1,1/),(/0,iblk_xxhp,0,0/),
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        else
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxc12,-idxr12,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
     &         (/1,1,iblk_xxhp,1/),(/0,0,iblk_xxhp,0/),
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        end if
      end if
      if(r12op.eq.2.or.r12op.ge.3)then
        do while(associated(fl_t_r_pnt%next))
          fl_t_r_pnt => fl_t_r_pnt%next
        enddo
      
        ! find xp|pp and xx|pp blocks of R12
        occ = 0
        occ(IPART,1) = 1
        occ(IEXTR,1) = 1
        occ(IPART,2) = 2
        iblk_pxpp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op)
        occ = 0
        occ(IEXTR,1) = 2
        occ(IPART,2) = 2
        iblk_xxpp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op)

        if (.not.bar) then
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxr12,idxcpp12,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
     &         (/1,iblk_pxpp,1,1/),(/0,iblk_pxpp,0,0/),
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        else
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxcpp12,-idxr12,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
     &         (/1,1,iblk_pxpp,1/),(/0,0,iblk_pxpp,0/),
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        end if
        do while(associated(fl_t_r_pnt%next))
          fl_t_r_pnt => fl_t_r_pnt%next
        enddo
        if (.not.bar) then
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxr12,idxcpp12,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
     &         (/1,iblk_xxpp,1,1/),(/0,iblk_xxpp,0,0/),
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        else
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxcpp12,-idxr12,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
     &         (/1,1,iblk_xxpp,1/),(/0,0,iblk_xxpp,0/),
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        end if
      end if
      ! if BAR: no non-linear terms for projection (currently ...)
      if (r12op.ge.4.and..not.bar) then
        do while(associated(fl_t_r_pnt%next))
          fl_t_r_pnt => fl_t_r_pnt%next
        enddo
      
        ! find xp|pp and xx|pp blocks of R12
        occ = 0
        occ(IPART,1) = 1
        occ(IEXTR,1) = 1
        occ(IPART,2) = 2
        iblk_pxpp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op)
        occ = 0
        occ(IEXTR,1) = 2
        occ(IPART,2) = 2
        iblk_xxpp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op)

        call expand_op_product2(fl_t_r_pnt,idxsop,
     &       1.0d0,5,4,
     &       (/idxsop,idxr12,idxc12,idxc12,idxsop/),
     &       (/1      ,2     ,3     ,4     ,1     /),
     &       (/1,iblk_pxpp,1,1,1/),(/0,iblk_pxpp,0,0,0/),
     &       (/2,3,2,4/),2,
     &       0,0,
     &       0,0,
     &       op_info)

        do while(associated(fl_t_r_pnt%next))
          fl_t_r_pnt => fl_t_r_pnt%next
        enddo
        call expand_op_product2(fl_t_r_pnt,idxsop,
     &       1d0,5,4,
     &       (/idxsop,idxr12,idxc12,idxc12,idxsop/),
     &       (/1      ,2     ,3       ,4   ,1     /),
     &       (/1,iblk_xxpp,1,1,1/),(/0,iblk_xxpp,0,0,0/),
     &       (/2,3,2,4/),2,
     &       0,0,
     &       0,0,
     &       op_info)
      endif

      end
