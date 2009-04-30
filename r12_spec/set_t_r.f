      subroutine set_t_r(flist_t_r,bar,set_r_r,
     &     idxsop,idxtop,
     &     idxr12,idxr12x,idxc12,idxcpp12,
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
     &     bar, set_r_r
      integer, intent(in) ::
     &     idxsop, idxtop, idxr12, idxr12x, idxc12, idxcpp12,
     &     r12op,r12fix

      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     fl_t_r_pnt

      integer ::
     &     iblk_pphp, iblk_pxhp, iblk_xxhp, iblk_pxpp, iblk_xxpp,
     &     occ(ngastp,2)

      integer, external ::
     &     iblk_occ
      
      fl_t_r_pnt => flist_t_r
      call new_formula_item(fl_t_r_pnt,command_set_target_init,idxsop)
      fl_t_r_pnt => fl_t_r_pnt%next
c replaced by set_primitive_formula, because operator block versions
c of top and sop should match
c      call expand_op_product(fl_t_r_pnt,idxsop,
c     &       1d0,1,idxtop,-1,-1,
c     &       0,0,.false.,op_info)
      call set_primitive_formula(fl_t_r_pnt,idxtop,
     &     1d0,idxsop,.false.,op_info)
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
        iblk_pxhp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op,1)
        occ = 0
        occ(IEXTR,1) = 2
        occ(IHOLE,2) = 1
        occ(IPART,2) = 1
        iblk_xxhp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op,1)
      

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
     &         1.0d0,4,3,
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
      if((r12op.eq.1.or.r12op.ge.3).and.set_r_r)then
        do while(associated(fl_t_r_pnt%next))
          fl_t_r_pnt => fl_t_r_pnt%next
        enddo

        ! find hx|hp R12X
        occ = 0
        occ(IHOLE,2) = 1
        occ(IPART,2) = 1
        occ(IHOLE,1) = 2
        iblk_pphp = iblk_occ(occ,.false.,op_info%op_arr(idxr12x)%op,1)
c dbg
        print *,'iblk_pphp : ',iblk_pphp,idxr12x
        print *,'set to 1'
        iblk_pphp = 1
c dbg

        if (.not.bar) then
c dbg
          print *,'call'
c dbg
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxr12x,idxc12,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
c     &         (/1,iblk_pphp,1,1/),(/0,iblk_pphp,0,0/),
     &         -1,-1,
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        else
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1d0,4,3,
     &         (/idxsop,idxc12,-idxr12x,idxsop/),
     &         (/1      ,2     ,3       ,1     /),
     &         (/1,1,iblk_pphp,1/),(/0,0,iblk_pphp,0/),
     &         (/2,3/),1,
     &         0,0,
     &         0,0,
     &         op_info)
        end if

c        
c        ! find pp|hp, xp|hp and xx|hp blocks of R12X/R12
c        occ = 0
c        occ(IPART,1) = 2
c        occ(IHOLE,2) = 1
c        occ(IPART,2) = 1
c        iblk_pphp = iblk_occ(occ,.false.,op_info%op_arr(idxr12x)%op)
c        occ = 0
c        occ(IPART,1) = 1
c        occ(IEXTR,1) = 1
c        occ(IHOLE,2) = 1
c        occ(IPART,2) = 1
c        iblk_pxhp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op)
c        occ = 0
c        occ(IEXTR,1) = 2
c        occ(IHOLE,2) = 1
c        occ(IPART,2) = 1
c        iblk_xxhp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op)
c      
c
cc        if (.not.bar) then
cc          call expand_op_product2(fl_t_r_pnt,idxsop,
cc     &         1d0,4,3,
cc     &         (/idxsop,idxr12x,idxr12,idxsop/),
cc     &         (/1      ,2     ,3       ,1     /),
cc     &         (/1,iblk_pphp,1,1/),(/0,iblk_pphp,0,0/),
cc     &         (/2,3/),1,
cc     &         0,0,
cc     &         0,0,
cc     &         op_info)
cc        else
cc          call expand_op_product2(fl_t_r_pnt,idxsop,
cc     &         1d0,4,3,
cc     &         (/idxsop,-idxr12,-idxr12x,idxsop/),
cc     &         (/1      ,2     ,3       ,1     /),
cc     &         (/1,1,iblk_pphp,1/),(/0,0,iblk_pphp,0/),
cc     &         (/2,3/),1,
cc     &         0,0,
cc     &         0,0,
cc     &         op_info)
cc        end if
c        do while(associated(fl_t_r_pnt%next))
c          fl_t_r_pnt => fl_t_r_pnt%next
c        enddo
c        if (.not.bar) then
c          call expand_op_product2(fl_t_r_pnt,idxsop,
c     &         1d0,4,3,
c     &         (/idxsop,idxr12,idxr12,idxsop/),
c     &         (/1      ,2     ,3       ,1     /),
c     &         (/1,iblk_pxhp,1,1/),(/0,iblk_pxhp,0,0/),
c     &         (/2,3/),1,
c     &         0,0,
c     &         0,0,
c     &         op_info)
c        else
c          call expand_op_product2(fl_t_r_pnt,idxsop,
c     &         1d0,4,3,
c     &         (/idxsop,-idxr12,-idxr12,idxsop/),
c     &         (/1      ,2     ,3       ,1     /),
c     &         (/1,1,iblk_pxhp,1/),(/0,0,iblk_pxhp,0/),
c     &         (/2,3/),1,
c     &         0,0,
c     &         0,0,
c     &         op_info)
c        end if
c        do while(associated(fl_t_r_pnt%next))
c          fl_t_r_pnt => fl_t_r_pnt%next
c        enddo
c        if (.not.bar) then
c          call expand_op_product2(fl_t_r_pnt,idxsop,
c     &         1.0d0,4,3,
c     &         (/idxsop,idxr12,idxr12,idxsop/),
c     &         (/1      ,2     ,3       ,1     /),
c     &         (/1,iblk_xxhp,1,1/),(/0,iblk_xxhp,0,0/),
c     &         (/2,3/),1,
c     &         0,0,
c     &         0,0,
c     &         op_info)
c        else
c          call expand_op_product2(fl_t_r_pnt,idxsop,
c     &         1d0,4,3,
c     &         (/idxsop,-idxr12,-idxr12,idxsop/),
c     &         (/1      ,2     ,3       ,1     /),
c     &         (/1,1,iblk_xxhp,1/),(/0,0,iblk_xxhp,0/),
c     &         (/2,3/),1,
c     &         0,0,
c     &         0,0,
c     &         op_info)
c        end if
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
        iblk_pxpp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op,1)
        occ = 0
        occ(IEXTR,1) = 2
        occ(IPART,2) = 2
        iblk_xxpp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op,1)

        if (.not.bar) then
          call expand_op_product2(fl_t_r_pnt,idxsop,
     &         1.0d0,4,3,
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
     &         1.0d0,4,3,
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
        iblk_pxpp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op,1)
        occ = 0
        occ(IEXTR,1) = 2
        occ(IPART,2) = 2
        iblk_xxpp = iblk_occ(occ,.false.,op_info%op_arr(idxr12)%op,1)

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
