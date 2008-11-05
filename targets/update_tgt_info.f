*----------------------------------------------------------------------*
      subroutine update_tgt_info(tgt_info)
*----------------------------------------------------------------------*
*     setup or update direct acces pointer array for linked list
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_target_info.h'

      type(target_info) ::
     &     tgt_info
      integer ::
     &     itgt

      ! treat target array
      if (associated(tgt_info%array)) deallocate(tgt_info%array)

      if (tgt_info%ntargets.gt.0) then
        allocate(tgt_info%array(tgt_info%ntargets))
        call tgt_list2arr(tgt_info%list,
     &       tgt_info%array,tgt_info%ntargets)
      end if

      ! update last_mod table
      if (associated(tgt_info%last_mod)) deallocate(tgt_info%last_mod)
      allocate(tgt_info%last_mod(tgt_info%ntargets))
      do itgt = 1, tgt_info%ntargets
        tgt_info%last_mod(itgt) = tgt_info%array(itgt)%tgt%last_mod
      end do

      ! treat action array
      if (associated(tgt_info%act_array))
     &     deallocate(tgt_info%act_array)

      if (tgt_info%nactions.gt.0) then
        allocate(tgt_info%array(tgt_info%nactions))
        call act_list2arr(tgt_info%act_list,
     &       tgt_info%act_array,tgt_info%nactions)
      end if
      

      return
      end

