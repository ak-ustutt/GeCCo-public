*----------------------------------------------------------------------*
      integer function idx_next_target(tgt_info)
*----------------------------------------------------------------------*
*     return the index of the next target to be evaluated
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'

      integer, parameter ::
     &     maxcycles = 1 000 000
      
      type(target_info), intent(in) ::
     &     tgt_info

      integer ::
     &     icycle, itgt
      type(target), pointer ::
     &     tgt

      icycle = 1
      idx_next_target = -1
      do itgt = 1, tgt_info%ntargets
        tgt => tgt_info%array(itgt)%tgt
        if (tgt%required) then
c dbg
c        print *,'next required target: ',trim(tgt%name)
c        print *,'    #### : ',tgt%my_idx
c dbg
          idx_next_target = next_target_rec(tgt,tgt_info)
c dbg
c          print *,'idx = ',idx_next_target
c dbg
          if (idx_next_target.gt.0) exit
        end if
      end do

      return

      contains

      recursive integer function next_target_rec(tgt_in,tgt_info)

      type(target), intent(in) ::
     &     tgt_in
      type(target_info), intent(in) ::
     &     tgt_info

      type(target), pointer ::
     &     tgt
      integer ::
     &     idep, idx

      next_target_rec = -1
      
      icycle = icycle+1
      if (icycle.gt.maxcycles)
     &     call quit(1,'idx_next_target','infinite loop?')

      ! any dependencies?
      if (tgt_in%n_depends_on.gt.0) then
        ! check whether any target we depend on needs to be re-made
        do idep = 1, tgt_in%n_depends_on
          tgt => tgt_info%array(tgt_in%idx_depends_on(idep))%tgt
          idx = next_target_rec(tgt,tgt_info)
          if (idx.gt.0) then
            next_target_rec = idx
            return
          end if
        end do
        ! check whether any target we depend on is more recent
        do idep = 1, tgt_in%n_depends_on
          idx = tgt_in%idx_depends_on(idep)
          if (tgt_info%last_mod(idx).gt.tgt_in%last_mod) then
            next_target_rec = tgt_in%my_idx
            return
          end if
        end do
      else
        ! no dependency? make only, if last_mod==-1
c dbg
c        print *,'no more dependency: ',trim(tgt_in%name)
c        print *,'last mod: ',tgt_in%last_mod
c dbg
        if (tgt_in%last_mod.lt.0)
     &       next_target_rec = tgt_in%my_idx
      end if

c dbg
c      print *,'returning: ',next_target_rec
c dbg
      return
      end function next_target_rec

      end
