*----------------------------------------------------------------------*
      integer function idx_next_target(tgt_info)
*----------------------------------------------------------------------*
*     return the index of the next target to be evaluated
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'
      include 'stdunit.h'

      integer, parameter ::
     &     maxcycles = 100 000 000,
     &     maxlevel  = 50,
     &     ntest = 000
      
      type(target_info), intent(in) ::
     &     tgt_info

      integer ::
     &     icycle, itgt, level
      type(target), pointer ::
     &     tgt

      level = 0
      icycle = 1
      idx_next_target = -1
      do itgt = 1, tgt_info%ntargets
        tgt => tgt_info%array(itgt)%tgt
        if (tgt%required) then

          if (ntest.ge.100) then
            write(luout,*) 'next required target: ',trim(tgt%name)
          end if

          idx_next_target = next_target_rec(tgt,tgt_info)

          if (ntest.ge.100) then
            if (idx_next_target.gt.0) then
              write(luout,*) 'need to make next: ',
     &             trim(tgt_info%array(idx_next_target)%tgt%name)
            end if
          end if

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

      if (ntest.ge.1000) then
        write(luout,*) icycle+1, level, 'checking dependecies of ',
     &       trim(tgt_in%name)
      end if
      
      icycle = icycle+1
      level  = level+1
      if (level.gt.maxlevel)
     &     call quit(1,'idx_next_target',
     &          'Too many nested dependencies ... '//
     &          'I suspect a circular dependency ... quitting')
      if (icycle.gt.maxcycles)
     &     call quit(1,'idx_next_target','infinite loop?')

      ! any dependencies?
      if (tgt_in%n_depends_on.gt.0) then
        ! check whether any target we depend on needs to be re-made
        do idep = 1, tgt_in%n_depends_on
          tgt => tgt_info%array(tgt_in%idx_depends_on(idep))%tgt
          idx = next_target_rec(tgt,tgt_info)
          if (idx.gt.0) then
            if (ntest.ge.1000) then
              write(luout,*) 'need to remake: ',
     &             trim(tgt_info%array(idx)%tgt%name)
            end if
            next_target_rec = idx
            level = level-1
            return
          end if
        end do
        ! check whether any target we depend on is more recent
        do idep = 1, tgt_in%n_depends_on
          idx = tgt_in%idx_depends_on(idep)
          if (tgt_info%last_mod(idx).gt.tgt_in%last_mod) then
            if (ntest.ge.1000) then
              write(luout,*) trim(tgt_info%array(idx)%tgt%name),
     &             ' is more recent: remake ',
     &             trim(tgt_in%name)
            end if
            next_target_rec = tgt_in%my_idx
            level = level-1
            return
          end if
        end do
      else
        ! no dependency? make only, if last_mod==-1
c dbg
c        print *,'no more dependency: ',trim(tgt_in%name)
c        print *,'last mod: ',tgt_in%last_mod
c dbg
        if (ntest.ge.1000) then
          write(luout,*) trim(tgt_in%name),': all depencies resolved'
          if (tgt_in%last_mod.lt.0)
     &         write(luout,*)
     &         ' target is still untouched, so I will make it!'
        end if
        if (tgt_in%last_mod.lt.0)
     &       next_target_rec = tgt_in%my_idx
      end if

c dbg
c      print *,'returning: ',next_target_rec
c dbg
      level = level-1
      return
      end function next_target_rec

      end
