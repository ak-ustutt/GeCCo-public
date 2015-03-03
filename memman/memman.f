*----------------------------------------------------------------------*
      module memman
*----------------------------------------------------------------------*
*     manager for linear real(8) and integer slices
*     started march 2007 by andreas
*
*     use the wrapper routines
*       mem_init
*       mem_alloc_int
*       mem_alloc_real
*       mem_register
*       mem_dealloc
*       mem_setmark
*       mem_flushmark
*       mem_clean
*----------------------------------------------------------------------*

      implicit none

      include 'def_filinf.h'
*----------------------------------------------------------------------*
*     parameters
*----------------------------------------------------------------------*
      integer, parameter ::
     &     mem_maxname = 32

      integer, parameter ::
     &     membuffer_maxmax_slots = 128*1024,
     &     membuffer_maxname = 8,
     &     membuffer_stat_max = 500

*----------------------------------------------------------------------*
*     types
*----------------------------------------------------------------------*
      type mem_slice
        character ::
     &       name*(mem_maxname)
        integer ::
     &       type, len
        real(8), pointer ::
     &       xmem(:)
        integer, pointer ::
     &       imem(:)
        type(mem_slice), pointer ::
     &       prev, next
      end type mem_slice

      type mem_section
        character ::
     &       name*(mem_maxname)
        type(mem_slice), pointer ::
     &       head, tail
        type(mem_section), pointer ::
     &       prev, next
      end type mem_section

      type section_stack
        type(mem_section), pointer ::
     &     section
        type(section_stack), pointer ::
     &     prev, next
      end type section_stack

      type mem_slice_array
        type(mem_slice), pointer ::
     &     mem_slc
      end type mem_slice_array

      type membuffer

      character ::
     &       name*(mem_maxname)
      type(filinf), pointer ::
     &     ffbuf
      integer ::
     &     id,                  ! a unique ID
     &     max_buf_len,         ! maximum length of buffer (8byte words)
     &     cur_buf_len,         ! current length
     &     max_buf_slots,       ! maximum number of slots
     &     cur_buf_slots

      integer ::
     &     idxlast,
     &     last10(10)           ! the 10 last accesses
      integer, pointer ::
     &     idx4id(:)            ! an array for finding slots
      integer, pointer ::
     &     slot_info(:)         ! contains: [ID, status]

      type(mem_slice_array), pointer ::
     &     slot(:)

      type(membuffer), pointer ::
     &     prev, next

      end type membuffer

      type(mem_section), pointer ::
     &     mem_root, mem_cursection, mem_tail

      type(mem_slice), pointer ::
     &     mem_curslice

      type(section_stack), pointer ::
     &     mem_stack_head, mem_stack_tail

      type(membuffer), pointer ::
     &     mem_buf_head, mem_buf_tail, mem_buf_pnt

      integer, parameter ::
     &     mtyp_int = 1,
     &     mtyp_rl8 = 2,
     &     mtyp_reg = 3   ! register only (in r8 words)

      real(8), parameter ::
     &     over_warn = 0.01d0,
     &     over_err  = 0.10d0

      ! should be made private
      integer, parameter ::
     &     npad = 1,
     &     ipad = 1234567890
      real(8), parameter ::
     &     xpad = 1234567890.0d0

      integer ::
     &     mem_free, mem_total, irat,
     &     max_blk, max_mem
      character ::
     &     name_max*(2*mem_maxname+1)

*----------------------------------------------------------------------*
      contains
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
      subroutine memman_init(mem_free_init)
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'
      include 'ifc_baserout.h'

      integer, intent(in) ::
     &     mem_free_init
      integer ::
     &     istat

      allocate(mem_root,stat=istat)
      if (istat.ne.0) call memman_error('init',istat)

      mem_root%name = 'static'
      
      nullify(mem_root%head)
      nullify(mem_root%tail)
      nullify(mem_root%prev)
      nullify(mem_root%next)
      nullify(mem_buf_head)
      nullify(mem_buf_tail)

      mem_free = mem_free_init
      mem_total = mem_free_init
      mem_cursection => mem_root  ! pointer to current section
      mem_tail => mem_root        ! pointer to last section
      nullify(mem_curslice)

      if (iprlvl.gt.0)
     &     write(lulog,'(x,a,i12,a,f8.2,a)')
     &     'Memory set to ',mem_free_init,
     &     ' r8-words = (',dble(mem_free_init)/1024d0/128d0,' Mb)'

      irat = zirat()
      if (iprlvl.ge.2)
     &     write(lulog,*) 'real-word/integer-word ratio = ',irat

      max_blk = 0
      max_mem = 0
      name_max = ' '

      end subroutine

*----------------------------------------------------------------------*
      subroutine memman_clean()
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'
      include 'ifc_baserout.h'

      integer ::
     &     ifree

      ! free all sections
      do
        if (associated(mem_tail,mem_root)) exit
        ifree = memman_remsection()
      end do

      ! free all slices in root section
      do
        if (.not.associated(mem_curslice)) exit
        ifree = memman_dealloc()
      end do

      ! and free the memory root
      deallocate(mem_root)

      end subroutine

      subroutine memman_initslice()

      implicit none
      include 'stdunit.h'

      integer ::
     &     istat

      allocate(mem_cursection%head,stat=istat)
      if (istat.ne.0) call memman_error('init_slice',istat)
      mem_cursection%tail => mem_cursection%head
      mem_curslice => mem_cursection%head

      nullify(mem_curslice%xmem)
      nullify(mem_curslice%imem)
      nullify(mem_curslice%prev)
      nullify(mem_curslice%next)

      end subroutine

*----------------------------------------------------------------------*
      integer function memman_alloc(type,nalloc,name,xpnt,ipnt)
*----------------------------------------------------------------------*
*
*     allocate a memory slice in current section
*     return free memory left after allocate
*
*     if type == mtyp_reg, register only memory consumption
*     (e.g. user allocate of more compicated type)
*
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'

      integer, pointer, optional ::
     &     ipnt(:)
      real(8), pointer, optional ::
     &     xpnt(:)
      integer, intent(in) ::
     &     type, nalloc
      character, intent(in) ::
     &     name*(*)

      integer ::
     &     mem_reg, len
      real(8) ::
     &     over
      integer ::
     &     istat

      if (.not.associated(mem_curslice)) then
        call memman_initslice()
      else
        allocate(mem_curslice%next,stat=istat)
        if (istat.ne.0) call memman_error('allocate(ini)',istat)
        mem_curslice%next%prev => mem_curslice
        nullify(mem_curslice%next%next)
        mem_curslice => mem_curslice%next
        mem_cursection%tail => mem_curslice
      end if

      len = len_trim(name)
      if (len.gt.mem_maxname)
     &     call quit(1,'memman_alloc',
     &     'identifier too long "'//trim(name)//'"')
      if (nalloc.lt.0)
     &     call quit(1,'memman_alloc',
     &     'negative length in allocation ('//trim(name)//')')
      mem_curslice%name = name
      mem_curslice%type = type
      mem_curslice%len = nalloc

      select case(type)
      case(mtyp_int)
        mem_reg = (nalloc+2*npad)/irat
        if (mod(nalloc,irat).ne.0) mem_reg = mem_reg+1
      case(mtyp_rl8)
        mem_reg = nalloc+2*npad
      case(mtyp_reg)
        mem_reg = nalloc+2*npad
      case default
        call quit(1,'memman_alloc','illegal type')
      end select

      mem_free = mem_free - mem_reg

      if (mem_free.lt.0) then
        over = -dble(mem_free)/dble(mem_total)
        if (over.gt.over_err) then
          write(lulog,'(x,a,e8.1,a)') 'ERROR: memory exceeded by ',
     &         over*100d0,' %'
          ! print memory map here
          write(lulog,'(x,2a)') 'trying to allocate slice: ',
     &         trim(mem_curslice%name)
          write(lulog,'(x,a,i25)') 'size of requested slice:  ',mem_reg
          write(lulog,'(x,a,i25)') 'avaible free space:       ',
     &         mem_free+mem_reg
          write(lulog,'(x,a,i25)') 'memmax had been set to:   ',
     &         mem_total
          call memman_map(lulog,.true.)
          call quit(0,'memman_alloc','memory exceeded')
        else if (over.gt.over_warn) then
          call warn('memman_alloc','memory exeeded')
          write(lulog,'(x,a,e8.1,a)') 'WARNING: memory exceeded by ',
     &         over*100d0,' %'
          write(lulog,'(x,2a)') 'trying to allocate slice: ',
     &         trim(mem_curslice%name)
        end if
            
      end if

      select case(type)
      case(mtyp_int)
        allocate(mem_curslice%imem(1-npad:nalloc+npad),stat=istat)
        if (istat.ne.0) call memman_error('allocate(int)',istat)
        if (npad.gt.0) then
          mem_curslice%imem(1-npad:0) = ipad
          mem_curslice%imem(nalloc+1:nalloc+npad) = ipad
        end if
        if (.not.present(ipnt))
     &       call quit(1,'memman_alloc','ipnt not present')
        ipnt => mem_curslice%imem(1:nalloc)
      case(mtyp_rl8)
        allocate(mem_curslice%xmem(1-npad:nalloc+npad),stat=istat)
        if (istat.ne.0) 
     &       call memman_error('allocate(rl8) '//trim(name),istat)
        if (npad.gt.0) then
          mem_curslice%xmem(1-npad:0) = xpad
          mem_curslice%xmem(nalloc+1:nalloc+npad) = xpad
        end if
        if (.not.present(xpnt))
     &       call quit(1,'memman_alloc','xpnt not present')
        xpnt => mem_curslice%xmem(1:nalloc)
      end select

      mem_reg = nalloc
      if (type.eq.mtyp_int) mem_reg = nalloc/irat + mod(nalloc,irat)
      ! add here statistic on largest block
      if (mem_reg.gt.max_blk) then
        max_blk = mem_reg
        name_max = trim(mem_cursection%name)//'.'
     &           //trim(mem_curslice%name)
      end if
      if (mem_total-mem_free.gt.max_mem) then
        max_mem = mem_total-mem_free
      end if

      memman_alloc = mem_free

      end function

*----------------------------------------------------------------------*
      integer function memman_dealloc(name,section)
*----------------------------------------------------------------------*
*
*     deallocate a memory slice in current section, the latest one
*     with name 'name'
*     if name is not given, the latest slice is freed
*     if section is given, use section instead of current section
*     the mem_curslice pointer remains at the last slice in the current
*     section
*
*     return the free memory after dealloc
*
*     if type == mtyp_reg, register only memory consumption
*     (e.g. user allocate of more compicated type)
*
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'ifc_baserout.h'

      character, intent(in), optional ::
     &     name*(*)
      type(mem_section), pointer, optional ::
     &     section

      type(mem_slice), pointer ::
     &     slice, slice_next
      type(mem_section), pointer ::
     &     cursection
      logical ::
     &     on_last_slice, in_last_section
      integer ::
     &     nalloc, mem_reg

      if (present(section)) then
        if (.not.associated(section))
     &       call quit(1,'memman_dealloc','invalid section pointer')
        cursection => section
      else
        cursection => mem_cursection
      end if
      if (.not.associated(cursection%tail)) then
        call quit(1,'memman_dealloc',
     &       'nothing to deallocate in section "'//
     &       trim(mem_cursection%name)//'"')
      end if

      ! point to last slice in section
      slice => cursection%tail

      if (present(name)) then
        ! search for node
        do
          if (trim(slice%name).eq.trim(name)) exit
          if (.not.associated(slice%prev)) then
            call memman_map(lulog,.true.)
            call quit(1,'memman_dealloc',
     &         'could not find a node named: "'//trim(name)//'"')
          end if
          slice => slice%prev
        end do
      end if

      ! flag whether we are on last slice
      on_last_slice = associated(slice,cursection%tail)
c      in_last_section = associated(cursection,mem_cursection)
      in_last_section = associated(cursection,mem_tail)

      nalloc = slice%len
      if (nalloc.lt.0) then
        call memman_map(lulog,.true.)
        write(lulog,*) 'nalloc = ',nalloc,' ?'
        call quit(1,'memman_dealloc','fishy length')
      end if

      ! deallocate and register freed memory
      select case(slice%type)
      case(mtyp_int)
        if (.not.associated(slice%imem))
     &       call quit(1,'memman_dealloc','slice was not allocated: '//
     &       trim(slice%name))        
        ! check paddings
        if (npad.gt.0) then
          if (.not.cmpiarr(ipad,slice%imem(1-npad:0),npad).or.
     &        .not.cmpiarr(ipad,slice%imem(nalloc+1:nalloc+npad),npad)) 
     &    then 
            call memman_map(lulog,.true.)
            call quit(1,'memman_dealloc','range error for '//
     &                                        trim(slice%name))
          end if
        end if
        deallocate(slice%imem)
        mem_reg = (nalloc+2*npad)/irat
        if (mod(nalloc,irat).ne.0) mem_reg = mem_reg+1
      case(mtyp_rl8)
        if (.not.associated(slice%xmem))
     &       call quit(1,'memman_dealloc','slice was not allocated: '//
     &       trim(slice%name))
        ! check paddings
        if (npad.gt.0) then
          if (.not.cmpxarr(xpad,slice%xmem(1-npad:0),npad).or.
     &        .not.cmpxarr(xpad,slice%xmem(nalloc+1:nalloc+npad),npad)) 
     &    then
            call memman_map(lulog,.true.)
            call quit(1,'memman_dealloc','range error for '//
     &                                    trim(slice%name))
          end if
        end if
        deallocate(slice%xmem)
        mem_reg = nalloc+2*npad
      case(mtyp_reg)
        mem_reg = nalloc+2*npad
      case default
        call memman_map(lulog,.true.)
        call quit(1,'memman_dealloc','illegal type')
      end select

      mem_free = mem_free + mem_reg

      ! remove slice from list
      if (on_last_slice) then
        if (.not.associated(slice%prev)) then
          ! a) this was the last member of subsection
          !    subsection is empty now
          deallocate(cursection%head)
          if (in_last_section) nullify(mem_curslice)
          nullify(cursection%head)
          nullify(cursection%tail)
        else
          ! b) remove slice and move mem_curslice and tail pointers
          slice => slice%prev
          deallocate(slice%next)
          nullify(slice%next)
          cursection%tail => slice
          if (on_last_slice) mem_curslice => slice
        end if
      else
        if (.not.associated(slice%prev)) then
          ! c) remove first slice in section
          slice_next => slice%next
          deallocate(cursection%head)
          cursection%head => slice_next
          nullify(cursection%head%prev)
        else
          ! d) remove a middle slice
          slice%prev%next => slice%next
          slice%next%prev => slice%prev
          deallocate(slice)
        end if
      end if

      memman_dealloc = mem_free

      end function

*----------------------------------------------------------------------*
      integer function memman_addsection(name)
*----------------------------------------------------------------------*
*     add a new section to memory manager
*----------------------------------------------------------------------*
      implicit none

      character, intent(in) ::
     &     name*(*)

      integer ::
     &     len, istat

      if (.not.associated(mem_tail))
     &     call quit(1,'memman_addsection','memman not initialized?')

      allocate(mem_tail%next,stat=istat)
      if (istat.ne.0) call memman_error('addsection',istat)
      mem_tail%next%prev => mem_tail
      nullify(mem_tail%next%next)
      mem_tail => mem_tail%next
      ! the new section is also the current section
      mem_cursection => mem_tail

      len = len_trim(name)
      if (len.gt.mem_maxname)
     &     call quit(1,'memman_addsection',
     &     'identifier too long "'//trim(name)//'"')
      mem_tail%name = name

      nullify(mem_cursection%head)
      nullify(mem_cursection%tail)
      nullify(mem_curslice)

      memman_addsection = mem_free

      return
      end function

*----------------------------------------------------------------------*
      integer function memman_remsection(name)
*----------------------------------------------------------------------*
*     remove last section (or section named 'name') from memory manager
*     free all memory slices allocated in this section
*     the mem_cursection pointer is afterwards associated with the
*     last section in list
*     then mem_curslice pointer is afterwards associated with the last
*     slice in the last section in the list
*----------------------------------------------------------------------*
      implicit none

      character, intent(in), optional ::
     &     name*(*)

      type(mem_section), pointer ::
     &     section, nextnext
      integer ::
     &     len, ifree

      if (.not.associated(mem_cursection))
     &     call quit(1,'memman_remsection','memman not initialized?')

      ! default: latest section
      section => mem_cursection
      if (present(name)) then
        ! else: start at last section and look for name
        section => mem_tail
        ! search for node
        do
          if (trim(section%name).eq.trim(name)) exit
          if (.not.associated(section%prev))
     &         call quit(1,'memman_remsection',
     &         'could not find a node name: "'//trim(name)//'"')
          section => section%prev
        end do
      end if

      if (associated(section,mem_root))
     &     call quit(1,'memman_remsection','cannot remove root section')

      ! deallocate all slices in this section
      do
        ! if everthing is deallocated, section%head points to NULL
        if (.not.associated(section%head)) exit
        ifree = memman_dealloc(section=section)
      end do

      ! current section to remove
      if (associated(section,mem_cursection)) then
        if (associated(section,mem_tail))
     &       mem_tail => section%prev

        nextnext => mem_cursection%next ! remember pointer to next el.
        mem_cursection => mem_cursection%prev
        deallocate(mem_cursection%next)
        mem_cursection%next => nextnext
        if (associated(mem_cursection%tail)) then
          mem_curslice => mem_cursection%tail
        else
          nullify(mem_curslice)
        end if
      ! last section to remove?
      else if (associated(section,mem_tail)) then
        mem_tail => mem_tail%prev
        deallocate(mem_tail%next)
        nullify(mem_tail%next)
      else
        ! middle section to remove
        section%prev%next => section%next
        section%next%prev => section%prev
        deallocate(section)
      end if

      memman_remsection = mem_free

      return
      end function

*----------------------------------------------------------------------*
      integer function memman_set_cursection(name)
*----------------------------------------------------------------------*
*     set mem_cursection pointer to the latest section with name <name>
*     if name is not given, we set mem_cursetion to mem_tail
*----------------------------------------------------------------------*

      implicit none

      character, intent(in), optional ::
     &     name*(*)

      type(mem_section), pointer ::
     &     section

      if (.not.associated(mem_tail))
     &    call quit(1,'memman_set_cursection','memman not initialized?')

      mem_cursection =>  mem_tail
      
      if (present(name)) then
        do
          if (trim(mem_cursection%name).eq.trim(name)) exit
          if (.not.associated(mem_cursection%prev))
     &         call quit(1,'memman_set_cursection',
     &         'could not find a node name: "'//trim(name)//'"')
          mem_cursection => mem_cursection%prev
        end do
      end if

      if (associated(mem_cursection%tail)) then
        mem_curslice => mem_cursection%tail
      else
        nullify(mem_curslice)
      end if

      memman_set_cursection = mem_free

      return
      end function

*----------------------------------------------------------------------*
      subroutine memman_section_stack(push_pop)
*----------------------------------------------------------------------*
*     manage a stack with section pointers
*     push:  save mem_cursection pointer on stack
*     pop:   get mem_cursection pointer from stack
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     push_pop

      if (push_pop.gt.0) then
        ! ---------------
        !       push
        ! ---------------
        if (.not.associated(mem_stack_head)) then
          allocate(mem_stack_head)
          mem_stack_tail => mem_stack_head
          nullify(mem_stack_head%prev)
          nullify(mem_stack_head%next)
        else
          allocate(mem_stack_tail%next)
          mem_stack_tail%next%prev => mem_stack_tail
          mem_stack_tail => mem_stack_tail%next
          nullify(mem_stack_tail%next)
        end if
        mem_stack_tail%section => mem_cursection
      else
        ! ---------------
        !       pop
        ! ---------------
        if (.not.associated(mem_stack_tail))
     &       call quit(1,'memman_section_stack','nothing to pop')
        mem_cursection => mem_stack_tail%section
        if (associated(mem_cursection%tail)) then
          mem_curslice => mem_cursection%tail
        else
          nullify(mem_curslice)
        end if
        if (associated(mem_stack_tail,mem_stack_head)) then
          deallocate(mem_stack_head)
        else
          mem_stack_tail => mem_stack_tail%prev
          deallocate(mem_stack_tail%next)
          nullify(mem_stack_tail%next)
        end if

      end if

      return
      end subroutine
*----------------------------------------------------------------------*
      subroutine memman_map(lulog,check)
*----------------------------------------------------------------------*
*     print a memory map
*----------------------------------------------------------------------*

      implicit none
      include 'ifc_baserout.h'

      integer, intent(in) ::
     &     lulog
      logical,intent(in) ::
     &     check

      type(mem_section), pointer ::
     &     cursec
      type(mem_slice), pointer ::
     &     curslc

      logical ::
     &     patchk1, patchk2
      integer ::
     &     memsum, nalloc, ierr
      character ::
     &     namscr*(mem_maxname)

      cursec => mem_root

      write(lulog,'(/x,"+",76("-"),"+",/x,"|",33x,a,33x,"|",'//
     &     '/x,"+",76("-"),"+")') 'memory map'
      memsum = 0
      ierr = 0
      do
        write(lulog,'(3x,a)') trim(cursec%name) 
        if (associated(cursec%head)) then
          curslc => cursec%head
          do
            nalloc = curslc%len
            select case(curslc%type)
            case (mtyp_int)
              memsum = memsum + nalloc/irat + 2*npad
              if (mod(nalloc,irat).ne.0) memsum = memsum+1
              if (check) then
                patchk1 = cmpiarr(ipad,curslc%imem(1-npad:0),npad)
                patchk2 = cmpiarr(ipad,
     &                     curslc%imem(nalloc+1:nalloc+npad),npad)
              end if
            case (mtyp_rl8)
              memsum = memsum + nalloc + 2*npad
              if (check) then
                patchk1 = cmpxarr(xpad,curslc%xmem(1-npad:0),npad)
                patchk2 = cmpxarr(xpad,
     &                     curslc%xmem(nalloc+1:nalloc+npad),npad)
              end if
            case default
              memsum = memsum + nalloc + 2*npad
            end select
            namscr(1:mem_maxname) = ' '
            namscr = curslc%name
            if ((curslc%type.eq.1.or.curslc%type.eq.2).and.check) then
              if (patchk1.and.patchk2) then
                write(lulog,'(6x,a,x,i2,x,i10,x,i10,x,l,2x,l)')
     &               namscr(1:mem_maxname),curslc%type,
     &               curslc%len,memsum,patchk1,patchk2
              else
                write(lulog,'(3x,"!",2x,a,x,i2,x,i10,'//
     &               'x,i10,x,l,2x,l,x,"!")')
     &               namscr(1:mem_maxname),curslc%type,
     &               curslc%len,memsum,patchk1,patchk2
                ierr = ierr+1
             end if
            else
              write(lulog,'(6x,a,x,i2,x,i10,x,i10,x,"N/A")')
     &               namscr(1:mem_maxname),curslc%type,
     &               curslc%len,memsum
            end if
            if (.not.associated(curslc%next)) exit
            curslc => curslc%next

          end do
        end if

        if (.not.associated(cursec%next)) exit
        cursec => cursec%next

      end do

      if (ierr.gt.0)
     &     write(lulog,*) '!! range errors detected (see above) !!'
c     &     call quit(1,'memman_map',
c     &     'range errors detected (see above)')

      end subroutine

*----------------------------------------------------------------------*
      subroutine memman_check(lulog,label)
*----------------------------------------------------------------------*
*     print a memory map
*----------------------------------------------------------------------*

      implicit none
      include 'ifc_baserout.h'

      integer, intent(in) ::
     &     lulog
      character, intent(in) ::
     &     label*(*)

      type(mem_section), pointer ::
     &     cursec
      type(mem_slice), pointer ::
     &     curslc

      logical ::
     &     patchk1, patchk2, ok
      integer ::
     &     nalloc

      cursec => mem_root

      main_loop: do       
        if (associated(cursec%head)) then
          curslc => cursec%head
          do
            nalloc = curslc%len
            select case(curslc%type)
            case (mtyp_int)
              patchk1 = cmpiarr(ipad,curslc%imem(1-npad:0),npad)
              patchk2 = cmpiarr(ipad,
     &                   curslc%imem(nalloc+1:nalloc+npad),npad)
              ok = patchk1.and.patchk2
            case (mtyp_rl8)
              patchk1 = cmpxarr(xpad,curslc%xmem(1-npad:0),npad)
              patchk2 = cmpxarr(xpad,
     &                   curslc%xmem(nalloc+1:nalloc+npad),npad)
              ok = patchk1.and.patchk2
            case default
              ok = .true.
            end select

            if (.not.ok) exit main_loop

            if (.not.associated(curslc%next)) exit
            curslc => curslc%next

          end do
        end if

        if (.not.associated(cursec%next)) exit
        cursec => cursec%next

      end do main_loop

      if (.not.ok) then
        write(lulog,*) 'Errors detected at check-point: ',trim(label)
        call memman_map(lulog,.true.)
        call quit(1,'memman_check','Check failed!')
      end if

      end subroutine

      subroutine memman_stat(lulog)

      implicit none

      integer, intent(in) ::
     &     lulog

      write(lulog,'(x,"+",76("-"),"+")')
      if (max_mem.lt.1024*1024*1024*64) then
      write(lulog,'(3x,a,i14,a,f10.2,a)')
     &     'Maximum allocated memory: ',max_mem,
     &     ' real(8)-words (',dble(max_mem)/1024d0/128d0,' Mb)'
      write(lulog,'(3x,a,i14,a,f10.2,a)')
     &     'Largest memory block:     ',max_blk,
     &     ' real(8)-words (',dble(max_blk)/1024d0/128d0,' Mb)'
      else
      write(lulog,'(3x,a,i14,a,f10.2,a)')
     &     'Maximum allocated memory: ',max_mem,
     &     ' real(8)-words (',dble(max_mem)/1024d0/1024d0/128d0,' Gb)'
      write(lulog,'(3x,a,i14,a,f10.2,a)')
     &     'Largest memory block:     ',max_blk,
     &     ' real(8)-words (',dble(max_blk)/1024d0/1024d0/128d0,' Gb)'
      end if
      write(lulog,'(3x,a,a)')
     &     'Name of largest block:    ',trim(name_max)
      write(lulog,'(x,"+",76("-"),"+")')

      end subroutine

*----------------------------------------------------------------------*
*     some routines for managing buffers with restricted sizes follow
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
      subroutine memman_init_vbuffer(ffbuf,name,max_len,max_slots)
*----------------------------------------------------------------------*
*     initialize
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
c      include 'def_filinf.h'

      type(filinf), intent(inout), target ::
     &     ffbuf
      character, intent(in) ::
     &     name*(*)
      integer, intent(in) ::
     &     max_len,max_slots

      integer ::
     &     len, ifree, idx, actual_len, mem_buf_id

      len = len_trim(name)
      if (len.gt.membuffer_maxname)
     &     call quit(1,'memman_init_vbuffer',
     &     'identifier too long "'//trim(name)//'"')
      if (max_len.le.0.or.max_slots.le.0)
     &     call quit(1,'memman_init_vbuffer',
     &     'negative dimensions encountered')
      if (max_len.gt.mem_free.or.max_slots.gt.membuffer_maxmax_slots)
     &     call quit(1,'memman_init_vbuffer',
     &     'too much memory or too many slots requested')
      ! note: we need to implement a mechanism to reserve memory space

      ifree = memman_addsection(name)

      if (.not.associated(mem_buf_pnt)) then
        allocate(mem_buf_head)
        mem_buf_tail => mem_buf_head
        mem_buf_pnt  => mem_buf_head
        nullify(mem_buf_pnt%prev)
        nullify(mem_buf_pnt%next)
        mem_buf_id = 1
      else
        mem_buf_id = mem_buf_pnt%id
        do while (associated(mem_buf_pnt%next))
          mem_buf_pnt => mem_buf_pnt%next
          mem_buf_id = mem_buf_pnt%id
        end do
        mem_buf_id = mem_buf_id+1
        allocate(mem_buf_pnt%next)
        mem_buf_pnt%next%prev => mem_buf_pnt
        mem_buf_pnt => mem_buf_pnt%next
        nullify(mem_buf_pnt%next)
      end if

      mem_buf_pnt%ffbuf => ffbuf
      mem_buf_pnt%name = name
      mem_buf_pnt%id = mem_buf_id
      ffbuf%buf_id = mem_buf_id
c      actual_len = (max_len-1)/irat+1
      mem_buf_pnt%max_buf_len = max_len
      mem_buf_pnt%cur_buf_len = 0
      mem_buf_pnt%max_buf_slots = max_slots
      mem_buf_pnt%cur_buf_slots = 0
      mem_buf_pnt%idxlast = 0
      mem_buf_pnt%last10(1:10) = 0
      
      allocate(mem_buf_pnt%slot(max_slots))
      ifree = memman_alloc(mtyp_reg,max_slots*8,'slots')
      ifree = memman_alloc(mtyp_int,max_slots*2,'slot_info',
     &     ipnt=mem_buf_pnt%slot_info)
      
      do idx = 1, max_slots
        mem_buf_pnt%slot_info(2*idx-1) =  0   ! ID
        mem_buf_pnt%slot_info(2*idx  ) =  0   ! priority
      end do

      ifree = memman_alloc(mtyp_int,max_slots,'idx4id',
     &     ipnt=mem_buf_pnt%idx4id)
      mem_buf_pnt%idx4id(1:max_slots) = 0

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine memman_clean_vbuffer(name)
*----------------------------------------------------------------------*
*     clean up
*     still missing: a sync mechanism for the file
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
c      include 'def_membuffer.h'

      character, intent(in) ::
     &     name*(*)

      integer ::
     &     ifree
      type(membuffer), pointer ::
     &     mem_buf_rem

      if (.not.associated(mem_buf_tail))
     &     call quit(1,'memman_clean_vbuffer','buffer list seems empty')
      mem_buf_pnt => mem_buf_tail
      do while(trim(mem_buf_pnt%name).ne.trim(name))
        if (.not.associated(mem_buf_pnt%prev))
     &       call quit(1,'memman_clean_vbuffer',
     &       'did not find: '//trim(name))
        mem_buf_pnt => mem_buf_pnt%prev
      end do

      ifree = memman_remsection(mem_buf_pnt%name)
      deallocate(mem_buf_pnt%slot)

      nullify(mem_buf_rem)
      ! remove current link from list
      if (associated(mem_buf_pnt%prev)) then
        mem_buf_pnt%prev%next => mem_buf_pnt%next
        if (associated(mem_buf_pnt,mem_buf_tail))
     &       mem_buf_head => mem_buf_pnt%prev
        mem_buf_rem => mem_buf_pnt%prev
      end if

      if (associated(mem_buf_pnt%next)) then
        mem_buf_pnt%next%prev => mem_buf_pnt%prev
        if (associated(mem_buf_pnt,mem_buf_head))
     &       mem_buf_head => mem_buf_pnt%next
        if (.not.associated(mem_buf_rem))
     &       mem_buf_rem => mem_buf_pnt%next
      end if

      deallocate(mem_buf_pnt)

      mem_buf_pnt => mem_buf_rem
      if (.not.associated(mem_buf_pnt)) then
        nullify(mem_buf_head)
        nullify(mem_buf_tail)
      end if

      return
      end subroutine
*----------------------------------------------------------------------*
      subroutine memman_vbuffer_get_int(id_buf,idxst,idxnd,ibuf)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     id_buf, idxst, idxnd
      integer, intent(out) ::
     &     ibuf(*)

      integer ::
     &     un, lenr, irecst, irecnd, ioffrec1, idxrecl, nread, irec,
     &     ioffbf, idx, buf_id
      integer, pointer ::
     &     ibufpnt(:)
*----------------------------------------------------------------------*

      ! find correct buffer
      if (mem_buf_pnt%id.ne.id_buf) then
        mem_buf_pnt => mem_buf_tail
        do while(mem_buf_pnt%id.ne.id_buf)
          if (.not.associated(mem_buf_pnt%prev))
     &         call quit(1,'memman_vbuffer_get_int','unknown buffer')
          mem_buf_pnt => mem_buf_pnt%prev
        end do
      end if
      
      ! do not use structure elements directly
      un = mem_buf_pnt%ffbuf%unit
      lenr = mem_buf_pnt%ffbuf%reclen*irat
      buf_id = mem_buf_pnt%ffbuf%buf_id

      ! first and last record to read from
      irecst = (idxst-1)/lenr+1
      irecnd = (idxnd-1)/lenr+1
      
      ! offset in first record
      ioffrec1 = idxst-1 - (irecst-1)*lenr
      ! last index in last record
      idxrecl = idxnd - (irecnd-1)*lenr

      if (irecst.eq.irecnd) then
        ! special case -- only one record to read:
        nread = idxrecl-ioffrec1
        call memman_idx_bufblk(buf_id,irecst,mtyp_int,ipnt=ibufpnt)
        if (.not.associated(ibufpnt)) then
          call memman_new_bufblk(buf_id,irecst,lenr,
     &         mtyp_int,ipnt=ibufpnt)
          read(un,rec=irecst) ibufpnt(1:lenr)
        end if
        ibuf(1:nread) = ibufpnt(ioffrec1+1:ioffrec1+nread)
      else
        ! first record
        nread = lenr-ioffrec1
        call memman_idx_bufblk(buf_id,irecst,mtyp_int,ipnt=ibufpnt)
        if (.not.(associated(ibufpnt))) then
          call memman_new_bufblk(buf_id,irecst,lenr,
     &         mtyp_int,ipnt=ibufpnt)
          read(un,rec=irecst) ibufpnt(1:lenr)
        end if
        ibuf(1:nread) = ibufpnt(ioffrec1+1:ioffrec1+nread)
        ! 2nd to (last-1)st record
        ioffbf = nread
        do irec = irecst+1, irecnd-1
          call memman_idx_bufblk(buf_id,irec,mtyp_int,ipnt=ibufpnt)
          if (.not.(associated(ibufpnt))) then
            call memman_new_bufblk(buf_id,irec,lenr,
     &           mtyp_int,ipnt=ibufpnt)
            read(un,rec=irec) ibufpnt(1:lenr)
          end if
          ibuf(ioffbf+1:ioffbf+lenr) = ibufpnt(1:lenr)
          ioffbf = ioffbf+lenr
        end do
        ! last record
        call memman_idx_bufblk(buf_id,irecnd,mtyp_int,ipnt=ibufpnt)
        if (.not.(associated(ibufpnt))) then
          call memman_new_bufblk(buf_id,irecnd,lenr,
     &         mtyp_int,ipnt=ibufpnt)
          read(un,rec=irecnd) ibufpnt(1:lenr)
        end if
        ibuf(ioffbf+1:ioffbf+idxrecl) = ibufpnt(1:idxrecl)
      end if

      return
      end subroutine
*----------------------------------------------------------------------*
      subroutine memman_vbuffer_put_int(id_buf,idxst,idxnd,ibuf)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     id_buf, idxst, idxnd
      integer, intent(in) ::
     &     ibuf(*)

      integer ::
     &     un, lenr, irecst, irecnd, ioffrec1, idxrecl, nwrite, irec,
     &     ioffbf, idx, buf_id
      integer, pointer ::
     &     ibufpnt(:)
*----------------------------------------------------------------------*

      ! find correct buffer
      if (mem_buf_pnt%id.ne.id_buf) then
        mem_buf_pnt => mem_buf_tail
        do while(mem_buf_pnt%id.ne.id_buf)
          if (.not.associated(mem_buf_pnt%prev))
     &         call quit(1,'memman_vbuffer_put_int','unknown buffer')
          mem_buf_pnt => mem_buf_pnt%prev
        end do
      end if
      
      ! do not use structure elements directly
      un = mem_buf_pnt%ffbuf%unit
      lenr = mem_buf_pnt%ffbuf%reclen*irat
      buf_id = mem_buf_pnt%ffbuf%buf_id

      ! first and last record to read from
      irecst = (idxst-1)/lenr+1
      irecnd = (idxnd-1)/lenr+1
      
      ! offset in first record
      ioffrec1 = idxst-1 - (irecst-1)*lenr
      ! last index in last record
      idxrecl = idxnd - (irecnd-1)*lenr

      if (irecst.eq.irecnd) then
        nwrite = idxrecl-ioffrec1
        call memman_idx_bufblk(buf_id,irecst,mtyp_int,ipnt=ibufpnt,
     &       modify=.true.)
        if (.not.associated(ibufpnt)) then
          call memman_new_bufblk(buf_id,irecst,lenr,
     &         mtyp_int,ipnt=ibufpnt,modify=.true.)
          if (idxrecl.ne.lenr.or.ioffrec1.ne.1)
     &         read(un,rec=irecst,err=10) ibufpnt(1:lenr)
          goto 20
 10       ibufpnt(1:ioffrec1) = 0
          ibufpnt(idxrecl+1:lenr) = 0
 20       continue
        end if
        ibufpnt(ioffrec1+1:ioffrec1+nwrite) = ibuf(1:nwrite)
      else
        ! first record
        nwrite = lenr-ioffrec1
        call memman_idx_bufblk(buf_id,irecst,mtyp_int,ipnt=ibufpnt,
     &       modify=.true.)
        if (.not.(associated(ibufpnt))) then
          call memman_new_bufblk(buf_id,irecst,lenr,
     &         mtyp_int,ipnt=ibufpnt,modify=.true.)
          if (ioffrec1.ne.1)
     &         read(un,rec=irecst,err=11) ibufpnt(1:ioffrec1)
          goto 21
 11       ibufpnt(1:ioffrec1) = 0
 21       continue
        end if
        ibufpnt(ioffrec1+1:ioffrec1+nwrite) = ibuf(1:nwrite)
        ! 2nd to (last-1)st record
        ioffbf = nwrite
        do irec = irecst+1, irecnd-1
          call memman_idx_bufblk(buf_id,irec,mtyp_int,ipnt=ibufpnt,
     &         modify=.true.)
          if (.not.(associated(ibufpnt))) then
            call memman_new_bufblk(buf_id,irec,lenr,
     &           mtyp_int,ipnt=ibufpnt,modify=.true.)
          end if
          ibufpnt(1:lenr) = ibuf(ioffbf+1:ioffbf+lenr)
          ioffbf = ioffbf+lenr
        end do
        ! last record
        call memman_idx_bufblk(buf_id,irecnd,mtyp_int,ipnt=ibufpnt,
     &       modify=.true.)
        if (.not.(associated(ibufpnt))) then
          call memman_new_bufblk(buf_id,irecnd,lenr,
     &         mtyp_int,ipnt=ibufpnt,modify=.true.)
          if (idxrecl.ne.lenr)
     &         read(un,rec=irecnd,err=12) ibufpnt(1:lenr)
          goto 22
 12       ibufpnt(idxrecl+1:lenr) = 0
 22       continue
        end if
        ibufpnt(1:idxrecl) = ibuf(ioffbf+1:ioffbf+idxrecl)
        
      end if

      return
      end subroutine
*----------------------------------------------------------------------*
      subroutine memman_new_bufblk(id_buf,id_slot,length,type,ipnt,xpnt,
     &     modify)
*----------------------------------------------------------------------*
*     allocate a new slot in the buffer (and flush other buffer(s),
*     if not enough memory is left
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'ifc_baserout.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     id_buf, id_slot, length, type
      
      integer, pointer, optional ::
     &     ipnt(:)
      real(8), pointer, optional ::
     &     xpnt(:)

      logical, optional ::
     &     modify

      logical ::
     &     free_slc
      integer ::
     &     actual_len, n_usage_min, n_usage_max,
     &     idx_slot, idx, i_usage, ifree, id_cur, istat
      integer, pointer ::
     &     slot_info(:), last10(:)
      type(mem_slice), pointer ::
     &     curslice

      character ::
     &     name_slot*8

      if (ntest.ge.100) then
        write(lulog,*) '-----------------'
        write(lulog,*) 'memman_new_buffer'
        write(lulog,*) '-----------------'
        write(lulog,*) ' ID buffer = ',id_buf
        write(lulog,*) ' ID slot   = ',id_slot
        write(lulog,*) ' length = ',length
      end if

      ! find correct buffer
      if (mem_buf_pnt%id.ne.id_buf) then
        mem_buf_pnt => mem_buf_tail
        do while(mem_buf_pnt%id.ne.id_buf)
          if (.not.associated(mem_buf_pnt%prev))
     &         call quit(1,'memman_new_bufblk','unknown buffer')
          mem_buf_pnt => mem_buf_pnt%prev
        end do
      end if

      if (id_slot.gt.mem_buf_pnt%max_buf_slots)
     &     call quit(1,'memman_new_bufblk',
     &     'requested slot is out of range')


      select case(type)
      case(mtyp_int) 
        actual_len = (length-1)/irat+1
      case(mtyp_rl8) 
        actual_len = length
      case default
        call quit(1,'memman_new_bufblk','illegal type encountered')
      end select

      if (ntest.ge.100) then
        write(lulog,*) ' length in 8b-words  = ',actual_len
        write(lulog,*) ' currently in use    = ',mem_buf_pnt%cur_buf_len
        write(lulog,*) ' currently available = ',
     &       mem_buf_pnt%max_buf_len-mem_buf_pnt%cur_buf_len
      end if

      if (actual_len.gt.mem_buf_pnt%max_buf_len)
     &     call quit(1,'memman_new_bufblk',
     &     'requested array is too large')

c      if (memman_idx_mem_buf_pnt(mem_buf_pnt,id_buf).ge.0)
c     &     call quit(1,'memman_new_mem_buf_pnt','ID is already in use!')

      slot_info => mem_buf_pnt%slot_info
      last10 => mem_buf_pnt%last10
      if (mem_buf_pnt%cur_buf_len+actual_len.gt.
     &    mem_buf_pnt%max_buf_len.or.
     &    mem_buf_pnt%cur_buf_slots.eq.mem_buf_pnt%max_buf_slots) then
        if (ntest.ge.100) then
          write(lulog,*) 'looking for some buffers to free:'
          write(lulog,*) ' last10 array: '
          write(lulog,'(x,10i5)') mem_buf_pnt%last10(1:10)
        end if
        ! 1) look for buffers, not accessed the last 10 times
        n_usage_min = huge(n_usage_max)
        n_usage_max = 0
        do idx = 1, mem_buf_pnt%max_buf_slots
          id_cur = slot_info(2*idx-1)
          istat  = slot_info(2*idx)
          if (id_cur.gt.0.and.abs(istat).lt.membuffer_stat_max) then
            n_usage_min = min(n_usage_min,abs(istat))
            n_usage_max = max(n_usage_max,abs(istat))
            if (imltlist(idx,last10,10,1).eq.0) then
              if (ntest.ge.100) then
                write(lulog,*) 'freeing slot: ',idx
              end if
              curslice => mem_buf_pnt%slot(idx)%mem_slc
              free_slc = length.ne.curslice%len
              if (istat.lt.0)
     &             call memman_write_buffer_slot(id_buf,idx,id_cur)
              call memman_release_buffer_slot(id_buf,idx,free_slc)
              idx_slot = idx
              if (mem_buf_pnt%cur_buf_len+actual_len.le.
     &            mem_buf_pnt%max_buf_len)
     &             exit
            end if
          end if
        end do
        if (ntest.ge.100) then
          write(lulog,*) 'free memory after 1st round: ',
     &         mem_buf_pnt%max_buf_len-mem_buf_pnt%cur_buf_len
        end if
        ! 2) remove mem_buf_pnts according to number of accesses
        if (mem_buf_pnt%cur_buf_len+actual_len.gt.
     &      mem_buf_pnt%max_buf_len) then
          if (n_usage_min.gt.n_usage_max)
     &         call quit(1,'memman_new_bufblk','inconsistency')
          usage_loop: do i_usage = n_usage_min, n_usage_max
            do idx = 1, mem_buf_pnt%max_buf_slots
              id_cur = slot_info(2*idx-1)
              istat  = slot_info(2*idx)
              if (id_cur.gt.0.and.
     &            abs(istat).eq.i_usage) then
                if (ntest.ge.100) then
                  write(lulog,*) 'freeing slot: ',idx
                end if
                curslice => mem_buf_pnt%slot(idx)%mem_slc
                free_slc = length.ne.curslice%len
                if (istat.lt.0)
     &               call memman_write_buffer_slot(id_buf,idx,id_cur)
                call memman_release_buffer_slot(id_buf,idx,free_slc)
                idx_slot = idx
                if (mem_buf_pnt%cur_buf_len+actual_len.le.
     &              mem_buf_pnt%max_buf_len)
     &               exit usage_loop
              end if
            end do
          end do usage_loop
          if (ntest.ge.100) then
            write(lulog,*) 'free memory after 2nd round: ',
     &           mem_buf_pnt%max_buf_len-mem_buf_pnt%cur_buf_len
          end if
        end if
      else
        ! signal that memory must be allocated
        free_slc = .true.
        ! look for gaps
        do idx = 1, mem_buf_pnt%max_buf_slots
          id_cur = slot_info(2*idx-1)
          if (id_cur.eq.0) then
            idx_slot = idx
            exit
          end if
        end do
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'next free slot: ',idx_slot
        write(lulog,*) 'current buffer length: ',mem_buf_pnt%cur_buf_len
      end if

      ! buffer usage:
      mem_buf_pnt%cur_buf_slots = mem_buf_pnt%cur_buf_slots+1
      mem_buf_pnt%cur_buf_len = mem_buf_pnt%cur_buf_len+actual_len

      ! set up slot entries
      mem_buf_pnt%slot_info(2*idx_slot-1) = id_slot
      mem_buf_pnt%slot_info(2*idx_slot) = 1
      if (present(modify)) then
        if (modify) mem_buf_pnt%slot_info(2*idx_slot) = -1
      end if
      mem_buf_pnt%idxlast = mod(mem_buf_pnt%idxlast,10)+1
      mem_buf_pnt%last10(mem_buf_pnt%idxlast) = idx_slot
c      mem_buf_pnt%slot(idx_slot)%length = actual_len
      if (free_slc) then
        ! re-allocate memory for buffer slot
        write(name_slot,'("sl",i6.6)') idx_slot
        call memman_section_stack(+1) !push
        ifree = memman_set_cursection(mem_buf_pnt%name)
        if (type.eq.mtyp_int) then
          ifree = memman_alloc(mtyp_int,length,name_slot,
     &         ipnt=ipnt)
        else
          ifree = memman_alloc(mtyp_rl8,length,name_slot,
     &         xpnt=xpnt)
        end if
        mem_buf_pnt%slot(idx_slot)%mem_slc => mem_curslice
        call memman_section_stack(-1) !pop
      else
        ! re-use old memory
        mem_buf_pnt%slot(idx_slot)%mem_slc => curslice
        if (type.eq.mtyp_int) then
          ipnt => curslice%imem
        else
          xpnt => curslice%xmem
        end if
      end if

      if (ntest.ge.500) then
        write(lulog,*) 'current memory map'
        call memman_map(lulog,.true.)
      end if
      if (ntest.ge.100) then
        write(lulog,*) 'buffer length on exit: ',mem_buf_pnt%cur_buf_len
      end if

      mem_buf_pnt%idx4id(id_slot)=idx_slot

      return
      end subroutine
*----------------------------------------------------------------------*
      subroutine memman_write_buffer_slot(id_buf,idx_slot,idx_rec)
*----------------------------------------------------------------------*
*     internal routine: write a buffer slot to disc
*----------------------------------------------------------------------*
      implicit none
      
      integer, intent(in) ::
     &     id_buf, idx_slot, idx_rec

      integer ::
     &     unit, type, length
      type(mem_slice), pointer ::
     &     curslice

      if (mem_buf_pnt%id.ne.id_buf)
     &     call quit(1,'memman_write_buffer_slot',
     &     'mem_buf_pnt must point to the correct buffer on input')

      mem_buf_pnt%slot_info(2*idx_slot) =
     &     abs(mem_buf_pnt%slot_info(2*idx_slot))
      curslice => mem_buf_pnt%slot(idx_slot)%mem_slc
      unit = mem_buf_pnt%ffbuf%unit
      if (unit.le.0)
     &     call quit(1,'memman_write_buffer_slot',
     &     'file not open: '//trim(mem_buf_pnt%ffbuf%name))
      type = curslice%type
      length = curslice%len
      if (type.eq.mtyp_int) then  
        write(unit,rec=idx_rec)
     &       curslice%imem(1:length)
      else if (type.eq.mtyp_rl8) then  
        write(unit,rec=idx_rec)
     &       curslice%xmem(1:length)
      else
        call quit(1,'memman_write_buffer_slot','illegal type')
      end if

      return
      end subroutine
*----------------------------------------------------------------------*
      subroutine memman_release_buffer_slot(id_buf,idx_slot,free_slc,
     &     modify)
*----------------------------------------------------------------------*
*     internal routine: remove a buffer slot
*----------------------------------------------------------------------*
      implicit none
      
      integer, intent(in) ::
     &     id_buf, idx_slot
      logical, intent(in) ::
     &     free_slc
      logical, optional ::
     &     modify

      integer ::
     &     idx, ifree, lenbuf, id_slot, type
      type(mem_slice), pointer ::
     &     curslice

      character ::
     &     name_slot*8

      if (mem_buf_pnt%id.ne.id_buf)
     &     call quit(1,'memman_release_buffer_slot',
     &     'mem_buf_pnt must point to the correct buffer on input')

      curslice => mem_buf_pnt%slot(idx_slot)%mem_slc
      ! save these two entries for post-processing (see below)
      lenbuf=curslice%len
      type = curslice%type
      if (type.eq.mtyp_int)  lenbuf = (lenbuf-1)/irat+1 
      id_slot = mem_buf_pnt%slot_info(2*idx_slot-1)
      if (free_slc) then
        ! remove allocated memory
        write(name_slot,'("sl",i6.6)') idx_slot
        ! switch to correct memory section and deallocate
        call memman_section_stack(+1) !push
        ifree = memman_set_cursection(mem_buf_pnt%name)
        ifree = memman_dealloc(name_slot)
        call memman_section_stack(-1) !pop
      end if
      ! update buffer usage
      mem_buf_pnt%cur_buf_len = mem_buf_pnt%cur_buf_len-lenbuf
      mem_buf_pnt%cur_buf_slots = mem_buf_pnt%cur_buf_slots-1
      ! reset slot entries
      mem_buf_pnt%slot_info(2*idx_slot-1) = 0
      mem_buf_pnt%slot_info(2*idx_slot) = 0
      nullify(mem_buf_pnt%slot(idx_slot)%mem_slc)

      ! remove entry from array of last 10 accesses to buffer
      do idx = 1, 10
        if (mem_buf_pnt%last10(idx).eq.idx_slot)
     &       mem_buf_pnt%last10(idx)=0
      end do

      ! remove entry from index array
      mem_buf_pnt%idx4id(id_slot) = 0

      return
      end subroutine
*----------------------------------------------------------------------*
      subroutine memman_idx_bufblk(id_buf,id_slot,type,ipnt,xpnt,modify)
*----------------------------------------------------------------------*
*     return a pointer to the slot associated with ID id_buf
*     NULL, if nothing found
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     id_buf, id_slot, type

      integer, pointer, optional ::
     &     ipnt(:)
      real(8), pointer, optional ::
     &     xpnt(:)
      logical, optional ::
     &     modify

      integer ::
     &     idx_slot, idx, istat, isgn
      integer, pointer ::
     &     slot_info(:)
      type(mem_slice), pointer ::
     &     curslice

      if (ntest.ge.100) then
        write(lulog,*) '-----------------'
        write(lulog,*) 'memman_idx_buffer'
        write(lulog,*) '-----------------'
        write(lulog,*) ' ID buffer: ',id_buf
        write(lulog,*) ' ID slot  : ',id_slot
      end if

      ! find correct buffer
      if (mem_buf_pnt%id.ne.id_buf) then
        mem_buf_pnt => mem_buf_tail
        do while(mem_buf_pnt%id.ne.id_buf)
          if (.not.associated(mem_buf_pnt%prev))
     &         call quit(1,'memman_idx_bufblk','unknown buffer')
          mem_buf_pnt => mem_buf_pnt%prev
        end do
      end if

      if (id_slot.gt.mem_buf_pnt%max_buf_slots) then
        write(lulog,*) ' ID buffer: ',id_buf
        write(lulog,*) ' ID slot  : ',id_slot
        write(lulog,*) ' ID slot(max): ',mem_buf_pnt%max_buf_slots
        write(lulog,*) ' modify set? ',present(modify)
        call print_vbuffer_int(lulog,id_buf,.true.)
        call memman_map(lulog,.true.)
        call quit(1,'memman_idx_bufblk',
     &     'requested slot is out of range')
      end if

      if (present(ipnt)) nullify(ipnt)
      if (present(xpnt)) nullify(xpnt)
      
      idx_slot = mem_buf_pnt%idx4id(id_slot)
      slot_info => mem_buf_pnt%slot_info
      if (idx_slot.gt.0) then
        istat = slot_info(2*idx_slot)
        isgn = 1
        if (istat.lt.0) isgn = -1
        slot_info(2*idx_slot) =
     &       isgn*(min(membuffer_stat_max,abs(istat)+1))
        if (present(modify).and.isgn.eq.1) then
          if (modify) slot_info(2*idx_slot) = -slot_info(2*idx_slot)
        end if
        mem_buf_pnt%idxlast = mod(mem_buf_pnt%idxlast,10)+1
        mem_buf_pnt%last10(mem_buf_pnt%idxlast) = idx_slot
        curslice => mem_buf_pnt%slot(idx_slot)%mem_slc
        if (curslice%type.eq.mtyp_int) then
          ipnt => curslice%imem
        else
          xpnt => curslice%xmem
        end if
      end if

      if (ntest.ge.100) then
        if (idx_slot.le.0) then          
          write(lulog,*) 'nothing found'
        else
          write(lulog,*) 'idx_slot = ',idx_slot
        end if
        write(lulog,*) 'current last10 array:'
        write(lulog,'(x,10i5)') mem_buf_pnt%last10(1:10)
      end if

      return
      end subroutine

      subroutine print_vbuffer_int(lulog,id_buf,short)
      implicit none

      integer, intent(in) ::
     &     lulog, id_buf
      logical, intent(in) ::
     &     short

      integer ::
     &     idx, irecmax, maxlen, length, unit 
      integer, pointer ::
     &     ibuf(:)

      ! find correct buffer
      if (mem_buf_pnt%id.ne.id_buf) then
        mem_buf_pnt => mem_buf_tail
        do while(mem_buf_pnt%id.ne.id_buf)
          if (.not.associated(mem_buf_pnt%prev))
     &         call quit(1,'print_vbuffer_int','unknown buffer')
          mem_buf_pnt => mem_buf_pnt%prev
        end do
      end if

      write(lulog,*) 'Info on buffer: ',trim(mem_buf_pnt%name)

      write(lulog,*) 'in-core length (current/max): ',
     &     mem_buf_pnt%cur_buf_len,
     &     mem_buf_pnt%max_buf_len

      write(lulog,*) 'last10 array: ',mem_buf_pnt%last10(1:10)
      write(lulog,*) 'idx4id array: ',
     &     mem_buf_pnt%idx4id(1:mem_buf_pnt%max_buf_slots)

      write(lulog,*) 'contents in-core:'
      irecmax = 1
      maxlen = 1
      do idx = 1, mem_buf_pnt%max_buf_slots
        write(lulog,*) 'slot #',idx
        write(lulog,*) 'ID, status: ',
     &       mem_buf_pnt%slot_info(idx*2-1),
     &       mem_buf_pnt%slot_info(idx*2)
        if (short) cycle
        irecmax = max(irecmax,mem_buf_pnt%slot_info(idx*2-1))
        if (mem_buf_pnt%slot_info(idx*2-1).gt.0) then
          length = mem_buf_pnt%slot(idx)%mem_slc%len
          maxlen = max(maxlen,length)
          call wrtimat2(mem_buf_pnt%slot(idx)%mem_slc%imem(1:length),
     &         1,length,1,length)
        end if
      end do

      write(lulog,*) 'Contents on disc:'
      unit = mem_buf_pnt%ffbuf%unit
      allocate(ibuf(maxlen))
      do idx = 1, irecmax
        read(unit,rec=idx,err=10) ibuf(1:maxlen)
        write(lulog,*) idx,'-> on disk'
        if (short) cycle
        call wrtimat2(ibuf,1,maxlen,1,maxlen)
        cycle
 10     write(lulog,*) idx,'-> empty record'
      end do
      deallocate(ibuf)

      return
      end subroutine

      subroutine memman_error(msg,stat)

      implicit none
      include 'stdunit.h'

      integer, intent(in) ::
     &     stat
      character(len=*), intent(in) ::
     &     msg

      write(lulog,*) 'ERROR in memory manager: '
      write(lulog,*) 'status: ',stat,' location: ',trim(msg)

      call memman_map(lulog,.false.)

      call quit(0,'memman','error in memory manager')

      end subroutine

      end module
