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

      integer, parameter ::
     &     mem_maxname = 16

      type mem_slice
        character ::
     &       name*(mem_maxname)
        integer ::
     &       type, len
        real(8), allocatable ::
     &       xmem(:)
        integer, allocatable ::
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

      type(mem_section), pointer ::
     &     mem_root, mem_cursection

      type(mem_slice), pointer ::
     &     mem_curslice

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

      allocate(mem_root)
      mem_root%name = 'static'
      
      nullify(mem_root%head)
      nullify(mem_root%tail)
      nullify(mem_root%prev)
      nullify(mem_root%next)

      mem_free = mem_free_init
      mem_total = mem_free_init
      mem_cursection => mem_root
      nullify(mem_curslice)

      if (iprlvl.gt.0)
     &     write(luout,'(x,a,i12,a,f8.2,a)')
     &     'Memory set to ',mem_free_init,
     &     ' r8-words = (',dble(mem_free_init)/1024d0/128d0,' Mb)'

      irat = zirat()
      if (iprlvl.ge.2)
     &     write(luout,*) 'real-word/integer-word ratio = ',irat

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
        if (associated(mem_cursection,mem_root)) exit
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

      allocate(mem_cursection%head)
      mem_cursection%tail => mem_cursection%head
      mem_curslice => mem_cursection%head

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

      integer, pointer, intent(out), optional ::
     &     ipnt(:)
      real(8), pointer, intent(out), optional ::
     &     xpnt(:)
      integer, intent(in) ::
     &     type, nalloc
      character, intent(in) ::
     &     name*(*)

      integer ::
     &     mem_reg, len
      real(8) ::
     &     over

      if (.not.associated(mem_curslice)) then
        call memman_initslice()
      else
        allocate(mem_curslice%next)
        mem_curslice%next%prev => mem_curslice
        nullify(mem_curslice%next%next)
        mem_curslice => mem_curslice%next
        mem_cursection%tail => mem_curslice
      end if

      len = len_trim(name)
      if (len.gt.mem_maxname)
     &     call quit(1,'memman_alloc',
     &     'identifier too long "'//trim(name)//'"')
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
        call quit(1,'mem_alloc','illegal type')
      end select

      mem_free = mem_free - mem_reg

      if (mem_free.lt.0) then
        over = -dble(mem_free)/dble(mem_total)
        if (over.gt.over_err) then
          write(luout,'(x,a,e8.1,a)') 'ERROR: memory exceeded by ',
     &         over*100d0,' %'
          ! print memory map here
          call quit(0,'memman','memory exceeded')
        else if (over.gt.over_warn) then
          write(luout,'(x,a,e8.1,a)') 'WARNING: memory exceeded by ',
     &         over*100d0,' %'
        end if
            
      end if

      select case(type)
      case(mtyp_int)
        allocate(mem_curslice%imem(1-npad:nalloc+npad))
        if (npad.gt.0) then
          mem_curslice%imem(1-npad:0) = ipad
          mem_curslice%imem(nalloc+1:nalloc+npad) = ipad
        end if
        if (.not.present(ipnt))
     &       call quit(1,'mem_alloc','ipnt not present')
        ipnt => mem_curslice%imem(1:nalloc)
      case(mtyp_rl8)
        allocate(mem_curslice%xmem(1-npad:nalloc+npad))
        if (npad.gt.0) then
          mem_curslice%xmem(1-npad:0) = xpad
          mem_curslice%xmem(nalloc+1:nalloc+npad) = xpad
        end if
        if (.not.present(xpnt))
     &       call quit(1,'mem_alloc','xpnt not present')
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
*     allocate a memory slice in current section, the latest one
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
      type(mem_section), pointer, intent(in), optional ::
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

c dbg
c      print *,'deallocating'
c dbg
      ! point to last slice in section
      slice => cursection%tail

      if (present(name)) then
        ! search for node
        do
c dbg
c          print *,'s: present: ',trim(slice%name),' ',trim(name),
c     &         trim(slice%name).eq.trim(name),associated(slice%prev),
c     &         len_trim(slice%name),len_trim(name)
c dbg
          if (trim(slice%name).eq.trim(name)) exit
          if (.not.associated(slice%prev))
     &         call quit(1,'memman_dealloc',
     &         'could not find a node named: "'//trim(name)//'"')
          slice => slice%prev
        end do
c dbg
c        print *,'node found'
c dbg
      end if

      ! flag whether we are on last slice
      on_last_slice = associated(slice,cursection%tail)
      in_last_section = associated(cursection,mem_cursection)

      nalloc = slice%len
      if (nalloc.lt.1) then
        write(luout,*) 'nalloc = ',nalloc,' ?'
        call quit(1,'memman_dealloc','fishy length')
      end if

      ! deallocate and register freed memory
      select case(slice%type)
      case(mtyp_int)
        if (.not.allocated(slice%imem))
     &       call quit(1,'memman_dealloc','slice was not allocated: '//
     &       trim(slice%name))        
        ! check paddings
        if (npad.gt.0) then
          if (.not.cmpiarr(ipad,slice%imem(1-npad),npad).or.
     &        .not.cmpiarr(ipad,slice%imem(nalloc+1),npad)) then
            call memman_map(luout,.true.)
            call quit(1,'memman','range error for '//trim(slice%name))
          end if
        end if
        deallocate(slice%imem)
        mem_reg = (nalloc+2*npad)/irat
        if (mod(nalloc,irat).ne.0) mem_reg = mem_reg+1
      case(mtyp_rl8)
        if (.not.allocated(slice%xmem))
     &       call quit(1,'memman_dealloc','slice was not allocated: '//
     &       trim(slice%name))
        ! check paddings
        if (npad.gt.0) then
          if (.not.cmpxarr(xpad,slice%xmem(1-npad),npad).or.
     &        .not.cmpxarr(xpad,slice%xmem(nalloc+1),npad)) then
            call memman_map(luout,.true.)
            call quit(1,'memman','range error for '//trim(slice%name))
          end if
        end if
        deallocate(slice%xmem)
        mem_reg = nalloc+2*npad
      case(mtyp_reg)
        mem_reg = nalloc+2*npad
      case default
        call quit(1,'mem_dealloc','illegal type')
      end select

      mem_free = mem_free + mem_reg
c dbg
c        print *,'successfully deallocated'
c dbg

      ! remove slice from list
      if (on_last_slice) then
        if (.not.associated(slice%prev)) then
          ! a) this was the last member of subsection
          !    subsection is empty now
c dbg
c          print *,'section a)'
c dbg
          deallocate(cursection%head)
          if (in_last_section) nullify(mem_curslice)
          nullify(cursection%head)
          nullify(cursection%tail)
        else
c dbg
c          print *,'section b)'
c dbg
          ! b) remove slice and move mem_curslice and tail pointers
          slice => slice%prev
          deallocate(slice%next)
          nullify(slice%next)
          cursection%tail => slice
          if (in_last_section) mem_curslice => slice
        end if
      else
        if (.not.associated(slice%prev)) then
c dbg
c          print *,'section c)'
c dbg
          ! c) remove first slice in section
          slice_next => slice%next
          deallocate(cursection%head)
          cursection%head => slice_next
          nullify(cursection%head%prev)
        else
c dbg
c          print *,'section d)'
c dbg
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
     &     len

      if (.not.associated(mem_cursection))
     &     call quit(1,'memman_addsection','memman not initialized?')

      allocate(mem_cursection%next)
      mem_cursection%next%prev => mem_cursection
      nullify(mem_cursection%next%next)
      mem_cursection => mem_cursection%next

      len = len_trim(name)
      if (len.gt.mem_maxname)
     &     call quit(1,'memman_addsection',
     &     'identifier too long "'//trim(name)//'"')
      mem_cursection%name = name

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
     &     section
      integer ::
     &     len, ifree

      if (.not.associated(mem_cursection))
     &     call quit(1,'memman_remsection','memman not initialized?')

      ! start at last section
      section => mem_cursection
      if (present(name)) then
        ! search for node
        do
c dbg
c          print *,'s: present: ',trim(section%name),' ',trim(name),
c     &        trim(section%name).eq.trim(name),associated(section%prev),
c     &         len_trim(section%name),len_trim(name)
c dbg
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

      ! last section to remove?
      if (associated(section,mem_cursection)) then
        mem_cursection => mem_cursection%prev
        deallocate(mem_cursection%next)
        nullify(mem_cursection%next)
        if (associated(mem_cursection%tail)) then
          mem_curslice => mem_cursection%tail
        else
          nullify(mem_curslice)
        end if
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
      subroutine memman_map(luout,check)
*----------------------------------------------------------------------*

      implicit none
      include 'ifc_baserout.h'

      integer, intent(in) ::
     &     luout
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

      write(luout,'(/x,"+",76("-"),"+",/x,"|",33x,a,33x,"|",'//
     &     '/x,"+",76("-"),"+")') 'memory map'
      memsum = 0
      ierr = 0
      do
        write(luout,'(3x,a)') trim(cursec%name) 
        if (associated(cursec%head)) then
          curslc => cursec%head
          do
            nalloc = curslc%len
            select case(curslc%type)
            case (mtyp_int)
              memsum = memsum + nalloc/irat + 2*npad
              if (mod(nalloc,irat).ne.0) memsum = memsum+1
              if (check) then
                patchk1 = cmpiarr(ipad,curslc%imem(1-npad),npad)
                patchk2 = cmpiarr(ipad,curslc%imem(nalloc+1),npad)
              end if
            case (mtyp_rl8)
              memsum = memsum + nalloc + 2*npad
              if (check) then
                patchk1 = cmpxarr(xpad,curslc%xmem(1-npad),npad)
                patchk2 = cmpxarr(xpad,curslc%xmem(nalloc+1),npad)
              end if
            case default
              memsum = memsum + nalloc + 2*npad
            end select
            namscr(1:mem_maxname) = ' '
            namscr = curslc%name
            if (curslc%type.eq.1.or.curslc%type.eq.2) then
              if (patchk1.and.patchk2) then
                write(luout,'(6x,a,x,i2,x,i10,x,i10,x,l,2x,l)')
     &               namscr(1:mem_maxname),curslc%type,
     &               curslc%len,memsum,patchk1,patchk2
              else
                write(luout,'(3x,"!",2x,a,x,i2,x,i10,'//
     &               'x,i10,x,l,2x,l,x,"!")')
     &               namscr(1:mem_maxname),curslc%type,
     &               curslc%len,memsum,patchk1,patchk2
                ierr = ierr+1
             end if
            else
              write(luout,'(6x,a,x,i2,x,i10,x,i10,x,"N/A")')
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
     &     call quit(1,'memman_map',
     &     'range errors detected (see above)')

      end subroutine

      subroutine memman_stat(luout)

      implicit none

      integer, intent(in) ::
     &     luout

      write(luout,'(x,"+",76("-"),"+")')
      write(luout,'(3x,a,i10,a,f8.2,a)')
     &     'Maximum allocated memory: ',max_mem,
     &     ' real(8)-words (',dble(max_mem)/1024d0/128d0,' Mb)'
      write(luout,'(3x,a,i10,a,f8.2,a)')
     &     'Largest memory block:     ',max_blk,
     &     ' real(8)-words (',dble(max_mem)/1024d0/128d0,' Mb)'
      write(luout,'(3x,a,a)')
     &     'Name of largest block:    ',trim(name_max)
      write(luout,'(x,"+",76("-"),"+")')

      end subroutine

      end module
