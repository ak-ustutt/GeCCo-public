*----------------------------------------------------------------------*
      program GeCCo
*----------------------------------------------------------------------*
*     started by andreas in march 2007 
*     from a pilot version within LUCIA
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'warnings.h'
      include 'ioparam.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'opdim.h'
      include 'event_counter.h'
      include 'ifc_memman.h'
      include 'ifc_input.h'

      logical ::
     &     l_infile,l_logfile,l_exit,l_molpro,one_more, do_stat
      character(256) ::
     &     name_infile, name_logfile, host
      character(32) ::
     &     env_type
      character(24) ::
     &     date
      integer ::
     &     idum, memmax, ifree, len
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0
      type(filinf) ::
     &     ffinput, fflog, ffwarn, ffstat
      type(orbinf) ::
     &     orb_info

      ! a few settings
      luout = 6      ! output unit
      lulog = 6      ! logfile (stdout by default)

c      iprlvl = 3     ! print level
      iprlvl = 10    ! print level

      call hostname(host)
      call datum(date)

      ! process arguments to GeCCo
      call arg_inp(l_exit,l_infile,l_logfile,l_molpro,
     &             name_infile,name_logfile)
      if (l_exit) goto 2308

      if (l_molpro) iprlvl = 0

      ! init the file-handler
      call fh_init(iprlvl)

      ! open the logfile
      if (l_logfile) then
        call file_init(fflog,trim(name_logfile),ftyp_sq_frm,idum)
        call file_open(fflog)
        lulog = fflog%unit
      end if

      write(lulog,'(x,"run starts at ",a,"   host: ",a)')
     &     trim(date),trim(host)
      if (lulog.ne.luout.and..not.l_molpro) ! a bit less verbose inside molpro
     &  write(luout,'(x,"run starts at ",a,"   host: ",a)')
     &     trim(date),trim(host)
      
      ! give information about compilation date etc.
      if (.not.l_molpro) call printversion(lulog)

      ! set internal counter to 0
      event_time = 0
      ! initialize timing routine
      call init_time()
      call atim_csw(cpu0,sys0,wall0)

      call printheader(lulog)
      if (luout.ne.lulog.and..not.l_molpro) call printheader(luout)

      ! warnings
      nwarn = 0
      call file_init(ffwarn,'WARNINGS',ftyp_sq_frm,idum)
      call file_open(ffwarn)
      luwarn = ffwarn%unit

      ! find out, which environment we are using, and where
      ! to get our input data from
      call where_am_I(l_infile,name_infile,env_type)

      ! assign handle to input file
      call file_init(ffinput,name_infile,ftyp_sq_frm,idum)
      ! get all possible information from environment:
      !  number of orbitals, symmetry etc.
      call read_env(env_type,orb_info)

      ! read and parse input file
      call read_input(ffinput)

      ! post-process input up to first "calculate" block
      ! (data resides in module parse_input)
      call process_input(one_more,orb_info)
      ! one_more is ignored, as we might have cases 
      ! (export/import stuff) where no "calculate" block is
      ! specified

      call get_argument_value('general','memmax',ival=memmax)
      call mem_init(memmax)

      ifree = mem_register(4000,'input')

      ! open statistics file, if requested
      call get_argument_value('general','statistics',lval=do_stat)
      if (do_stat) then
        call file_init(ffstat,'STATISTICS',ftyp_sq_frm,idum)
        call file_open(ffstat)
        lustat = ffstat%unit
      else
        lustat = -1
      end if


      ! loop over calculations
      do
      
        if (.not.one_more) exit

        call do_calc(orb_info,env_type,name_infile)

        ! post-process input up to next "calculate" block
        ! (data resides in module parse_input)
        call process_input(one_more,orb_info)

      end do

      call mem_clean

      if (nwarn.gt.0) then
        call file_close_keep(ffwarn)
        write(lulog,'(1x,a,i4,a)')
     &     'There were ',nwarn,' warnings, see file '//trim(ffwarn%name)
        if (lulog.ne.luout) write(luout,'(1x,a,i4,a)')
     &     'There were ',nwarn,' warnings, see file '//trim(ffwarn%name)
      else
        call file_close_delete(ffwarn)
      end if

      if (lustat.gt.0) call file_close_keep(ffstat)

      call atim_csw(cpu,sys,wall)
      call prtim(lulog,'total time in GeCCo run',
     &     cpu-cpu0,sys-sys0,wall-wall0)
      if (lulog.ne.luout.and..not.l_molpro)
     &   call prtim(luout,'total time in GeCCo run',
     &     cpu-cpu0,sys-sys0,wall-wall0)

 2308 call datum(date)
      write(lulog,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)
      if (lulog.ne.luout.and..not.l_molpro)
     &   write(luout,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)

      if (l_logfile) call file_close_keep(fflog)

      if (.not.l_molpro) stop '+++ GeCCo run finished +++'
      end


