*----------------------------------------------------------------------*
      program GeCCo
*----------------------------------------------------------------------*
*     started by andreas in march 2007 
*     from a pilot version within LUCIA
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'opdim.h'
      include 'ifc_memman.h'
      include 'ifc_input.h'

      logical ::
     &     l_infile,l_exit,one_more
      character ::
     &     name_infile*256, env_type*32
      integer ::
     &     idum, memmax, ifree
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0
      type(filinf) ::
     &     ffinput
      type(orbinf) ::
     &     orb_info

      ! a few settings
      luout = 6      ! output unit
c      iprlvl = 1     ! print level
      iprlvl = 10    ! print level
      ivale = 3      ! V is 3
      iextr = 4      ! X is 4

      ! give information about compilation date etc.
      call printversion()

      ! initialize timing routine
      call init_time()
      call atim_csw(cpu0,sys0,wall0)
      
      call printheader()

      ! process arguments to GeCCo
      call arg_inp(l_exit,l_infile,name_infile)
      if (l_exit) goto 2308

      ! init the file-handler
      call fh_init(iprlvl)

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
      ! one more is ignored, as we might have cases 
      ! (export/import stuff) where no "calculate" block is
      ! specified

      call get_argument_value('general','memmax',ival=memmax)
      call mem_init(memmax)

      ifree = mem_register(4000,'input')
      
c      call test_memman()

      ! loop over calculations
      do
      
        if (.not.one_more) exit

        call do_calc(orb_info,env_type)

        ! post-process input up to next "calculate" block
        ! (data resides in module parse_input)
        call process_input(one_more,orb_info)

      end do

      call mem_clean

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'total time in GeCCo run',
     &     cpu-cpu0,sys-sys0,wall-wall0)

 2308 stop '+++ GeCCo run finished +++'
      end


