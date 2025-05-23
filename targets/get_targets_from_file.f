*----------------------------------------------------------------------*
      subroutine get_targets_from_file(tgt_info, file_name)
*----------------------------------------------------------------------*
* Get targets from the file file_name 
*
* The format of such file is somehow boring and is described in the
* following. It is supposed to be generated by an interface.
*
*
*Target list generated by: <source info>
*Name <target name>
*Required <F or T>
*Dependencies <number of dependencies>
*<dependency 1>
*...
*<dependency n>
*Joined <number of joined targets>
*<joined 1>
*...
*<joined n>
*Rules <number of rules>
*<rule 1> <number of arguments>
*<argument 1> <type of argument 1> <number of values>
*<value 1>
*...
*<value n>
*<argument 2> <type of argument 2> <number of values>
*...
*<rule 1> <number of arguments>
*...
*<last value of last argument of last rule>
*---
*Name <target_name>
*...
*---
*
*
* If some error occured in the target generation, this can be passed
* to gecco if the first line is
*
* Error when creating list: <error description>
* 
* case when GeCCo stops and report the error.
*
*
* yuri, oct, 2014
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'def_filinf.h'
      include 'ifc_input.h'

      include 'ifc_targets.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      character(*) ::
     &     file_name

      type(filinf) ::
     &     fftgt
      integer ::
     &     idum, lutgt, ios_var, i1, i2, i3, n1, n2, n3

      character(50) ::
     &     thing, file_origin

      character(len_target_name) ::
     &     tgt_name, tgt_aux
      logical ::
     &     tgt_req
      character(len_command_name) ::
     &     rule_name
      character(len_target_name) ::
     &     arg_name
      integer, parameter ::
     &     max_len_val = 256, max_len_type = 7, len_val_number = 7
      character(max_len_val) ::
     &     args
      character(max_len_type) ::
     &     arg_type
      character(len_val_number) ::
     &     val_number
      integer ::
     &     print_level
      character(30), allocatable, dimension(:) ::
     &     label_list
      logical, allocatable, dimension(:) ::
     &     log_list
      integer, allocatable, dimension(:) ::
     &     int_list
      real(8), allocatable, dimension(:) ::
     &     real_list
      integer, parameter ::
     &     len_str = 2048
      character(len_str) ::
     &     str

      call get_argument_value('general','print',
     &     ival=print_level)

      call file_init(fftgt,trim(file_name),ftyp_sq_frm,idum)
      call file_open(fftgt)
      lutgt = fftgt%unit

      read(lutgt,FMT='(A25,A50)') thing, file_origin

      if (trim(thing).EQ."Error when creating list:") then
       read(lutgt,FMT='(A256)',iostat=ios_var) args
       do while(ios_var.eq.0)
        write(lulog,fmt='(A)') trim(args)
        if (lulog.ne.luout) write(luout,fmt='(A)') trim(args)
        read(lutgt,FMT='(A256)',iostat=ios_var) args
       end do
       call quit(0,'get_targets_from_file','Error signal from '//
     &      trim(file_name)//' ::'//trim(file_origin))
       end if

      if (trim(thing).NE."Target list generated by:")
     &     call quit(1,'get_targets_from_file',
     &     'Strange first line: '//
     &     'Target files have a well defined format.')

      if (print_level.ge.3)
     &     write(lulog,fmt='(" Reading target file ",A,":",A)')
     &     trim(file_name),trim(file_origin)
      
      read(lutgt,*,iostat=ios_var) thing, tgt_name
      do while(ios_var.eq.0)

       if (trim(thing).ne.'Name') call quit(1,'get_targets_from_file',
     &      'Bad formatted file. Expected "Name"')

       read(lutgt,*,iostat=ios_var) thing, tgt_req
       if (ios_var.ne.0) call quit(1,'get_targets_from_file',
     &      'Error in file reading: looking for required information.')
       if (trim(thing).ne.'Required')
     &      call quit(1,'get_targets_from_file',
     &      'Bad formatted file. Expected "Required"')

       call add_target2(trim(tgt_name),tgt_req,tgt_info)

       read(lutgt,*,iostat=ios_var) thing, n1
       if (ios_var.ne.0) call quit(1,'get_targets_from_file',
     &      'Error in file reading: '//
     &      'looking for number of dependencies.')
       if (trim(thing).ne.'Dependencies')
     &      call quit(1,'get_targets_from_file',
     &      'Bad formatted file. Expected "Dependencies"')

       do i1 = 1,n1
        read(lutgt,*,iostat=ios_var) tgt_aux
        if (ios_var.ne.0) call quit(1,'get_targets_from_file',
     &       'Error in file reading: looking for dependency.')
        call set_dependency(trim(tgt_name),trim(tgt_aux),tgt_info)
       end do

       read(lutgt,*,iostat=ios_var) thing, n1
       if (ios_var.ne.0) call quit(1,'get_targets_from_file',
     &      'Error in file reading: '//
     &      'looking for number of joined targets.')
       if (trim(thing).ne.'Joined')
     &      call quit(1,'get_targets_from_file',
     &      'Bad formatted file. Expected "Joined"')

       do i1 = 1,n1
        read(lutgt,*,iostat=ios_var) tgt_aux
        if (ios_var.ne.0) call quit(1,'get_targets_from_file',
     &       'Error in file reading: looking for joined.')
        call set_dependency(trim(tgt_name),trim(tgt_aux),tgt_info)
       end do

       read(lutgt,*,iostat=ios_var) thing, n1
       if (ios_var.ne.0) call quit(1,'get_targets_from_file',
     &      'Error in file reading: looking for number of rules.')
       if (trim(thing).ne.'Rules')
     &      call quit(1,'get_targets_from_file',
     &      'Bad formatted file. Expected "Rules"')

       do i1 = 1,n1
        
        read(lutgt,*,iostat=ios_var) rule_name, n2
        if (ios_var.ne.0) call quit(1,'get_targets_from_file',
     &       'Error in file reading: looking for number of arguments.')

        call set_rule2(trim(tgt_name),trim(rule_name),tgt_info)

        do i2 = 1,n2

         read(lutgt,*,iostat=ios_var) arg_name, arg_type, n3
         if (ios_var.ne.0) call quit(1,'get_targets_from_file',
     &        'Error in file reading: looking for argument type and '//
     &        'number of values.')

         select case( trim(arg_type))
         case('label')
          allocate(label_list(n3))
         case('log')
          allocate(log_list(n3))
         case('int')
          allocate(int_list(n3))
         case('real')
          allocate(real_list(n3))
         end select

         do i3 = 1,n3

          select case( trim(arg_type))
          case('label')
           read(lutgt,'(A)',iostat=ios_var) label_list(i3)
          case('log')
           read(lutgt,*,iostat=ios_var) log_list(i3)
          case('int')
           read(lutgt,*,iostat=ios_var) int_list(i3)
          case('occ')
           call quit(1,'get_targets_from_file',
     &          'Weird! "occ" not yet implemented for interfaces!')
          case('restr')
           call quit(1,'get_targets_from_file',
     &          'Weird! "restr" not yet implemented for interfaces!')
          case('real')
           read(lutgt,*,iostat=ios_var) real_list(i3)
          case('str')
           call get_whole_line()
          case default
           call quit(1,'get_targets_from_file',
     &          'Unknown type of argument: ',trim(arg_type))
          end select

         end do

         select case( trim(arg_type))
         case('label')
          call set_arg(trim(tgt_name), trim(rule_name), trim(arg_name),
     &         n3, tgt_info, val_label = label_list)
          deallocate(label_list)
         case('log')
          call set_arg(trim(tgt_name), trim(rule_name), trim(arg_name),
     &         n3, tgt_info, val_log = log_list)
          deallocate(log_list)
         case('int')
          call set_arg(trim(tgt_name), trim(rule_name), trim(arg_name),
     &         n3, tgt_info, val_int = int_list)
          deallocate(int_list)
         case('occ')
          call quit(1,'get_targets_from_file',
     &         'Weird! "occ" not yet implemented for interfaces!')
         case('restr')
          call quit(1,'get_targets_from_file',
     &         'Weird! "restr" not yet implemented for interfaces!')
         case('real')
          call set_arg(trim(tgt_name), trim(rule_name), trim(arg_name),
     &         n3, tgt_info, val_rl8 = real_list)
          deallocate(real_list)
         case('str')
          call set_arg(trim(tgt_name), trim(rule_name), trim(arg_name),
     &         n3, tgt_info, val_str = str)
         end select
        end do
        
       end do

       read(lutgt,*,iostat=ios_var) thing
       if (ios_var.ne.0) call quit(1,'get_targets_from_file',
     &      'Error in file reading: looking for "---"')
       if (trim(thing).ne.'---') call quit(1,'get_targets_from_file',
     &      'Bad formatted file. Expected "---", get '//trim(thing))
       
       read(lutgt,*,iostat=ios_var) thing, tgt_name
      end do
      
      call file_close_keep(fftgt)

      contains

      subroutine get_whole_line()
      implicit none

      character(1) :: char
      integer :: ilett

      str = ''
      ilett = 1
      do
       read(lutgt,fmt='(a1)',advance='no', iostat = ios_var) char
       if (ios_var.ne.0) exit
       str(ilett:ilett) = char
       ilett = ilett+1
       if(ilett.gt.len_str) then
        write(lulog,fmt='("get_targets_from_file: Warning: "'//
     &       '"string is larger than len_str.")')
        exit
       end if
      end do

      end subroutine get_whole_line

      end
