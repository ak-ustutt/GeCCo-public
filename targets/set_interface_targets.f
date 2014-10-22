*----------------------------------------------------------------------*
      subroutine set_interface_targets(tgt_info,name_infile,
     &     name_orbinfo)
*----------------------------------------------------------------------*
*
* Set targets form script interfaces: python iterface
* name_infile is the input file name and name_orbinfo
* is the name of the file with orbital informations, generated
* by interfaces/put_orbinfo
*
* yuri, oct, 2014
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'
      include 'ifc_input.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      character(*), intent(in) ::
     &     name_infile, name_orbinfo

      integer ::
     &     nkey, ikey, narg, iarg, len, pos1, pos2
      character(len=256) ::
     &     file_name
      character(len=256), parameter ::
     &     tgt_sufix = ".tgt_list"
      logical ::
     &     file_exists
      
      call verify_python_module()

      nkey = is_keyword_set('calculate.interfaces')
      do ikey = 1, nkey
       narg = is_argument_set('calculate.interfaces','file',
     &      keycount=ikey)
       do iarg = 1, narg
        file_name = ""
        call get_argument_value('calculate.interfaces','file',
     &       keycount=ikey,argcount=iarg,str=file_name)

        inquire(file=trim(file_name)//trim(tgt_sufix),
     &       exist=file_exists)
        
        if (file_exists)
     &       call system("rm "//trim(file_name)//trim(tgt_sufix))

        call system("python "//trim(file_name)//" "//
     &       trim(name_infile)//" "//trim(name_orbinfo))
        
        inquire(file=trim(file_name)//trim(tgt_sufix),
     &       exist=file_exists)
        
        if (.not.file_exists)
     &       call quit(1,'set_interface_targets',
     &       "Target file not generated!")
        
        call get_targets_from_file(tgt_info,
     &       trim(file_name)//trim(tgt_sufix))
        
       end do
      end do

      return

      contains

      ! Check if the python module is accessible:
      ! is in the work directory or in the $PYTHONPATH
      subroutine verify_python_module()
      
      implicit none

      character(len=256) ::
     &     py_path
      character(len=256), parameter ::
     &     py_module = "/gecco_interface.py"

      inquire(file=trim(py_module),
     &     exist=file_exists)
      
      if (.not.file_exists) then

       py_path = ""
       call get_environment_variable( "PYTHONPATH", value=py_path,
     &      length = len)
       
       pos2 = 0
       do 
        pos1 = pos2
        pos2 = index( py_path(pos2+1:), ":")
        if (pos2.eq.0) then
         pos2 = len_trim( py_path) + 1
        else
         pos2 = pos2 + pos1
        end if
        
        inquire(file=py_path(pos1+1:pos2-1)//trim(py_module),
     &       exist=file_exists)
        if (file_exists) exit
        
        if (pos2.eq.len_trim( py_path)+1) exit
        
       end do

      end if

      if (.not.file_exists)
     &     call quit(1,'set_interface_targets',
     &     "It seems that the module for Python interface is "//
     &     "not acessible. Add the proper directory to the "//
     &     "PYTHONPATH environment variable.")

      end
      end

