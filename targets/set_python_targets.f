*----------------------------------------------------------------------*
      subroutine set_python_targets(tgt_info,python_file,name_infile,
     &     name_orbinfo)
*----------------------------------------------------------------------*
*
* Set targets to use python script as target files.
* python_file is the python script.
* name_infile is the input file name and name_orbinfo
* is the name of the file with orbital informations, generated
* by interfaces/put_orbinfo
*
* pradipta, yuri, nov, 2014
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      character(*), intent(in) ::
     &     python_file, name_infile, name_orbinfo

      character(len=256), parameter ::
     &     tgt_sufix = ".tgt_list"
      logical ::
     &     file_exists
      integer ::
     &     pos

      pos = index(python_file,'/',back=.true.)
      pos = pos+1

      inquire(file=trim(python_file(pos:))//trim(tgt_sufix),
     &     exist=file_exists)
      if (file_exists)
     &     call system("rm "//trim(python_file(pos:))//trim(tgt_sufix))

      call system("python "//trim(python_file)//" "//
     &     trim(name_infile)//" "//trim(name_orbinfo))
c dbg
c       write(lulog,*) 'called python to process: ',trim(python_file)
c dbg
     
      inquire(file=trim(python_file(pos:))//trim(tgt_sufix),
     &     exist=file_exists)
      if (.not.file_exists)
     &     call quit(1,'set_python_targets',
     &     "Target file not generated!")
  
      call get_targets_from_file(tgt_info,
     &     trim(python_file(pos:))//trim(tgt_sufix))

      return

      contains

      ! Check if the python module is accessible:
      ! It is in the work directory or in the $PYTHONPATH
      subroutine verify_python_module()
      
      implicit none

      character(len=256) ::
     &     py_path
      character(len=256), parameter ::
     &     py_module = "/gecco_interface.py"


      integer ::
     &     pos1, pos2


      inquire(file=trim(py_module),
     &     exist=file_exists)
      
      if (.not.file_exists) then

       py_path = ""
       call get_environment_variable( "PYTHONPATH", value=py_path)
       
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
     &     call quit(1,'verify_python_module',
     &     "It seems that the module for Python interface is "//
     &     "not acessible. Add the proper directory to the "//
     &     "PYTHONPATH environment variable.")

      end subroutine verify_python_module



      end
