*----------------------------------------------------------------------*
      subroutine set_interface_targets(tgt_info,name_infile,
     &     name_orbinfo)
*----------------------------------------------------------------------*
*
* Set targets form script interfaces
*
      implicit none

      include 'mdef_target_info.h'
      include 'ifc_input.h'


      type(target_info), intent(inout) ::
     &     tgt_info
      character(*), intent(in) ::
     &     name_infile, name_orbinfo

      integer ::
     &     nkey, ikey, narg, iarg
      character(len=256) ::
     &     file_name
      character(len=256), parameter ::
     &     tgt_sufix = ".tgt_list"
      logical ::
     &     file_exists
      

      nkey = is_keyword_set('calculate.interfaces')
      do ikey = 1, nkey
       narg = is_argument_set('calculate.interfaces','file',
     &      keycount=ikey)
       do iarg = 1, narg
        file_name = ""
        call get_argument_value('calculate.interfaces','file',
     &       keycount=ikey,argcount=iarg,str=file_name)

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
      end

