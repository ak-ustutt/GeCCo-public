*----------------------------------------------------------------------*
      subroutine set_interface_targets(tgt_info,name_infile,
     &     name_orbinfo)
*----------------------------------------------------------------------*
*
* Set targets from script interfaces: python iterface
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
     &     nkey, ikey, narg, iarg
      character(len=256) ::
     &     file_name
      
      nkey = is_keyword_set('calculate.interfaces')
      do ikey = 1, nkey
       narg = is_argument_set('calculate.interfaces','file',
     &      keycount=ikey)
       do iarg = 1, narg
        file_name = ""
        call get_argument_value('calculate.interfaces','file',
     &       keycount=ikey,argcount=iarg,str=file_name)

        call set_python_targets(tgt_info,file_name,name_infile,
     &       name_orbinfo)
        
       end do
      end do

      return

      end

