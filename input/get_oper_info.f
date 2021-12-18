*----------------------------------------------------------------------*
      subroutine get_oper_info(env_type,prop_type,labelMD,psym,trplt)
*----------------------------------------------------------------------*
*     Get informations about the one-electron operator prop_type
*
*     labelMD is a concatenation of the possible labels for this
*          labelMD =  <labelM>:<labelD>
*          where <labelM> is the nomenclature of Molpro
*                <labelD> is the nomenclature of Molpro
*
*     psym is the parity
*     trplt is true for triplet (spin dependent) operators
*
*----------------------------------------------------------------------*
      implicit none

      include 'def_filinf.h'
      include 'stdunit.h'

      character(len=*), intent(in) ::
     &     env_type
      character(len=*), intent(in) ::
     &     prop_type
      character(len=*), intent(out) ::
     &     labelMD
      integer, intent(out) ::
     &     psym
      logical, intent(out) ::
     &     trplt

      integer ::
     &     len
      character(len=256) ::
     &     gecco_path
      type(filinf) ::
     &     ffopinfo
      character(len=100) ::
     &     line
      character(len=15) ::
     &     labelM, labelD
      logical ::
     &     found

      call get_environment_variable( "GECCO_DIR", value=gecco_path,
     &     length = len)

      call file_init(ffopinfo,
     &     trim(gecco_path)//'/input/one_el_op_convention.txt'
     &     ,ftyp_sq_frm,0)

      call file_open(ffopinfo)

      found = .false.
      do while (.not.found)
        read( ffopinfo%unit, '(a100)', end=1, err=2) line
        if (line(1:1) .NE. '#') then
          backspace( ffopinfo%unit)
          read( ffopinfo%unit, *, end=1, err=2)
     &         labelM, labelD, psym, trplt
          if (trim(prop_type) .EQ. trim(labelM) .OR.
     &         trim(prop_type) .EQ. trim(labelD)) then
            found = .true.
            labelMD = trim(labelM)//':'//trim(labelD)
          end if
        end if
      end do

      call file_close_keep(ffopinfo)

      select case(trim(env_type))
      case ('d','D')
        
        if (labelD.EQ.'-') call quit(1,'get_oper_info',
     &       'One electron integrals not available for '
     &       //trim(env_type))
        
      case ('m','M')
        
        if (labelM.EQ.'-') call quit(1,'get_oper_info',
     &       'One electron integrals not available for '
     &       //trim(env_type))
        
      case default
        call quit(1,'get_oper_info',
     &       'One electron integrals not available for '
     &       //trim(env_type))
        
      end select
      
      return

 1    write(lulog, 7) trim(prop_type)
      call quit(0,'get_oper_info',
     &     'One electron operator label not found: end of file.')
c
 2    write(lulog, 7) trim(prop_type)
      call quit(0,'get_oper_info',
     &     'One electron operator label not found: error.') 
c
 7    format(/' GeCCo does not know this one electron operator:', A )

      end subroutine
