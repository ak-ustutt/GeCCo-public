*---------------------------------------------------------------------------*
      subroutine set_keywords()
*---------------------------------------------------------------------------*
*
*     set the keywords to interpret the gecco input file
*     using the info in $GECCO_DIR/keyword_registry.xml
*       (sloppy read of XML data)
*
*     reactivated from old pre-2016 version
*
*     andreas, sept 2021
*
*---------------------------------------------------------------------------*

c      use parse_input
      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'
      include 'def_filinf.h'
      include 'def_word_list.h'

      integer, parameter ::
     &     ntest = 00
      integer, parameter ::
     &     maxlength = 256

      logical ::
     &     log_array(256)
      integer ::
     &     int_array(256)
      real(8) ::
     &     rl8_array(256)
      character(len=1) ::
     &     str_array(256)

      integer ::
     &     iprint, ival, type, length, length2, idx, i
      logical ::
     &     name_found, type_found, value_found, length_found
      
      character(len=256) ::
     &     keyword_file_name, context, word, name, word_cut
      character(len=1) ::
     &     sep
      type(filinf) ::
     &     keyword_file
      type(word_list) ::
     &     keyword_defs

      logical, external ::
     &    file_exists, advance_word_list_entry

      iprint = max(ntest,iprlvl)

      ! get keyword file
      keyword_file_name = get_keyword_file()

      write(lulog,*) 'reading keywords from '//trim(keyword_file_name)

      ! open keyword file
      call file_init(keyword_file,keyword_file_name,ftyp_sq_frm,0)

      if (.not.file_exists(keyword_file)) then
        call quit(0,'set_keywords','file not found: "'//
     &              trim(keyword_file_name)//'"')     
      end if

      ! read file contents to word list
      call init_word_list(keyword_defs)
      call lex_file(keyword_defs,keyword_file,
     &                    '','')
c     &                    'gecco keywords','end gecco keywords')
      if (ntest.ge.100)
     &     call print_word_list(lulog,keyword_defs)

      call keyword_init()

      ! we start with no context
      context = ""

      ! make sure that "current" pointer is at the head of the list
      keyword_defs%current => keyword_defs%head

      ! loop over contents
      main_loop: do

        if (ntest.ge.100)
     &      write(lulog,*) 'current context: '//trim(context)

        call get_word_list_entry(word,sep,keyword_defs)

        if (ntest.ge.100)
     &      write(lulog,*) 'top split level: '//trim(word)

        !-------------------------------------------------------------
        ! "keyword"
        !-------------------------------------------------------------
        if (index(trim(word),"<keyword").gt.0) then
           ! get next entry
           if (.not.advance_word_list_entry(keyword_defs,' ')) then
             call quit(1,'set_keywords','expected "name" but file ends')
           end if
           call get_word_list_entry(word,sep,keyword_defs)
           if (trim(word)=="name") then
              if (.not.advance_word_list_entry(keyword_defs,' ')) then
                 call quit(1,'set_keywords',
     &                       'expected keyword name but file ends')
              end if
              call get_word_list_entry(word,sep,keyword_defs)
              ! add keyword
              if (trim(context)=="") then
                 call keyword_add(trim(word))
                ! update context
                context = trim(word)
              else
                 call keyword_add(trim(word),context=trim(context))
                ! update context
                context = trim(context)//"."//trim(word)
              end if
           else
              call quit(1,'set_keywords',
     &          'expected "name" but found '//trim(word))
           end if
        !-------------------------------------------------------------
        ! "argument"
        !-------------------------------------------------------------
        else if (index(trim(word),"<argument").gt.0) then
           name_found = .false.
           length_found = .false.
           type_found = .false.
           value_found = .false.
           arg_loop: do
              if (.not.advance_word_list_entry(keyword_defs,' ')) then
                call quit(1,'set_keywords',
     &                      'scanning "argument" but file ends')
              end if
              call get_word_list_entry(word,sep,keyword_defs)
              if (ntest.ge.100)
     &            write(lulog,*) 'at top of argument loop: '//trim(word)
              if (trim(word)=="name") then
                 name_found = .true.
                 if (.not.advance_word_list_entry(keyword_defs,' '))then
                    call quit(1,'set_keywords',
     &                      'scanning "argument name" but file ends')
                 end if
                 call get_word_list_entry(word,sep,keyword_defs)
                 name = trim(word)
                if (ntest.ge.100)
     &              write(lulog,*) 'argument: '//trim(word)
              else if (trim(word)=="length") then
                 length_found = .true.
                 if (.not.advance_word_list_entry(keyword_defs,' '))then
                    call quit(1,'set_keywords',
     &                      'scanning "argument length" but file ends')
                 end if
                 call get_word_list_entry(word,sep,keyword_defs)
                 read(word,*,err=100) length
                 if (length.gt.maxlength)
     &              call quit(1,'set_keyword','length > 256 detected')
              else if (trim(word)=="type") then
                 type_found = .true.
                 if (.not.advance_word_list_entry(keyword_defs,' '))then
                    call quit(1,'set_keywords',
     &                      'scanning "argument type" but file ends')
                 end if
                 call get_word_list_entry(word,sep,keyword_defs)
                 read(word,*,err=100) type
              else if (trim(word)=="value") then
                 value_found = .true.
c                 if (.not.advance_word_list_entry(keyword_defs,' '))then
c                    call quit(1,'set_keywords',
c     &                      'scanning "argument value" but file ends')
c                 end if
                 if (.not.type_found.or..not.length_found) then
                    call quit(1,'set_keywords',
     &                   'found "argument value" but type and length '//
     &                   'are not yet defined: '//trim(context) )
                 end if
                 if (.not.advance_word_list_entry(keyword_defs,' '))then
                    call quit(1,'set_keywords',
     &                      'scanning "argument value entries"'//
     &                      ' but file ends')
                 end if
                 call get_word_list_entry(word,sep,keyword_defs)
                 if (type==vtyp_log.or.type==vtyp_int.or.type==vtyp_rl8)
     &           then
                    log_array = .false.
                    int_array = 0
                    rl8_array = 0d0
                    word_cut = word
                    do ival = 1, length
                       length2 = len_trim(word)
                       ! here: cut and interpret word
                       if (ival.lt.length) then
                          idx = index(word,',')
                          word_cut = word(1:idx)
                          word = word(idx+1:length2)
                       end if
                       select case(type)
                        case (vtyp_log)
                          if (trim(word_cut)=="false") then
                            log_array(ival) = .false.
                           else if (trim(word_cut)=="true") then
                            log_array(ival) = .true.
                           else
                             call quit(1,'set_keyword',
     &                     'expected true or false, found:'
     &                       //trim(word_cut))
                           end if
                        case (vtyp_int)
                           read(word_cut,*,err=101) int_array(ival)
                        case (vtyp_rl8)
                           read(word_cut,*,err=101) rl8_array(ival)
                       end select
                    end do
                 else if (type==vtyp_str) then
                    str_array = " "
                    length2 = len_trim(word)
                    do ival = 1, min(length2,length)
                       str_array(ival) = word(ival:ival)
                    end do
                 else
                    call quit(1,'set_keywords',
     &                      'unknown "argument type"'//
     &                      ' (use 1, 2, 4, 8)')
                 end if
              ! comment begins (and does not directly end)
              else if (index(trim(word),'<!').gt.0.and.
     &                 index(trim(word),'-->').eq.0) then
                 comment_loop: do
                    if (.not.advance_word_list_entry(keyword_defs,' '))
     &              then
                       call quit(1,'set_keywords',
     &                      'scanning comment line'//
     &                      ' but file ends')
                    end if
                    call get_word_list_entry(word,sep,keyword_defs)                   
                    if (index(word,'-->').gt.1) exit comment_loop
                 end do comment_loop
c              else
c                 call quit(1,'set_keywords',
c     &                      'unexpected after "argument": '//
c     &                      trim(word))
              end if
              if (sep == "E") exit arg_loop
           end do arg_loop
           if (.not.name_found.or..not.type_found.or..not.length_found)
     &       call quit(1,'set_keywords','incomplete argument statement')
           ! add argument


           if (value_found) then
             select case(type)
              case (vtyp_log)
                 call argument_add(trim(name),context=trim(context),
     &                             type=vtyp_log,
     &                             len=length,ldef=log_array(1:length))
              case (vtyp_int)
                 call argument_add(trim(name),context=trim(context),
     &                             type=vtyp_int,
     &                             len=length,idef=int_array(1:length))
              case (vtyp_rl8)
                 call argument_add(trim(name),context=trim(context),
     &                             type=vtyp_rl8,
     &                             len=length,xdef=rl8_array(1:length))
              case (vtyp_str)
                 call argument_add(trim(name),context=trim(context),
     &                             type=vtyp_str,
     &                             len=length,cdef=str_array(1:length))
             end select
           else
             call argument_add(trim(name),context=trim(context),
     &                         type=type,
     &                         len=length)
           end if
        !-------------------------------------------------------------
        ! "</keyword>"
        !-------------------------------------------------------------
        else if (index(trim(word),"</keyword").gt.0) then
          ! update context
          idx = index(context,".",back=.true.)
          if (idx.ge.0) then
             context = context(1:idx-1)
          else
             context = ""
          end if
        ! ignore anything else
c        else
c          ! some info...
c          call quit(1,'set_keywords',
c     &   'unexpected content (expected "keyword", "argument" or "end"')
        end if

         if (.not.advance_word_list_entry(keyword_defs,' '))
     &             exit main_loop

      end do main_loop

      if (trim(context).ne."") then
         call quit(1,'set_keywords','final context non-empty: '//
     &              trim(context))
      end if

      call clean_word_list(keyword_defs)

      if (iprint.ge.50)
     &     call show_keywords(lulog)

      return

 100  call quit(1,'set_keywords','troubly reading from string "'
     &             //trim(word)//'"')
 101  call quit(1,'set_keywords','troubly reading from string "'
     &             //trim(word_cut)//'"')

      contains
 
*----------------------------------------------------------------------*
!>    returns name of the keyword_file
!!    reads an environment variable and appends a specified relative path
!!    variable and path are named in par_keyreg.h
*----------------------------------------------------------------------*
      function get_keyword_file() result(file_name)
*----------------------------------------------------------------------*
      implicit none
      include "par_keyreg.h"
      integer,parameter::
     &     ntest= 00
      character(len=16),parameter ::
     &     i_am="get_keyword_file"
      character(len=256)::
     &     path_name,file_name
      integer ::
     &     length,file_loc_len

      file_loc_len=len(keyreg_file_loc)

      call get_environment_variable( keyreg_variable, value=path_name,
     &     length = length)

      if (length.EQ.0)
     &     call quit(0,i_am,
     &     "Please, set the "//keyreg_variable
     &     //" environment variable.")

      if (length .gt.(256-file_loc_len)  )
     &     call quit(0,i_am,
     &     "GECCO_DIR to long, cannot set keyword_registry")

      file_name=trim(path_name)//keyreg_file_loc

      end function

      end
