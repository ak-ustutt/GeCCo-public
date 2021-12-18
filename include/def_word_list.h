      integer, parameter ::
     &     maxlen_word = 512

      type word_list_entry
        character(len=maxlen_word) ::
     &       word
        character(len=1) ::
     &       sep
        integer ::
     &       line, col
        type(word_list_entry), pointer ::
     &       next, down, up
      end type word_list_entry

      type word_list
        type(word_list_entry), pointer ::
     &     head, tail, current
      end type word_list
