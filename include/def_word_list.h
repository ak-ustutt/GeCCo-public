      integer, parameter ::
     &     maxlen_word = 128

      type word_list_entry
        character(len=maxlen_word) ::
     &       word
        character(len=1) ::
     &       sep
        type(word_list_entry), pointer ::
     &       next
      end type word_list_entry

      type word_list
        type(word_list_entry), pointer ::
     &     head, tail, current
      end type word_list
