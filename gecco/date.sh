#!/bin/sh
# eine kleine TURBOMOLE-Leihgabe

# probably very inelegant way to keep strings limited:
echo $2 | split -b60
optstr=`cat xaa`
rm xa*

cmp=$1
opt1=$optstr
opt2=""
opt3=""
host=`hostname`
set `date`
cat > date.h <<EOF
      character, parameter :: vers*80 = 
     ,  "Version compiled $3 $2 $6 at $4 on "//
     ,  "$host "
      character, parameter :: cmp*40 = 
     ,  "$cmp "
      character, parameter :: opt1*60 =
     ,  "$opt1 "
      character, parameter :: opt2*60 =
     ,  "$opt2 "
      character, parameter :: opt3*60 =
     ,  "$opt3 "
EOF

