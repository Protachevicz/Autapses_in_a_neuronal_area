#!/bin/bash

for i in {1..11}; do ### 1..9
  icc -fast -o z$i.z IandF_exc_ini_atraso_bash.c -lm
 ##  gcc -O3 -o z$i.z IandF_exc_ini_atraso_bash.c -lm
##	icc -fast -o z$i.z IandF_exc_ini_atraso_bash.c -lm
  echo 
done

for i in {1..11}; do 
   echo $i | ./z$i.z & 
   echo 
   sleep 1
   rm z$i.z 
done
