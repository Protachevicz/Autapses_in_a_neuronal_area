> POK_CV_F0.dat

for i in {1..11}; do
   let "j=$i-1" 
   cat  POK_CV_F$j.dat POK_CV_F_$i.dat > POK_CV_F$i.dat
   echo 
   rm POK_CV_F$j.dat
   rm POK_CV_F_$i.dat
done
