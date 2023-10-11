% for the not banded matrix

D = [ 0.00000
     0.00000
     0.00000
     0.00000
     0.00000
     0.00000
    -0.29700
    -0.09277
    -0.27683
    -0.00000
    -0.29700
     0.09277
    -0.58792
    -0.09208
    -0.58737
    -0.00000
    -0.58792
     0.09208
    -0.88746
    -0.08979
    -0.88839
    -0.00000
    -0.88746
     0.08979];

Kmat = [  1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   0.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00  16.81   0.00   1.87   0.00   0.00   0.00  -5.14  -0.23  -4.20  -3.04   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00   0.00  16.81   0.00 -10.27   0.00   0.00   0.23   0.93  -3.04  -4.20   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00   1.87   0.00  33.63   0.00   1.87   0.00  -4.20   3.04 -10.27   0.00  -4.20  -3.04   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00   0.00 -10.27   0.00  33.63   0.00 -10.27   3.04  -4.20   0.00   1.87  -3.04  -4.20   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   1.87   0.00  16.81   0.00   0.00   0.00  -4.20   3.04  -5.14   0.23   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00 -10.27   0.00  16.81   0.00   0.00   3.04  -4.20  -0.23   0.93   0.00   0.00   0.00   0.00   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00  -5.14   0.23  -4.20   3.04   0.00   0.00  16.81   0.00   1.87   0.00   0.00   0.00  -5.14  -0.23  -4.20  -3.04   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00  -0.23   0.93   3.04  -4.20   0.00   0.00   0.00  16.81   0.00 -10.27   0.00   0.00   0.23   0.93  -3.04  -4.20   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00  -4.20  -3.04 -10.27   0.00  -4.20   3.04   1.87   0.00  33.63   0.00   1.87   0.00  -4.20   3.04 -10.27   0.00  -4.20  -3.04
          0.00   0.00   0.00   0.00   0.00   0.00  -3.04  -4.20   0.00   1.87   3.04  -4.20   0.00 -10.27   0.00  33.63   0.00 -10.27   3.04  -4.20   0.00   1.87  -3.04  -4.20
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -4.20  -3.04  -5.14  -0.23   0.00   0.00   1.87   0.00  16.81   0.00   0.00   0.00  -4.20   3.04  -5.14   0.23
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -3.04  -4.20   0.23   0.93   0.00   0.00   0.00 -10.27   0.00  16.81   0.00   0.00   3.04  -4.20  -0.23   0.93
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -5.14   0.23  -4.20   3.04   0.00   0.00   8.41  -3.04   0.93  -0.23   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -0.23   0.93   3.04  -4.20   0.00   0.00  -3.04   8.41   0.23  -5.14   0.00   0.00
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -4.20  -3.04 -10.27   0.00  -4.20   3.04   0.93   0.23  16.81   0.00   0.93  -0.23
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -3.04  -4.20   0.00   1.87   3.04  -4.20  -0.23  -5.14   0.00  16.81   0.23  -5.14
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -4.20  -3.04  -5.14  -0.23   0.00   0.00   0.93   0.23   8.41   3.04
          0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -3.04  -4.20   0.23   0.93   0.00   0.00  -0.23  -5.14   3.04   8.41];

Dfirst = [  0.00000
          0.00000
          0.00000
          0.00000
          0.00000
          0.00000
         -4.69658
         -1.46701
         -9.52926
          0.95320
         -5.18271
          1.10889
         -5.29638
         -2.77334
        -10.01048
          0.70110
         -5.16308
          2.41829
          1.54712
         -0.84067
          5.68102
          0.21061
          3.03006
          0.95401];

Dfirst_matlab = [ 0
                0         
                0         
                0         
                0         
                0   
          -0.0002    
           0.0053   
          -0.0096         
                0   
          -0.0002  
          -0.0053    
           0.0000    
           0.0057   
          -0.0177
                0    
           0.0000 
          -0.0057
          -2.5491
           0.0024
          -5.0952
                0   
          -2.5491   
          -0.0024];

C_matlab = 9.0616;

D_actual = D'*Kmat;
C_actual = D'*Kmat*D;
