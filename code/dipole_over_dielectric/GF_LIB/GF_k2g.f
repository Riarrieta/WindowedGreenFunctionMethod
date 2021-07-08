      program GF_k2g        ! Calculation of the half-space
                            ! electromagnetic Green's tensor
                            ! by the short-dsitance expansion.
                            ! This program reads the expansion
                            !   coefficients K** from previously
                            !   created files and then constructs both
                            !   GR and GT
c-------------------------------------------------------------------------
c               Authors:
c                  George Y. Panasyuk, Vadim A. Markel, John C. Schotland
c                  University of Pennsylvania
c               E-mail address for questions and comments:
c                  vmarkel@mail.med.upenn.edu
c               For updates and additional info check
c                  URL: http://whale.seas.upenn.edu/CODES/
c----------------------------------------------------------------------------
      implicit none

      integer*4 idummy, kl, nl, key, nord, ianord, n, nn
      logical*4 ldummy, dicrepancy
      character*4 hdummy, output_prec

      real*8 rdummy
      real*8 z1, z2, z_p, rho, k1, lambda
      real*8 lambda_a_xx, lambda_a_yy, lambda_a_zz, lambda_a_xz
      real*8 lambda_t_xx, lambda_t_yy, lambda_t_zz, lambda_t_xz
      real*8 ru, rz, twopi
      real*8 eps_1
      real*8 L, Lhat, Lhat_nn


      real*8 Kxx_r(0:7), Kyy_r(0:7), Kzz_r(0:7), Kxz_r(0:7)
      real*8 Kxx_i(0:7), Kyy_i(0:7), Kzz_i(0:7), Kxz_i(0:7)

      complex*16 cu, cz
      complex*16 Kxx(0:7), Kyy(0:7), Kzz(0:7), Kxz(0:7)
      complex*16 GRxx, GRyy, GRzz, GRxz
c
c------ Mathematical constants ------------------------------------
c
      ru = 1.0d0
      rz = 0.0d0
      cu = (0.0d0, 1.0d0)
      cz = (0.0d0, 0.0d0)
      twopi = 4.0d0*dasin(ru)
c
c---------- Reading input parameters -----------------------------------
c
      open(unit=70, file='GF.par', status='old', err=101)
      read(70,*) hdummy
      read(70,*) hdummy

      read(70,*, err=101) idummy

      read(70,*,err=101) nord
      if(nord .lt. -7) then 
         write(*,*) 'Error. Incorrect parameter: NORD<-7'
         write(*,*) 'The program GF_k2g allows -7<=NORD<=7'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      else if(nord .lt. 0) then 
         write(*,*) 'You have given me NORD<0'
         write(*,*) 'This means that I will take the coefficients K'
         write(*,*) '   ... from K**_asub.dat for n=0,2,3'
         write(*,*) '   ... and from K**_tsub.dat for n>3,'
         write(*,*) '   ... up to n=|NORD|'
         write(*,*) 'Input file consistency is required'
         write(*,*) 'Continue?'
         pause
      else if (nord .eq. 1) then
         write(*,*) 'The first order correction is identicallyz zero'
         write(*,*) '    (see paper)'
         write(*,*) 'Therefore, there iz no difference between NORD=0'
         write(*,*) '    and NORD=1'
         write(*,*) 'I will set NORD=0 and continue'
         write(*,*) 'The result is not affected. Goodby.'
         nord = 0
      else if(nord .gt. 7) then 
         write(*,*) 'Error. Incorrect parameters: NORD>7'
         write(*,*) 'GF_k2g only allows -7<=NORD<=7'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,2,err=101) output_prec
      if(output_prec .eq. 'S') goto 3
      if(output_prec .eq. 's') goto 3
      if(output_prec .eq. 'D') goto 3
      if(output_prec .eq. 'd') goto 3
         write(*,*) 'Incorrect parameter OUTPUT_PREC'
         write(*,*) 'Only chars S, s, D, d are allowed'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
 3    continue

      read(70,4,err=101) hdummy

      read(70,*,err=101) rdummy, rdummy, idummy
       
      read(70,*,err=101) z1, z2, rho
      if(z1.le.rz) then
         write(*,*) 'Error. Incorrect parameter: z1<=0; z1=', z1
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if
      if(z2.le.rz) then
         write(*,*) 'Error. Incorrect parameter: z2<=0; z2=', z2
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if
      if(rho.le.rz) then
         write(*,*) 'Error. Incorrect parameter: rho<=0; rho=', rho
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*,err=101) rdummy, rdummy
      read(70,*,err=101) rdummy
      read(70,*,err=101) rdummy
      read(70,*,err=101) rdummy

      read(70,*,err=101) eps_1
      if(eps_1.le.rz) then
         write(*,*) 'Error. Incorrect parameter: eps_1<=0; eps_1=', 
     $                                                     eps_1
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*,err=101) idummy
      read(70,*,err=101) idummy
      read(70,*,err=101) idummy      
      read(70,*,err=101) rdummy, rdummy
      
      close(unit=70)
c
c------- End reading parameters -----------------------------
c


c
c------- Input-dependent constants --------------------------
c
      ianord = iabs(nord)
      z_p    = z2 + z1
      L      = dsqrt(rho*rho + z_p*z_p)
c
c-------- End input-dependent constants block ------------------------
c
      if(nord .lt. -3) then

         open(71, file='Kxx_asub', status='unknown', err=102)
         open(72, file='Kyy_asub', status='unknown', err=102)
         open(73, file='Kzz_asub', status='unknown', err=102)
         open(74, file='Kxz_asub', status='unknown', err=102)

         open(75, file='Kxx_tsub', status='unknown', err=102)
         open(76, file='Kyy_tsub', status='unknown', err=102)
         open(77, file='Kzz_tsub', status='unknown', err=102)
         open(78, file='Kxz_tsub', status='unknown', err=102)

      else if(nord .le. 3) then

         open(71, file='Kxx_asub', status='unknown', err=102)
         open(72, file='Kyy_asub', status='unknown', err=102)
         open(73, file='Kzz_asub', status='unknown', err=102)
         open(74, file='Kxz_asub', status='unknown', err=102)

      else if(nord .le. 7) then

         open(75, file='Kxx_tsub', status='unknown', err=102)
         open(76, file='Kyy_tsub', status='unknown', err=102)
         open(77, file='Kzz_tsub', status='unknown', err=102)
         open(78, file='Kxz_tsub', status='unknown', err=102)

      else

         write(*,*) 'Unexpected error 1; exiting'
         write(*,*) 'Please report this error to code developers'
         stop

      end if 

      open(81, file='GRxx_k2g', status='unknown', err=102)
      open(82, file='GRyy_k2g', status='unknown', err=102)
      open(83, file='GRzz_k2g', status='unknown', err=102)
      open(84, file='GRxz_k2g', status='unknown', err=102)
      
      Kxx(2) = cz
      Kyy(2) = cz
      Kzz(2) = cz
      Kxz(2) = cz
c
c------- Start loop over wavelengths ---------------------------------
c

      kl = 0
 1    continue
      kl = kl+1

      if(nord .lt. -3) then

         read(71,*,err=102, end=103) lambda_a_xx, 
     $                             Kxx_r(0), Kxx_i(0), 
     $                             Kxx_r(2), Kxx_i(2),
     $                             Kxx_r(3), Kxx_i(3)
         read(72,*,err=102, end=103) lambda_a_yy, 
     $                             Kyy_r(0), Kyy_i(0), 
     $                             Kyy_r(2), Kyy_i(2),
     $                             Kyy_r(3), Kyy_i(3)
         read(73,*,err=102, end=103) lambda_a_zz, 
     $                             Kzz_r(0), Kzz_i(0), 
     $                             Kzz_r(2), Kzz_i(2),
     $                             Kzz_r(3), Kzz_i(3)
         read(74,*,err=102, end=103) lambda_a_xz, 
     $                             Kxz_r(0), Kxz_i(0), 
     $                             Kxz_r(2), Kxz_i(2),
     $                             Kxz_r(3), Kxz_i(3)

         read(75,*,err=102, end=103) lambda_t_xx, 
     $                             rdummy, rdummy,
     $                             rdummy, rdummy,
     $                             rdummy, rdummy,
     $                             Kxx_r(4), Kxx_i(4), 
     $                             Kxx_r(5), Kxx_i(5),
     $                             Kxx_r(6), Kxx_i(6),
     $                             Kxx_r(7), Kxx_i(7)
         read(76,*,err=102, end=103) lambda_t_yy, 
     $                             rdummy, rdummy,
     $                             rdummy, rdummy,
     $                             rdummy, rdummy,
     $                             Kyy_r(4), Kyy_i(4), 
     $                             Kyy_r(5), Kyy_i(5),
     $                             Kyy_r(6), Kyy_i(6),
     $                             Kyy_r(7), Kyy_i(7)
         read(77,*,err=102, end=103) lambda_t_zz, 
     $                             rdummy, rdummy,
     $                             rdummy, rdummy,
     $                             rdummy, rdummy,
     $                             Kzz_r(4), Kzz_i(4), 
     $                             Kzz_r(5), Kzz_i(5),
     $                             Kzz_r(6), Kzz_i(6),
     $                             Kzz_r(7), Kzz_i(7)
         read(78,*,err=102, end=103) lambda_t_xz, 
     $                             rdummy, rdummy,
     $                             rdummy, rdummy,
     $                             rdummy, rdummy,
     $                             Kxz_r(4), Kxz_i(4), 
     $                             Kxz_r(5), Kxz_i(5),
     $                             Kxz_r(6), Kxz_i(6),
     $                             Kxz_r(7), Kxz_i(7)

         if     (lambda_a_xx .ne. lambda_a_yy) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kyy_a'
         else if(lambda_a_xx .ne. lambda_a_zz) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kzz_a'
         else if(lambda_a_xx .ne. lambda_a_xz) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kxz_a'
         else if(lambda_a_xx .ne. lambda_t_xx) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kxx_t'
         else if(lambda_a_xx .ne. lambda_t_xx) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kxx_t'
         else if(lambda_a_xx .ne. lambda_t_yy) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kyy_t'
         else if(lambda_a_xx .ne. lambda_t_zz) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kzz_t'
         else if(lambda_a_xx .ne. lambda_t_xz) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kxz_t'
         end if

         lambda = lambda_a_xx

      else if(nord .le. 3) then
         
         read(71,*,err=102, end=103) lambda_a_xx, 
     $                             Kxx_r(0), Kxx_i(0), 
     $                             Kxx_r(2), Kxx_i(2),
     $                             Kxx_r(3), Kxx_i(3)
         read(72,*,err=102, end=103) lambda_a_yy, 
     $                             Kyy_r(0), Kyy_i(0), 
     $                             Kyy_r(2), Kyy_i(2),
     $                             Kyy_r(3), Kyy_i(3)
         read(73,*,err=102, end=103) lambda_a_zz, 
     $                             Kzz_r(0), Kzz_i(0), 
     $                             Kzz_r(2), Kzz_i(2),
     $                             Kzz_r(3), Kzz_i(3)
         read(74,*,err=102, end=103) lambda_a_xz, 
     $                             Kxz_r(0), Kxz_i(0), 
     $                             Kxz_r(2), Kxz_i(2),
     $                             Kxz_r(3), Kxz_i(3)

         if     (lambda_a_xx .ne. lambda_a_yy) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kyy_a'
         else if(lambda_a_xx .ne. lambda_a_zz) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kzz_a'
         else if(lambda_a_xx .ne. lambda_a_xz) then
            write(*,*) 'Wavelength discrepancy in files Kxx_a,Kxz_a'
         end if

         lambda = lambda_a_xx

      else if(nord .le. 7) then

         read(75,*,err=102, end=103) lambda_t_xx, 
     $                             Kxx_r(0), Kxx_i(0), 
     $                             Kxx_r(2), Kxx_i(2),
     $                             Kxx_r(3), Kxx_i(3),
     $                             Kxx_r(4), Kxx_i(4), 
     $                             Kxx_r(5), Kxx_i(5),
     $                             Kxx_r(6), Kxx_i(6),
     $                             Kxx_r(7), Kxx_i(7)
         read(76,*,err=102, end=103) lambda_t_yy, 
     $                             Kyy_r(0), Kyy_i(0), 
     $                             Kyy_r(2), Kyy_i(2),
     $                             Kyy_r(3), Kyy_i(3),
     $                             Kyy_r(4), Kyy_i(4), 
     $                             Kyy_r(5), Kyy_i(5),
     $                             Kyy_r(6), Kyy_i(6),
     $                             Kyy_r(7), Kyy_i(7)
         read(77,*,err=102, end=103) lambda_t_zz, 
     $                             Kzz_r(0), Kzz_i(0), 
     $                             Kzz_r(2), Kzz_i(2),
     $                             Kzz_r(3), Kzz_i(3),
     $                             Kzz_r(4), Kzz_i(4), 
     $                             Kzz_r(5), Kzz_i(5),
     $                             Kzz_r(6), Kzz_i(6),
     $                             Kzz_r(7), Kzz_i(7)
         read(78,*,err=102, end=103) lambda_t_xz, 
     $                             Kxz_r(0), Kxz_i(0), 
     $                             Kxz_r(2), Kxz_i(2),
     $                             Kxz_r(3), Kxz_i(3),
     $                             Kxz_r(4), Kxz_i(4), 
     $                             Kxz_r(5), Kxz_i(5),
     $                             Kxz_r(6), Kxz_i(6),
     $                             Kxz_r(7), Kxz_i(7)
                                   
         if     (lambda_t_xx .ne. lambda_t_yy) then
            write(*,*) 'Wavelength discrepancy in files Kxx_t,Kyy_t'
         else if(lambda_t_xx .ne. lambda_t_zz) then
            write(*,*) 'Wavelength discrepancy in files Kxx_t,Kzz_t'
         else if(lambda_t_xx .ne. lambda_t_xz) then
            write(*,*) 'Wavelength discrepancy in files Kxx_t,Kxz_t'
         end if

         lambda = lambda_t_xx

      else

         write(*,*) 'Unexpected error 2; exiting'
         write(*,*) 'Please report this error to code developers'
         stop

      end if 
      
      write(*,*) kl, ':   lambda=', sngl(lambda), '[nm]'
      do 20, n=0,ianord
         Kxx(n) = Kxx_r(n) + cu*Kxx_i(n)
         Kyy(n) = Kyy_r(n) + cu*Kyy_i(n)
         Kzz(n) = Kzz_r(n) + cu*Kzz_i(n)
         Kxz(n) = Kxz_r(n) + cu*Kxz_i(n)
 20   continue
                                               
      k1     = twopi*dsqrt(eps_1)/lambda
      Lhat   = L*k1

      GRxx = cz
      GRyy = cz
      GRzz = cz
      GRxz = cz
      do 40, n=0,ianord
         nn = n-3
         Lhat_nn = Lhat**nn
         GRxx = GRxx + Kxx(n)*Lhat_nn
         GRyy = GRyy + Kyy(n)*Lhat_nn
         GRzz = GRzz + Kzz(n)*Lhat_nn
         GRxz = GRxz + Kxz(n)*Lhat_nn
 40   continue


      if((output_prec.eq.'S') .or. (output_prec.eq.'s')) then 

         write(81,6,err=102) sngl(lambda),
     $      sngl(dreal(GRxx)), sngl(dimag(GRxx))
         write(82,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRyy)), sngl(dimag(GRyy))
         write(83,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRzz)), sngl(dimag(GRzz))
         write(84,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRxz)), sngl(dimag(GRxz))

      else if((output_prec.eq.'D') .or. (output_prec.eq.'d')) then

         write(81,10,err=102) sngl(lambda),
     $                       dreal(GRxx), dimag(GRxx)
         write(82,10,err=102) sngl(lambda), 
     $                       dreal(GRyy), dimag(GRyy)
         write(83,10,err=102) sngl(lambda), 
     $                       dreal(GRzz), dimag(GRzz)
         write(84,10,err=102) sngl(lambda), 
     $                       dreal(GRxz), dimag(GRxz)

      else

          write(*,*) 'Unexpected error 3; exiting'
          write(*,*) 'Please report this error to code developers'
          stop

      end if


      goto 1
 103  continue
      nl = kl-1
      write(*,*) 'Tital lines read:', nl

      close(unit=81)
      close(unit=82)
      close(unit=83)
      close(unit=84)

      if(nord .lt. -3) then

         close(unit=71)
         close(unit=72)
         close(unit=73)
         close(unit=74)

         close(unit=75)
         close(unit=76)
         close(unit=77)
         close(unit=78)

      else if(nord .le. 3) then

         close(unit=71)
         close(unit=72)
         close(unit=73)
         close(unit=74)

      else if(nord .le. 7) then

         close(unit=75)
         close(unit=76)
         close(unit=77)
         close(unit=78)

      else

         write(*,*) 'Unexpected error 4; exiting'
         write(*,*) 'Please report this error to code developers'
         stop

      end if 


      stop

 101  continue
      write(*,*) 'Error opening or reading input file GF.par; exiting'
      stop

 102  continue
      write(*,*) 'Error opening or writing to an output file; exiting'
      stop

 2    format(a1)
 4    format(a3)
 6    format(3G16.7E3)
 8    format(15G16.7E3)
 10   format(3G23.14E3)
 12   format(15G23.14E3)

      end
