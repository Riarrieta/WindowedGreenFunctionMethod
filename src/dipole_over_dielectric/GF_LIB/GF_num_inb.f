      program GF_num_inb ! Numerical calculation of the half-space
                         ! electromagnetic Green's tensor
                         ! by discretization of the Sommerfeld integrals
                         !
                         ! This program uses intrinsic Bessel functions
c-----------------------------------------------------------------------
c               Authors:
c                  George Y. Panasyuk, Vadim A. Markel, John C. Schotland
c                  University of Pennsylvania
c               E-mail address for questions and comments:
c                  vmarkel@mail.med.upenn.edu
c               For updates and additional info check
c                  URL: http://whale.seas.upenn.edu/CODES/
c--------------------------------------------------------------------------
      implicit none         

      integer*4 key, nl, kl, idummy
      logical*4 logscale
      character*4 hdummy, output_prec, lambda_scale

      integer*8 N1, N2, N3, k

      real*8 ru, rz, two, three, four, half, pi, twopi, power
      real*8 lambda, lambda_max, lambda_min, dlambda, rlambda, lambda_p
      real*8 w, J0, J1, J0Hat, J1Hat
      real*8 y, y2, y3, y4, dy1, dy2, dy3, y_max
      real*8 t, gd, delta
      real*8 Rhat, Zhat_m, Zhat_p
      real*8 z1, z2, z_p, z_m, rho
      real*8 Rti, k1
      real*8 rpe, div, r_p, r_m
      real*8 xh, xh_2, zh_p, zh_p_2, zh_m, zh_m_2
      real*8 RJ01p, RJ0p, RJ1p
      real*8 DBESJ0, DBESJ1
      real*8 eps_2_r, eps_2_i, eps_1, eps0

      complex*16 cu, cz, eps, eps_2, eps1, kzeps
      complex*16 R, Rpr
      complex*16 c2, c3, cpe
      complex*16 J01p, J0p, J1p
      complex*16 s1xx,s1yy,s1xz,s1zz,s2xx,s2yy,s2xz,s2zz
      complex*16 sum1xz,sum2xz
      complex*16 sum1xx_a,sum1xx_b,sum1xx_c
      complex*16 sum1yy_a,sum1yy_b,sum1yy_c
      complex*16 sum1zz_a,sum1zz_b,sum1zz_c
      complex*16 sum2xx_a,sum2xx_b,sum2xx_c
      complex*16 sum2yy_a,sum2yy_b,sum2yy_c
      complex*16 sum2zz_a,sum2zz_b,sum2zz_c
      complex*16 Gxx1a,Gxx1b,Gxx2a1,Gxx2a2,Gxx2ah,Gxx2bh,Gxx2b1,Gxx2b2
      complex*16 Gyy1a,Gyy1b,Gyy2a1,Gyy2a2,Gyy2ah,Gyy2bh,Gyy2b1,Gyy2b2
      complex*16 Gzz1,Gzz2b1,Gzz2b2,Gzz2
      complex*16 Gxz1,Gxz2a1,Gxz2a2,Gxz2
      complex*16 GTxx,GTyy,GTzz,GTxz,GTzx
      complex*16 GFxx,GFyy,GFzz,GFxz,GFzx
      complex*16 GRxx,GRyy,GRzz,GRxz,GRzx

c
c------ Mathematical constants ------------------------------------
c
      rz    = 0.0d0
      ru    = 1.0d0
      cz    = (0.0d0,0.0d0)
      cu    = (0.0d0,1.0d0)
      two   = 2.0d0
      three = 3.0d0
      four  = 4.0d0
      half  = 0.5d0
      pi    = two*dasin(1.0d0)
      twopi = two*pi
      delta = 1.d-100
c
c---------- Reading input parameters -----------------------------------
c
      open(unit=70, file='GF.par', status='old', err=101)
      read(70,*) hdummy
      read(70,*) hdummy

      read(70,*) key
      if(key .eq. 0) then
         write(*,*) 'Assuming a metallic substrate'
         write(*,*) 'Drude formula will be used for eps_2'
      else if (key .eq. 1) then
         write(*,*) 'Assuming that the substrate has a'
         write(*,*) '   ... wavelength-independent permittivity'
         write(*,*) '   ... eps_2 = eps_2_r + i*eps_2_i'
         write(*,*) '   ... where eps_2_r, eps_2_i are given in GF.par'
         write(*,*) 'If you dont like this, use KEY=2'
         write(*,*) '   ... and program your own formula for eps_2'
      else if (key .eq. 2) then
         write(*,*) 'KEY=2 is used for user-defined permittivity'
         write(*,*) 'This option requires additional programming'
         write(*,*) 'But this has not been done so far'
         write(*,*) 'Go to the line "if(key .eq. 2)" below'
         write(*,*) ' ...and add your code directly after that line'
         write(*,*) 'Then comment this message out'
         write(*,*) 'For now, I am exiting. Goodby'
         stop
      else 
         write(*,*) 'Error. Incorrect parameter: KEY=', key
         write(*,*) 'KEY is allowed to take only the values 0,1,2'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*) idummy

      read(70,2) output_prec
      if(output_prec .eq. 'S') goto 3
      if(output_prec .eq. 's') goto 3
      if(output_prec .eq. 'D') goto 3
      if(output_prec .eq. 'd') goto 3
         write(*,*) 'Incorrect parameter OUTPUT_PREC'
         write(*,*) 'Only chars S, s, D, d are allowed'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
 3    continue

      read(70,4) lambda_scale
      if(lambda_scale .eq. 'LOG') goto 5
      if(lambda_scale .eq. 'log') goto 5
      if(lambda_scale .eq. 'LIN') goto 5
      if(lambda_scale .eq. 'lin') goto 5
         write(*,*) 'Incorrect parameter LAMBDA_SCALE'
         write(*,*) 'Only "LOG", "log", "LIN", "lin" are allowed'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
 5    continue

      read(70,*) lambda_min, lambda_max, nl
      write(*,*) 'lambda_min,lambda_max:', sngl(lambda_min), 
     $                                     sngl(lambda_max), '  [nm]'
      write(*,*) 'nl:', nl
      if(lambda_min.lt.rz) then
         write(*,*) 'Error. Incorrect parameters: lambda_min<0'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      else if(lambda_min.gt.lambda_max) then
         write(*,*) 'Error. Incorrect parameters: lambda_min>lambda_max'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      else if((lambda_min.eq.lambda_max) .and. (nl.ne.1) ) then
         write(*,*) 'Warning. lambda_min=lambda_max but nl=/=1'
         write(*,*) 'Setting nl=1 and continuing'
         nl=1
      else if((lambda_min .lt. lambda_max) .and. (nl.lt.2)) then
         write(*,*) 'Warning. lambda_min<lambda_max but nl<2'
         write(*,*) 'Setting nl=2 and continuing'
         nl=2
      end if
       
      read(70,*) z1, z2, rho
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

      read(70,*) eps_2_r, eps_2_i
      if(eps_2_i.lt.rz) then
         write(*,*) 'Error. Incorrect parameter: eps_2_i<0; eps_2_i=', 
     $                                                     eps_2_i
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*) lambda_p
      if(lambda_p.le.rz) then
         write(*,*) 'Error. Incorrect parameter: lambda_p<=0; 
     $                                           lambda_p=', lambda_p
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*) gd
      if(gd.lt.rz) then
         write(*,*) 'Error. Incorrect parameter: gd<0; gd=', gd
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      else if(gd. lt. 0.0002) then
         write(*,*) 'Warning. Unrealistically small value: gd=', gd
         write(*,*) 'This can cause slow numerical convergence'
         write(*,*) 'Check convergence by increasing N1, N2, N3'
         write(*,*) '   ... by the factor of 10'
         write(*,*) 'Or better yet, use GD>=0.002 (the value for Ag)'
      end if

      read(70,*) eps0

      read(70,*) eps_1
      if(eps_1.le.rz) then
         write(*,*) 'Error. Incorrect parameter: eps_1<=0; eps_1=', 
     $                                                     eps_1
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*) N1
      read(70,*) N2
      read(70,*) N3
      if(N1.lt.10) then
         write(*,*) 'Error. Incorrect parameter: N1<10; N1=', N1
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if
      if(N2.lt.10) then
         write(*,*) 'Error. Incorrect parameter: N2<10; N2=', N2
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if
      if(N3.lt.10) then
         write(*,*) 'Error. Incorrect parameter: N2<10; N2=', N2
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if
      if((N1 .ne. N2) .or. (N1 .ne. N3)) then
         write(*,*) 'You are using different values for N1,N2,N3'
         write(*,*) 'Are you sure you know what you are doing?'
      end if
      
      read(70,*) y_max, div
      if(y_max.le.rz) then
         write(*,*) 'Error. Incorrect parameter: y_max<=0; y_max=', 
     $                                                     y_max
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if
      if(div.le.rz) then
         write(*,*) 'Error. Incorrect parameter: div<=0; div=', div
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if
      
      close(unit=70)
c
c------- End reading parameters -----------------------------
c


c
c------- Input-dependent constants --------------------------
c
      if(nl .gt. 1) then
         dlambda = (lambda_max - lambda_min)/dfloat(nl-1)
      else
         dlambda = rz
      end if

      z_p = z2 + z1
      z_m = z2 - z1
      r_p = dsqrt(rho*rho + z_p*z_p)
      r_m = dsqrt(rho*rho + z_m*z_m)
      zh_p   = z_p/r_p
      zh_p_2 = zh_p*zh_p
      zh_m   = z_m/r_m   ! used only for G^F
      zh_m_2 = zh_m*zh_m
      xh     = rho/r_m      ! used only for G^F
      xh_2   = xh*xh

      N1  = 2*(N1/2)
      N2  = 2*(N2/2)
      N3  = 2*(N3/2)
      dy1 = ru/dfloat(N1)
      dy2 = div/dfloat(N2)
      dy3 = (y_max-div) / dfloat(N3)

      logscale = .false.
      if(lambda_scale .eq. 'LOG') logscale=.true.
      if(lambda_scale .eq. 'log') logscale=.true.
      if(logscale) then
         if(lambda_min .eq. rz) then
            write(*,*) 'If using logscale, lambda_min=0 is not allowed'
            write(*,*) 'Correct GF.par and try again; exiting'
            stop
         end if
         rlambda = lambda_max/lambda_min
      end if
c
c-------- End input-dependent constants block ------------------------
c
      open(71, file='GRxx_num', status='unknown', err=102)
      open(72, file='GRyy_num', status='unknown', err=102)
      open(73, file='GRzz_num', status='unknown', err=102)
      open(74, file='GRxz_num', status='unknown', err=102)
      open(75, file='GRzx_num', status='unknown', err=102)

      open(76, file='GTxx_num', status='unknown', err=102)
      open(77, file='GTyy_num', status='unknown', err=102)
      open(78, file='GTzz_num', status='unknown', err=102)
      open(79, file='GTxz_num', status='unknown', err=102)
      open(80, file='GTzx_num', status='unknown', err=102)
c
c------- Start loop over wavelengths ---------------------------------
c
      do 1, kl=1,nl
      if(logscale) then
         power = dfloat(kl-1)/dfloat(nl-1)
         lambda = lambda_min * rlambda**power
      else
         lambda = lambda_min + dlambda*dfloat(kl-1) 
      end if
      write(*,*) kl, ':   lambda=', sngl(lambda), '[nm]'

      k1 = twopi*dsqrt(eps_1)/lambda

      if(key .eq. 0) then
         t = lambda_p/lambda
         eps_2 = eps0 - ru/t/(t + cu*gd)
      else if(key .eq. 1) then
         eps_2 = eps_2_r + cu*eps_2_i
      else if(key .eq. 2) then
         ! Enter your own formula for eps_2 here
         continue
         ! Do not edit past this line
      else
         write(*,*) 'Unexpected error 1; exiting'
         write(*,*) 'Please report this error to code developers'
         stop
      end if 
      
      eps_2 = dreal(eps_2) + cu*max(delta,dimag(eps_2))
      eps  = eps_2/eps_1
      eps1 = eps - ru

      Zhat_p = k1*z_p       ! Used for G^R
      Zhat_m = k1*z_m       ! Used for G^F
      Rhat   = k1*rho     
      Rti    = dsqrt(Rhat*Rhat + Zhat_m*Zhat_m)   ! Argument for G^F

      J0Hat = DBESJ0(Rhat)
      J1Hat = DBESJ1(Rhat)/Rhat
c
c-------- Integrals over the finite interval (0,1) ---------
c                 dy = dy1 = 1/N1   
c
      s2xx = -J1Hat 
      s2yy =  rz    
      s2zz =  J0Hat 
      s2xz =  rz    

      kzeps = cdsqrt(eps1 + ru)
      R   = (ru-kzeps) /(ru+kzeps)
      Rpr = (kzeps-eps)/(kzeps+eps)
      cpe = half*cdexp(cu*Zhat_p)
      s1xx = s2xx + cpe*R    
      s1yy = s2yy + cpe*Rpr  
      s1zz = s2zz             
      s1xz = s2xz             

      sum1xx_a = cz
      sum1yy_a = cz
      sum1zz_a = cz
      sum1xx_b = cz
      sum1yy_b = cz
      sum1xz   = cz
      do 11,k=1,N1-1,2
         y     = k*dy1
         y2    = y*y
         y3    = ru-y2
         y4    = y*y3
         kzeps = cdsqrt(eps1 + y2)
         R     = (y-kzeps)/(y+kzeps)
         Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
         cpe   = cdexp(cu*y*Zhat_p)
         w     = Rhat*dsqrt(y3)
            J0 = DBESJ0(w)
            J1 = DBESJ1(w)/w
         J01p  = (J0-J1)*cpe
         J1p   = J1*cpe
         J0p   = J0*cpe
         sum1xx_a = sum1xx_a + J1p  * R
         sum1yy_a = sum1yy_a + J1p  * Rpr * y2
         sum1zz_a = sum1zz_a + J0p  * Rpr * y3
         sum1xx_b = sum1xx_b + J01p * Rpr * y2
         sum1yy_b = sum1yy_b + J01p * R
         sum1xz   = sum1xz   + J1p  * Rpr * y4
 11   continue

      sum2xx_a = cz
      sum2yy_a = cz
      sum2zz_a = cz
      sum2xx_b = cz
      sum2yy_b = cz
      sum2xz   = cz
      do 12, k=2,N1-2,2
         y     = k*dy1
         y2    = y*y
         y3    = ru-y2
         y4    = y*y3
         kzeps = cdsqrt(eps1 + y2)
         R     = (y-kzeps)/(y+kzeps)
         Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
         cpe   = cdexp(cu*y*Zhat_p)
         w     = Rhat*dsqrt(y3)
            J0 = DBESJ0(w)
            J1 = DBESJ1(w)/w
         J01p  = (J0-J1)*cpe
         J1p   = J1*cpe
         J0p   = J0*cpe
         sum2xx_a = sum2xx_a + J1p  * R             
         sum2yy_a = sum2yy_a + J1p  * Rpr * y2
         sum2zz_a = sum2zz_a + J0p  * Rpr * y3
         sum2xx_b = sum2xx_b + J01p * Rpr * y2
         sum2yy_b = sum2yy_b + J01p * R
         sum2xz   = sum2xz   + J1p  * Rpr * y4
 12   continue

      Gxx1a = dy1*(two*sum2xx_a + four*sum1xx_a + s1xx)/three
      Gyy1a = dy1*(two*sum2yy_a + four*sum1yy_a + s1yy)/three
      Gzz1  = dy1*(two*sum2zz_a + four*sum1zz_a + s1zz)/three
      Gxz1  = dy1*(two*sum2xz   + four*sum1xz   + s1xz)/three
        
      s2xx = rz
      s2yy = J1Hat-J0Hat

      kzeps = cdsqrt(eps1 + ru)
      Rpr   = (kzeps-eps)/(kzeps+eps)
      R     = (ru-kzeps)/(ru+kzeps)
      cpe   = half*cdexp(cu*Zhat_p)
      s1xx  = s2xx + cpe*Rpr
      s1yy  = s2yy + cpe*R

      Gxx1b = dy1*(two*sum2xx_b + four*sum1xx_b + s1xx)/three
      Gyy1b = dy1*(two*sum2yy_b + four*sum1yy_b + s1yy)/three
c
c-------- Integrals over the infinite interval (0,Inf) ---------
c-------- The integlas is broken in two parts ------------------
c-------- a) First, compute over the interval (0,div)
c                 dy = dy2 = div/N2
c
      s2xx = -J1Hat
      s2yy =  rz
      s2xz =  rz

      y     = div
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      R     = (y-kzeps)/(y+kzeps)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      w     = Rhat*dsqrt(ru+y2)
         J1 = DBESJ1(w)/w
      s1xx  = s2xx + rpe*J1*R
      s1yy  = s2yy + y2*rpe*J1*Rpr
      s1xz  = s2xz + y*(y2+ru)*rpe*J1*Rpr

      sum1xx_a = cz
      sum1yy_a = cz
      sum1xx_b = cz
      sum1yy_b = cz
      sum1zz_b = cz
      sum1xz   = cz
      do 21, k = 1, N2-1, 2
         y      = k*dy2
         y2     = y*y
         y3     = ru+y2
         y4     = y*y3
         rpe    = dexp(-y*Zhat_p)
         kzeps  = cdsqrt(y2-eps1)
         R      = (y-kzeps)/(y+kzeps)
         Rpr    = (kzeps-eps*y)/(kzeps+eps*y)
         w      = Rhat*dsqrt(y3)
             J0 = DBESJ0(w)
             J1 = DBESJ1(w)/w
         RJ01p  = (J0-J1)*rpe
         RJ1p   = J1*rpe
         RJ0p   = J0*rpe
         sum1xx_a = sum1xx_a + RJ1p  * R
         sum1yy_a = sum1yy_a + RJ1p  * Rpr * y2
         sum1xx_b = sum1xx_b + RJ01p * Rpr * y2
         sum1yy_b = sum1yy_b + RJ01p * R
         sum1zz_b = sum1zz_b + RJ0p  * Rpr * y3
         sum1xz   = sum1xz   + RJ1p  * Rpr * y4 
 21   continue

      sum2xx_a = cz
      sum2yy_a = cz
      sum2xx_b = cz
      sum2yy_b = cz
      sum2zz_b = cz
      sum2xz   = cz
      do 22, k = 2, N2-2, 2
         y      = k*dy2
         y2     = y*y
         y3     = ru+y2
         y4     = y*y3
         rpe    = dexp(-y*Zhat_p)
         kzeps  = cdsqrt(y2-eps1)
         R      = (y-kzeps)/(y+kzeps)
         Rpr    = (kzeps-eps*y)/(kzeps+eps*y)
         w      = Rhat*dsqrt(y3)
             J0 = DBESJ0(w)
             J1 = DBESJ1(w)/w
         RJ01p  = (J0-J1)*rpe
         RJ1p   = J1*rpe
         RJ0p   = J0*rpe
         sum2xx_a = sum2xx_a + RJ1p  * R
         sum2yy_a = sum2yy_a + RJ1p  * Rpr * y2
         sum2xx_b = sum2xx_b + RJ01p * Rpr * y2
         sum2yy_b = sum2yy_b + RJ01p * R
         sum2zz_b = sum2zz_b + RJ0p  * Rpr * y3
         sum2xz   = sum2xz   + RJ1p  * Rpr * y4
 22   continue

      Gxx2a1 = dy2*(two*sum2xx_a + four*sum1xx_a + s1xx)/three
      Gyy2a1 = dy2*(two*sum2yy_a + four*sum1yy_a + s1yy)/three
      Gxz2a1 = dy2*(two*sum2xz   + four*sum1xz   + s1xz)/three
c
c-------- Integrals over the infinite interval (0,Inf) ---------
c-------- The integlas is broken in two parts ------------------
c-------- b) Second part, compute over the interval (div,y_max)
c                 dy = dy3 = (y_max-div)/N3
c
      y     = div
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      R     = (y-kzeps)/(y+kzeps)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      w     = Rhat*dsqrt(ru+y2)
         J1 = DBESJ1(w)/w
      s2xx  = rpe*J1*R
      s2yy  = y2*rpe*J1*Rpr
      s2xz  = y*(y2+ru)*rpe*J1*Rpr

      y     = y_max
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      R     = (y-kzeps)/(y+kzeps)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      w     = Rhat*dsqrt(ru+y2)
         J1 = DBESJ1(w)/w
      s1xx  = s2xx + rpe*J1*R
      s1yy  = s2yy + y2*rpe*J1*Rpr
      s1xz  = s2xz + y*(y2+ru)*rpe*J1*Rpr

      sum1xx_a = cz
      sum1yy_a = cz
      sum1xz   = cz
      sum1xx_c = cz
      sum1yy_c = cz
      sum1zz_c = cz
      do 31, k = 1, N3-1, 2
         y      = div + k*dy3
         y2     = y*y
         y3     = ru+y2
         y4     = y*y3
         rpe    = dexp(-y*Zhat_p)
         kzeps  = cdsqrt(y2-eps1)
         R      = (y-kzeps)/(y+kzeps)
         Rpr    = (kzeps-eps*y)/(kzeps+eps*y)
         w      = Rhat*dsqrt(y3)
             J0 = DBESJ0(w)
             J1 = DBESJ1(w)/w
         RJ01p  = (J0-J1)*rpe
         RJ1p   = J1*rpe
         RJ0p   = J0*rpe
         sum1xx_a = sum1xx_a + RJ1p  * R
         sum1yy_a = sum1yy_a + RJ1p  * Rpr * y2
         sum1xz   = sum1xz   + RJ1p  * Rpr * y4
         sum1xx_c = sum1xx_c + RJ01p * Rpr * y2
         sum1yy_c = sum1yy_c + RJ01p * R
         sum1zz_c = sum1zz_c + RJ0p  * Rpr * y3
 31   continue
        
      sum2xx_a = cz
      sum2yy_a = cz
      sum2xz   = cz
      sum2xx_c = cz
      sum2yy_c = cz
      sum2zz_c = cz
      do 32, k = 2, N3-2, 2
         y      = div + k*dy3
         y2     = y*y
         y3     = ru+y2
         y4     = y*y3
         rpe    = dexp(-y*Zhat_p)
         kzeps  = cdsqrt(y2-eps1)
         R      = (y-kzeps)/(y+kzeps)
         Rpr    = (kzeps-eps*y)/(kzeps+eps*y)
         w      = Rhat*dsqrt(y3)
             J0 = DBESJ0(w)
             J1 = DBESJ1(w)/w
         RJ01p  = (J0-J1)*rpe
         RJ1p   = J1*rpe
         RJ0p   = J0*rpe
         sum2xx_a = sum2xx_a + RJ1p  * R
         sum2yy_a = sum2yy_a + RJ1p  * Rpr * y2
         sum2xz   = sum2xz   + RJ1p  * Rpr * y4
         sum2xx_c = sum2xx_c + RJ01p * Rpr * y2
         sum2yy_c = sum2yy_c + RJ01p * R  
         sum2zz_c = sum2zz_c + RJ0p  * Rpr * y3
 32   continue

      Gxx2a2 = dy3*(two*sum2xx_a + four*sum1xx_a + s1xx)/three
      Gyy2a2 = dy3*(two*sum2yy_a + four*sum1yy_a + s1yy)/three
      Gxz2a2 = dy3*(two*sum2xz   + four*sum1xz   + s1xz)/three
      Gxx2ah = Gxx2a1 + Gxx2a2
      Gyy2ah = Gyy2a1 + Gyy2a2
      Gxz2   = Gxz2a1 + Gxz2a2

      s2xx = rz
      s2yy = J1Hat - J0Hat
      s2zz = J0Hat

      y     = div
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      R     = (y-kzeps)/(y+kzeps)
      w     = Rhat*dsqrt(ru+y2)
         J0 = DBESJ0(w)
         J1 = DBESJ1(w)/w
      s1xx  = s2xx + y2*rpe*Rpr*(J0-J1)
      s1yy  = s2yy + rpe*R*(J0-J1)
      s1zz  = s2zz + (y2+ru)*rpe*Rpr*J0

      Gxx2b1 = dy2*(two*sum2xx_b + four*sum1xx_b + s1xx)/three
      Gyy2b1 = dy2*(two*sum2yy_b + four*sum1yy_b + s1yy)/three
      Gzz2b1 = dy2*(two*sum2zz_b + four*sum1zz_b + s1zz)/three

      y     = div
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      R     = (y-kzeps)/(y+kzeps)
      w     = Rhat*dsqrt(ru+y2)
         J0 = DBESJ0(w)
         J1 = DBESJ1(w)/w
      s2xx  = y2*rpe*Rpr*(J0-J1)
      s2yy  = rpe*R*(J0-J1)
      s2zz  = (y2+ru)*rpe*Rpr*J0

      y     = y_max
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      R     = (y-kzeps)/(y+kzeps)
      w     = Rhat*dsqrt(ru+y2)
         J0 = DBESJ0(w)
         J1 = DBESJ1(w)/w
      s1xx  = s2xx + y2*rpe*Rpr*(J0-J1)
      s1yy  = s2yy + rpe*R*(J0-J1)
      s1zz  = s2zz + (y2+ru)*rpe*Rpr*J0

      Gxx2b2 = dy3*(two*sum2xx_c + four*sum1xx_c + s1xx)/three
      Gyy2b2 = dy3*(two*sum2yy_c + four*sum1yy_c + s1yy)/three
      Gzz2b2 = dy3*(two*sum2zz_c + four*sum1zz_c + s1zz)/three

      Gxx2bh = Gxx2b1 + Gxx2b2
      Gyy2bh = Gyy2b1 + Gyy2b2
      Gzz2   = Gzz2b1 + Gzz2b2

c     Dimensionless G^R/k1**3
      GRxx = cu*(Gxx1a + Gxx1b) + Gxx2ah - Gxx2bh
      GRyy = cu*(Gyy1a + Gyy1b) + Gyy2bh - Gyy2ah
      GRzz = -cu*Gzz1 - Gzz2
      GRxz = -Rhat*(Gxz1 + Gxz2)
      GRzx = -GRxz

c     Dimensionless G^F/k1**3
      c2 = cdexp(cu*Rti)/Rti 
      c3 = (cu - ru/Rti)/Rti
      GFxx = c2*(ru - xh_2 + c3*(ru - three*xh_2))  
      GFyy = c2*(ru + c3)                           
      GFzz = c2*(ru - zh_m_2 + c3*(ru-three*zh_m_2))
      GFxz = c2*xh*zh_m*(three*c3 + ru)
      GFzx = GFxz

c     Dimensionless G^T/k1**3, where G^T=G^F+G^R
      GTxx =  GRxx + GFxx
      GTyy =  GRyy + GFyy     
      GTzz =  GRzz + GFzz     
      GTxz =  GRxz + GFxz 
      GTzx =  GRzx + GFzx

      if((output_prec.eq.'S') .or. (output_prec.eq.'s')) then 

         write(71,6,err=102) sngl(lambda),
     $      sngl(dreal(GRxx)), sngl(dimag(GRxx))
         write(72,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRyy)), sngl(dimag(GRyy))
         write(73,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRzz)), sngl(dimag(GRzz))
         write(74,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRxz)), sngl(dimag(GRxz))
         write(75,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRzx)), sngl(dimag(GRzx))


         write(76,6,err=102) sngl(lambda), 
     $      sngl(dreal(GTxx)), sngl(dimag(GTxx))
         write(77,6,err=102) sngl(lambda), 
     $      sngl(dreal(GTyy)), sngl(dimag(GTyy))
         write(78,6,err=102) sngl(lambda), 
     $      sngl(dreal(GTzz)), sngl(dimag(GTzz))
         write(79,6,err=102) sngl(lambda), 
     $      sngl(dreal(GTxz)), sngl(dimag(GTxz))
         write(80,6,err=102) sngl(lambda), 
     $      sngl(dreal(GTzx)), sngl(dimag(GTzx))

      else if((output_prec.eq.'D') .or. (output_prec.eq.'d')) then

         write(71,8,err=102) sngl(lambda), dreal(GRxx), dimag(GRxx)
         write(72,8,err=102) sngl(lambda), dreal(GRyy), dimag(GRyy)
         write(73,8,err=102) sngl(lambda), dreal(GRzz), dimag(GRzz)
         write(74,8,err=102) sngl(lambda), dreal(GRxz), dimag(GRxz)
         write(75,8,err=102) sngl(lambda), dreal(GRzx), dimag(GRzx)

         write(76,8,err=102) sngl(lambda), dreal(GTxx), dimag(GTxx)
         write(77,8,err=102) sngl(lambda), dreal(GTyy), dimag(GTyy)
         write(78,8,err=102) sngl(lambda), dreal(GTzz), dimag(GTzz)
         write(79,8,err=102) sngl(lambda), dreal(GTxz), dimag(GTxz)
         write(80,8,err=102) sngl(lambda), dreal(GTzx), dimag(GTzx)

      else

          write(*,*) 'Unexpected error 2; exiting'
          write(*,*) 'Please report this error to code developers'
          stop

      end if

c$$$  These lines to be used for debugging only
c$$$       write(*,*) '-----------------------------------------'
c$$$       write(*,*) 'XX-I:', dimag(GRxx) 
c$$$       write(*,*) 'XX-R:', dreal(GRxx) 
c$$$       write(*,*) 'YY-I:', dimag(GRyy) 
c$$$       write(*,*) 'YY-R:', dreal(GRyy) 
c$$$       write(*,*) 'ZZ-I:', dimag(GRzz) 
c$$$       write(*,*) 'ZZ-R:', dreal(GRzz) 
c$$$       write(*,*) 'XZ-I:', dimag(GRxz) 
c$$$       write(*,*) 'XZ-R:', dreal(GRxz) 

 1    continue
      close(unit=71)
      close(unit=72)
      close(unit=73)
      close(unit=74)
      close(unit=75)
      close(unit=76)
      close(unit=77)
      close(unit=78)
      close(unit=79)
      close(unit=80)

      stop

 101  continue
      write(*,*) 'Error opening input file GF.par; exiting'
      stop

 102  continue
      write(*,*) 'Error opening an output file or writing data; exiting'
      stop

 2    format(a1)
 4    format(a3)
 6    format(3G16.7E3)
 8    format(3G23.14E3)

      end
