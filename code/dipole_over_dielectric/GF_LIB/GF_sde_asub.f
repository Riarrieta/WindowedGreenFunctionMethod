      program GF_sde_asub   ! Calculation of the half-space
                            ! electromagnetic Green's tensor
                            ! by the short-dsitance expansion
                            ! for the case of a metal substrate.
                            !
                            ! ONLY the reflected part, GR, is computed
c------------------------------------------------------------------
c               Authors:
c                  George Y. Panasyuk, Vadim A. Markel, John C. Schotland
c                  University of Pennsylvania
c               E-mail address for questions and comments:
c                  vmarkel@mail.med.upenn.edu
c               For updates and additional info check
c                  URL: http://whale.seas.upenn.edu/CODES/
c------------------------------------------------------------------
      implicit none

      integer*4 idummy, kl, nl, key, nord
      logical*4 logscale
      character*4 hdummy, output_prec, lambda_scale

      real*8 rdummy
      real*8 t, gd, ru, rz, two, three, twopi
      real*8 z1, z2, z_p, eps0, dlambda, rlambda, power
      real*8 lambda_p, lambda, lambda_min, lambda_max, k1, rho
      real*8 zh_p, zh_p_2, eps_1, eps_2_r, eps_2_i
      real*8 L, Lhat, Lhat_2, Lhat_3, Zhat, Rhat, RZLhat, LZRhat

      complex*16 cu, cz
      complex*16 eps, eps_m, eps_p, eps_s, eps_c
      complex*16 c1, c2, c3, c4, c5, c6, c7, Lam
      complex*16 eps_2, s
      complex*16 GRxx(0:3), GRyy(0:3), GRzz(0:3), GRxz(0:3)

      complex*16 Kxx0, Kyy0, Kzz0, Kxz0
      complex*16 Kxx2, Kyy2, Kzz2, Kxz2
      complex*16 Kxx3, Kyy3, Kzz3, Kxz3

      complex*16 GRxx0, GRyy0, GRzz0, GRxz0
      complex*16 GRxx2, GRyy2, GRzz2, GRxz2
      complex*16 GRxx3, GRyy3, GRzz3, GRxz3
c
c------ Mathematical constants ------------------------------------
c
      rz = 0.0d0
      ru = 1.0d0
      cu = (0.0d0, 1.0d0)
      cz = (0.0d0, 0.0d0)

      two   = 2.0d0
      three = 3.0d0
      twopi = 4.0d0*dasin(ru)
c
c---------- Reading input parameters -----------------------------------
c
      open(unit=70, file='GF.par', status='old', err=101)
      read(70,*) hdummy
      read(70,*) hdummy

      read(70,*, err=101) key
      if(key .eq. 0) then
         write(*,*) 'Assuming a metallic substrate;'
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
         write(*,*) 'key is allowed to take only the values 0,1,2'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*,err=101) nord
      if(nord .lt. 0) then
         write(*,*) 'Warning. You have requested NORD<0.'
         write(*,*) 'This option is used only in GF_k2g'
         write(*,*) 'For GF_sde_asub, the sign of NORD does not matter'
         write(*,*) 'I will set NORD=|NORD| and continue.'
         nord = -nord
         write(*,*) 'NORD=', nord
         write(*,*) '--------------------------------------------'
      end if
      if (nord .eq. 1) then
         write(*,*) 'The first order correction is identicallyz zero'
         write(*,*) '    (see paper)'
         write(*,*) 'Therefore, there is no difference between NORD=0'
         write(*,*) '    and NORD=1'
         write(*,*) 'I will set NORD=0 and continue'
         write(*,*) 'The result is not affected. Goodby.'
         nord = 0
      else if(nord .gt. 3) then 
         write(*,*) 'Warning. You have requested NORD>3'
         write(*,*) 'GF_sde_asub only allows NORD<=3'
         write(*,*) 'Note that GF_sde_tsub will allow |NORD|<=7'
         write(*,*) 'I will now set NORD=3 and continue'
         write(*,*) 'If this is not OK, either use GF_sde_tsub'
         write(*,*) '   ... or edit GF.par and try again'
         write(*,*) 'Continue with NORD=3?'
         pause
         nord = 3
         write(*,*) 'NORD=', nord
      else if(nord .gt. 7) then 
         write(*,*) 'Error. Incorrect parameter: NORD>7'
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

      read(70,4,err=101) lambda_scale
      if(lambda_scale .eq. 'LOG') goto 5
      if(lambda_scale .eq. 'log') goto 5
      if(lambda_scale .eq. 'LIN') goto 5
      if(lambda_scale .eq. 'lin') goto 5
         write(*,*) 'Incorrect parameter LAMBDA_SCALE'
         write(*,*) 'Only "LOG", "log", "LIN", "lin" are allowed'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
 5    continue

      read(70,*,err=101) lambda_min, lambda_max, nl
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

      read(70,*) eps_2_r, eps_2_i
      if(eps_2_i.lt.rz) then
         write(*,*) 'Error. Incorrect parameter: eps_2_i<0; eps_2_i=', 
     $                                                     eps_2_i
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*,err=101) lambda_p
      if(lambda_p.le.rz) then
         write(*,*) 'Error. Incorrect parameter: lambda_p<=0; 
     $                                           lambda_p=', lambda_p
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*,err=101) gd
      if(gd.lt.rz) then
         write(*,*) 'Error. Incorrect parameter: gd<0; gd=', gd
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      else if(gd. lt. 0.001) then
         write(*,*) 'Warning. Unrealistically small value: gd=', gd
         write(*,*) 'This can cause slow numerical convergence'
         write(*,*) 'Check convergence by increasing N1, N2, N3'
         write(*,*) '   ... by the factor of 10'
         write(*,*) 'Or better yet, use GD>=0.002 (the value for Ag)'
      end if

      read(70,*,err=101) eps0

      read(70,*,err=101) eps_1
      if(eps_1.le.rz) then
         write(*,*) 'Error. Incorrect parameter: eps_1<=0; eps_1=', 
     $                                                     eps_1
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      end if

      read(70,*) idummy
      read(70,*) idummy
      read(70,*) idummy
      
      read(70,*) rdummy, rdummy
      
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

      z_p    = z2 + z1
      L      = dsqrt(rho*rho + z_p*z_p)
      zh_p   = z_p/L
      zh_p_2 = zh_p*zh_p

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
      open(71, file='GRxx_asub', status='unknown', err=102)
      open(72, file='GRyy_asub', status='unknown', err=102)
      open(73, file='GRzz_asub', status='unknown', err=102)
      open(74, file='GRxz_asub', status='unknown', err=102)

      open(76, file='Kxx_asub', status='unknown', err=102)
      open(77, file='Kyy_asub', status='unknown', err=102)
      open(78, file='Kzz_asub', status='unknown', err=102)
      open(79, file='Kxz_asub', status='unknown', err=102)
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
      Lhat   = L*k1
      Lhat_2 = Lhat*Lhat
      Lhat_3 = Lhat_2*Lhat
      Zhat   = z_p*k1
      Rhat   = rho*k1
      LZRhat = (Lhat - Zhat)/Rhat
      RZLhat = Rhat*Zhat/Lhat_2

      if(key.eq.0) then
         t = lambda_p/lambda
         eps_2 = (eps0 - ru/t/(t + cu*gd))
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

      eps   = eps_2/eps_1
      eps_m = eps - ru
      eps_p = eps + ru
      eps_s = eps*eps
      eps_c = eps_s*eps
      s     = eps_m/eps_p

      c1 = cdsqrt(eps_p)
      c2 = eps_p*eps_p
      c3 = cdsqrt(eps)
      c4 = (ru + c1)/(eps + c3*c1)
      c5 = three*(c3 + ru)
      c6 = (ru - three*(c3-eps) + two*eps**2) / c5
      c7 = ru + c3 + two*eps*(ru+c3-eps*c3) - eps_c
      Lam = (eps/eps_m/c1)*cdlog(c4)

      Kxx3 = (cu*eps/c2)*(c6 + Lam)
      Kyy3 = Kxx3
      Kxz3 = rz
      Kzz3 = -(Two*cu/c2)*(c7/c5 + Lam*eps_s)

      Kxx0 = s*(three*zh_p_2 - Two)
      Kyy0 = s
      Kzz0 = s*(three*zh_p_2 - ru)
      Kxz0 = s*three*RZLhat

      Kxx2 = s*zh_p*((s - ru)/(ru + zh_p) + zh_p) / two
      Kyy2 = s*(s + zh_p)/Two/(ru + zh_p)
      Kzz2 = s*(zh_p_2 + s + Two)/Two
      Kxz2 = s*(eps*LZRhat/eps_p + RZLhat/two)

      GRxx0 = Kxx0/Lhat_3
      GRyy0 = Kyy0/Lhat_3
      GRzz0 = Kzz0/Lhat_3
      GRxz0 = Kxz0/Lhat_3
 
      GRxx2 = Kxx2/Lhat
      GRyy2 = Kyy2/Lhat
      GRzz2 = Kzz2/Lhat
      GRxz2 = Kxz2/Lhat

      GRxx3 = Kxx3
      GRyy3 = Kyy3
      GRzz3 = Kzz3
      GRxz3 = Kxz3

      GRxx(0) = GRxx0
      GRxz(0) = GRxz0
      GRyy(0) = GRyy0
      GRzz(0) = GRzz0

      GRxx(2) = GRxx(0) + GRxx2
      GRyy(2) = GRyy(0) + GRyy2
      GRzz(2) = GRzz(0) + GRzz2
      GRxz(2) = GRxz(0) + GRxz2

      GRxx(3) = GRxx(2) + GRxx3 
      GRyy(3) = GRyy(2) + GRyy3 
      GRzz(3) = GRzz(2) + GRzz3 
      GRxz(3) = GRxz(2) + GRxz3 

      if((output_prec.eq.'S') .or. (output_prec.eq.'s')) then 

         write(71,6,err=102) sngl(lambda),
     $      sngl(dreal(GRxx(nord))), sngl(dimag(GRxx(nord)))
         write(72,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRyy(nord))), sngl(dimag(GRyy(nord)))
         write(73,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRzz(nord))), sngl(dimag(GRzz(nord)))
         write(74,6,err=102) sngl(lambda), 
     $      sngl(dreal(GRxz(nord))), sngl(dimag(GRxz(nord)))

         write(76,8,err=102) sngl(lambda),
     $      sngl(dreal(Kxx0)), sngl(dimag(Kxx0)),
     $      sngl(dreal(Kxx2)), sngl(dimag(Kxx2)),
     $      sngl(dreal(Kxx3)), sngl(dimag(Kxx3))
         write(77,8,err=102) sngl(lambda),
     $      sngl(dreal(Kyy0)), sngl(dimag(Kyy0)),
     $      sngl(dreal(Kyy2)), sngl(dimag(Kyy2)),
     $      sngl(dreal(Kyy3)), sngl(dimag(Kyy3))
         write(78,8,err=102) sngl(lambda),
     $      sngl(dreal(Kzz0)), sngl(dimag(Kzz0)),
     $      sngl(dreal(Kzz2)), sngl(dimag(Kzz2)),
     $      sngl(dreal(Kzz3)), sngl(dimag(Kzz3))
         write(79,8,err=102) sngl(lambda),
     $      sngl(dreal(Kxz0)), sngl(dimag(Kxz0)),
     $      sngl(dreal(Kxz2)), sngl(dimag(Kxz2)),
     $      sngl(dreal(Kxz3)), sngl(dimag(Kxz3))

      else if((output_prec.eq.'D') .or. (output_prec.eq.'d')) then

         write(71,10,err=102) sngl(lambda),
     $                       dreal(GRxx(nord)), dimag(GRxx(nord))
         write(72,10,err=102) sngl(lambda), 
     $                       dreal(GRyy(nord)), dimag(GRyy(nord))
         write(73,10,err=102) sngl(lambda), 
     $                       dreal(GRzz(nord)), dimag(GRzz(nord))
         write(74,10,err=102) sngl(lambda), 
     $                       dreal(GRxz(nord)), dimag(GRxz(nord))

         write(76,12,err=102) sngl(lambda),
     $      dreal(Kxx0), dimag(Kxx0),
     $      dreal(Kxx2), dimag(Kxx2),
     $      dreal(Kxx3), dimag(Kxx3)
         write(77,12,err=102) sngl(lambda),
     $      dreal(Kyy0), dimag(Kyy0),
     $      dreal(Kyy2), dimag(Kyy2),
     $      dreal(Kyy3), dimag(Kyy3)
         write(78,12,err=102) sngl(lambda),
     $      dreal(Kzz0), dimag(Kzz0),
     $      dreal(Kzz2), dimag(Kzz2),
     $      dreal(Kzz3), dimag(Kzz3)
         write(79,12,err=102) sngl(lambda),
     $      dreal(Kxz0), dimag(Kxz0),
     $      dreal(Kxz2), dimag(Kxz2),
     $      dreal(Kxz3), dimag(Kxz3)

      else

          write(*,*) 'Unexpected error 2; exiting'
          write(*,*) 'Please report this error to code developers'
          stop

      end if

c$$$  These lines to be used for debugging only       
c$$$      write(*,*) '---------------------------------------------'
c$$$      write(*,*) 'XX:', dreal(GRxx(nord)), dimag(GRxx(nord))
c$$$      write(*,*) 'YY:', dreal(GRyy(nord)), dimag(GRyy(nord))
c$$$      write(*,*) 'ZZ:', dreal(GRzz(nord)), dimag(GRzz(nord))
c$$$      write(*,*) 'XZ:', dreal(GRxz(nord)), dimag(GRxz(nord))

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
 8    format(7G16.7E3)
 10   format(G16.7E3,2G23.14E3)
 12   format(G16.7E3,6G23.14E3)

      end
