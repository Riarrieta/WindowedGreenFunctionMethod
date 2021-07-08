      program GF_sde_tdss  ! Calculation of the half-space
                           ! electromagnetic Green's tensor
                           ! by the short-dsitance expansion
                           ! for the case of a transparent, 
                           ! (weakly) DISPERSIVE substrate.
                           !
                           ! ONLY the reflected part, GR, is computed
c
c--- Additional programming is required to use this code
c
c===== WARNING: Garbage in - garbage out ======================================
c
c--- It is the user's responsibility to make sure that the 
c--- permittivity of the substrate, eps_2,  is only WEAKLY DISPERSIVE,
c--- so that the assumption that eps_2 is purely real (which is used below)
c--- remains physically reasonable.
c
c--------------------------------------------------------------------------------
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
      real*8 p1, p2, p3, p4 
      real*8 ru, rz, pi, twopi

      real*8 k1, rho, z1, z2, z_p, s
      real*8 lambda, lambda_min, lambda_max, dlambda, rlambda, power
      real*8 eps, eps_1, eps_2, da, db
      real*8 eps_m,eps_p, ep2,ep3,ep4,ep5

      real*8 zh, zh2, zh3, zh4, rhoh, rhoh2, rhoh3, rhoh4
      real*8 Lam, sf, tt, p5,p6,p7, bi1,bi2,bi3,bi4
      real*8 v,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14, v_p
      real*8 pi1,pi2,pi3,pi4,pi5,pi6,pi7,pi8,pi9,pi10,pi11,pi12,pi13
      real*8 d1,d2,d3,d4,d5,d6,d7
      real*8 L,Lhat,Lhat2,Lhat3,Lhat4
      real*8 s3, s4, s5, s6  

      complex*16 cu, cz

      complex*16 Kxx0,Kyy0,Kzz0,Kxz0, Kxx2,Kyy2,Kzz2,Kxz2
      complex*16 Kxx3,Kyy3,Kzz3,Kxz3, Kxx4,Kyy4,Kzz4,Kxz4
      complex*16 Kxx5,Kyy5,Kzz5,Kxz5, Kxx6,Kyy6,Kzz6,Kxz6
      complex*16 Kxx7,Kyy7,Kzz7,Kxz7
      complex*16 GRxx(0:7), GRyy(0:7), GRzz(0:7), GRxz(0:7), GRzx(0:7)
c
c------ Mathematical constants ------------------------------------
c
      rz = 0.0d0
      ru = 1.0d0
      cu = (0.0d0, 1.0d0)
      cz = (0.0d0, 0.0d0)

      pi = 2*dasin(1.0d0)
      twopi = 2.0d0*pi

      sf = 64.
      tt = 32.
      bi1 = 20160.
      bi2 = 3360.
      bi3 = 840.
      bi4 = 105.
c
c---------- Reading input parameters -----------------------------------
c
      open(unit=70, file='GF.par', status='old', err=101)
      read(70,*) hdummy
      read(70,*) hdummy

      read(70,*, err=101) key
      if((key .eq. 0) .or. (key .eq. 1)) then
         write(*,*) 'GF_sde_tdss ignores KEY=/=2 options'
         write(*,*) 'Even though you have given me KEY=', key, ','
         write(*,*) '   ... I will set KEY=2 and assume'
         write(*,*) '   ... a DISPERSIVE transparent substrate'
         key = 2
      end if
      if (key .eq. 2) then
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

      read(70,*,err=101) nord
      if(nord .lt. 0) then
         write(*,*) 'Warning. You have requested NORD<0.'
         write(*,*) 'This option is used only in GF_k2g'
         write(*,*) 'For GF_sde_tnds, the sign of NORD does not matter'
         write(*,*) 'I will set NORD=|NORD| and continue.'
         nord = -nord
         write(*,*) 'NORD=', nord
         write(*,*) '--------------------------------------------'
      end if

      if(nord .gt. 7) then 
         write(*,*) 'Error. Incorrect parameter: NORD>7'
         write(*,*) 'Correct the error in GF.par and try again'
         stop
      else if (nord .eq. 1) then
         write(*,*) 'The first order correction is identicallyz zero'
         write(*,*) '    (see paper)'
         write(*,*) 'Therefore, there is no difference between NORD=0'
         write(*,*) '    and NORD=1'
         write(*,*) 'I will set NORD=0 and continue'
         write(*,*) 'The result is not affected. Goodby.'
         nord = 0
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

      read(70,*) rdummy, rdummy
      read(70,*) rdummy
      read(70,*) rdummy
      read(70,*) rdummy

      read(70,*) eps_1
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
      z_p   = z2 + z1
      L     = dsqrt(rho*rho + z_p*z_p)
      zh    = z_p/L
      zh2   = zh*zh
      zh3   = zh2*zh
      zh4   = zh2*zh2
      rhoh  = rho/L
      rhoh2 = rhoh*rhoh
      rhoh3 = rhoh2*rhoh
      rhoh4 = rhoh2*rhoh2

      if(nl .gt. 1) then
         dlambda = (lambda_max - lambda_min)/dfloat(nl-1)
      else
         dlambda = rz
      end if
      logscale = .false.
      if(lambda_scale .eq. 'LOG') logscale=.true.
      if(lambda_scale .eq. 'log') logscale=.true.
      if(logscale) then
         write(*,*) 'Using logscale'
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
      open(71, file='GRxx_tsub', status='unknown', err=102)
      open(72, file='GRyy_tsub', status='unknown', err=102)
      open(73, file='GRzz_tsub', status='unknown', err=102)
      open(74, file='GRxz_tsub', status='unknown', err=102)

      open(76, file='Kxx_tsub', status='unknown', err=102)
      open(77, file='Kyy_tsub', status='unknown', err=102)
      open(78, file='Kzz_tsub', status='unknown', err=102)
      open(79, file='Kxz_tsub', status='unknown', err=102)

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
      Lhat  = L*k1
      Lhat2 = Lhat*Lhat
      Lhat3 = Lhat2*Lhat
      Lhat4 = Lhat3*Lhat

      if(key. eq. 2) then
         ! Enter your own formula for eps_2 here
         ! For examle: 
         !  eps_2 = 2.5d0 + 0.001*(400./lambda)**2
         continue
         ! Do not edit past this line
      else
         write(*,*) 'Unexpected error 1; exiting'
         write(*,*) 'Please report this error to code developers'
         stop
      end if 

      eps =  eps_2/eps_1

      eps_m = eps - ru
      eps_p = eps + ru

      s = (ru-eps)/(ru+eps)
      ep2 = eps*eps
      ep3 = ep2*eps
      ep4 = ep3*eps
      ep5 = ep4*eps

      d1 = dsqrt(eps_p)
      d2 = dsqrt(eps)
      d3 = (ru + d1)/(eps + d2*d1)
      d4 = (eps + ru)**2
      d5 = 3.*(d2 + ru)
      d6 = (ru - 3.*(d2 - eps) + 2.*ep2)/d5
      d7 = ru + d2 + 2.*eps*(ru - d2*eps_m) - ep3

      Lam = (eps/eps_m/d1)*dlog(d3)

      Kxx3 =  cu * (eps/d4)*(d6 + Lam)
      Kyy3 = Kxx3
      Kzz3 = -cu * (2./d4)*(d7/d5 + Lam*ep2)
      Kxz3 = cz

      p1 = 2.*ep3 + 5.*ep2 + ru
      p2 = ep2 + 6.*eps + ru
      p3 = eps_m*(2.*ep4+7.*ep3+7.*ep2-eps+ru)
      p4 = 4.*ep5+17.*ep4+22.*ep3-14.*ep2+2.*eps+ru
      p5 = 4.*ep5+19.*ep4+34.*ep3+22.*ep2+14.*eps+3.
      p6 = eps_m*(ep3+5.*ep2+9.*eps+ru)
      p7 = ep4+6.*ep3+18.*ep2+6.*eps+ru

      v = d2
      v2 = v*v
      v3 = v2*v
      v4 = v3*v
      v5 = v4*v
      v6 = v5*v
      v7 = v6*v
      v8 = v7*v
      v9 = v8*v
      v10 = v9*v
      v11 = v10*v
      v12 = v11*v
      v13 = v12*v
      v14 = v13*v
      v_p = v + ru

        pi1 = 4.*v10+4.*v9+4.*v8+8.*v7 -
     $        14.*v6+2.*v5-21.*v4+11.*v3-9.*v2 -
     $        2.*v-2.
        pi2 = v2*(8.*v8+8.*v7+26.*v6+18.*v5 +
     $        15.*v4-33.*v3-v2+2.*v+2.)
        pi3 = 8.*v10+8.*v9+30.*v8+14.*v7 +
     $        29.*v6-19.*v5-3.*v4-18.*v3-18.*v2
     $        -8.*v-8.
        pi4 = 2.*v10+2.*v9+6.*v8+5.*v7+2.*v6 -
     $        10.*v5+3.*v3+3.*v2+v+ru
        pi5 = 2.*v10+2.*v9+11.*v8+15.*v7+7.*v6 -
     $        10.*v5-10.*v4-12.*v3-12.*v2-4.*v-4.
        pi6 = 16.*v14+16.*v13+16.*v12+16.*v11 -
     $        52.*v10-70.*v9-30.*v8-120.*v7+154.*v6-
     $        34.*v5+183.*v4-77.*v3+63.*v2+12.*v+12.
        pi7 = 32.*v14+32.*v13+80.*v12+80.*v11-
     $        24.*v10-8.*v9-246.*v8-150.*v7-177.*v6+
     $        191.*v5-33.*v4-38.*v3-38.*v2-8.*v-8.
        pi8 = eps*(48.*v12+48.*v11+200.*v10+200.*v9
     $        +298.*v8+250.*v7+7.*v6-425.*v5-137.*v4
     $        +10.*v3+10.*v2+8.*v+8.)
        pi9 = eps*(32.*v12+32.*v11+96.*v10+96.*v9+
     $        40.*v8+64.*v7-154.*v6-50.*v5-187.*v4+
     $        37.*v3-75.*v2-18.*v-18.)
        pi10 = 48.*v14+48.*v13+232.*v12+232.*v11+
     $         482.*v10+338.*v9+443.*v8-85.*v7+11.*v6
     $         -382.*v5-382.*v4-344.*v3-344.*v2-96.*v-96
        pi11 = 8.*v14+8.*v13+18.*v12+18.*v11-
     $         14.*v10-11.*v9-73.*v8-50.*v7-43.*v6+
     $         67.*v5-3.*v4-12.*v3-12.*v2-3.*v-3.
        pi12 = 8.*v14+8.*v13+32.*v12+32.*v11+
     $         42.*v10+38.*v9-17.*v8-85.*v7-29.*v6+
     $         18.*v5+18.*v4+16.*v3+16.*v2+4.*v+4.
        pi13 = 8.*v14+8.*v13+46.*v12+46.*v11+
     $         133.*v10+157.*v9+109.*v8-50.*v7-50.*v6-
     $         136.*v5-136.*v4-96.*v3-96.*v2-24.*v-24.

      s3 = ru/eps_p**3
      s4 = s3/eps_p
      s5 = s4*pi*eps_m
      s6 = s3*pi*eps_m
      
c     O(1/r/r/r) and O(1/r) contributions in refl part:

      Kxx0 = -s*(3.*zh2 - 2.)
      Kyy0 = -s 
      Kzz0 = -s*(3.*zh2 - ru)
      Kxz0 = -s*3.*rhoh*zh
 
      Kxx2 =  s*zh*((ru + s)/(ru + zh) - zh)/2. 
      Kyy2 =  s*(s - zh)/2./(ru + zh)
      Kzz2 = -s*((zh2 - s)/2. + ru)
      Kxz2 = -s*(eps*(L - z_p)/eps_p/rho + zh*rhoh/2.)

c     higher orders in h-expansion in refl part:

      Kxx4 = -cu*s6*zh*p1/16.
      Kyy4 =  Kxx4
      Kzz4 = -cu*s6*zh*eps*p2/8.
      Kxz4 =  cu*s6*rhoh*eps*p2/16.

      da = 4.*(pi1/v_p - 15.*eps*Lam)
      db = pi2/v_p + 45.*ep2*Lam
      Kxx5 = cu*s3*(da*zh2 - db*rhoh2)/120.

      db = pi3/v_p + 15.*ep2*Lam
      Kyy5 = cu*s3*(da*zh2 - db*rhoh2)/120.

      da = 2.*(pi4/v_p + 15.*ep2*Lam)
      db = pi5/v_p - 15.*ep3*Lam
      Kzz5 = cu*s3*(da*zh2 - db*rhoh2)/30.

      Kxz5 = -cu*s3*zh*rhoh*da/15./2.

      Kxx6 = cu*s5*zh*(rhoh2*p4/4. - p3*zh2/3.)/sf
      Kyy6 = cu*s5*zh*(rhoh2*p5/4. - p3*zh2/3.)/sf
      Kzz6 = cu*s5*eps*zh*(rhoh2*p7/2. - p6*zh2/3.)/tt
      Kxz6 = cu*s5*eps*rhoh*(p6*zh2 - rhoh2*p7/4.)/sf

      Kxx7 = cu*s4*(8.*zh4*(pi6/v_p+bi4*eps*Lam)+12.*zh2*rhoh2*
     $  (315.*ep2*Lam-pi7/v_p)+rhoh4*(pi8/v_p+525.*ep3*Lam))/bi1
      Kyy7 = cu*s4*(8.*zh4*(pi6/v_p+bi4*eps*Lam)+12.*zh2*rhoh2*
     $  (bi4*ep2*Lam-pi9/v_p)+rhoh4*(pi10/v_p+bi4*ep3*Lam))/bi1
      Kzz7 = cu*s4*(8.*zh4*(pi11/v_p-bi4*ep2*Lam)/3.-8.*zh2*
     $       rhoh2*(bi4*ep3*Lam+pi12/v_p) 
     $      + rhoh4*(pi13/v_p - bi4*ep4*Lam))/bi2
      Kxz7 = cu*s4*(-4.*zh3*rhoh*(pi11/v_p - bi4*ep2*Lam)/3. +
     $  zh*rhoh3*(bi4*ep3*Lam + pi12/v_p))/bi3

      GRxx(0) = Kxx0/Lhat3
      GRyy(0) = Kyy0/Lhat3
      GRzz(0) = Kzz0/Lhat3
      GRxz(0) = Kxz0/Lhat3

      GRxx(1) = GRxx(0)
      GRyy(1) = GRyy(0)
      GRzz(1) = GRzz(0)
      GRxz(1) = GRxz(0)

      GRxx(2) = GRxx(1) + Kxx2/Lhat
      GRyy(2) = GRyy(1) + Kyy2/Lhat
      GRzz(2) = GRzz(1) + Kzz2/Lhat
      GRxz(2) = GRxz(1) + Kxz2/Lhat

      GRxx(3) = GRxx(2) + Kxx3
      GRyy(3) = GRyy(2) + Kyy3
      GRzz(3) = GRzz(2) + Kzz3
      GRxz(3) = GRxz(2) + Kxz3

      GRxx(4) = GRxx(3) + Kxx4*Lhat
      GRyy(4) = GRyy(3) + Kyy4*Lhat
      GRzz(4) = GRzz(3) + Kzz4*Lhat
      GRxz(4) = GRxz(3) + Kxz4*Lhat

      GRxx(5) = GRxx(4) + Kxx5*Lhat2
      GRyy(5) = GRyy(4) + Kyy5*Lhat2
      GRzz(5) = GRzz(4) + Kzz5*Lhat2
      GRxz(5) = GRxz(4) + Kxz5*Lhat2

      GRxx(6) = GRxx(5) + Kxx6*Lhat3
      GRyy(6) = GRyy(5) + Kyy6*Lhat3
      GRzz(6) = GRzz(5) + Kzz6*Lhat3
      GRxz(6) = GRxz(5) + Kxz6*Lhat3

      GRxx(7) = GRxx(6) + Kxx7*Lhat4
      GRyy(7) = GRyy(6) + Kyy7*Lhat4
      GRzz(7) = GRzz(6) + Kzz7*Lhat4
      GRxz(7) = GRxz(6) + Kxz7*Lhat4

      GRzx(0) = -GRxz(0)
      GRzx(1) = -GRxz(1)
      GRzx(2) = -GRxz(2)
      GRzx(3) = -GRxz(3)
      GRzx(4) = -GRxz(4)
      GRzx(5) = -GRxz(5)
      GRzx(6) = -GRxz(6)
      GRzx(7) = -GRxz(7)

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
     $      sngl(dreal(Kxx3)), sngl(dimag(Kxx3)),
     $      sngl(dreal(Kxx4)), sngl(dimag(Kxx4)),
     $      sngl(dreal(Kxx5)), sngl(dimag(Kxx5)),
     $      sngl(dreal(Kxx6)), sngl(dimag(Kxx6)),
     $      sngl(dreal(Kxx7)), sngl(dimag(Kxx7))
         write(77,8,err=102) sngl(lambda),
     $      sngl(dreal(Kyy0)), sngl(dimag(Kyy0)),
     $      sngl(dreal(Kyy2)), sngl(dimag(Kyy2)),
     $      sngl(dreal(Kyy3)), sngl(dimag(Kyy3)),
     $      sngl(dreal(Kyy4)), sngl(dimag(Kyy4)),
     $      sngl(dreal(Kyy5)), sngl(dimag(Kyy5)),
     $      sngl(dreal(Kyy6)), sngl(dimag(Kyy6)),
     $      sngl(dreal(Kyy7)), sngl(dimag(Kyy7))
         write(78,8,err=102) sngl(lambda),
     $      sngl(dreal(Kzz0)), sngl(dimag(Kzz0)),
     $      sngl(dreal(Kzz2)), sngl(dimag(Kzz2)),
     $      sngl(dreal(Kzz3)), sngl(dimag(Kzz3)),
     $      sngl(dreal(Kzz4)), sngl(dimag(Kzz4)),
     $      sngl(dreal(Kzz5)), sngl(dimag(Kzz5)),
     $      sngl(dreal(Kzz6)), sngl(dimag(Kzz6)),
     $      sngl(dreal(Kzz7)), sngl(dimag(Kzz7))
         write(79,8,err=102) sngl(lambda),
     $      sngl(dreal(Kxz0)), sngl(dimag(Kxz0)),
     $      sngl(dreal(Kxz2)), sngl(dimag(Kxz2)),
     $      sngl(dreal(Kxz3)), sngl(dimag(Kxz3)),
     $      sngl(dreal(Kxz4)), sngl(dimag(Kxz4)),
     $      sngl(dreal(Kxz5)), sngl(dimag(Kxz5)),
     $      sngl(dreal(Kxz6)), sngl(dimag(Kxz6)),
     $      sngl(dreal(Kxz7)), sngl(dimag(Kxz7))

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
     $      dreal(Kxx3), dimag(Kxx3),
     $      dreal(Kxx4), dimag(Kxx4),
     $      dreal(Kxx5), dimag(Kxx5),
     $      dreal(Kxx6), dimag(Kxx6),
     $      dreal(Kxx7), dimag(Kxx7)
         write(77,12,err=102) sngl(lambda),
     $      dreal(Kyy0), dimag(Kyy0),
     $      dreal(Kyy2), dimag(Kyy2),
     $      dreal(Kyy3), dimag(Kyy3),
     $      dreal(Kyy4), dimag(Kyy4),
     $      dreal(Kyy3), dimag(Kyy5),
     $      dreal(Kyy6), dimag(Kyy6),
     $      dreal(Kyy7), dimag(Kyy7)
         write(78,12,err=102) sngl(lambda),
     $      dreal(Kzz0), dimag(Kzz0),
     $      dreal(Kzz2), dimag(Kzz2),
     $      dreal(Kzz3), dimag(Kzz3),
     $      dreal(Kzz4), dimag(Kzz4),
     $      dreal(Kzz5), dimag(Kzz5),
     $      dreal(Kzz6), dimag(Kzz6),
     $      dreal(Kzz7), dimag(Kzz7)
         write(79,12,err=102) sngl(lambda),
     $      dreal(Kxz0), dimag(Kxz0),
     $      dreal(Kxz2), dimag(Kxz2),
     $      dreal(Kxz3), dimag(Kxz3),
     $      dreal(Kxz4), dimag(Kxz4),
     $      dreal(Kxz5), dimag(Kxz5),
     $      dreal(Kxz6), dimag(Kxz6),
     $      dreal(Kxz7), dimag(Kxz7)

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

      close(unit=76)
      close(unit=77)
      close(unit=78)
      close(unit=79)

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
 8    format(15G16.7E3)
 10   format(G16.7E3,  2G23.14E3)
 12   format(G16.7E3, 14G23.14E3)

      end
