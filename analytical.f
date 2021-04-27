!This program computes the probability of finding a particle for the
!Classical Harmonic oscillator (analitically), together with the
!density probability of the Quantum Harmonic Oscillator in order to
!compare each order and test the Correspondence Principle      

        implicit none
        real*16 Pi, r, A, t, TT, ex, exx, psi0, psi1, P, dx, x, y, yy
        real*16 Z, psi, Prob, H1, Hn, fat, alpha, alphan, areaq, AA
        real*16 xx, areac
        real*16 sin, cos, ABS, exp
        integer m, k, l, ii, NP, EM

        dimension areac(0:30)
        dimension areaq(0:30)
        dimension Z(0:90) !Hermite Polynomaisl
        external H1 !Hermite 1
        external Hn !Hermite n
        external fat !Factorial
        dimension psi(0:30) !Wave function
        dimension Prob(0:30) !Quantum Probability

        !Parameters
        Pi = acos(-1.d0) !Yeah, FORTRAN <3
        NP = 1000 !Number of points
        EM = 30 !Maximum Energy
        r = (1.d0/sqrt(Pi)) !Quantum normalization factor
        AA = sqrt(1.d0 + (2.d0*EM)) + 2.d0 !Classical amplitude
        t = Pi/NP !Parameter of variation in time

        !Resetting everything
        do l = 0, 30
               psi(l) = 0.d0
               Prob(l) = 0.d0
               areac(l) = 0.d0
               areaq(l) = 0.d0
        end do
        do ii = 0, 90
               Z(ii) = 0.d0
        end do

        !Beggining of the program
        do m = 0, 1000
             x = AA*cos(m*t) !Particle position
             !Classical probability
             if (m.LT.1000) then
               do k = 0, EM
                  A = sqrt(1.d0 + (2.d0*k))
                  xx = A*(cos(m*t)+cos(m*t+t))*0.5d0 !Particle position
                  dx = (-1.d0)*A*sin(m*t) !Particle velocity
                  P = (ABS((1.d0)/dx))/2.d0 !Classical probability
                  write(600+k,*) xx, P
                  areac(k) = areac(k) + P/NP
                end do
             endif
             !Quantum density probability 
             ex = (exp(((-1.d0)*(x**2))/2.d0)) !Exponencitial adjust
             alpha = r
             !Hermite Polynomials computations
             !n = 0
             Z(0) = 1.d0
             psi0 = Z(0)*ex
             Prob(0) = ((psi0)**2)*alpha
             write(200,*) x, Prob(0)
             areaq(0) = areaq(0) + Prob(0)/NP
             !n = 1
             Z(1) = H1(x)
             psi1 = Z(1)*ex
             Prob(1) = ((psi1)**2)*(alpha)/(2.d0)
             write(201,*) x, Prob(1)
             areaq(1) = areaq(1) + Prob(1)/NP
             !Other ns
             do k = 2, EM
                alphan = r/((2**(k))*fat(k))
                exx = exp(((-1.d0)*(x**2))/2.d0)
                Z(k) = Hn(k, x, Z(k-2), Z(k-1))
                psi(k) = Z(k)*exx
                Prob(k) = ((psi(k))**2)*alphan
                write(200+k,*) x, Prob(k) !QHO: analytic result
                areaq(k) = areaq(k) + Prob(k)/NP
             end do
        end do 

        do m = 0, EM
            write(19,*) areaq(m)
            write(18,*) areac(m)
        end do
 
        stop
        end

        !Hermite Polynomials
        !First not constant
        function H1(x) 
                real*16 x, H1
                H1 = 2.d0*x
        return
        end

        !Recorrence expression
        function Hn(n, x, Hn_2, Hn_1)
                real*16 x, Hn_1, Hn_2, Hn
                integer n
                Hn = 2.d0*x*Hn_1 - 2.d0*(n-1)*Hn_2
        return
        end

        !Factorial - to normalization
        function fat(k)
                real*16 fat
                integer k, i
                fat = 1.d0
                if((k.EQ.0).AND.(k.EQ.1))then
                        return
                end if
                do i = 1, k
                        fat = fat*i
                end do
         return
         end
