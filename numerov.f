!This program computes the probability density for the Quantum Harmonic
!Oscillator using Numerov's Method

       implicit none
       real*16 x, dx, w, k, m, h, psi, E, ynp1, area
       real*16 fn, gn, gamma, y0, y1, f0, f1, V, mp
       integer i, j, n, np, l

       external mp !Verifies the maximum number of points to use for each energy level
       external fn
       external gn
       external ynp1 !New point
       external E !Discretized energy

       !Parameters
       Parameter (n = 30) !Maximum level of energy
       dimension area(0:n)

       k = 1.d0 !k = w²m
       m = 1.d0 !mass
       w = 1.d0 !w² = k/m
       h = 1.d0 !hbar
       dx = 1.d-3 !Step in x
       gamma = (2.d0*m)/(h**2) !Normalization factor

       do i = 0, n
       x = 2.d0*dx
       np = (mp(i))/dx !Number of points to use in each energy level
       area(i) = 0.d0
       !Beginning of the program: I see if the number n is odd or even to start the guesses
       if(mod(i, 2).EQ.0)then
          write(*,*)'Even', i
          y0 = E(i) !Random number => removed in the normalization
          f0 = fn(x, dx, gn(x, E(i), gamma))
          f1 = fn((x + dx), dx, gn(x, E(i), gamma))
          y1 = ((12.d0 - 10.d0*f0)*y0)/(2.d0*f1)
          write(700 + i, *) 0.d0, (y0**2)
          write(700 + i, *) dx, (y1**2)
          area(i) = area(i) + ((y0**2)*dx) + ((y1**2)*dx)
       end if
       if(mod(i, 2).NE.0)then
          write(*,*)'Odd', i
          y0 = 0.d0
          y1 = E(i) !Random number => removed in the normalization
          f0 = fn(x, dx, gn(x, E(i), gamma))
          f1 = fn((x + dx), dx, gn(x, E(i), gamma))
          write(700 + i, *) 0.d0, (y0**2)
          write(700 + i, *) dx, (y1**2)
          area(i) = area(i) + ((y0**2)*dx) + ((y1**2)*dx)
       end if  
       !Computing the points of the functions
       do j = 0, np
          x = x + dx
          psi = ynp1(x, y1, y0, f1, fn(x, dx, gn(x, E(i), gamma)), f0)
          y0 = y1
          f0 = f1
          f1 = fn(x, dx, gn(x, E(i), gamma))
          y1 = psi
          write(700 + i, *) x, (psi**2)
          area(i) = area(i) + ((psi**2)*dx)
       end do
       !Values of the area (already multiplied by 2)
       write(666,*) i, (area(i)*2.d0) !File 666 has the values of area
       end do

       !Renormalizing the eigenfunctions
       do i = 0, n
           rewind(700+i)
       end do

       do i = 0, n
          np = mp(i)/dx
          do j = 0, np
             read(700+i,*) x, psi
             write(800+i,*) (-1.d0*x), psi/(area(i)*2.d0)
          end do
       end do

       do i = 0, n
           rewind(700+i)
       end do

       do i = 0, n
          np = mp(i)/dx
          do j = 0, np
             read(700+i,*) x, psi
             write(800+i,*) x, psi/(area(i)*2.d0)
          end do
       end do

       stop
       end      

       !Eigenvalues of the energy
       function E(n)
             real*16 E
             integer n
             E = (0.5d0 + n)
       return
       end

       !Maximum number of points to use for each n
       function mp(n)
             real*16 mp
             integer n
             mp = sqrt(2.d0*(0.5d0+n)) + 1.5d0
       return
       end

       !Parameter gn
       function gn (x, E, gamma)
             real*16 gn, V, E, gamma, k, x
             k = 1.d0
             V = (k*(x**2))/2.d0
             gn = gamma*(E - V)
       return
       end

       !Parameter fn
       function fn(x, dx, gn)
            real*16 fn, x, dx, gn
            integer n
            fn = 1.d0 + (gn*(dx**2))/12.d0
       return
       end

       !Next points in the eigenfunctions
       function ynp1(x, yn, ynm1, fn, fnp1, fnm1)
            real*16 ynp1, x, yn, ynm1, fn, fnp1, fnm1
            ynp1 = ((12.d0 - 10.d0*fn)*yn - (fnm1)*(ynm1))/(fnp1)
       return
       end
