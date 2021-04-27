!This program computes eigenvalues for the energy of the Quantum Harmonic
!Oscillator using Numerov's Method      

       implicit none
       real*16 x, dx, w, k, m, h, psi, E, ynp1, a, b, c, d, g, f
       real*16 fn, gn, gamma, y0, y1, f0, f1, V, mp, dE
       integer i, j, n, np, l, t, mit

       external fn
       external gn
       external ynp1 !New point

       !Parameters
       Parameter (mit = 3100) !Maximum number of iterations
       Parameter (np = 10000) !Maximum number of points to use
       Dimension a(0:mit)
       Dimension b(0:mit)

       k = 1.d0 !k = w²m
       m = 1.d0 !mass
       w = 1.d0 !w² = k/m
       h = 1.d0 !hbar
       dx = 1.d-3 !Step in x
       dE = 0.1d-1 !Step for eigenvalues of energy
       gamma = (2.d0*m)/(h**2) !Normalization factor

       !Even numbers
       E = 0.d0
       do t = 0, mit
          E = E + dE
          x = 2.d0*dx
          y0 = E !Random number => removed in the normalization
          f0 = fn(x, dx, gn(x, E, gamma))
          f1 = fn((x + dx), dx, gn(x, E, gamma))
          y1 = ((12.d0 - 10.d0*f0)*y0)/(2.d0*f1)
          !Eigenfunctions
          do j = 0, np
             x = x + dx
             psi = ynp1(x, y1, y0, f1, fn(x, dx, gn(x, E, gamma)), f0)
             y0 = y1
             f0 = f1
             f1 = fn(x, dx, gn(x, E, gamma))
             y1 = psi
             if(j.EQ.np)then
                b(t) = psi
             end if
          end do
       end do

       !Odd numbers
       E = 0.d0
       do t = 0, mit
          E = E + dE
          x = 2.d0*dx
          y0 = 0.d0
          y1 = E !Random number => removed in the normalization
          f0 = fn(x, dx, gn(x, E, gamma))
          f1 = fn((x+dx), dx, gn(x, E, gamma))
          !Eigenfunctions
          do j = 0, np
             x = x + dx
             psi = ynp1(x, y1, y0, f1, fn(x, dx, gn(x, E, gamma)), f0)
             y0 = y1
             f0 = f1
             f1 = fn(x, dx, gn(x, E, gamma))
             y1 = psi
             if(j.EQ.np)then
                a(t) = psi
             end if
          end do
       end do

       !Computing the eingenvalues of the energy
       do t = 1, mit
          c = a(t)
          d = a(t - 1)
          g = b(t)
          f = b(t - 1)
          if(((c.GE.0).AND.(d.LE.0)).OR.((c.LE.0).AND.(d.GE.0)))then
             open(42, file = 'energies/energies.dat')
             write(42,*) t*dE
          end if
          if(((g.GE.0).AND.(f.LE.0)).OR.((g.LE.0).AND.(f.GE.0)))then
             open(42, file = 'energies/energies.dat')
             write(42,*) t*dE
          end if
       end do

       stop
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
