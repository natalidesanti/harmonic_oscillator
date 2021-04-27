!This program computes the Hermite polynomials until order 5
!The polynomials are saved in files whith the following correspondence:
! [10, 11, 12, 13, 14, 15] = [0, 1, 2, 3, 4, 5]	
	
	external H1

	do i = -4000, 4000
	   x = i*1.d-3
	   write(10,*) x, 1.d0
	   write(11,*) x, H1(x)
	end do

	rewind(10)
	rewind(11)

	do j = 0, 3
	   do k = 1, 8000
	      read(10+j,*) x, H
	      read(11+j,*) xx, HH
	      Hn = 2.d0*x*HH - 2.d0*(j+2)*H
	      write(12+j,*) x, Hn
	   end do
	rewind(10+j)
	rewind(11+j)
	rewind(12+j)
	end do
	
	stop
	end

	function H1 (x)
		H1 = 2.d0*x
	return
	end
