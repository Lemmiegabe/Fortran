c DYLAN: The purpose of this program is to see if a disk much lighter
c        than a cluster whose angular momentum vector is misaligned 
c        from that of the cluster will partialy align with the 
c        disk. (Doing Your Lovely Amazing Numerics)

      implicit none 
c All the variables used in the program
      integer n, l, nh, j, st, stsh, k 
      parameter (n = 100, l = 1000000)
      real*8 v(3,n), x(3,n), dt, t, R2, a(3,n), r, c, phi, pi
      real*8 m(n), w, da, dx(3), f1, f2, xh, yh, xf, yf, zf, vxf
      real f(1000000)
      parameter (pi = 3.14159265359)
      common /var /v, x, dt, t, m
      common /coeff /R2, a, c, phi
      common /string /r, da, dx, w, f
      common /mem /st, nh
      

      real*8 sm_device
      ! Max simulation time
      real*8 T_max
      parameter (T_max = 10.0)

    
c If the device is not null device then print you can not open          
      if (sm_device('X11 -bg BLACK -g 1000x1000') .lt. 0) then
         print *,'Can''t open output device.'
         stop
      endif
 
      call sm_graphics
c Equivalent to Define name value. Makes SM understand TeX strings.
      call sm_defvar("TeX_strings","1")
      call sm_ptype(1.,1)
            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      stsh = 10000

      call init

      do st = 0, 1000000000

         if ((st / stsh) * stsh.eq.st) then
            call show
         endif

         call step
          
c         if(t.gt.T_max) then 
c            exit
c         endif
         

      enddo

      STOP

      END
          

            
 
   
     

      
      
       




      subroutine init
      implicit none 
c   init: The intial conditions for the many body sim 
      integer n, l, st, nh, j, k
      parameter (n = 100)
      parameter(l = 1000000)
      real*8 v(3,n),x(3,n),dt,t,R2,a(3,n),r,c,phi,pi
      real*8  m(n),w,da,dx(3),va,b,d(3),e(3),lz,lx
      real   f(1000000)
      parameter (pi = 3.14159265359)
      common /var /v, x, dt, t, m
      common /coeff /R2, a, c, phi
      common /string /r, da, dx, w, f
      common /mem /st, nh 

      write (*,*) 'dt=?'
      read (*,*)  dt
     
      
      t = 0
      
      nh = 0
      
      m(1) = 2000 ! Mass of the black hole
      

      do  k = 2, n

         r = 2 * (rand(0)) ** (1 / 3) ! random radius within cluster_radius
         c=(-1 + 2 * rand(0))     ! random cos 
         phi=2 * pi * rand(0)     ! random angle 

        
         
       

         m(k) = 1.0 ! Mass of cluster particles 
         x(1,k) = r * sqrt(1 - c ** 2) * cos(phi)
         x(2,k) = r * sqrt(1 - c ** 2) * sin(phi)
         x(3,k) = r * c
         r = sqrt(x(1,k) ** 2 + x(2,k) ** 2 + x(3,k) ** 2)
         c = (-1 + 2 * rand(0))
         phi = 2 * pi * rand(0)
         va = (2 * (m(1) / r) ** (1 / 2)) * (rand(0)) ** (1 / 3)
         v(1,k) = va * sqrt(1 - c ** 2) * cos(phi)
         v(2,k) = va * sqrt(1 - c ** 2) * sin(phi)
         v(3,k) = va * c
          
        
      enddo
      go to 99
      do k = n / 2 + 1, n

         r = sqrt(rand(0))  ! Random radius within disk_radius
         phi = 2 * pi * rand(0)               ! Random angle
         x(1,k) = r * cos(phi)
         x(2,k) = r * sin(phi)
         x(3,k) = 0.0
         r = sqrt(x(1,k) ** 2 + x(2,k) ** 2)
         va = (2 * (m(1) / r) ** (1 / 2)) * (rand(0)) ** (1 / 3)
         v(1,k) = va * cos(phi)
         v(2,k) = va * sin(phi)
         v(3,k) = 0.0
         m(k) = 1.0  ! Mass of disk particles
      enddo
 99   continue

      do j = 1, 3
         x(j,1) = 0
         v(j,1) = 0
      enddo
   
      

      b = 0
      do k = 1, n
         b = b + m(k)
      enddo

      do j = 1, 3
         e(j) = 0
      enddo
      do j = 1, 3
         d(j) = 0
      enddo
 
      do k = 1, n
         do j = 1, 3
            e(j) = e(j) + m(k) * x(j,k)
            d(j) = d(j) + m(k) * v(j,k) 
         enddo
      enddo

      do 2 k = 1, n
      do 2 j = 1, 3
         v(j,k) = v(j,k) - d(j) / b
         x(j,k) = x(j,k) - e(j) / b

           
 2       continue

         do k = 2, n
            lz = (x(1, k) * v(2, K) - x(2, k) * v(1, k))! z-angular momentum cluster    
            if (lz.le.0) then
               v(1,K) = -v(1,k)
               v(2,k) = -v(2,k)
               v(3,k) = -v(3,k)
            endif
         enddo
         go to 98
         do k = n / 2 + 1, n
            lx = x(2,k) * v(3,k) - x(3,k) * v(2,k)! x-angular momentum disk 
            if (lx.le.0) then 
               v(1,K) = -v(1,k)
               v(2,k) = -v(2,k)
               v(3,k) = -v(3,k)
            endif
         enddo
 98      continue

    
      

      return
     
     
      end

cccccccccccccccccccccccccccc step cccccccccccccccccccccccccccccccccc

      
     
      subroutine step
      implicit none  

      integer n, l, st, sth, nh
      parameter (n=100)
      parameter(l=1000000)
      real*8 v(3,n), x(3,n), dt, t, R2, a(3,n), r, c, phi, pi
      real*8 m(n), w, da, dx(3), b, xf, yf, zf
      real f(1000000), g(l), h(l)
      parameter (pi = 3.14159265359)
      common /var /v, x, dt, t, m
      common /coeff /R2, a, phi, c
      common /string /r, w, da, dx, f
      common /mem /st, nh
      integer j, k, k1, k2
     
      t=t+dt

      do 1 k =  1, n
      do 1 j = 1, 3
         
         x(j,k) = x(j,k) + dt * v(j,k)
         
1     continue


    
      do k = 1, n
         do j = 1, 3
            a(j,k) = 0
         enddo
      enddo
     
      
      do 2 k1 = 1, n - 1
      do 2 k2 = k1 + 1, n
        

         do j = 1, 3
            dx(j) = x(j,k1) - x(j,k2)
         enddo
         b = .1
         R2 = (dx(1)) ** 2 + (dx(2)) ** 2 + (dx(3)) ** 2
         w = 1/((R2 + b ** 2) * sqrt(R2 + b ** 2))
         
         do j = 1, 3
            da = -w * dx(j)
            a(j,k1) = a(j,k1) + m(k2) * da
            a(j,k2) = a(j,k2) - m(k1) * da
         enddo 
         
2     continue
     
     
      do 3 k = 1, n
      do 3 j = 1, 3
         
         v(j,k) = v(j,k) + dt * a(j,k)

3     continue
     
   
    
      
      
      return
      end
       


cccccccccccccccccccc show cccccccccccccccccccccccccccccc

      subroutine show
      implicit none 
      integer n 
      parameter (n = 100)
      integer l, st, sth, nh
      parameter(l = 1000000)
      real*8 v(3,n), x(3,n), dt, t, R2, a(3,n), r, c, phi, P(3), E, pi
      real*8 m(n), w, da, dx(3), Ke, U, Lx, Ly, Lz, Xcm(3), b, d(3), q, s, z
      real f(l), g(l), h(l), i(l), f2, g2, h1, h2, i2, xh, yh, o(l)
      real xf, yf, zf, vf, rs, i1, o1, o2
      parameter (pi = 3.14159265359)
      common /var /v, x, dt, t, m 
      common /coeff /R2, a, c, phi
      common /string /r, w, da, dx, f
      common /mem /st, nh
      integer k, j, k1, k2
     
      
      ! graph 
      call sm_limits(-3.,3.,-3.,3.)
      call sm_ctype('white')
      call sm_box(1,2,0,0)
      call sm_lweight(2.)
      call sm_ctype('yellow')
      call sm_gflush 
      ! x-z plane cluster
      s = .3
      do k = 1, n
         z = -1.5 + s * x(1,k)
         q = -1.5 + s * x(3,k)
         call sm_relocate(real(z),real(q))
         call sm_dot          
      enddo
      ! x-y plane cluster
      do k =1, n
         z = -1.5 + s * x(1,k)
         q = 1.5 + s * x(2,k)
         call sm_relocate(real(z),real(q))
         call sm_dot     
      enddo
      ! z-y plane cluster
      do k = 1, n
         z = 1.5 + s * x(3,k)
         q = 1.5 + s * x(2,k)
         call sm_relocate(real(z),real(q))
         call sm_dot   
      enddo
      
      go to 100
      ! x-z plane disk
      call sm_ctype('blue')
      do k = n / 2 + 1, n
         z = -1.5 + s * x(1,k)
         q = -1.5 + s * x(3,k)
         call sm_relocate(real(z),real(q))
         call sm_dot
      enddo 
      ! x-y plane disk
      do k = n / 2 + 1, n
         z = -1.5 + s * x(1,k)
         q = 1.5 + s * x(2,k)
         call sm_relocate(real(z),real(q))
         call sm_dot
      enddo
      ! z-y plane disk
      do k = n / 2 + 1, n
         z = 1.5 + s * x(3,k)
         q = 1.5 + s * x(2,k)
         call sm_relocate(real(z),real(q))
         call sm_dot
      enddo
 100  continue 
      call sm_ctype('white')
      call sm_gflush
      call sm_relocate(-3.,0.)
      call sm_draw(3.,0.)
      call sm_relocate(0.,-3.)
      call sm_draw(0.,3.)
      call sm_relocate(0.,-.75)
      call sm_draw(3.,-.75)
      call sm_relocate(0.,-1.5)
      call sm_draw(3.,-1.5)
      call sm_relocate(0.,-2.25)
      call sm_draw(3.,-2.25)
      call sm_gflush
     
      
   
     
      
       
      U = 0
      do 1 k1 = 1, n - 1 
         do 1 k2 = k1 + 1, n
            do j = 1, 3
               dx(j) = x(j,k1) - x(j,k2)
            enddo

            R2 = dx(1) ** 2 + dx(2) ** 2 + dx(3) ** 2
            b = .1


           
           
          
               
            U = U - m(k1) * m(k2) / sqrt(R2 + b ** 2)
                      
          
1          continue 

        
           b = 0
           do k = 1, n
              b = b + m(k)
           enddo

           do j = 1, 3
              d(j) = 0
           enddo
 
           do k = 1, n
              do j = 1, 3
                 d(j) = d(j) + m(k) * x(j,k)
              enddo
           enddo   
      
           do j = 1, 3 
              Xcm(j) = 0
              P(j) = 0
           enddo
     
           do j = 1, 3
              Xcm(j) = Xcm(j) + d(j) / b
           enddo
           Lx = 0
           Ly = 0
           Lz = 0
           Ke = 0
           do k = 1, n
   
              Lx = Lx + m(k) * ((x(2,k) * v(3,k)) - (x(3,k) * v(2,k)))
              Ly = Ly + m(k) * ((x(3,k) * v(1,k)) - (x(1,k) * v(3,k)))
              Lz = Lz + m(k) * ((x(1,k) * v(2,k)) - (x(2,k) * v(1,k)))
           
              do j = 1, 3
                 P(j) = P(j) + m(k) * v(j,k)
                 Ke = Ke + m(k) * (v(j,k)) ** 2

              enddo
           enddo
           
      nh = nh + 1


     
      Ke = Ke / 2
      E = Ke + U
      xf = 0.
      yf = 0.
      zf = 0.
      vf = 0.
      rs = 0.
         
      f(nh) = -E
      g(nh) = Lx ** 2 + Ly ** 2 + Lz ** 2
      do k = 2, n
         xf = xf + x(1,k) ** 2
         yf = yf + x(2,k) ** 2
         zf = zf + x(3,k) ** 2
      enddo
      h(nh) = ((xf + yf - 2 * zf) / (n - 1)) / ((xf + yf + zf) / (n - 1))
      do k = 2, n
         vf = vf + (sqrt(v(1,k) ** 2 + v(2,k) ** 2 + v(3,k) ** 2)) ** 4
      enddo
       do k = 2, n
         rs = rs +(sqrt(x(1,k) ** 2 + x(2,k) ** 2 + x(3,k) ** 2)) ** 2
      enddo
      i(nh) = (vf / (n - 1))
      o(nh) = (rs / (n - 1))
      
      f2 = 1
      g2 = 1
      h1 = 1
      h2 = 1
      i2 = 1
      i1 = 1
      o1 = 1
      o2 = 1
      do k = 1, nh
         if (f(k).gt.f2) then
            f2 = f(k)
         endif
      enddo
      do k = 1, nh
         if (o(k).gt.o2) then
            o2 = o(k)
         endif
      enddo
       do k = 1, nh
         if (g(k).gt.g2) then
            g2 = g(k)
         endif
      enddo
       do k = 1, nh
         if (h(k).gt.h2) then
            h2 = h(k)
            if (h(k).lt.h1) then
               h1 = h(k)
            endif
         endif
      enddo
       do k = 1, nh
         if (i(k).gt.i2) then
            i2 = i(k)
         endif
      enddo
      call sm_ctype('yellow')
      call sm_relocate(0.,-.375)
      do k = 2, nh
         xh = 3. * (k - 1) / (nh - 1)
         yh = -.375 * f(k) / f2
         call sm_draw(xh,yh)
      enddo
      call sm_relocate(0.,-1.125)
      do k = 2, nh
         xh = 3. * (k - 1) / (nh - 1)
         yh = -1.125 * g(k) / g2
         call sm_draw(xh,yh)
      enddo
      call sm_relocate(0.,-1.875)
      call sm_ctype('white')
      do k = 2, nh
         xh = 3. * (k - 1) / (nh - 1)
         yh =-1.5 + .65 * (h(k) - h1) / (h2) 
         call sm_draw(xh,yh)
      enddo
      call sm_ctype('red')
      call sm_relocate(0.,-1.875)
       do k = 2, nh
         xh = 3. * (k - 1) / (nh - 1)
         yh =-1.5 + .005 * (o(k) - o1) / (o2 - o1)
         
         call sm_draw(xh,yh)
      enddo
      call sm_ctype('blue')
      call sm_relocate(0.,-2.625)
      do k = 2, nh
         xh = 3. * (k - 1) / (nh - 1)
         yh = -2.25 + .75 * ((i(k) - i1) * (-3. + 2.25) / (i2 - i1))
         call sm_draw(xh,yh)
      enddo
      
      call sm_gflush
      call sm_erase

      
      write(*,*) 't, E           =', t, E
      write(*,*) 'Xc, Yc, Zc     =', Xcm(1), Xcm(2), Xcm(3)
      write(*,*) 'Px, Py, Pz     =', P(1), P(2), P(3)
      write(*,*) 'Lx, Ly, Lz     =', Lx, Ly, Lz

      return
      end

      