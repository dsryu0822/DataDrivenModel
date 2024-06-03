      program lyapunov_lorenz
      implicit double precision (a-h,o-z)
      parameter ( N = 3)
      dimension iseed(4)
      dimension x(N),y(N),dlambda(N),stretch(N),dsort(N)
      double precision lambda_old(N),lambda_new(N)
      dimension dx(N,N),v(N,N),u(N,N)
      external flush

      common/parameters/sigma,r,b

      open(18,file='d.d')
c      open(16,file='traj.dat')

   
      sigma =10.d0              !lorenz eq's parameter
      b = 8.0d0/3.d0

      ipre=100/h
      max_iterate=500
      h=0.01d0
      h2=0.01d0

      pi2 = 8.d0*datan(1.d0)
             
      iseed(1) = 18
      iseed(2) = 68
      iseed(3) = 22
      iseed(4) = 35

      DO Nvar =1,1
         r= 10.d0+0.01*dble(Nvar-1) ! lorenz parameter
c         r=15.d0
         do i = 1, N
            x(i) = ran(iseed)
            do j = 1, N
               dx(i,j) = ran(iseed)
            enddo
         enddo
              
C-------Initial transients---------------------

         Do iterate = 1, ipre

            call gram_schmidt(dx,v,u)
c             v12 = 0.d0
c             v13 = 0.d0
c             v23 = 0.d0
c         do j = 1, N
c             v12 = v12 + v(1,j)*v(2,j)
c             v13 = v13 + v(1,j)*v(3,j)
c             v23 = v23 + v(2,j)*v(3,j)
c         enddo
c         print *, 'v12,v13,v23=', v12,v13,v23
c         pause


            Do i = 1, N
               do j = 1, N
                  y(j) = u(i,j)
               enddo
               call rk4_tangent(h2,x,y)
               do j = 1, N
                  dx(i,j) = y(j)
               enddo
            Enddo
            call rk4_map(h2,x)
         Enddo  
c      pause
C----------Accumulate exponent----------------------
         do i = 1, N
            lambda_old(i) = 0.d0
            lambda_new(i) = 0.d0
            stretch(i) = 0.d0
         enddo
         DO iterate = 1, max_iterate
c            if(mod(iterate,1) .eq. 0.d0)then
c               write(16,160) iterate*h,x(1),x(2),x(3)
c 160           format(4(1x,e20.12))
c            endif
            call gram_schmidt(dx,v,u)            
            time = iterate*h
            dlambda_max = 0.d0

            Do i = 1, N
               v_modu = 0.d0
               do j = 1, N
                  v_modu = v_modu + v(i,j)**2
               enddo
               v_modu = dsqrt(v_modu)
               stretch(i) = stretch(i) + dlog(v_modu)
               lambda_new(i) = stretch(i)/time
               dlambda(i) = dabs(lambda_new(i) - lambda_old(i))
               if(dlambda(i) .ge. dlambda_max)then
                  dlambda_max = dlambda(i)
               endif
               lambda_old(i) = lambda_new(i)
            Enddo
            Do i = 1, N
               do j = 1, N
                  y(j) = u(i,j)
               enddo
               call rk4_tangent(h,x,y)
               do j = 1, N
                  dx(i,j) = y(j)
               enddo
            Enddo
            call rk4_map(h,x)
         ENDDO

         write(18,180) r, (lambda_new(i),i=1,N)
         call flush(18)
 180     format(4(1x,e20.12))
c         print*,'lorenz is','r',' ',(lambda_new(i),i=1,N)

      ENDDO
      
      stop
      end

c******************************************************************
c
c******************************************************************
      subroutine rk4_tangent(h,z,y)
      implicit double precision (a-h,o-z)
      parameter(N = 3)
      dimension z(N), zt(N), dzdt(N), dzt(N), dzm(N)
      dimension y(N), dydt(N), yt(N), dyt(N), dym(N)

      hh = h*0.5d0
      h6 = h/6.d0
         call derivs_tangent(z,y,dydt)
         call derivs_map(z,dzdt)
      do i = 1, N
         yt(i) = y(i) + hh*dydt(i)
         zt(i) = z(i) + hh*dzdt(i)
      enddo
         call derivs_tangent(zt,yt,dyt)
         call derivs_map(zt,dzt)
      do i = 1, N
         yt(i) = y(i) + hh*dyt(i)
         zt(i) = z(i) + hh*dzt(i)
      enddo
         call derivs_tangent(zt,yt,dym)
         call derivs_map(zt,dzm)
      do i = 1, N
         yt(i) = y(i) + h*dym(i)
         dym(i) = dyt(i) + dym(i)
         zt(i) = z(i) + h*dzm(i)
      enddo
         call derivs_tangent(zt,yt,dyt)
      do i = 1, N
         y(i) = y(i) + h6*(dydt(i) + dyt(i) + 2.d0*dym(i))
      enddo

      return
      end

c*********************************************************************
c x --> f(x)
c********************************************************************
      subroutine rk4_map(h,z)
      implicit double precision (a-h,o-z)
      parameter(N = 3)
      dimension z(N), dzdt(N)
      dimension zt(N), dzt(N), dzm(N)

      hh = h*0.5d0
      h6 = h/6.d0
         call derivs_map(z,dzdt)
      do i = 1, N
         zt(i) = z(i) + hh*dzdt(i)
      enddo
         call derivs_map(zt,dzt)
      do i = 1, N
         zt(i) = z(i) + hh*dzt(i)
      enddo
         call derivs_map(zt,dzm)
      do i = 1, N
         zt(i) = z(i) + h*dzm(i)
         dzm(i) = dzt(i) + dzm(i)
      enddo
         call derivs_map(zt,dzt)
      do i = 1, N
         z(i) = z(i) + h6*(dzdt(i) + dzt(i) + 2.d0*dzm(i)) 
      enddo

      return
      end


c**********************************************************
C: x(N) = velocity field; y(t) = tangent velocity field   *
c**********************************************************

      subroutine derivs_tangent(x,y,dydt)
      implicit double precision (a-h,o-z)
      parameter (N = 3)
      dimension x(N),y(N),dydt(N),df(N,N)
     
      common/parameters/sigma,r,b

      dydt(1)=sigma*(y(2)-y(1))
      dydt(2)=(r-x(3))*y(1)-y(2)-x(1)*y(3)
      dydt(3)=x(2)*y(1)+x(1)*y(2)-b*y(3)

      return
      end

c******************************************************************
c x(N) = velocity field;                                          *
c******************************************************************
      subroutine derivs_map(x,dxdt)
      implicit double precision (a-h,o-z)
      parameter (N = 3)
      dimension x(N),dxdt(N)

      common/parameters/sigma,r,b

      dxdt(1)=sigma*(x(2)-x(1))
      dxdt(2)=r*x(1)-x(2)-x(1)*x(3)
      dxdt(3)=-b*x(3)+x(1)*x(2)

      return
      end
      
c*******************************************************************
      subroutine gram_schmidt(dx,v,u)
      implicit double precision (a-h,o-z)
      parameter (N = 3)
      dimension dx(N,N),v(N,N),u(N,N)
      dimension product(N-1)
 
      do j = 1, N
         v(1,j) = dx(1,j)
      enddo
      v_modu1 = 0.d0
      do j = 1, N
         v_modu1 = v_modu1 + v(1,j)**2
      enddo
      v_modu1 = dsqrt(v_modu1)
      do j = 1, N
         u(1,j) = v(1,j)/v_modu1
      enddo
      DO i = 2, N
         do j = 1, N
            v(i,j) = dx(i,j)
         enddo
         do ip = 1, i-1
            product(ip) = 0.d0
            do j = 1, N
               product(ip) = product(ip) + dx(i,j)*u(ip,j)
            enddo
         enddo
         do ip = 1, i - 1
            do j = 1, N
               v(i,j) = v(i,j) - product(ip)*u(ip,j)
            enddo
         enddo
         v_modu = 0.d0
         do j = 1, N
            v_modu = v_modu + v(i,j)**2
         enddo
         v_modu = dsqrt(v_modu)
         do j = 1, N
            u(i,j) = v(i,j)/v_modu
         enddo
      ENDDO

      return
      end

c*********************************************************************
      double precision function ran( iseed )
      integer            iseed( 4 )
      integer            m1, m2, m3, m4
      parameter          ( m1 = 494, m2 = 322, m3 = 2508, m4 = 2549 )
      double precision   one
      parameter          ( one = 1.0d+0 )
      integer            ipw2
      double precision   r
      parameter          ( ipw2 = 4096, r = one / ipw2 )
      integer            it1, it2, it3, it4
      intrinsic          dble, mod
      it4 = iseed( 4 )*m4
      it3 = it4 / ipw2
      it4 = it4 - ipw2*it3
      it3 = it3 + iseed( 3 )*m4 + iseed( 4 )*m3
      it2 = it3 / ipw2
      it3 = it3 - ipw2*it2
      it2 = it2 + iseed( 2 )*m4 + iseed( 3 )*m3 + iseed( 4 )*m2
      it1 = it2 / ipw2
      it2 = it2 - ipw2*it1
      it1 = it1 + iseed( 1 )*m4 + iseed( 2 )*m3 + iseed( 3 )*m2 +
     &  iseed( 4 )*m1
      it1 = mod( it1, ipw2 )
      iseed( 1 ) = it1
      iseed( 2 ) = it2
      iseed( 3 ) = it3
      iseed( 4 ) = it4
      ran = r*( dble( it1 )+r*( dble( it2 )+r*( dble( it3 )+r*
     &  ( dble( it4 ) ) ) ) )
      return
      end
C********************************************************************
