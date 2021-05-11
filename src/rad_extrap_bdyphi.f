c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the scalar field phi of the boundary CFT
c----------------------------------------------------------------------

        subroutine bdyphi_radextrap(bdyphi,
     &                  leadordcoeff_phi1,
     &                  xpbdy,ypbdy,zpbdy,
     &                  chrbdy,numbdypoints,
     &                  bdy_extrap_order,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        integer bdy_extrap_order
        real*8 chrbdy(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        integer i,j,k,is,ie,js,je,ks,ke,lind
        integer ipa,jpa,kpa
        integer ipb,jpb,kpb
        integer ipc,jpc,kpc
        integer ipd,jpd,kpd
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0
        real*8 xpa,xpb,xpc,xpd
        real*8 ypa,ypb,ypc,ypd
        real*8 zpa,zpb,zpc,zpd

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 leadordcoeff_phi1(Nx,Ny,Nz)
        real*8 leadordcoeff_phi1_p1

        real*8 leadordcoeff_phi1_pa
        real*8 leadordcoeff_phi1_pb
        real*8 leadordcoeff_phi1_pc
        real*8 leadordcoeff_phi1_pd

        real*8 leadordcoeff_phi1_p2
        real*8 leadordcoeff_phi1_p3
        real*8 leadordcoeff_phi1_p4
        real*8 bdyphi(numbdypoints)

        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)

        real*8 xp1,yp1,zp1
        real*8 xp2,yp2,zp2
        real*8 xp3,yp3,zp3
        real*8 xp4,yp4,zp4
        real*8 xex,yex,zex,rhoex,chiex,xiex
        real*8 rhop1,chip1,xip1
        real*8 rhop2,chip2,xip2
        real*8 maxxyzp1

        real*8 bilinear_interp
        real*8 firstord_extrap
        real*8 secondord_extrap
        real*8 thirdord_extrap
        real*8 fourthord_extrap
!----------------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)


        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke

           if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

              xp1=x(i)
              yp1=y(j)
              zp1=z(k)
              rhop1=sqrt(xp1**2+yp1**2+zp1**2)
              chip1=(1/PI)*acos(xp1/rhop1)
              if (zp1.lt.0) then
                xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
              else
                xip1=(1/(2*PI))*atan2(zp1,yp1)
              end if
  
              leadordcoeff_phi1_p1=leadordcoeff_phi1(i,j,k)

              xex=xpbdy(lind)
              yex=ypbdy(lind)
              zex=zpbdy(lind)
              rhoex=1
              chiex=(1/PI)*acos(xex/rhoex)
              if (zex.lt.0) then
                xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
              else
                xiex=(1/(2*PI))*atan2(zex,yex)
              end if



            if (bdy_extrap_order.eq.1) then

                        bdyphi(lind)=
                        firstord_bdyphi_radextrap(bdyphi,
     &                  leadordcoeff_phi1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)




!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.1


           end if !closes condition on chrbdy(i,j,k).ne.ex




          end do
         end do
        end do


        return
        end
c--------------------------------------------------------------------------------------


c----------------------------------------------------------------------
c bilinear interpolation on a plane using values 
c Tpa at (xpa,ypa), Tpb at (xpb,ypb), Tpc at (xpc,ypc), Tpd at (xpd,ypd)  
c to obtain the value at (xp0,yp0) 
c----------------------------------------------------------------------
        real*8 function bilinear_interp(Tpa,Tpb,Tpc,Tpd,
     &                                  xpa,xpb,xpc,xpd,
     &                                  ypa,ypb,ypc,ypd,
     &                                  xp0,yp0)
        implicit none
        real*8 Tpa,Tpb,Tpc,Tpd
        real*8 xpa,xpb,xpc,xpd
        real*8 ypa,ypb,ypc,ypd
        real*8 xp0,yp0

        !--------------------------------------------------------------

        bilinear_interp=(1/((xpc-xpa)*(ypc-ypa)))*
     &    (Tpa*(xpc-xp0)*(ypc-yp0)+
     &     Tpb*(xpc-xp0)*(yp0-ypa)+
     &     Tpc*(xp0-xpa)*(yp0-ypa)+
     &     Tpd*(xp0-xpa)*(ypc-yp0))

        return
        end
c--------------------------------------------------------------------------------------