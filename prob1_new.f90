
!This is the final code for Tully's 1st problem.
!Author: Maitrayee Ghosh
!======================================================

program cal_c
 implicit none

 integer :: mom_index, i_traj, i_surf, &
            total_traj, tran_on_s1, tran_on_s2, &
             refl_on_s1, refl_on_s2 
            

 real(kind=8) :: p, rand_num, zeta, direction

 !integer, parameter :: p_divisions = 20
 real(kind=8), parameter :: mass=2.0d3, h_bar = 1.0d0, dt = 0.01d0, &
                   A = 0.01d0, B= 1.6d0 , C= 0.005d0, D= 1.0d0, &
                   p_start = 5.0d0, p_end = 10.0d0, &
                   delta_p = 0.5d0

 real(kind=8) :: E_plus, E_minus, &
                V11, V22, V12, &
                x_old, x_new,&
                x_start = -10.0d0, x_end = 10.0d0, &
                v_old, v_adj, v_new, &
                delta_E, KE, d12, & 
                b12, b21, &
                k1_x, k2_x, k3_x, k4_x, &
                k1_v, k2_v, k3_v, k4_v


 complex(kind=8) :: C1_old, C2_old, &
                    C1_new,C2_new , &
                    conj_c1_old, conj_c2_old, &
                    k1_c1, k1_c2, &
                    k2_c1, k2_c2, &
                    k3_c1, k3_c2, &
                    k4_c1, k4_c2
 
 complex(kind=8), parameter :: im_i = cmplx(0.0d0, 1.0d0)
integer :: i


open(unit=40, file='switch1.out')
open(unit=42, file='probability1.out')

 do mom_index = 1, 11 !1, p_divisions+1             !!!!!!!!!!!!!!!!!!!!!!!!!!!


   total_traj = 0
   tran_on_s1 = 0
   tran_on_s2 = 0
   refl_on_s1 = 0
   refl_on_s2 = 0

   write(40,*)'momentum: ', mom_index
 ! do i = -10, 10
  ! x_old = dble(i)
   ! d12 = get_d12(x_old)
    ! call get_E(x_old, E_plus, E_minus)
    !write(44,*) x_old, d12 !, E_plus, E_minus
   !enddo
  ! stop

   do i_traj = 1, 2000 

     write(40,*) '--> traj: ', i_traj

     ! Initialize the trajectory
     C1_new = cmplx(1.0d0, 0.0d0)
     C2_new = cmplx(0.0d0, 0.0d0)     

     x_new = -10.0d0
     p = p_start + (mom_index-1)*delta_p
     write(40,*)'p = ',p 
      !p = dble(mom_index)
     v_new = p/mass
     i_surf = 1

     do while ( x_new <= x_end .AND. x_new >= x_start)

        v_old = v_new
        x_old = x_new
        C1_old = C1_new
        C2_old = C2_new

         !write(40,*)'   x = ', x_old
         !1. Check if hopping is possible
         !--------------------------------
  
         !Generate random number, rand_num
         rand_num = rand(0) 
         conj_c1_old = conjg(c1_old)
         conj_c2_old = conjg(c2_old)

         d12 = get_d12(x_old)

         b21 = 2.0d0*Real(conj_c2_old*c1_old*v_old*d12)            

                 !V12 = 0 for adiabatic basis

         b12 = -2.0d0*Real(conj_c1_old*c2_old*v_old*d12)
                 !V12 = 0 for adiabatic basis


         !Hopping factor (zeta)
         
         if (i_surf == 1) then         !!For state 1 to state 2 hop 
             zeta = dt*b21/(c1_old*conj_c1_old)
              
         elseif (i_surf == 2) then     !!For state 2 to state 1 hop
             zeta = dt* b12 /(c2_old* conj_c2_old)

         endif


         !Check if hop can occur 

         call get_E(x_old, E_plus, E_minus)

         delta_E = E_plus - E_minus
         KE = 0.5d0*mass*V_old*V_old

         V_adj = V_old

         if ( abs(v_old) > 1.0d-12 ) then   ! just as a safety precaution
            direction = v_old/(abs(v_old))
         endif

 
            if ((i_surf == 1) .AND. (zeta > rand_num) .AND. (KE-delta_E >= 0) )then
           
                i_surf = 2
                v_adj = direction*sqrt( 2.0d0 * (KE - delta_E) / mass )
                write(40,*) '    hop 1 ---> 2 at x = ', x_old,'and v = ',v_old


            elseif ((i_surf == 2) .AND. (zeta > rand_num) ) then
 
                i_surf = 1
                v_adj = direction*sqrt( 2.0d0 * (KE + delta_E) / mass )
                write(40,*) '    hop  2 ---> 1 at x = ', x_old,'and v = ',v_old

            endif


        v_old = v_adj





        !2. RK4
         

        !Evaluate x_new, v_new, c1 and c2 using RK4
  

         ! get k1 
         ! (ki_j = ki corresponding to Cj)
         call get_k(x_old, v_old, C1_old, C2_old, k1_x, k1_v, k1_c1, k1_c2)



         !Calculate k2

         call get_k(x_old+k1_x*dt/2.0d0, v_old+k1_v*dt/2.0d0, &
                    C1_old + k1_c1*dt/2.0d0, C2_old+k1_c2*dt/2.0d0,&
                    k2_x, k2_v, k2_c1, k2_c2)
      
 
         !get k3
         call get_k(x_old+k2_x*dt/2.0d0, v_old+k2_v*dt/2.0d0, &
                    C1_old + k2_c1*dt/2.0d0, C2_old+k2_c2*dt/2.0d0,&
                    k3_x, k3_v, k3_c1, k3_c2)


         !get k4
         call get_k(x_old+k3_x*dt, v_old+k3_v*dt, &
                     C1_old + k3_c1*dt, C2_old+k3_c2*dt,&
                    k4_x, k4_v, k4_c1, k4_c2)

         !Evaluate x_new, v_new, c1_new and c2_new
 
         x_new = x_old + (dt/6.0d0)*(k1_x + 2.0d0*k2_x + 2.0d0*k3_x + k4_x)
         v_new = v_old + (dt/6.0d0)*(k1_v + 2.0d0*k2_v + 2.0d0*k3_v + k4_v)

         c1_new = c1_old + (dt/6.0d0)*(k1_c1 + 2.0d0*k2_c1 + 2.0d0*k3_c1 + k4_c1)
         c2_new = c2_old + (dt/6.0d0)*(k1_c2 + 2.0d0*k2_c2 + 2.0d0*k3_c2 + k4_c2)


       end do 

       total_traj = total_traj + 1

       if (x_new > x_end ) then  ! => transmitted

           if (i_surf == 1 ) then
               tran_on_s1 = tran_on_s1 + 1
           else
               tran_on_s2 = tran_on_s2 + 1
           endif

       else  ! => reflected

           if (i_surf == 1 ) then
               refl_on_s1 = refl_on_s1 + 1
           else
               refl_on_s2 = refl_on_s2 + 1
           endif

       end if
      

   end do 

   write(42, 400) p, dble(tran_on_s1)/dble(total_traj), &
                     dble(tran_on_s2)/dble(total_traj), &
                     dble(refl_on_s1)/dble(total_traj), &
                     dble(refl_on_s2)/dble(total_traj)

end do 

400 format(f10.3, 2x, E16.5, 2x, E16.5, 2x, E16.5, 2x, E16.5)

 close(unit=40)
 close(unit=42)

contains

 !Define the subroutine get_V
 !---------------------------------
 subroutine get_V(x_old,V11, V22, V12)
   real(kind=8), intent(in) :: x_old
   real(kind=8), intent(out) :: V11, V22, V12

      if (X_old >= 0) then
          V11 = A*(1.0d0 - exp(-B*X_old))
      else 
          V11 = -A*(1.0d0 - exp(B*X_old))
      endif


   V22 = -V11

   V12 = C*exp(-D*X_old*X_old)


 end subroutine get_V


 !Define the subroutine get_dV
 !-----------------------------
 subroutine get_dV(x_old, dV11, dV22, dV12)
   real(kind=8),intent(in) :: x_old
   real(kind=8),intent(out) :: dV11, dV22, dV12

   if (x_old >= 0) then
       dV11 = A*B*exp(-B*x_old)
   else
       dV11 = A*B*exp(B*x_old)
   endif

   dV22 = -dV11

   dV12 = -2.0d0*C*D*x_old*exp(-D*x_old*x_old)


  end subroutine get_dV
 


 !Define the subroutine get_E
 !------------------------------------

 subroutine get_E(x_old, E_plus, E_minus)
   real(kind=8):: x_old, V11, V22, V12
   real(kind=8),intent(out) :: E_plus, E_minus

   call get_V(x_old, V11, V22, V12)
   E_plus = 0.5d0*((V11+V22) + sqrt((V11-V22)*(V11-V22)  & 
            + 4.0d0*V12*V12))

   E_minus = 0.5d0*((V11+V22) - sqrt((V11-V22)*(V11-V22) &
             + 4.0d0*V12*V12))
   

  end subroutine get_E  
  

 !Define the function get_dE
 !------------------------------------

 real(kind=8) function get_dE(i_surf, x_old) result(dE)
   integer :: i_surf
   real(kind=8) :: x_old, dV11, dV22, dV12 
   
    call get_V(x_old, V11, V22, V12)
    call get_dV(x_old, dV11, dV22, dV12)

    if(i_surf == 1) then 
       dE = 0.5d0*((dV11+dV22)- &
            ((V11-V22)*(dV11-dV22)+4.0d0*V12*dV12)/ &
             sqrt((V11-V22)**2+4.0d0*V12*V12))

   else
       dE = 0.5d0*((dV11+dV22)+ &
            ((V11-V22)*(dV11-dV22)+4.0d0*V12*dV12)/ &
             sqrt((V11-V22)**2+4.0d0*V12*V12))

   endif


  end function get_dE





 
 !Define the function get_d12
 !------------------------------------
 real(kind=8) function get_d12(x) result(d12)
 real(kind=8) :: x, V11, V22, V12, &
                   dV11, dV22, dV12, &
                   z, u, du_over_dz, dz_over_dx

   call get_V(x, V11, V22, V12)
   call get_dV(x, dV11, dV22, dV12)

   z = (V11 - V22)/(2.0d0*V12)
   u = z - sqrt(z*z + 1.0d0)
   du_over_dz = 1.0d0 - z/ sqrt(z*z+1.0d0)
   dz_over_dx = (0.5d0/(V12*V12))  &
         *( V12*(dV11 -dV22) - (V11-V22)*dV12 )
          
   d12 = du_over_dz * dz_over_dx / (1.0d0 + u*u)
   


 end function get_d12



 !Define the subroutine get_k
 !----------------------------------
 subroutine get_k(x, v, C1, C2, k1_x, k1_v, k1_c1, k1_c2)
   
   real(kind=8), intent(in) :: x, v
   complex(kind=8), intent(in) :: C1, C2
   real(kind=8), intent(out) :: k1_x, k1_v
   complex(kind=8), intent(out) :: k1_c1, k1_c2
   real(kind=8) :: d12, E_plus, E_minus, dE

   d12 = get_d12(x)
   call get_E(x, E_plus, E_minus)
   k1_x = v
   dE = get_dE(i_surf, x)
   k1_v = -dE/mass 
   k1_c1 = -im_i *(c1*E_minus - c2*im_i*h_bar*v*d12)/h_bar
   k1_c2 = -im_i *(c1*im_i*h_bar*v*d12 + c2*E_plus)/h_bar
 end subroutine get_k
   
 end program cal_c
