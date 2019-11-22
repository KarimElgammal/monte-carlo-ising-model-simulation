!ising_mpi_reduce.f90
program isingmodel_mpi
   implicit none
   include 'mpif.h'
   
   real :: kT, paccept,ptrial
   integer :: Estart,Eflip, iEtot,iEtotloc
   integer :: col, row, east,west,north,south
   integer :: nMC
   integer, parameter :: numL = 8 ! the number of row or column
   integer, parameter :: J = 1    ! define a ferromagnetic system
   integer, dimension(:,:), allocatable :: spin

        !mpi processes local index
   integer, parameter :: ni=numL, nj=numL
   integer :: ierr, myid, nprocs, i1, i2, j1, j2, i1p, i2m, j1p, j2m, &
              i1n, i2n, ninom, njnom, niproc, njproc, nitot
   integer :: status(mpi_status_size), row_type
   integer :: myparity,print_id
   integer,parameter :: H = 1, D = 0        ! tags of passing messages
   
   
     ! initialize MPI

   call mpi_init(ierr)
   call mpi_comm_rank(mpi_comm_world, myid, ierr)
   call mpi_comm_size(mpi_comm_world, nprocs, ierr)   
 
   ! domain decomposition

  ! nominal number of points per proc., without ghost cells,
  ! assume numbers divide evenly; niproc and njproc are the
  ! numbers of procs in the i and j directions.
   niproc = nprocs
   njproc = 1
   ninom  = ni/niproc
   njnom  = nj/niproc

  ! nominal starting and ending indices, without ghost cells
   i1n = myid*ninom + 1
   i2n = i1n + ninom - 1

  ! local starting and ending index, including 2 ghost cells
  ! in each direction (at beginning and end)
   i1  = i1n - 1
   i1p = i1 + 1
   i2  = i2n + 1
   i2m = i2 - 1
   j1  = 0
   j1p = j1 + 1
   j2  = nj + 1
   j2m = j2 - 1
   nitot = i2 - i1 + 1
   
      ! allocate arrays
   allocate( spin(i1:i2,j1:j2) )

      ! random number generator
   call random_seed()

      ! initialize spin configuration
   do col=j1p,j2m
      do row=i1p,i2m
         call random_number(paccept)
         if (paccept .lt. 0.5) then
            spin(row,col) = +1
         else 
            spin(row,col) = -1
         end if
      end do
   end do

   do print_id = 0,nprocs-1
      if(myid ==print_id) then
         write(*,*), 'initialized spin in Processor=',myid,'is' 
         write(*,*), spin(i1p:i2m,j1p:j2m)
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierr)   
   end do    
   
  ! Create derived type for single row of array.
  ! There are nj "blocks," each containing 1 element,
  ! with a stride of nitot between the blocks

   call mpi_type_vector(nj+2, 1, nitot, mpi_integer, row_type, ierr);
   call mpi_type_commit(row_type, ierr);
   
    !  iterate
    !  messages passing divided in two steps depending on myid even or odd 
   myparity = mod(myid,2)
   

   Temp_decreasing : do kT=4.0, 0.1, -0.05      !JS: Start at kT=4, and go to kT=0.1 in steps of -0.05


 !  transfor data to ghost cells via subroutine transfor_data

   call transfor_data(spin,i1,i1p,i2,i2m,j1,j1p,j2,j2m,row_type,nprocs,myid, &
        myparity,mpi_comm_world,ierr,status,H,D)

! loop MC trials
      MC_iteration: do nMC=1,1000
        
      
      do col=j1p,j2m
         do row=i1p,i2m

            east = col+1
            west = col-1
            north = row-1
            south = row+1

       ! try flip a spin at position (row,col)

       !call transfor_data(spin,i1,i1p,i2,i2m,j1,j1p,j2,j2m,row_type,nprocs,myid, &
       !     myparity,mpi_comm_world,ierr,status,H,D)
            
            
            Estart = -J*(spin(row,col)*spin(row,east) &
                 + spin(row,col)*spin(row,west) &
                 + spin(row,col)*spin(north,col) &
                 + spin(row,col)*spin(south,col) )

            Eflip =  -J*((-spin(row,col))*spin(row,east) &
                 + (-spin(row,col))*spin(row,west) &
                 + (-spin(row,col))*spin(north,col) &
                 + (-spin(row,col))*spin(south,col) )

      ! acceptance/rejectance for MC algorithm 

            if (Eflip .lt. Estart) then 
               spin(row,col) = -spin(row,col)
            else  
               call random_number(paccept)
               ptrial = EXP(-(Eflip-Estart)/kT)

               if (paccept .lt. ptrial) then
                  spin(row,col) = -spin(row,col)
               end if
            end if

         end do
      end do
      
      enddo MC_iteration
   enddo Temp_decreasing
! calculate enery in each decomposition domain
! mpi_reduce to summerize the energy in master processor and print the results

   call sum_energy(spin,iEtot,i1,i1p,i2,i2m,j1,j1p,j2,j2m)
   
   if(myid==0) then 
     write(*,*),"Temperature =", kT     
     write(*,*),"the total energy at this Temperature =", iEtot
   endif
   
   do print_id = 0,nprocs-1
      if(myid ==print_id) then
         write(*,*), 'spin configurations in Processor=',myid,'is' 
         write(*,*), spin(i1p:i2m,j1p:j2m)
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierr)   
   end do
   deallocate(spin)
   
                                       
   call mpi_finalize(ierr)   

contains 

   subroutine sum_energy(spin,iEtot,i1,i1p,i2,i2m,j1,j1p,j2,j2m)
      integer,intent(in) :: i1,i1p,i2,i2m,j1,j1p,j2,j2m
      integer,intent(in), dimension(i1:i2,j1:j2) :: spin
      integer,intent(out) :: iEtot
      integer :: iEtot_loc
      integer :: east,west,north,south,col,row
      iEtot = 0
      do col=j1p,j2m
         do row=i1p,i2m

            east = col+1
            west = col-1
            north = row-1
            south = row+1
            
            iEtot = iEtot -J*(spin(row,col)*spin(row,east) &
                 + spin(row,col)*spin(row,west) &
                 + spin(row,col)*spin(north,col) &
                 + spin(row,col)*spin(south,col) )
                       
         end do
      end do 
      
      if(nprocs >1) then
        iEtot_loc = iEtot
        call mpi_reduce(iEtot_loc, iEtot, 1, mpi_integer, &
             mpi_sum, 0, mpi_comm_world, ierr)  
      end if           
   end subroutine sum_energy  

subroutine transfor_data(spin,i1,i1p,i2,i2m,j1,j1p,j2,j2m,row_type,nprocs,myid, &
           myparity,mpi_comm_world,ierr,status,H,D)
   integer,intent(in) :: i1,i1p,i2,i2m,j1,j1p,j2,j2m 
   integer,intent(in) :: row_type,nprocs,myid,myparity,mpi_comm_world, & 
                         ierr,status(mpi_status_size)
   integer,parameter :: H=1, D=0         
   integer,intent(inout), dimension(i1:i2,j1:j2) :: spin 
   
        ! left and right boundary conditions
      spin(i1p:i2m,0)  = spin(i1p:i2m,j2m)
      spin(i1p:i2m,j2) = spin(i1p:i2m,1)
     
      if(nprocs == 1) then

        ! top and bottom boundary conditions

         spin(i1,:) = spin(i2m,:)
         spin(i2,:) = spin(i1p,:)

        ! corners

         spin(i1,j1) = spin(i2m,j2m)
         spin(i1,j2) = spin(i2m,j1p)
         spin(i2,j2) = spin(i1p,j1p)
         spin(i2,j1) = spin(i1p,j2m)

      elseif(myid == 0) then

        ! top and bottom boundary conditions

         call mpi_send(spin(i1p,j1), 1, row_type, &
              nprocs-1, D, mpi_comm_world, ierr)
         call mpi_recv(spin(i1,j1),  1, row_type, &
              nprocs-1, D, mpi_comm_world, status, ierr)
         call mpi_send(spin(i2m,j1), 1, row_type, &
              myid+1, H, mpi_comm_world, ierr)
         call mpi_recv(spin(i2,j1),  1, row_type, &
              myid+1, H, mpi_comm_world, status, ierr)

        ! corners

         call mpi_send(spin(i1p,j1p), 1, mpi_integer, &
              nprocs-1, 2, mpi_comm_world, ierr)
         call mpi_recv(spin(i1, j1 ), 1, mpi_integer, &
              nprocs-1, 3, mpi_comm_world, status, ierr)
         call mpi_send(spin(i1p,j2m), 1, mpi_integer, &
              nprocs-1, 4, mpi_comm_world, ierr)
         call mpi_recv(spin(i1, j2 ), 1, mpi_integer, &
              nprocs-1, 5, mpi_comm_world, status, ierr)
      elseif(myid < nprocs-1) then

        ! the processes in between (except P0 and the P-last ) only have top and bottom boundary conditions
         if(myparity.eq.1) then
            call mpi_send(spin(i2m,j1), 1, row_type, &
                 myid+1, D, mpi_comm_world, ierr)
            call mpi_recv(spin(i2,j1),  1, row_type, &
                 myid+1, D, mpi_comm_world, status, ierr)
            call mpi_send(spin(i1p,j1), 1, row_type, &
                 myid-1, H, mpi_comm_world, ierr)
            call mpi_recv(spin(i1,j1),  1, row_type, &
                 myid-1, H, mpi_comm_world, status, ierr)

         else                                                 !(myparity.eq.0)
            call mpi_send(spin(i1p,j1), 1, row_type, &
                 myid-1, D, mpi_comm_world, ierr)
            call mpi_recv(spin(i1,j1),  1, row_type, &
                 myid-1, D, mpi_comm_world, status, ierr)
            call mpi_send(spin(i2m,j1), 1, row_type, &
                 myid+1, H, mpi_comm_world, ierr)  
            call mpi_recv(spin(i2,j1),  1, row_type, &
                 myid+1, H, mpi_comm_world, status, ierr)

         endif                       
      else

        ! top and bottom boundary conditions

         call mpi_send(spin(i2m,j1), 1, row_type, &
              0, D, mpi_comm_world, ierr)
         call mpi_recv(spin(i1,j1),  1, row_type, &
              myid-1, H, mpi_comm_world, status, ierr)
         call mpi_send(spin(i1p,j1), 1, row_type, &
              myid-1, H, mpi_comm_world, ierr)
         call mpi_recv(spin(i2,j1),  1, row_type, &
              0, D, mpi_comm_world, status, ierr)
        ! corners

         call mpi_recv(spin(i2, j2 ), 1, mpi_integer, &
              0, 2, mpi_comm_world, status, ierr)
         call mpi_send(spin(i2m,j2m), 1, mpi_integer, &
              0, 3, mpi_comm_world, ierr)
         call mpi_recv(spin(i2, j1 ), 1, mpi_integer, &
              0, 4, mpi_comm_world, status, ierr)
         call mpi_send(spin(i2m,j1p), 1, mpi_integer, &
              0, 5, mpi_comm_world, ierr)
      endif   
end subroutine transfor_data                                                 
end program isingmodel_mpi