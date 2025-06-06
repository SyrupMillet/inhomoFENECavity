!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private

   !> Single config
   type(config), public :: cfg

   public :: geometry_init

contains


   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid

      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         real(WP), dimension(:), allocatable :: weights
         real(WP) :: center, delta, ks, dominator

         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))

         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)*Lx/real(nx,WP)
         end do

         ! refine the mesh near to the wall, use a step function for dy
         ! allocate(weights(ny))
         ! center = (real(ny,WP)-1.0_WP)/5.0_WP   ! center of step function
         ! delta = 1.0_WP/2.0_WP                  ! dy(1):dy(ny) = 1-delta:1+delta
         ! ks = 0.1_WP*128.0_WP/real(ny,WP)       ! control the variance
         ! do j=0,ny-1
         !    if (ks*center .gt. 0.0_WP) then
         !       dominator = tanh(ks*center)
         !    else
         !       dominator = 1.0_WP
         !    end if
         !    weights(j+1) = 1.0_WP + delta*tanh(ks*(real(j,WP)-center))/dominator
         ! end do
         ! weights = weights/sum(weights)
         ! y(1) = 0.0_WP
         ! do j=2,ny+1
         !    y(j) = y(j-1) + weights(j-1)*Ly
         ! end do
         ! deallocate(weights)

         do j=1,ny+1
            y(j)=real(j-1,WP)*Ly/real(ny,WP)
         end do

         do k=1,nz+1
            z(k)=real(k-1,WP)*Lz/real(nz,WP)-0.5_WP*Lz
         end do

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='box')

      end block create_grid

      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_cfg

      ! Create walls for x<0
      create_walls: block
         real(WP) :: Lx
         integer :: i,j,k

         call param_read('Lx',Lx)

         cfg%VF=1.0_WP
         do k = cfg%kmino_, cfg%kmaxo_
            do j = cfg%jmino_, cfg%jmaxo_
               do i = cfg%imino_, cfg%imaxo_
                  if (j.lt.cfg%jmin) then
                     cfg%VF(i,j,k)=0.0_WP
                  end if
               end do
            end do
         end do

      end block create_walls

   end subroutine geometry_init


end module geometry

