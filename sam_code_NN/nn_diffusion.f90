
! neural_net diffusion emulator

module nn_diffusion_mod


!----------------------------------------------------------------------
use netcdf
use vars
use grid
use params , only: fac_cond, fac_fus, tprmin, a_pr, cp
implicit none
private


!---------------------------------------------------------------------
!  ---- public interfaces ----

   public  nn_diffusion, nn_diffusion_init, check

!-----------------------------------------------------------------------
!   ---- version number ----

 character(len=128) :: version = '$Id: nn_diffusion.f90,v 1 2017/08 fms Exp $'
 character(len=128) :: tag = '$Name: fez $'

!-----------------------------------------------------------------------
!   ---- local/private data ----

    logical :: do_init=.true.

    integer :: n_in ! Input dim features
    integer :: n_h1 ! hidden dim
    integer :: n_h2 ! hidden dim
    integer :: n_h3 ! hidden dim
    integer :: n_h4 ! hidden dim
    integer :: n_out ! outputs dim
    integer :: nrf  ! number of vertical levels the NN uses
    integer :: nrfq ! number of vertical levels the NN uses
    integer :: it,jt
    integer :: o_var_dim ! number of output variables  (different types)

    real(4), allocatable, dimension(:,:)     :: r_w1
    real(4), allocatable, dimension(:,:)     :: r_w2
    real(4), allocatable, dimension(:,:)     :: r_w3
    real(4), allocatable, dimension(:,:)     :: r_w4
    real(4), allocatable, dimension(:,:)     :: r_w5
    real(4), allocatable, dimension(:)       :: r_b1
    real(4), allocatable, dimension(:)       :: r_b2
    real(4), allocatable, dimension(:)       :: r_b3
    real(4), allocatable, dimension(:)       :: r_b4
    real(4), allocatable, dimension(:)       :: r_b5
    real(4), allocatable, dimension(:)       :: xscale_mean
    real(4), allocatable, dimension(:)       :: xscale_stnd

    real(4), allocatable, dimension(:)       :: z1
    real(4), allocatable, dimension(:)       :: z2
    real(4), allocatable, dimension(:)       :: z3
    real(4), allocatable, dimension(:)       :: z4
    real(4), allocatable, dimension(:)       :: z5

    real(4), allocatable, dimension(:)       :: yscale_mean
    real(4), allocatable, dimension(:)       :: yscale_stnd


contains

!#######################################################################

   subroutine nn_diffusion_init

!-----------------------------------------------------------------------
!
!        initialization for nn diffusion
!
!-----------------------------------------------------------------------
integer  unit,io,ierr

! This will be the netCDF ID for the file and data variable.
integer :: ncid
integer :: in_dimid, h1_dimid, out_dimid, single_dimid
integer :: h2_dimid, h3_dimid, h4_dimid
integer :: r_w1_varid, r_w2_varid, r_b1_varid, r_b2_varid
integer :: r_w3_varid, r_w4_varid, r_b3_varid, r_b4_varid
integer :: r_w5_varid, r_b5_varid


integer :: xscale_mean_varid, xscale_stnd_varid
integer :: yscale_mean_varid, yscale_stnd_varid

character(len=256) :: nn_filename
 
      call task_rank_to_index(rank,it,jt)
 

!-------------allocate arrays and read data-------------------------



    if(masterproc)  write(*,*) rf_filename_tkh
    nn_filename = rf_filename_tkh

! Open the file. NF90_NOWRITE tells netCDF we want read-only access
! Get the varid or dimid for each variable or dimension based on its name.

      call check( nf90_open(     trim(nn_filename),NF90_NOWRITE,ncid ))

      call check( nf90_inq_dimid(ncid, 'N_in', in_dimid))
      call check( nf90_inquire_dimension(ncid, in_dimid, len=n_in))

      call check( nf90_inq_dimid(ncid, 'N_h1', h1_dimid))
      call check( nf90_inquire_dimension(ncid, h1_dimid, len=n_h1))
      call check( nf90_inq_dimid(ncid, 'N_h2', h2_dimid))
      call check( nf90_inquire_dimension(ncid, h2_dimid, len=n_h2))

      call check( nf90_inq_dimid(ncid, 'N_h3', h3_dimid))
      call check( nf90_inquire_dimension(ncid, h3_dimid, len=n_h3))

      call check( nf90_inq_dimid(ncid, 'N_h4', h4_dimid))
      call check( nf90_inquire_dimension(ncid, h4_dimid, len=n_h4)) 
      call check( nf90_inq_dimid(ncid, 'N_out', out_dimid))
      call check( nf90_inquire_dimension(ncid, out_dimid, len=n_out))

      call check( nf90_inq_dimid(ncid, 'N_out_dim', out_dimid))
      call check( nf90_inquire_dimension(ncid, out_dimid, len=o_var_dim))

     print *, 'size of features', n_in
     print *, 'size of outputs', n_out

     nrf = 15 ! Size in the vertical 
     nrfq = 14 !Size in the vertical  for advection


      call check( nf90_open(     trim(nn_filename),NF90_NOWRITE,ncid ))

      allocate(r_w1(n_in,n_h1))
      allocate(r_w2(n_h1,n_h2))
      allocate(r_w3(n_h2,n_h3))
      allocate(r_w4(n_h3,n_h4))
      allocate(r_w5(n_h4,n_out))

      allocate(r_b1(n_h1))
      allocate(r_b2(n_h2))
      allocate(r_b3(n_h3))
      allocate(r_b4(n_h4))
      allocate(r_b5(n_out))
      allocate(z1(n_h1))
      allocate(z2(n_h2))
      allocate(z3(n_h3))
      allocate(z4(n_h4))
      allocate(z5(n_out))

     
      allocate(xscale_mean(n_in))
      allocate(xscale_stnd(n_in))

      allocate(yscale_mean(o_var_dim))
      allocate(yscale_stnd(o_var_dim))

      call check( nf90_inq_varid(ncid, "w1", r_w1_varid))
      call check( nf90_get_var(ncid, r_w1_varid, r_w1))
      call check( nf90_inq_varid(ncid, "w2", r_w2_varid))
      call check( nf90_get_var(ncid, r_w2_varid, r_w2))

      call check( nf90_inq_varid(ncid, "w3", r_w3_varid))
      call check( nf90_get_var(ncid, r_w3_varid, r_w3))
      call check( nf90_inq_varid(ncid, "w4", r_w4_varid))
      call check( nf90_get_var(ncid, r_w4_varid, r_w4))

      call check( nf90_inq_varid(ncid, "w5", r_w5_varid))
      call check( nf90_get_var(ncid, r_w5_varid, r_w5))

      call check( nf90_inq_varid(ncid, "b1", r_b1_varid))
      call check( nf90_get_var(ncid, r_b1_varid, r_b1))
      call check( nf90_inq_varid(ncid, "b2", r_b2_varid))
      call check( nf90_get_var(ncid, r_b2_varid, r_b2))

      call check( nf90_inq_varid(ncid, "b3", r_b3_varid))
      call check( nf90_get_var(ncid, r_b3_varid, r_b3))
      call check( nf90_inq_varid(ncid, "b4", r_b4_varid))
      call check( nf90_get_var(ncid, r_b4_varid, r_b4))
      call check( nf90_inq_varid(ncid, "b5", r_b5_varid))
      call check( nf90_get_var(ncid, r_b5_varid, r_b5))

      call check( nf90_inq_varid(ncid,"fscale_mean",     xscale_mean_varid))
      call check( nf90_get_var(  ncid, xscale_mean_varid,xscale_mean      ))
      call check( nf90_inq_varid(ncid,"fscale_stnd",     xscale_stnd_varid))
      call check( nf90_get_var(  ncid, xscale_stnd_varid,xscale_stnd      ))

      call check( nf90_inq_varid(ncid,"oscale_mean",     yscale_mean_varid))
      call check( nf90_get_var(  ncid, yscale_mean_varid,yscale_mean      ))
      call check( nf90_inq_varid(ncid,"oscale_stnd",     yscale_stnd_varid))
      call check( nf90_get_var(  ncid, yscale_stnd_varid,yscale_stnd      ))

    
! Close the file
      call check( nf90_close(ncid))

      write(*, *) 'Finished reading NN diffusion file.'

      do_init=.false.
   end subroutine nn_diffusion_init


!#######################################################################

   subroutine nn_diffusion

!-----------------------------------------------------------------------
!
!  NN subgrid parameterization
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!   input:  tabs               absolute temperature 
!           q                  total non-precipitating water
!           v                  meridional wind
!           u                  zonal wind
!           surface wind 
!           sst                sea surface temperatutre
!   changes: t                 liquid static energy as temperature
!            q                 total non-precipitating water
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!---------------------- local data -------------------------------------
   real,   dimension(nrf)             :: omp,fac ! Do not predict surface adv flux
   real,   dimension(nzm)             :: qsat
   real(4), dimension(n_in)     :: features
   real(4), dimension(n_out)      :: outputs
   real    omn, rev_dz, rev_dz2
   integer  i, j, k,dd, out_dim_counter, out_var_control
   integer   out_counter, ic, jc
   real     dtndxvar, dtndyvar,lat_v

   if (do_init) call error_mesg('nn_diffusion_init has not been called.')
   if (.not. rf_uses_qp) then
    ! initialize precipitation 
    if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) precsfc(:,:)=0.
   end if
 
   rev_dz = 1/dz
   rev_dz2 = 1/(dz*dz)

   do j=1,ny
    jc = j + 1
    lat_v = real((dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
    do i=1,nx

        ! Initialize variables
        features = 0.
        outputs = 0.
        z1 = 0.
        z2 = 0.
        z3 = 0.
        z4 = 0.
        z5 = 0.        
        omp = 0.
        fac = 0.
        out_counter = 0
       ! Combine all features into one vector
        features(1:nrf) = real(t_i(i,j,1:nrf),4)


        if (rf_uses_rh) then
         ! generalized relative humidity is used as a feature
         do k=1,nzm
          omn = omegan(tabs(i,j,k))
          qsat(k) = omn*qsatw(tabs(i,j,k),pres(k))+(1.-omn)*qsati(tabs(i,j,k),pres(k))
         end do
         features(nrf+1:2*nrf) = real(q(i,j,1:nrf)/qsat(1:nrf),4)
        else
         ! non-precipitating water is used as a feature
         features(nrf+1:2*nrf) = real(q_i(i,j,1:nrf),4)
        endif
        dtndxvar = dtn/dx
        do k=1,nrf
         dtndyvar = 1/rho(k)/dtndxvar
         features(2*nrf + k) = real((0.5*(u(i,j,k)*dtndyvar + u(ic,j,k)*dtndyvar)),4)
         if (lat_v.GE.0) then
          features(3*nrf + k) = real((0.5*(v(i,j,k)*dtndyvar + v(i,jc,k)*dtndyvar)),4)
         else
          features(3*nrf + k) = -real((0.5*(v(i,j,k)*dtndyvar + v(i,jc,k)*dtndyvar)),4)
         end if
        end do
     
        dtndyvar = 1/rho(1)/dtndxvar
       features(4*nrf+1) = sqrt(real(0.5*(u(i,j,1)* dtndyvar +u(ic,j,1)* dtndyvar) ,4)**2 + real(0.5*(v(i,j,1)* dtndyvar + v(i,jc,1)*dtndyvar),4) **2)
       !Use SST here
       features(4*nrf+2) = sstxy(i,j)    

!Noralize features
       features = (features - xscale_mean) / xscale_stnd
        
! calculate predicted values using NN
! Apply trained regressor network to data using rectifier activation function
! forward prop to hiddelayer

        z1 = matmul( features, r_w1) + r_b1
        !print *, 'SHAPE of Z1', shape(z1)
! rectifier
        where (z1 .lt. 0.0)  z1 = 0.0
! forward prop to output layer
        z2 = matmul( z1,r_w2) + r_b2
        where (z2 .lt. 0.0)  z2 = 0.0

        z3 = matmul( z2,r_w3) + r_b3
        where (z3 .lt. 0.0)  z3 = 0.0

        z4 = matmul( z3,r_w4) + r_b4
        where (z4 .lt. 0.0)  z4 = 0.0

       outputs = matmul( z4,r_w5) + r_b5

        out_var_control =1
        outputs(1) = (outputs(1) * yscale_stnd(out_var_control))  +  yscale_mean(out_var_control) 
        out_var_control = out_var_control + 1
        outputs(2) = (outputs(2) * yscale_stnd(out_var_control))  +  yscale_mean(out_var_control)
        out_var_control = out_var_control + 1
        outputs(3:2 + nrf) = (outputs(3:2 + nrf) * yscale_stnd(out_var_control))  +  yscale_mean(out_var_control)    

        if (do_q_surf_fluxes_out_rf) then
         t(i,j,1) = t(i,j,1) + outputs(1)*dtn
         q(i,j,1) = q(i,j,1) + outputs(2)*dtn
         out_counter = 2
        end if

        if (j.eq.0.or.i.eq.0) then
         tkh_x_rf(i,j,1:nrf) = outputs(out_counter+1:out_counter+nrf)
         tkh_y_rf(i,j,1:nrf)  = outputs(out_counter+1:out_counter+nrf)
        else
         tkh_x_rf(i,j,1:nrf) = outputs(out_counter+1:out_counter+nrf)
         tkh_y_rf(i,j,1:nrf) = outputs(out_counter+1:out_counter+nrf)
         tkh_z_rf(i,j,1:nrf) = outputs(out_counter+1:out_counter+nrf)
       end if
        do k=1,nrf
         tkh_z_rf(i,j,k) = max(0.,tkh_z_rf(i,j,k))
        end do


       end do
     end do
    
   end subroutine nn_diffusion





!#######################################################################


!##############################################################################
  subroutine check(status)

    ! checks error status after each netcdf, prints out text message each time
    !   an error code is returned. 

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
       write(*, *) trim(nf90_strerror(status))
    end if
  end subroutine check

!#######################################################################
 subroutine error_mesg (message)
  character(len=*), intent(in) :: message

!  input:
!      message   message written to output   (character string)

    if(masterproc) print*, 'Neural network  module: ', message
    stop

 end subroutine error_mesg



!#######################################################################


end module nn_diffusion_mod

