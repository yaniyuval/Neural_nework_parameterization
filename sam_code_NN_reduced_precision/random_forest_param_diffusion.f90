
! random forest parameterization

module random_forest_param_diffusion


!----------------------------------------------------------------------
use netcdf
use vars
use grid


implicit none
private


!---------------------------------------------------------------------
!  ---- public interfaces ----

   public  random_forest_diffusion, random_forest_init_diffusion, check_diffusion


!-----------------------------------------------------------------------
!   ---- local/private data ----

    logical :: do_init=.true.
    integer :: n_features ! number of features
    integer :: max_n_nodes ! maximum number of nodes across trees
    integer :: n_trees  ! number of trees
    integer :: n_outputs  ! number of outputs
    integer :: nrf  ! number of vertical levels the random forest uses
    integer :: it,jt

    integer(4), allocatable, dimension(:,:)  :: children_left
    integer(4), allocatable, dimension(:,:)  :: children_right
    integer(4), allocatable, dimension(:,:)  :: split_feature
    real(4), allocatable, dimension(:,:)     :: threshold
    real(4), allocatable, dimension(:,:,:)   :: values_predicted
    integer(4), allocatable, dimension(:)    :: n_nodes

! random tree mod
!real(4), allocatable, dimension(:,:)     :: rand_tree



! do_init           has the module been initialized

! children_left     left child of node (n_trees, max_n_nodes)
! children_right    right child of node
! split_feature     feature used for splitting the node
! threshold         threshold value at the node
! values_predicted  prediction for a given node (n_trees, max_n_nodes, n_lev)
! n_nodes           number of nodes in a given tree (n_trees)


!-----------------------------------------------------------------------


contains

!#######################################################################

   subroutine random_forest_init_diffusion

!-----------------------------------------------------------------------
!
!        initialization for random forest
!
!-----------------------------------------------------------------------
integer  unit,io,ierr

! This will be the netCDF ID for the file and data variable.
integer :: ncid
integer :: n_nodes_varid, children_left_varid, children_right_varid
integer :: trees_dimid, features_dimid, outputs_dimid, nodes_dimid
integer :: split_feature_varid, threshold_varid, values_predicted_varid
!character(len=256) :: rf_filename_diffusion
      call task_rank_to_index(rank,it,jt)
  
!-------------allocate arrays and read data-------------------------


      
      if(masterproc)  write(*,*) rf_filename_diffusion
      !rf_filename_diffusion ='/glade/work/janniy/mldata/gcm_regressors/qobsFFTFTTTFF0FFTFTF48TFFFFFFFFF_F-NoSc_O-Stan_Ntr4277070_Nte225180_F_Tin_qin_uin_vin_usurf_O_Tout_qout_RF_NTr10_MinS10max_d27_maxzinf_nocos_te83_tr87.nc'
! Open the file. NF90_NOWRITE tells netCDF we want read-only access
! Get the varid or dimid for each variable or dimension based on its name.

      call check_diffusion( nf90_open(     trim(rf_filename_diffusion),NF90_NOWRITE,ncid ))

      call check_diffusion( nf90_inq_dimid(ncid, 'dim_trees', trees_dimid))
      call check_diffusion( nf90_inquire_dimension(ncid, trees_dimid, len=n_trees))

      call check_diffusion( nf90_inq_dimid(ncid, 'dim_features', features_dimid))
      call check_diffusion( nf90_inquire_dimension(ncid, features_dimid, len=n_features))

      call check_diffusion( nf90_inq_dimid(ncid, 'dim_outputs', outputs_dimid))
      call check_diffusion( nf90_inquire_dimension(ncid, outputs_dimid, len=n_outputs))

      call check_diffusion( nf90_inq_dimid(ncid, 'dim_nodes', nodes_dimid))
      call check_diffusion( nf90_inquire_dimension(ncid, nodes_dimid, len=max_n_nodes))

      allocate(n_nodes(n_trees))

      call check_diffusion( nf90_inq_varid(ncid, "n_nodes", n_nodes_varid))
      call check_diffusion( nf90_get_var(ncid, n_nodes_varid, n_nodes))

! basic checks of the expected data structure
    !  if(n_outputs /= n_features) call error_mesg_diffusion('n_outputs not equal to n_features.')

      if(max_n_nodes /= maxval(n_nodes)) call error_mesg_diffusion('max_n_nodes not equal to maximum of n_nodes')

!      if (rf_uses_qp) then
!       nrf = n_features/3
!      else
!       nrf = n_features/2
!      endif
      nrf = 48

      allocate(children_left(n_trees, max_n_nodes))
      allocate(children_right(n_trees, max_n_nodes))
      allocate(split_feature(n_trees, max_n_nodes))
      allocate(threshold(n_trees, max_n_nodes))
      allocate(values_predicted(n_trees, max_n_nodes, n_outputs))
! mod random tree
!allocate(rand_tree(nx, ny))
	


print *, 'size of features and  outputs in diffusion RF is:'
      print *, n_features,n_outputs

      call check_diffusion( nf90_inq_varid(ncid, "children_left", children_left_varid))
      call check_diffusion( nf90_get_var(ncid, children_left_varid, children_left))
      
      call check_diffusion( nf90_inq_varid(ncid, "children_right", children_right_varid))
      call check_diffusion( nf90_get_var(ncid, children_right_varid, children_right))

      call check_diffusion( nf90_inq_varid(ncid, "split_feature", split_feature_varid))
      call check_diffusion( nf90_get_var(ncid, split_feature_varid, split_feature))

      call check_diffusion( nf90_inq_varid(ncid, "threshold", threshold_varid))
      call check_diffusion( nf90_get_var(ncid, threshold_varid, threshold))

      call check_diffusion( nf90_inq_varid(ncid, "values_predicted", values_predicted_varid))
      call check_diffusion( nf90_get_var(ncid, values_predicted_varid, values_predicted))
    
    
      ! Close the file
      call check_diffusion( nf90_close(ncid))

      write(*, *) 'Finished reading regression diffusion file.'

      ! change from python to f90 convection for array indexing
      children_left = children_left + 1
      children_right = children_right + 1
      split_feature = split_feature + 1

      do_init=.false.
   end subroutine random_forest_init_diffusion


!#######################################################################

   subroutine random_forest_diffusion

!-----------------------------------------------------------------------
!
!       Random-forest subgrid parameterization
!
!-----------------------------------------------------------------------
!
!   input:  tabs               absolute temperature 
!           q                  total non-precipitating water
!           qp                 precipitating water
!
!   changes: t                 liquid static energy as temperature
!            q                 total non-precipitating water
!            qp                precipitating water
!
!-----------------------------------------------------------------------
!---------------------- local data -------------------------------------





   real,   dimension(nzm)             :: t_tendency, q_tendency, qp_tendency
   real,   dimension(nzm)             :: qsat
   real(4), dimension(n_features)     :: features
   real(4), dimension(n_outputs)      :: outputs

   integer  i, j, k, inode, itree, next_node
   real     omn

   if (do_init) call error_mesg_diffusion('random_forest_init has not been called.')

 !  if (.not. rf_uses_qp) then
 !   ! initialize precipitation 
 !   if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) precsfc(:,:)=0.
 !  end if


! mod random tree
!call random_number(rand_tree)


   do j=1,ny
    do i=1,nx

        ! Initialize variables
        features = 0.
        outputs = 0.
        t_tendency = 0.
        q_tendency = 0.
        qp_tendency = 0.

        ! Combine all features into one vector
        features(1:nrf) = real(tabs(i,j,1:nrf),4)


        if (rf_uses_rh) then
         ! generalized relative humidity is used as a feature
         do k=1,nzm
          omn = omegan(tabs(i,j,k))
          qsat(k) = omn*qsatw(tabs(i,j,k),pres(k))+(1.-omn)*qsati(tabs(i,j,k),pres(k))
         end do
         features(nrf+1:2*nrf) = real(q(i,j,1:nrf)/qsat(1:nrf),4)
        else
         ! non-precipitating water is used as a feature
         features(nrf+1:2*nrf) = real(q(i,j,1:nrf),4)
        endif
       features(2*nrf + 1 : 3*nrf) = real(u(i,j,1:input_ver_dim),4)
       features(3*nrf + 1 : 4*nrf) = real(v(i,j,1:input_ver_dim),4)



!       features(4*nrf+1) = sqrt(real(u(i,j,1),4)**2 + real(v(i,j,1),4) **2)

! mod feature y
!features(2*nrf) = real(abs(dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
   
 !       if (rf_uses_qp) then
 !        features(2*nrf+1:3*nrf) = real(qp(i,j,1:nrf),4)
 !       endif

        ! calculate predicted values using RF

        ! loop over trees in forest
        do itree=1,n_trees

! mod rand tree (also remove do loop and division by n_trees)
!itree = ceiling(rand_tree(i,j)*n_trees)

         ! set current node to root node
         inode = 1 

         ! loop as long as not at leaf node
         do while (children_left(itree,inode) .ne. children_right(itree,inode))
          if (features(split_feature(itree,inode)).le.threshold(itree,inode)) then
            next_node = children_left(itree,inode)
          else
            next_node = children_right(itree,inode)
          endif
          inode = next_node
          if (inode .gt. n_nodes(itree)) call error_mesg_diffusion('inconsistent tree')
         enddo

         outputs = outputs + values_predicted(itree, inode, :)
  
        enddo ! end loop over trees
        outputs = outputs/n_trees


        ! Separate out outputs into heating and moistening tendencies
        t_tendency(1:nrf) = outputs(1:nrf)
        q_tendency(1:nrf) = outputs(nrf+1:2*nrf)

        t(i,j,:) = t(i,j,:) + t_tendency*dtn
        q(i,j,:) = q(i,j,:) + q_tendency*dtn

        where (qn(i,j,:).gt.0.0)
         qn(i,j,:) = qn(i,j,:) + q_tendency*dtn
        end where

        where (qn(i,j,:).lt.0.0)
         qn(i,j,:) = 0.0
        end where

        ! For precipitation rate, note that output of prec_fall 
        ! includes a dtn/dz factor that is reversed
        ! by a dz/dt factor in statistics.f90 and write_fields2D.f90


      !  if (rf_uses_qp) then
      !   qp_tendency(1:nrf) = outputs(2*nrf+1:3*nrf)
      !   qp(i,j,:) = qp(i,j,:) + qp_tendency*dtn
      !   do k=1,nzm
      !      precsfc(i,j) = precsfc(i,j)-qp_tendency(k)*adz(k)*dz*rho(k)*(dtn/dz)
      !      prec_xy(i,j) = prec_xy(i,j)-qp_tendency(k)*adz(k)*dz*rho(k)*(dtn/dz)
      !   end do
      !  endif


       end do
     end do

   end subroutine random_forest_diffusion





!#######################################################################


!##############################################################################
  subroutine check_diffusion(status)

    ! checks error status after each netcdf, prints out text message each time
    !   an error code is returned. 

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
       write(*, *) trim(nf90_strerror(status))
    end if
  end subroutine check_diffusion

!##############################################################################
 subroutine error_mesg_diffusion (message)
  character(len=*), intent(in) :: message

!  input:
!      message   message written to output   (character string)

    if(masterproc) print*, 'Random forest module: ', message
    stop

 end subroutine error_mesg_diffusion



!#######################################################################

end module random_forest_param_diffusion

