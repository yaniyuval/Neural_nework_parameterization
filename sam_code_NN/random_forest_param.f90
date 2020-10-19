! random forest parameterization
!Adding surface wind as a variable. 
module random_forest_param


!----------------------------------------------------------------------
use netcdf
use vars
use grid


implicit none
private


!---------------------------------------------------------------------
!  ---- public interfaces ----

   public  random_forest, random_forest_init, check


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

    !Yani: Here I should allocate a vector for the list of sizes and names of outputs that I have (Only ones that I care about their outputs...)
    !Also I should allocate vector for the input size and name (can I just read strings ? )

    !flags - think if I want to work with flags - and how.... I need to think of features (gradient or values or both and wheather I want to have wind)... Could I read names of inputs and outputs instead ?  


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

   subroutine random_forest_init

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

      call task_rank_to_index(rank,it,jt)
  
!-------------allocate arrays and read data-------------------------


      
      if(masterproc)  write(*,*) rf_filename

! Open the file. NF90_NOWRITE tells netCDF we want read-only access
! Get the varid or dimid for each variable or dimension based on its name.

      call check( nf90_open(     trim(rf_filename),NF90_NOWRITE,ncid ))

      call check( nf90_inq_dimid(ncid, 'dim_trees', trees_dimid))
      call check( nf90_inquire_dimension(ncid, trees_dimid, len=n_trees))

      call check( nf90_inq_dimid(ncid, 'dim_features', features_dimid))
      call check( nf90_inquire_dimension(ncid, features_dimid, len=n_features))
      
     print *, 'size of features', n_features
        
      call check( nf90_inq_dimid(ncid, 'dim_outputs', outputs_dimid))
      call check( nf90_inquire_dimension(ncid, outputs_dimid, len=n_outputs))

      call check( nf90_inq_dimid(ncid, 'dim_nodes', nodes_dimid))
      call check( nf90_inquire_dimension(ncid, nodes_dimid, len=max_n_nodes))

      allocate(n_nodes(n_trees))

      call check( nf90_inq_varid(ncid, "n_nodes", n_nodes_varid))
      call check( nf90_get_var(ncid, n_nodes_varid, n_nodes))

! basic checks of the expected data structure
!      if(n_outputs /= n_features) call error_mesg('n_outputs not equal to n_features.')

      print *, 'size of features and  outputs is'
      print *, n_features,n_outputs 
        
      if(max_n_nodes /= maxval(n_nodes)) call error_mesg('max_n_nodes not equal to maximum of n_nodes')

      !if (rf_uses_qp) then
      ! nrf = n_features/3
      !else
      ! nrf = n_features/2
      !endif
      nrf = 48 ! Size in the vertical 
    
      allocate(children_left(n_trees, max_n_nodes))
      allocate(children_right(n_trees, max_n_nodes))
      allocate(split_feature(n_trees, max_n_nodes))
      allocate(threshold(n_trees, max_n_nodes))
      allocate(values_predicted(n_trees, max_n_nodes, n_outputs))
! mod random tree
!allocate(rand_tree(nx, ny))



      call check( nf90_inq_varid(ncid, "children_left", children_left_varid))
      call check( nf90_get_var(ncid, children_left_varid, children_left))
      
      call check( nf90_inq_varid(ncid, "children_right", children_right_varid))
      call check( nf90_get_var(ncid, children_right_varid, children_right))

      call check( nf90_inq_varid(ncid, "split_feature", split_feature_varid))
      call check( nf90_get_var(ncid, split_feature_varid, split_feature))

      call check( nf90_inq_varid(ncid, "threshold", threshold_varid))
      call check( nf90_get_var(ncid, threshold_varid, threshold))

      call check( nf90_inq_varid(ncid, "values_predicted", values_predicted_varid))
      call check( nf90_get_var(ncid, values_predicted_varid, values_predicted))
    
    
      ! Close the file
      call check( nf90_close(ncid))

      write(*, *) 'Finished reading regression file.'

      ! change from python to f90 convection for array indexing
      children_left = children_left + 1
      children_right = children_right + 1
      split_feature = split_feature + 1

      do_init=.false.
   end subroutine random_forest_init


!#######################################################################

   subroutine random_forest

!-----------------------------------------------------------------------
!
!       Random-forest subgrid parameterization
!
!-----------------------------------------------------------------------
!
!   input:  tabs               absolute temperature 
!           q                  total non-precipitating water
!           qp                 precipitating water
!           u                  zonal wind
!           v                  meridional wind      
!   changes: t                 liquid static energy as temperature
!            q                 total non-precipitating water
!            qp                precipitating water
!
!-----------------------------------------------------------------------
!---------------------- local data -------------------------------------





   real,   dimension(nzm)             :: t_tendency, q_tendency, qp_tendency
   real,   dimension(nzm)             ::  t_rad_tend ! In case we predict radiation tendency seperately. 
   real,   dimension(nzm)             :: qsat
   real(4), dimension(n_features)     :: features
   real(4), dimension(n_outputs)      :: outputs

   integer  i, j, k, inode, itree, next_node, dim_counter, out_dim_counter
   real     omn

   if (do_init) call error_mesg('random_forest_init has not been called.')

   if (.not. rf_uses_qp) then
    ! initialize precipitation 
    if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) precsfc(:,:)=0.
   end if


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
	dim_counter = 0
        ! Combine all features into one vector
        if (Tin_feature_rf) then
         features(dim_counter+1:dim_counter + input_ver_dim) = real(t_i(i,j,1:input_ver_dim),4)
	 dim_counter = dim_counter + input_ver_dim
	endif
	if (qin_feature_rf) then
        if (rf_uses_rh) then
         ! generalized relative humidity is used as a feature
         do k=1,nzm
          omn = omegan(tabs(i,j,k))
          qsat(k) = omn*qsatw(tabs(i,j,k),pres(k))+(1.-omn)*qsati(tabs(i,j,k),pres(k))
         end do
         features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim)/qsat(1:input_ver_dim),4)
         dim_counter = dim_counter + input_ver_dim
	else
         ! non-precipitating water is used as a feature
         features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim),4)
         dim_counter =  dim_counter + input_ver_dim
        endif
        endif
	if (rf_uses_qp) then
         features(dim_counter+1:dim_counter+input_ver_dim) = real(qp_i(i,j,1:input_ver_dim),4)
         dim_counter = dim_counter + input_ver_dim
        endif



!JY added gradient
	if (Tin_z_grad_feature_rf) then
	  features(dim_counter + 1) = real(t_i(i,j,1),4) 
	  dim_counter = dim_counter + 1
	  do k=1,nzm-1
	    features(dim_counter + 1) = real(t_i(i,j,k+1),4) - real(t_i(i,j,k),4) 
	    dim_counter = dim_counter + 1
	  end do
	endif 
	if (qin_z_grad_feature_rf) then
          features(dim_counter + 1) = real(q_i(i,j,1),4)              
          dim_counter = dim_counter + 1
          do k=1,nzm-1
            features(dim_counter + 1) = real(q_i(i,j,k+1),4) - real(q_i(i,j,k),4)                            
            dim_counter = dim_counter + 1
          end do
	endif
        if (do_wind_input_rf) then  !NOTE THAT THIS FLAG IS ONLY FOR HORIZONTAL WIND 
            features(dim_counter+1:dim_counter+input_ver_dim) = real(u_i(i,j,1:input_ver_dim),4)
            dim_counter = dim_counter + input_ver_dim
            features(dim_counter+1:dim_counter+input_ver_dim) = real(v_i(i,j,1:input_ver_dim),4)
            dim_counter = dim_counter + input_ver_dim
        endif
        if (do_vert_wind_input_rf) then !NOTE THAT THIS FLAG IS FOR VERTICAL WIND 
            features(dim_counter+1:dim_counter+input_ver_dim) = real(w_i(i,j,1:input_ver_dim),4) !Need to think which dimentions should go here since w has one more dim...Also I need to insert w from the beginning of the time step! before it was changed...  
            dim_counter = dim_counter + input_ver_dim
        endif
        if (do_surf_wind_rf) then
         features(dim_counter+1) = sqrt(real(u_i(i,j,1),4)**2 + real(v_i(i,j,1),4) **2)
         dim_counter = dim_counter + 1
	endif

! mod feature y
        if(do_yin_input) then
         features(dim_counter+1) = real(abs(dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
         dim_counter = dim_counter+1
        endif

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
          if (inode .gt. n_nodes(itree)) call error_mesg('inconsistent tree')
         enddo

         outputs = outputs + values_predicted(itree, inode, :)
  
        enddo ! end loop over trees
        outputs = outputs/n_trees


        ! Separate out outputs into heating and moistening tendencies
        t_tendency(1:nrf) = outputs(1:nrf)
        out_dim_counter = nrf
        q_tendency(1:nrf) = outputs(out_dim_counter+1:out_dim_counter + nrf)
        out_dim_counter = out_dim_counter + nrf
        if (rf_uses_qp) then
         qp_tendency(1:nrf) = outputs(out_dim_counter+1:out_dim_counter +nrf)
         qp(i,j,:) = qp(i,j,:) + qp_tendency*dtn
         out_dim_counter = out_dim_counter + nrf
         do k=1,nzm
            precsfc(i,j) = precsfc(i,j)-qp_tendency(k)*adz(k)*dz*rho(k)*(dtn/dz)
            prec_xy(i,j) = prec_xy(i,j)-qp_tendency(k)*adz(k)*dz*rho(k)*(dtn/dz)
         end do
        endif


        if (do_radiation_output_rf) then
         t_tendency(1:rad_lev_pred) = t_tendency(1:rad_lev_pred) + outputs(out_dim_counter+1:out_dim_counter+rad_lev_pred)
         out_dim_counter = out_dim_counter + rad_lev_pred !before I had simulations that included zeros in levels I didn't predict. To save space I don't save these zeros. 
        endif 

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

        do k=1, nzm
            precsfc(i,j) = precsfc(i,j)-q_tendency(k)*adz(k)*dz*rho(k)*(dtn/dz)
            prec_xy(i,j) = prec_xy(i,j)-q_tendency(k)*adz(k)*dz*rho(k)*(dtn/dz)
        end do
      !  if (do_q_surf_fluxes_out_rf) then
      !   precsfc(i,j) = precsfc(i,j) + outputs(out_dim_counter+1)*adz(1)*dz*rho(1)*(dtn/dz)
      !   prec_xy(i,j) = prec_xy(i,j) + outputs(out_dim_counter+1)*adz(1)*dz*rho(1)*(dtn/dz)        
      !  endif 


       end do
     end do

   end subroutine random_forest





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

!##############################################################################
 subroutine error_mesg (message)
  character(len=*), intent(in) :: message

!  input:
!      message   message written to output   (character string)

    if(masterproc) print*, 'Random forest module: ', message
    stop

 end subroutine error_mesg



!#######################################################################

end module random_forest_param

