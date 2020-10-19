
! random forest parameterization

module random_forest_param_tkh


!----------------------------------------------------------------------
use netcdf
use vars
use grid


implicit none
private


!---------------------------------------------------------------------
!  ---- public interfaces ----

   public  random_forest_tkh, random_forest_init_tkh, check_tkh


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

   subroutine random_forest_init_tkh

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
!character(len=256) :: rf_filename_tkh
      call task_rank_to_index(rank,it,jt)
  
!-------------allocate arrays and read data-------------------------


      
      if(masterproc)  write(*,*) rf_filename_tkh

!      rf_filename_tkh = '/glade/work/janniy/mldata/gcm_regressors/qobsFFTFTTTFF0FFTFTF48TFFFFFFFFF_F-NoSc_O-Stan_Ntr4277070_Nte225180_F_Tin_qin_uin_vin_usurf_O_Tout_qout_RF_NTr10_MinS10max_d27_maxzinf_nocos_te83_tr87.nc'
! Open the file. NF90_NOWRITE tells netCDF we want read-only access
! Get the varid or dimid for each variable or dimension based on its name.

      call check_tkh( nf90_open(     trim(rf_filename_tkh),NF90_NOWRITE,ncid ))

      call check_tkh( nf90_inq_dimid(ncid, 'dim_trees', trees_dimid))
      call check_tkh( nf90_inquire_dimension(ncid, trees_dimid, len=n_trees))

      call check_tkh( nf90_inq_dimid(ncid, 'dim_features', features_dimid))
      call check_tkh( nf90_inquire_dimension(ncid, features_dimid, len=n_features))

      call check_tkh( nf90_inq_dimid(ncid, 'dim_outputs', outputs_dimid))
      call check_tkh( nf90_inquire_dimension(ncid, outputs_dimid, len=n_outputs))

      call check_tkh( nf90_inq_dimid(ncid, 'dim_nodes', nodes_dimid))
      call check_tkh( nf90_inquire_dimension(ncid, nodes_dimid, len=max_n_nodes))

      allocate(n_nodes(n_trees))

      call check_tkh( nf90_inq_varid(ncid, "n_nodes", n_nodes_varid))
      call check_tkh( nf90_get_var(ncid, n_nodes_varid, n_nodes))

      if(max_n_nodes /= maxval(n_nodes)) call error_mesg_tkh('max_n_nodes not equal to maximum of n_nodes')

      nrf = diffusivity_levels_rf

      allocate(children_left(n_trees, max_n_nodes))
      allocate(children_right(n_trees, max_n_nodes))
      allocate(split_feature(n_trees, max_n_nodes))
      allocate(threshold(n_trees, max_n_nodes))
      allocate(values_predicted(n_trees, max_n_nodes, n_outputs))
! mod random tree
!allocate(rand_tree(nx, ny))
	


print *, 'size of features and  outputs in tkh RF is:'
      print *, n_features,n_outputs

      call check_tkh( nf90_inq_varid(ncid, "children_left", children_left_varid))
      call check_tkh( nf90_get_var(ncid, children_left_varid, children_left))
      
      call check_tkh( nf90_inq_varid(ncid, "children_right", children_right_varid))
      call check_tkh( nf90_get_var(ncid, children_right_varid, children_right))

      call check_tkh( nf90_inq_varid(ncid, "split_feature", split_feature_varid))
      call check_tkh( nf90_get_var(ncid, split_feature_varid, split_feature))

      call check_tkh( nf90_inq_varid(ncid, "threshold", threshold_varid))
      call check_tkh( nf90_get_var(ncid, threshold_varid, threshold))

      call check_tkh( nf90_inq_varid(ncid, "values_predicted", values_predicted_varid))
      call check_tkh( nf90_get_var(ncid, values_predicted_varid, values_predicted))
    
    
      ! Close the file
      call check_tkh( nf90_close(ncid))

      write(*, *) 'Finished reading regression tkh file.'

      ! change from python to f90 convection for array indexing
      children_left = children_left + 1
      children_right = children_right + 1
      split_feature = split_feature + 1

      do_init=.false.
   end subroutine random_forest_init_tkh


!#######################################################################

   subroutine random_forest_tkh

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

   integer  i, j, k, inode, itree, next_node, out_counter, ic, jc
   real     omn, dtndxvar, dtndyvar,lat_v

   if (do_init) call error_mesg_tkh('random_forest_init has not been called.')

 !  if (.not. rf_uses_qp) then
 !   ! initialize precipitation 
 !   if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) precsfc(:,:)=0.
 !  end if


! mod random tree
!call random_number(rand_tree)


   do j=1,ny !Modified to try to calculate tkh_x tkh_y at additional grid point
    jc = j + 1
    lat_v = real((dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
    do i=1,nx
        ic = i + 1
        ! Initialize variables
        features = 0.
        outputs = 0.
        t_tendency = 0.
        q_tendency =0.
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
         !features(2*nrf + k) = real((u(i,j,k))**2,4)
         !features(3*nrf + k) = real((v(i,j,k))**2,4)
!        ! if (i.eq.1.and.j.eq.1.and.k.eq.1) then
!        !  print *, 'v v_i:', v_i(i,j,k), v(i,j,k)*dtndyvar, dtndyvar, rho(k),v(i,j,k)
!        !  print *, 'u u_i:', u_i(i,j,k), u(i,j,k)*dtndyvar
!        ! end if 
        end do

!       features(2*nrf + 1 : 3*nrf) = real(u_i(i,j,1:nrf),4)
!       features(3*nrf + 1 : 4*nrf) = real(v_i(i,j,1:nrf),4)
       dtndyvar = 1/rho(1)/dtndxvar
       features(4*nrf+1) = sqrt(real(0.5*(u(i,j,1)* dtndyvar +u(ic,j,1)* dtndyvar) ,4)**2 + real(0.5*(v(i,j,1)* dtndyvar + v(i,jc,1)*dtndyvar),4) **2)
       ! features(4*nrf+1) = sqrt(real((u(i,j,1)) ,4)**2 + real((v(i,j,1)),4) **2)       
!features(4*nrf + 1 : 5*nrf) = real(w(i,j,1:input_ver_dim),4)
! mod feature y
      ! features(4*nrf+2) = real(abs(dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
   features(4*nrf+2) = abs(lat_v)
         
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
          if (inode .gt. n_nodes(itree)) call error_mesg_tkh('inconsistent tree')
         enddo
         
         outputs = outputs + values_predicted(itree, inode, :)
  
        enddo ! end loop over trees
        outputs = outputs/n_trees


        ! Separate out outputs into heating and moistening tendencies

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
       end do
     end do

   end subroutine random_forest_tkh





!#######################################################################


!##############################################################################
  subroutine check_tkh(status)

    ! checks error status after each netcdf, prints out text message each time
    !   an error code is returned. 

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
       write(*, *) trim(nf90_strerror(status))
    end if
  end subroutine check_tkh

!##############################################################################
 subroutine error_mesg_tkh (message)
  character(len=*), intent(in) :: message

!  input:
!      message   message written to output   (character string)

    if(masterproc) print*, 'Random forest module: ', message
    stop

 end subroutine error_mesg_tkh



!#######################################################################

end module random_forest_param_tkh

