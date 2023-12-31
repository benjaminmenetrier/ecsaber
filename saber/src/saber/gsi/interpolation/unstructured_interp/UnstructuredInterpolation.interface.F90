! (C) Copyright 2020- UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module saber_interpolator_unstrc_interface

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding
use kinds
use missing_values_mod
use saber_unstructured_interpolation_mod

private
! ------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Create Unstructured Interpolator from lonlat
!!
subroutine saber_unstrc_create_c(c_key_unstrc, c_comm, c_lonlat1, &
                           c_lonlat2, c_config) bind(c, name='saber_unstrc_create_f90')
implicit none

! Passed variables
integer(c_int), intent(inout)  :: c_key_unstrc !< bump interpolator
type(c_ptr), value, intent(in) :: c_comm       !< MPI Communicator
type(c_ptr), intent(in),value  :: c_lonlat1    !< source grid (atlas field)
type(c_ptr), intent(in),value  :: c_lonlat2    !< target grid (atlas field)
type(c_ptr), value, intent(in) :: c_config     !< Configuration

! local variables
type(saber_unstrc_interp), pointer :: unstrc_int
type(fckit_mpi_comm) :: f_comm
type(fckit_configuration) :: f_config
type(atlas_field) :: lonlat_field1, lonlat_field2
real(kind_real), pointer :: lonlat1(:,:), lonlat2(:,:)
real(kind_real), allocatable :: lon1(:), lat1(:)
real(kind_real), allocatable :: lon2(:), lat2(:)
integer :: ngrid1, ngrid2, nnearest
character(len=:), allocatable :: us_interp_type, infile
logical :: read_from_file

f_comm = fckit_mpi_comm(c_comm)
f_config = fckit_configuration(c_config)

! check for specification of weight types and number of
! neighbors in the configuration.  If absent, implement defaults

call saber_unstrc_interp_registry%init()
call saber_unstrc_interp_registry%add(c_key_unstrc)
call saber_unstrc_interp_registry%get(c_key_unstrc, unstrc_int)

! This may need to be revised for multiple MPI tasks
lonlat_field1 = atlas_field(c_lonlat1)
call lonlat_field1%data(lonlat1)
ngrid1 = size(lonlat1(1,:))

lonlat_field2 = atlas_field(c_lonlat2)
call lonlat_field2%data(lonlat2)
ngrid2 = size(lonlat2(1,:))

if (f_config%has("read_from_file")) then
    call f_config%get_or_die("read_from_file", read_from_file)
else
    read_from_file = .false.
endif

if (read_from_file) then
    call fckit_log%info("Reading interpolator from file")
    call f_config%get_or_die("infile", infile)
    call unstrc_int%create(f_comm, infile)

    if ((unstrc_int%ngrid_in /= ngrid1) .or. (unstrc_int%ngrid_out /= ngrid2)) &
        call abor1_ftn('unstructured interpolator input file is not consistent with input and output functionspaces')

else

  us_interp_type = 'barycent'
  if (f_config%has("us_interp_type")) then
      call f_config%get_or_die("interpolation_type", us_interp_type)
  endif

  nnearest = 4
  if (f_config%has("nnearest")) then
      call f_config%get_or_die("nnearest", nnearest)
  endif

  ! For now copy the atlas fields into an unstructured grid
  ! in the future we can make this more efficient
  allocate(lon1(ngrid1), lat1(ngrid1))
  lon1 = lonlat1(1,:)
  lat1 = lonlat1(2,:)

  allocate(lon2(ngrid2), lat2(ngrid2))
  lon2 = lonlat2(1,:)
  lat2 = lonlat2(2,:)

  call unstrc_int%create(f_comm, nnearest, us_interp_type, ngrid1, lat1, lon1, &
                                                           ngrid2, lat2, lon2)

  deallocate(lon1, lat1, lon2, lat2)
endif

! release pointers
call lonlat_field1%final()
call lonlat_field2%final()

end subroutine saber_unstrc_create_c

!-------------------------------------------------------------------------------
!> Write Unstructured Interpolator to a file
!!
subroutine saber_unstrc_write_c(c_key_unstrc, c_config) bind(c, name='saber_unstrc_write_f90')
implicit none

! Passed variables
integer(c_int), intent(inout)  :: c_key_unstrc !< unstructured interpolator
type(c_ptr), value, intent(in) :: c_config     !< Configuration

! local variables
type(saber_unstrc_interp), pointer :: unstrc_int
type(fckit_configuration) :: f_config
character(len=:), allocatable :: outfile

f_config = fckit_configuration(c_config)

! read name of output file from configuration
call f_config%get_or_die("outfile",outfile)

call saber_unstrc_interp_registry%get(c_key_unstrc, unstrc_int)
call unstrc_int%write(outfile)

end subroutine saber_unstrc_write_c

!-------------------------------------------------------------------------------
!> Apply Unstructured Interpolator
!!
subroutine saber_unstrc_apply_c(c_key_unstrc, c_infield, c_outfield) bind(c, name='saber_unstrc_apply_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_unstrc  !< key to unstructured interpolator
type(c_ptr), intent(in), value :: c_infield  !< input field
type(c_ptr), intent(in), value :: c_outfield !< output field

! Local variables
type(saber_unstrc_interp), pointer :: unstrc_int
type(atlas_field) :: infield, outfield
real(kind_real), pointer :: fin_r1(:), fout_r1(:)
real(kind_real), pointer :: fin_r2(:,:), fout_r2(:,:)
integer :: jlev, jlev1, jlev2

call saber_unstrc_interp_registry%get(c_key_unstrc, unstrc_int)

infield = atlas_field(c_infield)
outfield = atlas_field(c_outfield)

jlev1 = infield%levels()
jlev2 = outfield%levels()

if (jlev1 .ne. jlev2) call abor1_ftn("interpolator_unstrc_interface.unstrc_apply_c: number"// &
                                     " of levels for the two fields does not match.")

if (jlev1 == 0) then

  call infield%data(fin_r1)
  call outfield%data(fout_r1)

  call unstrc_int%apply(fin_r1, fout_r1)

else

  call infield%data(fin_r2)
  call outfield%data(fout_r2)

  do jlev = 1, jlev1
     call unstrc_int%apply(fin_r2(jlev,:), fout_r2(jlev,:))
  enddo

endif



! free up pointers and arrays
call infield%final()
call outfield%final()

end subroutine saber_unstrc_apply_c

!-------------------------------------------------------------------------------
!> Apply Unstructured Interpolator Adjoint
!!
subroutine saber_unstrc_apply_ad_c(c_key_unstrc, c_field2, c_field1) bind(c, name='saber_unstrc_apply_ad_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_unstrc  !< key to unstructured interpolator
type(c_ptr), intent(in), value :: c_field2  !< input field
type(c_ptr), intent(in), value :: c_field1  !< output field

! Local variables
type(saber_unstrc_interp), pointer :: unstrc_int
type(atlas_field) :: field_grid2, field_grid1
real(kind_real), pointer :: f2_r1(:),   f1_r1(:)
real(kind_real), pointer :: f2_r2(:,:), f1_r2(:,:)
integer :: jlev, jlev1, jlev2

call saber_unstrc_interp_registry%get(c_key_unstrc, unstrc_int)

field_grid2 = atlas_field(c_field2)
field_grid1 = atlas_field(c_field1)

jlev1 = field_grid1%levels()
jlev2 = field_grid2%levels()

if (jlev1 .ne. jlev2) call abor1_ftn("interpolator_unstrc_interface.unstrc_apply_ad_c: number"// &
                                     " of levels for the two fields does not match.")

if (jlev1 == 0) then

  call field_grid1%data(f1_r1)
  call field_grid2%data(f2_r1)

  call unstrc_int%apply_ad(f1_r1, f2_r1)

else

  call field_grid1%data(f1_r2)
  call field_grid2%data(f2_r2)

  do jlev = 1, jlev1
     call unstrc_int%apply_ad(f1_r2(jlev,:), f2_r2(jlev,:))
  enddo

endif

! free up pointers
call field_grid2%final()
call field_grid1%final()

end subroutine saber_unstrc_apply_ad_c

!------------------------------------------------------------------------------
!> Delete unstructured_interpolator
subroutine saber_unstrc_delete_c(c_key_unstrc) bind(c, name='saber_unstrc_delete_f90')

implicit none

! Passed variables
integer(c_int), intent(inout) :: c_key_unstrc

! Local variables
type(saber_unstrc_interp), pointer :: unstrc_int

call saber_unstrc_interp_registry%get(c_key_unstrc, unstrc_int)

call unstrc_int%delete()

call saber_unstrc_interp_registry%remove(c_key_unstrc)

end subroutine saber_unstrc_delete_c

! ------------------------------------------------------------------------------

end module saber_interpolator_unstrc_interface
