module co2
    !--------------------------------------------------------------------------
    ! TODO
    !--------------------------------------------------------------------------

    use constants
    use namelist_define
    use mapinfo_define
    use obs_define
    use mpi_module
    use map_utils
    use netcdf
    use wrf_tools

    implicit none

contains
    !--------------------------------------------------------------------------
    ! height_to_k
    !
    ! Index of closest model level corresponding to height.
    !
    !--------------------------------------------------------------------------
    function height_to_k(ph, height)
        implicit none

        real                                :: height_to_k
        real, dimension(:), intent(in)      :: ph
        real, intent(in)                    :: height

        integer                             :: k, kx, k1
        real, dimension(:), allocatable     :: ph_unstaggered
        real                                :: height_lower, height_upper
        real                                :: alpha

        kx = size(ph)
        allocate(ph_unstaggered(kx-1))

        ph_unstaggered = (ph(1:kx-1) + ph(2:kx)) / 2.

        ! Find closest vertical level below 'height'
        do k = 1, kx-1
            if (ph_unstaggered(k) < height) then
                k1 = k
            end if
        end do

        ! Linearly interpolate to find k at 'height'
        height_lower = ph_unstaggered(k1)
        height_upper = ph_unstaggered(k1+1)

        alpha = (height - height_lower)/(height_upper - height_lower)
        height_to_k = k1 + alpha

        if (height_to_k < 1) then
            ! TODO: Raise a warning
            height_to_k = 1
        end if

        ! write(*, *) height, height_lower, height_upper, k1, height_to_k

        return

    end function height_to_k

    !--------------------------------------------------------------------------
    ! get_co2_tower_obs
    !
    ! Read CO2 tower files for the current time into raw%co2_tower variables.
    ! 
    ! Input files:
    !   tower_%Y%m%d%H%M%S.dat (e.g. tower_20170503103800.dat)
    !
    !--------------------------------------------------------------------------
    subroutine get_co2_tower_obs(times, ph, proj)
        implicit none

        character(len=80), intent(in)               :: times
        real, dimension(:,:,:), intent(in)          :: ph
        type(proj_info), intent(in)                 :: proj

        integer                                     :: n
        character (len=80)                          :: input_file
        integer                                     :: ierr
        integer                                     :: nlines

        character (len=80)                          :: tower_name
        real                                        :: lat, lon, elev, co2
        integer                                     :: iwrf, jwrf, kwrf

        real                                        :: aio, ajo, ako
        integer                                     :: io, jo, ko

        input_file = 'tower_' // times(1:4) // times(6:7) // times(9:10) // &
                     times(12:13) // times(15:16) // times(18:19) // ".dat"
        open(10, file=trim(input_file), status='old', form='formatted', &
             iostat=ierr)
        if (ierr /= 0) then
            if (my_proc_id == 0) then
                write(*, *) trim(input_file), ' does not exist.'
                write(*, *) 'Tower CO2 concentration data will not be &
                            &assimilated!'
            end if
        end if

        if (my_proc_id == 0) then
            write(*, *) ''
            write(*, *) '---------------------------------------------------'
            write(*, *) '.... Getting tower CO2 data ....'
        end if

        nlines = 0
        do
            read(10, *, iostat=ierr)
            if (ierr /= 0) then
                exit
            end if
            nlines = nlines + 1
        end do
        raw%co2_tower%num = nlines

        allocate(raw%co2_tower%latitude(raw%co2_tower%num))
        allocate(raw%co2_tower%longitude(raw%co2_tower%num))
        allocate(raw%co2_tower%elevation(raw%co2_tower%num))
        allocate(raw%co2_tower%co2(raw%co2_tower%num))
        allocate(raw%co2_tower%ii(raw%co2_tower%num))
        allocate(raw%co2_tower%jj(raw%co2_tower%num))
        allocate(raw%co2_tower%kk(raw%co2_tower%num))

        rewind(10)
        do n = 1, raw%co2_tower%num
            read(10, '(a5, f10.4, f10.4, f7.1, i4, i4, i3, f9.4)', &
                 iostat=ierr) tower_name, lat, lon, elev, jwrf, iwrf, kwrf, co2

            raw%co2_tower%latitude(n) = lat
            raw%co2_tower%longitude(n) = lon
            raw%co2_tower%elevation(n) = elev
            raw%co2_tower%co2(n) = co2

            ! Convert from 0-based indexing to 1-based indexing
            iwrf = iwrf + 1
            jwrf = jwrf + 1
            kwrf = kwrf + 1

            call latlon_to_ij(proj, lat, lon, aio, ajo)
            io = nint(aio)
            jo = nint(ajo)

            ako = height_to_k(ph(io,jo,:), elev)
            ko = nint(ako)

            ! Real data
            ! raw%co2_tower%ii(n) = aio
            ! raw%co2_tower%jj(n) = ajo
            ! raw%co2_tower%kk(n) = ako

            ! OSSE
            raw%co2_tower%ii(n) = iwrf
            raw%co2_tower%jj(n) = jwrf
            raw%co2_tower%kk(n) = kwrf

            ! Diagnostics
            ! write(*, *) trim(tower_name)
            ! write(*, *) io, iwrf
            ! write(*, *) jo, jwrf
            ! write(*, *) ko, kwrf
            ! write(*, *) lat
            ! write(*, *) lon
            ! write(*, *) elev
            ! write(*, *) co2

        end do
        close(10)

        return

    end subroutine get_co2_tower_obs

    !--------------------------------------------------------------------------
    ! sort_co2_tower_data
    !
    ! Prepare CO2 tower data to be assimilated.
    !
    !--------------------------------------------------------------------------
    subroutine sort_co2_tower_data(wrf_file, ix, jx, kx, proj, instrument, &
                                   datathin, hroi, vroi, grid_id )
        implicit none

        real, parameter                      :: CO2_ERROR = 1
        character(len=10), intent(in)        :: wrf_file
        integer, intent(in)                  :: ix, jx, kx
        type(proj_info), intent(in)          :: proj
        character(len=8), intent(in)         :: instrument
        integer, intent(in)                  :: datathin, hroi, vroi, grid_id

        integer                              :: n

        do n = 1, raw%co2_tower%num
            if (raw%co2_tower%co2(n) > 0) then
                obs%num                 = obs%num + 1
                obs%dat     (obs%num  ) = raw%co2_tower%co2(n)
                obs%type    (obs%num  ) = 'co2tower  '
                obs%err     (obs%num  ) = CO2_ERROR
                obs%position(obs%num,1) = raw%co2_tower%ii(n)
                obs%position(obs%num,2) = raw%co2_tower%jj(n)
                obs%position(obs%num,3) = raw%co2_tower%kk(n)
                obs%roi     (obs%num,1) = hroi
                obs%roi     (obs%num,2) = vroi
            endif
        end do

        return

    end subroutine sort_co2_tower_data

    !--------------------------------------------------------------------------
    ! xb_to_co2tower
    !
    ! Translate background state vector to CO2 concentration value at tower
    ! locations.
    !
    ! The tower value is obtained through linear interpolation horizontally and
    ! vertically.
    !
    !--------------------------------------------------------------------------
    subroutine xb_to_co2tower(inputfile, xa, ix, jx, kx, nv, iob, xb)
        implicit none
        character(len=10), intent(in)           :: inputfile
        character(len=10)                       :: obstype
        integer, intent(in)                     :: ix, jx, kx, nv, iob
        real, dimension(3,3,kx+1,nv), intent(in) :: xa  ! 1st guest
        real, intent(out)                       :: xb

        real, dimension(ix,jx,kx)               :: co2
        real, dimension(kx)                     :: co2_vert
        integer                                 :: m
        integer                                 :: i_co2 = 0
        integer                                 :: obs_ii, obs_jj, obs_kk
        integer                                 :: i1, j1, k1
        real                                    :: dx, dxm, dy, dym, dz, dzm

        obstype = obs%type(iob)
        obs_ii = obs%position(iob,1)
        obs_jj = obs%position(iob,2)
        obs_kk = obs%position(iob,3)

        write(*,*) obs_ii, obs_jj, obs_kk

        i1 = int(obs_ii)
        j1 = int(obs_jj)
        k1 = int(obs_kk)

        dx = obs_ii - real(i1)
        dy = obs_jj - real(j1)
        dz = obs_kk - real(k1)
        dxm = real(i1+1) - obs_ii
        dym = real(j1+1) - obs_jj
        dzm = real(k1+1) - obs_kk

        do m = 1, nv
            if (trim(enkfvar(m)) == 'CO2') then
                i_co2 = m
            end if
        end do

        if (i_co2 > 0) then
            co2(i1:i1+1,j1:j1+1,1:kx) = xa(1:2,1:2,1:kx,i_co2)
        else
            call get_variable3d(inputfile, 'CO2', ix, jx, kx, 1, co2)
        end if

       co2_vert(:) = dym*(dx*co2(i1+1,j1,:) + dxm*co2(i1,j1,:)) + &
                     dy*(dx*co2(i1+1,j1+1,:) + dxm*co2(i1,j1+1,:))

        xb = dzm*co2_vert(k1) + dz*co2_vert(k1+1)

        return

    end subroutine xb_to_co2tower

end module co2
