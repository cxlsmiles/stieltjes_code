program stieltjes
use interpolation
implicit none

integer :: npt
integer, parameter :: QR_K = selected_real_kind (32)
real (kind=QR_K), allocatable, dimension(:) :: e, g

integer :: nmin = 10, nmax = 30 ! according to Mueller-Plathe and Diercksen Stieltjes is inaccurate for n>=15
!real (kind=QR_K), dimension(0:nmax) :: sk
!real (kind=QR_K), dimension(nmax,nmax) ::  e1, g1, e2, g2
real (kind=QR_K), allocatable, dimension(:) :: sk
real (kind=QR_K), allocatable, dimension(:,:) ::  e1, g1
real (kind=QR_K) :: shift1 = 0.1d0
integer :: imax1, imax2
integer :: inmax, ishift
integer :: i, j
integer :: exit_cycle
character(len=60) :: fname

! Tsveta
real (kind=QR_K) :: g_
real (kind=QR_K) :: temp, pi

pi = dacos(-1d0)

! reads and sorts energy and matrix elements
call Get_Cmdline(nmin, nmax)

read(*,*)npt
allocate(e(npt),g(npt))

do i = 1, npt
  read(*,*)e(i),g(i)
enddo
 
if(sum(g)==0d0) then
 write(*,*)"All matrix elements are zero, I stop"
 stop
endif

call sort2(npt,e,g)
print*,e(1)

ishift = 2
    shift1 = ishift * 0.1d0
    write(fname, '(A8,I1,A4)')"gamma.sh",ishift,".dat"
    open(237,file=fname)


    do inmax = nmin, nmax
       allocate(sk(0:inmax))
       allocate(e1(inmax,inmax), g1(inmax,inmax))       
       
       ! computes the nmax moments of the matrix elements
        open(unit=10,file='moments.txt')
        sk(:) = 0q0
        do i = 0, inmax
           do j = 1, npt
              sk(i) = sk(i) + e(j)**i*g(j)
           enddo
         write(10,*)i,abs(sk(i))
        enddo
        close(10)
    
       ! "images" the matrix elements for two different shift values, may be used to evaluate the highest accurate order
        shift1 = shift1 + abs(e(1))
        e(:) = e(:) + shift1
        call imaging(npt,e,g,inmax,imax1,e1,g1)
        e(:) = e(:) - shift1
       
        do i = 1, imax1
         do  j = 1, i-1
          write(100+i,'(2(f20.16,1X))')e1(j,i)-shift1,g1(j,i)
         enddo
        enddo
       
        call interp(e1(:,imax1),g1(:,imax1),imax1-1,shift1,g_)  
        g_ = 2.0d0*pi*g_
        shift1 = shift1 - abs(e(1))
        !write(*,'(1X,A9,f8.3)')"Shift 1: ",shift1-abs(e(1))
        !write(*,'(I2,A3,E23.15)')inmax," ",g_
       
        ! March 03, 2016 Tsveta
        write(237,'(I3,A3,E23.15)')inmax," ",g_
       
        deallocate(sk,e1,g1)
    end do
close(237)
deallocate(e,g)

end program stieltjes

subroutine Get_cmdline(nmin, nmax)
implicit none
integer,intent(inout)   :: nmin, nmax
character(len=100)      :: arg
integer                 :: i

i=0
do while (i < command_argument_count())
  i=i+1
  call get_command_argument(i, arg)
  
  select case (arg)
  case ('-h')
     write(*,*)'Call as e.g. stieltjes -i 15 -n 35 < INPUT_FILE'
     stop 1
  case ('-i')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)nmin  
  case ('-n')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)nmax  
  case default
    write(*,*)'Invalid command line argument!'
    stop 1
  end select

end do

end subroutine Get_cmdline
