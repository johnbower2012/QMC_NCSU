   program harmonicoscillator
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: i8=selected_int_kind(15)
   integer, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), allocatable, dimension(:) :: xold,xnew
   real(kind=r8) :: omega,dt,etrial,xtest
   integer(kind=i4) :: nwalk0,nwalk,nw,neq,nav,nstep,nb
   integer(kind=i4) :: i,j,k,n
   real(kind=r8) :: esum,e2sum,eblock,e2block,wsum,wblock
   real(kind=r8) :: eloc,sigma,weight,ene,ene2,err
   real(kind=r8) :: errblock,prob,psinew,psiold
   integer(kind=i4) :: iwt,idmc,accept,nblocks,nmoves
   real(kind=r8) :: rgauss,csi
   integer(kind=i8) :: irn
   real(kind=r8) :: alpha
   real(kind=r8) :: tau
   read(5,*) irn   ! number used to initialize the random number generator
   read(5,*) omega ! strenght of the harmonic potential
   read(5,*) alpha ! alpha parameter of the trial wave function, psi=exp(-alpha*x**2/2)
   read(5,*) idmc
   read(5,*) dt    ! time step
   read(5,*) etrial ! trial energy
   read(5,*) nwalk0 ! number of walkers
   read(5,*) neq   ! number of blocks for equilibration
   read(5,*) nav   ! number of blocks for average
   read(5,*) nstep ! number of steps per block
   write(6,'(''One dimensional harmonic oscillator'')')
   select case (idmc)
      case(1) 
         write(6,'(''Variational Monte Carlo'')')
      case(2)
         write(6,'(''Diffusion Monte Carlo'')')
      case default
         write(6,'(''Invalid idmc!!!'')')
         stop
   end select
   write(6,'(''Random number seed = '',t40,i20)') irn
   write(6,'(''Omega strenght of the potential = '',t50,f10.5)') omega
   write(6,'(''Alpha parameter of the trial wave function = '',t50,f10.5)') alpha
   write(6,'(''Time step = '',t50,f10.5)') dt
   write(6,'(''Trial energy = '',t50,f10.5)') etrial
   write(6,'(''Walkers = '',t50,i10)') nwalk0
   write(6,'(''Number of blocks for equilibration = '',t50,i10)') neq
   write(6,'(''Number of blocks for statistics = '',t50,i10)') nav
   write(6,'(''Number of steps for each block = '',t50,i10)') nstep
   write(6,*) 
   select case (idmc)
      case (1)
         write(6,'(''iteration    acceptance    block energy    error         total energy    error'')')
      case (2)
         write(6,'(''iteration    tau           energy block    error         total energy    error'')')
   end select
   sigma=sqrt(dt)
   allocate(xold(2*nwalk0),xnew(2*nwalk0))
! walkers initializations
! give a random position to each walker
   do i=1,nwalk0
      call ran(csi,irn)
      xold(i)=1.0_r8/omega*(0.5-csi)
   enddo
   nwalk=nwalk0
   esum=0.0_r8
   e2sum=0.0_r8
   nblocks=0
   wsum=0.0_r8
   nb=0
   tau=0.0_r8
   do i=1,nav+neq
      if (i.eq.neq+1) then
         esum=0.0_r8
         e2sum=0.0_r8
         nblocks=0
         wsum=0.0_r8
         nb=0
         write(6,*) 'Equilibration done!'
         write(6,*)
      endif
      nblocks=nblocks+1
      eblock=0.0_r8
      e2block=0.0_r8
      wblock=0.0_r8
      accept=0.0_r8
      nmoves=0
      do j=1,nstep
         nw=0
! loop over the walkers
         do n=1,nwalk
            select case (idmc)
               case (1) ! VMC step
                  psiold=exp(-0.5_r8*alpha*xold(n)**2)
                  call ran(csi,irn)
! do a step according to   R = R' + sigma*random (random=linear between -1/2 and 1/2)
                  xtest=xold(n)+sigma*(csi-0.5_r8)
                  psinew=exp(-0.5_r8*alpha*xtest**2)
                  prob=(psinew/psiold)**2
                  call ran(csi,irn)
                  if (csi.lt.prob) then ! the move is accepted!
                     accept=accept+1
                  else ! the move has been rejected
                     xtest=xold(n)
                  endif
                  weight=1.0_r8
                  eloc=-0.5_r8*((alpha*xtest)**2-alpha)+0.5_r8*omega**2*xtest**2
                  nmoves=nmoves+1
               case(2) ! DMC step
! do a step according to   R = R' + sigma**2*dpsi(R')/psi(R') + sigma*random (random=gaussian)
                  xtest=xold(n)-sigma**2*alpha*xold(n)+sigma*rgauss(irn)
                  eloc=-0.5_r8*((alpha*xtest)**2-alpha)+0.5_r8*omega**2*xtest**2
                  weight=exp(-dt*(eloc-etrial))
            end select
! update the sum of expectation values
            if (j.eq.nstep) then
               eblock=eblock+weight*eloc
               e2block=e2block+weight*eloc**2
               wblock=wblock+weight
            endif
            select case (idmc)
               case (1)
                  csi=0.0_r8
               case (2)
! do the population control
                  weight=weight*real(nwalk0)/real(nwalk)
                  call ran(csi,irn)
! iwt = [w + random] 
! note that iwt is an integer, and (weight+csi) will be automatically recasted with a
! truncation to the near integer number
            end select
            iwt=weight+csi
            do k=1,iwt
! save a number iwt of the n-th walker to the new array
               nw=nw+1
               if (nw.gt.2*nwalk0) then
                  write(6,*) 'too many walkers, stopping...'
                  stop
               endif
               xnew(nw)=xtest
            enddo
         enddo
         do n=1,nw
! copy the new stack into the old one
            xold(n)=xnew(n)
         enddo
         nwalk=nw
         nb=nb+1
         tau=tau+dt
      enddo
! compute all the expectation values and errors
      eblock=eblock/wblock
      e2block=e2block/wblock
      errblock=sqrt(abs(e2block-eblock**2)/wblock)
      esum=esum+eblock
      e2sum=e2sum+eblock**2
      ene=esum/nblocks
      ene2=e2sum/nblocks
      err=sqrt(abs(ene2-ene**2)/nblocks)
      wsum=wsum+wblock
! write the output
      select case (idmc)
         case (1) 
            write(6,'(i7,5f15.10)') nb,real(accept)/real(nmoves),eblock,errblock,ene,err
         case (2) 
            write(6,'(i7,5f15.10)') nb,tau,eblock,errblock,ene,err
      end select
    enddo
    end program harmonicoscillator

    function rgauss(irn)
    implicit none
    integer, parameter :: i8=selected_int_kind(15)
    integer, parameter :: r8=selected_real_kind(15,9)
    real(kind=r8) :: pi
    real(kind=r8) :: rgauss,x1,x2
    integer(kind=i8) :: irn
    call ran(x1,irn)
    call ran(x2,irn)
    pi=4.0_r8*atan(1.0_r8)
    rgauss=sqrt(-2.0_r8*log(x1))*cos(2.0_r8*pi*x2)
    return
    end function rgauss
    
    subroutine ran(rn,irn)
    implicit none
    integer, parameter :: i4=selected_int_kind(9)
    integer, parameter :: i8=selected_int_kind(15)
    integer, parameter :: r8=selected_real_kind(15,9)
    integer(kind=i8),  parameter :: mask24 = ishft(1_i8,24)-1
    integer(kind=i8),  parameter :: mask48 = ishft(1_i8,48_i8)-1_i8
    real(kind=r8),  parameter :: twom48=2.0_r8**(-48)
    integer(kind=i8),  parameter :: mult1 = 44485709377909_i8
    integer(kind=i8),  parameter :: m11 = iand(mult1,mask24)
    integer(kind=i8),  parameter :: m12 = iand(ishft(mult1,-24),mask24)
    integer(kind=i8),  parameter :: iadd1 = 96309754297_i8
    integer(kind=i8) :: irn
    real(kind=r8) :: rn
    integer(kind=i8) :: is1,is2
    is2 = iand(ishft(irn,-24),mask24)
    is1 = iand(irn,mask24)
    irn = iand(ishft(iand(is1*m12+is2*m11,mask24),24)+is1*m11+iadd1,mask48)
    rn = ior(irn,1_i8)*twom48
    return
    end subroutine ran
