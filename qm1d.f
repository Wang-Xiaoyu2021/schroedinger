      program qm1d
c     *****************************************************************
c     find eigenfunctions of one-dimensional Schroedinger:
c     discretize derivative and use LAPACK to solve eivenvalue problem
c     use uniform mesh: 0,dx,2*dx,2*dx,..., (NPTS-1)*dx
c     boundary conditions are coded into discretized derivative
c     to compile and run unter Unix:
c     > g77 qm1d.f && a.out
c     >>> no attempt has been made to make this code efficient <<<
c     >>> so don't use it for production! <<<         EK MPI-FKF XII'03
c     *****************************************************************
      implicit none
      integer :: fid
      character(100) :: tmp
      integer    NPTS,     NSTM
      parameter (NPTS=201, NSTM=20)
      integer    i,n,nst
      real*8     dx, v(NPTS), vmin
      real*8     diag(NPTS),subd(NPTS), ee(NPTS),ev(NPTS,NSTM), de
      integer    m, iwork(5*NPTS), ifail(NPTS), info
      real*8     abstol, work(5*NPTS), dlamch
c     -----------------------------------------------------------------
      fid = 40
c     -----------------------------------------------------------------
c --- specify units:
c     rewrite Schroedinger equation: v(x)=2m*V(x)/hbar^2; ee=2m*E/hbar^2
c
c     [    d^2       2m        ]           2m*E
c     [ - -----  + ------ V(x) ] phi(x) = ------  phi(x)
c     [   d x^2    hbar^2      ]          hbar^2
c
c --- discretization
      dx=10d0/dble(NPTS-1)
c --- define potential
      write(*,'("particle in a box")')
      nst=10
      do i=1,NPTS
        v(i)=0d0
      enddo
c --- discretized form of kinetic energy operator (incl. boundary cond.)
c     f''(x_i) \approx (f(x_{i-1})-2*f(x_i)+f(x_{i+1}))/dx^2
c     by dropping terms from outside the mesh, we have chosen boundary
c     conditions such that the phi vanishes outside the mesh
      do i=1,NPTS
        diag(i)=v(i)+2d0/dx**2
        subd(i)=-1d0/dx**2
      enddo
c --- call LAPACK to diagonalize tridiagonal matrix
      if(nst.lt.1) nst=1
      if(nst.gt.NSTM) then
        write(*,'(">>> increase NSTM to get more states")')
        nst=NSTM
      endif
      abstol=2d0*dlamch('s')
      call dstevx('v','i',NPTS,diag,subd, 0d0,1d0, 1,nst, abstol,
     &            m,ee,ev,NPTS,work,iwork,ifail,info)
      if(info.ne.0) stop '>>> diagonalization failed'
c --- write eigenenergies
      write(*,'("eigenenergies:")')
      do n=1,nst
        write(*,'(I4,F20.10)') n,ee(n)
      enddo
c --- output for gnuplot -----------------------------------------------
c     minimum of potential for plotting range
      vmin=1d10
      do i=1,NPTS
        if(v(i).lt.vmin) vmin=v(i)
      enddo
c     average spacing of energy levels (for adjusting scale of ev)
      de=(ee(nst)-ee(1))/dble(nst-1)
      !-----------------------------------------------------------------
      ! output
      !
      open(fid, file="qm1d.dat")
         write(fid, '(a15)', advance='no') "x"
         write(fid, '(a15)', advance='no') "v"
         do i=1, nst
            write(tmp, '(a,i0,a)') "E(",i,")"
            write(fid, '(a15)', advance='no') trim(tmp)
            write(tmp, '(a,i0,a)') "V(",i,")"
            write(fid, '(a15)', advance='no') trim(tmp)
         end do
         write(fid,*)
         do n=1, NPTS
            write(fid, '(f15.7)', advance='no') (n-1) * dx
            write(fid, '(f15.7)', advance='no') v(n)
            do i=1, nst
               write(fid, '(f15.7)', advance='no') ee(i)
               write(fid, '(f15.7)', advance='no') ev(n,i) + ee(i)
            end do
            write(fid,*)
         end do
      close(fid)
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      open(10,file='qm1d.gnu')
c --- set up plotting options
      write(10,'("set term post ''Helvetica'' 21")')
      write(10,'("set output ''qm1d.ps''")')
      write(10,'("set nokey")')
      write(10,'("#set data style lines")')
      write(10,'("set yrange [",F12.8,":",F12.8,"]")') Vmin,ee(m)+de
c --- write plot commands for potential and eigenfunctions
      de=0.15d0*sqrt(dble(NPTS))*de
      write(10,'("plot ''-'' 1, ",A1)') char(92)
      do n=1,nst-1
        write(10,'("''-'' u 1:(",F12.8,"+",F12.8,"*$2) 1,",A1)')
     .        ee(n),de,char(92)
      enddo
      write(10,'("''-'' u 1:(",F12.8,"+",F12.8,"*$2) 1")') ee(nst),de
c --- data to plot
      write(10,'("# potential")')
      do i=1,NPTS
        write(10,'(I6,E20.10)') i-1,v(i)
      enddo
      write(10,'("e")')
      do n=1,nst
        write(10,'("# eng=",E20.10)') ee(n)
        do i=1,NPTS
          write(10,'(I6,E20.10)') i-1,ev(i,n)
        enddo
        write(10,'(I6,E20.10)') 0,ev(1,n)
        write(10,'("e")')
      enddo
c     for plotting run
c     > gnuplot qm1d.gnu && gv qm1d.ps
      end

c --- LAPACK routines for diagonalization
c     download e.g. from www.netlib.org
c     if you have LAPACK/BLAS installed, link your libraries instead
      include 'dstevx.f' 
      include 'dstevx-blas.f'

