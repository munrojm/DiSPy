* Sample Fortran source code for reading the data file, PIR_data.txt
* This data file contains the 2011 version of the matrices for the
* physically irreducible representations of the 230 crystallographic space
* groups for rational k vectors in 3-dimensional space.
*****************************************************************************
* Routines
*  pir_data_read: read data from file, PIR_data.txt
*  pir_data_unpack_operatormatrix: get matrix for given operator in space group
*  pir_data_unpack_irmatrix_pointoperator: get point operator part of
*    IR matrix for given operator
*  pir_data_unpack_irmatrix_translation: get translation part of
*    IR matrix for given translation
*  pir_data_get_irmatrix:  get IR matrix for given space group operator
*  pir_data_get_transform2complex: get matrix that transforms PIR to
*    complex block-diagonal form
* Routines for testing data
*  pir_data_write: write data into file, PIR_data2.txt.  Compare PIR_data2.txt
*    with PIR_data.txt to determine if the data file was read correctly
*  pir_data_test_multiplication_table: check that the IR matrices have the
*    same multiplication table as the space group operators
* Additional routines used by the above routines
*  pir_data_constant: store 16-byte double precision numbers as 1-byte integers
*  pir_data_vadd: add two vectors with rational components
*  pir_data_vsub: subtract two vectors with rational components
*  pir_data_vmlt: multiply a vector of rational numbers by a matrix of integers
*  pir_data_opmlt: multiply two operator matrices
*  pir_data_dmatmlt: multiply two matrices of double precision numbers
*  pir_data_factor: remove greatest common factor from a set of integers
*****************************************************************************
      module pir_data_module
* storage area for data from table
* number of IRs
      integer ircount
      parameter (ircount=10294)
* for the ith IR:
* spacegroupnumber(i), space group number (1-230)
      integer spacegroupnumber(ircount)
* spacegroupsymbol(i), space group symbol
      character spacegroupsymbol(ircount)*10
* irsymbol(i), IR symbol (Miller-Love convention)
      character irsymbol(ircount)*8
* irdimension(i), dimension of IR
      integer irdimension(ircount)
* irtype(i), type of IR (1,2,3)
      integer irtype(ircount)
* kcount(i), number of k vectors in the star of k
      integer kcount(ircount)
* pmkcount(i), number of k vectors in the star of +/-k
      integer pmkcount(ircount)
* operatorcount(i), number of representative operators in space group
      integer operatorcount(ircount)
* kvectorpointer(i), pointer to k vectors in kvector
      integer kvectorpointer(ircount)
* kvector(m), k vectors stored sequentially beginning at
*   kvectorpointer(i).  Each k vector is stored as kvec(4,4), where
*   kvec(1:3,1) are the x,y,z components of the constant term,
*   kvec(1:3,2) are the x,y,z components of the alpha term,
*   kvec(1:3,3) are the x,y,z components of the beta term,
*   kvec(1:3,4) are the x,y,z components of the gamma term,
*   kvec(4,m) are the common denominators.
*   alpha,beta,gamma are free parameters for non-special k vectors
      integer*1 kvector(350000)
* kspecial(i), true if special k vector (k point of symmetry)
      logical kspecial(ircount)

* for the jth space group operator:
* operatormatrixpointer(j,i), pointer to operator matrix in
*   operatormatrix
      integer operatormatrixpointer(48,ircount)
* operatormatrix(m), operator matrix stored sequentially beginning
*   at operatormatrixpointer(j,i).  This matrix is augmented with
*   the point operation in A(1:3,1:3), the translation in A(1:3,4),
*   and a common denominator in A(4,4).
      integer*1 operatormatrix(2000000)
* irtranslation(k,j,i), kth component (x,y,z) of the translational
*   part of the IR matrix.  This is usually, but not always, the
*   same as the translational part of the operator.
*   irtranslation(4,j,i) is a common denominator.
      integer irtranslation(4,48,ircount)
* irmatrixpointer(j,i), pointer to IR matrix in irmatrix
      integer irmatrixpointer(48,ircount)
* irmatrix(m), IR matrix stored sequentially beginning at
*   irmatrixpointer(j,i).  The elements of this matrix are stored
*   as integers which are mapped onto double precision values using
*   the routine, pir_data_constant.
      integer*1 irmatrix(9000000)
      end module pir_data_module


      subroutine pir_data_read
* Read data file, PIR_data.txt

      use pir_data_module

      implicit none
      integer i,j,k,m,kp,op,ip,irnumber
      double precision irtemp(2304)

* initialize data
      spacegroupnumber=0
      spacegroupsymbol=' '
      irsymbol=' '
      irdimension=0
      pmkcount=0
      operatorcount=0
      kvectorpointer=0
      kspecial=.true.
      kvector=0
      operatormatrix=0
      irmatrix=0
* initialize pointers
      kp=0
      op=0
      ip=0
* open data file
      open(30,file='PIR_data.txt')
* skip title line
      read(30,*)
      read(30,*)
      read(30,*)
* read next IR
      do while(.true.)
* read IR number, space group number, space group symbol, IR label,
* IR dimension, number of k vectors in the star of +/-k,
* number of operators in space group
        read(30,*,end=1)irnumber,spacegroupnumber(irnumber),
     $       spacegroupsymbol(irnumber),irsymbol(irnumber),
     $       irdimension(irnumber),irtype(irnumber),kcount(irnumber),
     $       pmkcount(irnumber),operatorcount(irnumber)
* read each k vector in the star of +/-k
        kvectorpointer(irnumber)=kp+1
        read(30,*)kvector(kp+1:kp+16*pmkcount(irnumber))
* is k vector special?
        kspecial(irnumber)=.true.
        do k=1,3
        do m=1,3
          if(kvector(kp+4*k+m).ne.0)kspecial(irnumber)=.false.
        enddo
        enddo
        kp=kp+16*pmkcount(irnumber)
* do each operator in space group
        do i=1,operatorcount(irnumber)
* read operator in extended matrix form
          operatormatrixpointer(i,irnumber)=op+1
          read(30,*)operatormatrix(op+1:op+16)
          op=op+16
* if nonspecial k vector, read contribution of translation to IR
          if(.not.kspecial(irnumber))then
            read(30,*)irtranslation(1:4,i,irnumber)
          endif
* read IR matrix
          irmatrixpointer(i,irnumber)=ip+1
          read(30,*)irtemp(1:irdimension(irnumber)**2)
          do j=1,irdimension(irnumber)**2
            call pir_data_constant(2,irmatrix(ip+j),irtemp(j))
          enddo
          ip=ip+irdimension(irnumber)**2
* next operator
        enddo
* next IR
      enddo
 1    close(30)
      end
**************************************************************************
      subroutine pir_data_constant(choice,n,x)
* store 16-byte double precision numbers as 1-byte integers
* arguments:
*   choice (input): choice=1, input n, output x
*                   choice=2, input x, output n
*   n (input/output): 1-byte integer
*   x (input/output): 16-byte double precision
      implicit none
      integer choice,i
      integer*1 n
      double precision x,constants(17)
      data constants/0,1,-1,0.5,-0.5,0.25,-0.25,
     $     0.866025403784439,-0.866025403784439,  ! sqrt(3)/2
     $     0.707106781186548,-0.707106781186548,  ! sqrt(2)/2
     $     0.433012701892219,-0.433012701892219,  ! sqrt(3)/4
     $     0.683012701892219,-0.683012701892219,  ! cos(15)/sqrt(2)
     $     0.183012701892219,-0.183012701892219/  ! sin(15)/sqrt(2)
* choice=1: input integer and output double precision
      if(choice.eq.1)then
        x=constants(n)
* choice=2: input double precision, output integer
      else
        do i=1,17
          if(dabs(x-constants(i)).lt.1d-4)then
            n=i
            return
          endif
        enddo
        stop 'Error in pir_data_constant: value of x not found.'
      endif
      end
*****************************************************************************
      subroutine pir_data_write
* check whether the data base was read correctly bu writing the data to
* PIR_data2.txt and then comparing it with PIR_data.txt.

      use pir_data_module

      implicit none
      integer i,j,k,m,n,kp,op,ip,i1,i2,nd,irnumber
      double precision xmat(48,48),x
      character lineout*10000

* open file
      open(20,file='PIR_data2.txt')
* title line
      write(20,'(a)')'Physically Irreducible Representations of the '
     $     //'230 Crystallographic Space Groups'
      write(20,'(a)')'2011 Version'
      write(20,'(a)')'Harold T. Stokes and Branton J. Campbell, 2013'
* do each IR
      do irnumber=1,ircount
* write IR number, space group number, space group symbol, IR label,
* IR dimension, number of k vectors in the star of +/-k,
* number of operators in space group
        write(20,'(i5,i4,a,5i3)')irnumber,spacegroupnumber(irnumber),
     $       ' "'//spacegroupsymbol(irnumber)//'" "'
     $       //irsymbol(irnumber)//'"',
     $       irdimension(irnumber),irtype(irnumber),kcount(irnumber),
     $       pmkcount(irnumber),operatorcount(irnumber)
* write each k vector in the star of +/-k
        kp=kvectorpointer(irnumber)-1
        write(20,'(24i3)')kvector(kp+1:kp+16*pmkcount(irnumber))
* do each operator in space group
        do i=1,operatorcount(irnumber)
* write operator in extended matrix form
          op=operatormatrixpointer(i,irnumber)-1
          write(20,'(20i4)')operatormatrix(op+1:op+16)
* if nonspecial k vector, write contribution of translation to IR
          if(.not.kspecial(irnumber))then
            write(20,'(4i3)')irtranslation(1:4,i,irnumber)
          endif
* write IR matrix
          ip=irmatrixpointer(i,irnumber)-1
          n=0
          do j=1,irdimension(irnumber)
          do k=1,irdimension(irnumber)
            n=n+1
            call pir_data_constant(1,irmatrix(ip+n),xmat(j,k))
          enddo
          enddo
        lineout=' '
        m=0
        nd=irdimension(irnumber)
        do j=1,nd
        do k=1,nd
          x=xmat(j,k)
          if(dabs(x).lt.1d-4)then
            lineout(m+1:m+2)=' 0'
            m=m+2
          else if(dabs(x-1).lt.1d-4)then
            lineout(m+1:m+2)=' 1'
            m=m+2
          else if(dabs(x+1).lt.1d-4)then
            lineout(m+1:m+3)=' -1'
            m=m+3
          else
            write(lineout(m+1:m+10),'(f10.5)')x
            m=m+1
            do n=m+1,m+9
              if(lineout(n:n).ne.' ')then
                lineout(m+1:m+1)=lineout(n:n)
                m=m+1
              endif
            enddo
            lineout(m+1:m+10)=' '
          endif
        enddo
        enddo
        i1=1
        i2=i1+79
        if(i2.gt.m)i2=m
        do while(i1.lt.m)
          do j=i2+1,i1,-1
            if(lineout(j:j).eq.' ')then
              write(20,'(a)')lineout(i1:j-1)
              i1=j
              i2=i1+79
              if(i2.gt.m)i2=m
              exit
            endif
          enddo
        enddo
* next operator
        enddo
* next IR
      enddo
 1    close(20)
      end
*****************************************************************************
      subroutine pir_data_unpack_operatormatrix(irnumber,opnumber,
     $     operatormatrix_out)
* get operators and IR matrices for a selected IR
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   opnumber (input): selected operator number
*   operatormatrix_out(k,j) (output): operator matrix

      use pir_data_module

      implicit none
      integer i,j,n,op,irnumber,opnumber,operatormatrix_out(4,4)

* check for valid selection
      if(irnumber.lt.1.or.irnumber.gt.ircount)then
        write(6,'(a)')'Error in pir_data_unpack_operatormatrix: invalid'
     $       //' value of irnumber'
        stop
      endif
      if(opnumber.lt.1.or.opnumber.gt.operatorcount(irnumber))then
        write(6,'(a)')'Error in pir_data_unpack_operatormatrix: invalid'
     $       //' value of opnumber'
        stop
      endif
* pointer to matrix
      op=operatormatrixpointer(opnumber,irnumber)
* unpack matrix
      do i=1,4
      do j=1,4
        operatormatrix_out(i,j)=operatormatrix(op)
        op=op+1
      enddo
      enddo
      end
*****************************************************************************
      subroutine pir_data_unpack_irmatrix_pointoperator(irnumber,
     $     opnumber,irmatrix_out,nr)
* get point operator part of IR matrix for a selected IR and operator
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   opnumber (input): selected operator number
*   irmatrix_out(k,j) (output): IR matrix
*   nr (input): number of rows in the array, irmatrix_out

      use pir_data_module

      implicit none
      integer nr
      integer i,j,n,ip,irnumber,opnumber
      double precision irmatrix_out(nr,nr)

* check for valid selection
      if(irnumber.lt.1.or.irnumber.gt.ircount)then
        write(6,'(a)')'Error in irr_data_unpack_irmatrix_pointoperator:'
     $       //' invalid value of irnumber'
        stop
      endif
      if(opnumber.lt.1.or.opnumber.gt.operatorcount(irnumber))then
        write(6,'(a)')'Error in irr_data_unpack_irmatrix_pointoperator:'
     $       //' invalid value of opnumber'
        stop
      endif
      if(nr.lt.irdimension(irnumber))then
        write(6,'(a)')'Error in irr_data_unpack_irmatrix_pointoperator:'
     $       //' invalid value of nr'
        stop
      endif
* pointer to matrix
      ip=irmatrixpointer(opnumber,irnumber)
* unpack matrix
      do i=1,irdimension(irnumber)
      do j=1,irdimension(irnumber)
        call pir_data_constant(1,irmatrix(ip),irmatrix_out(i,j))
        ip=ip+1
      enddo
      enddo
      end
*****************************************************************************
      subroutine pir_data_get_irmatrix_translation(irnumber,
     $     kvectorparameter,tvector,irmatrix_out,nr)
* get translation part of IR matrix for a selected IR and translation
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   kvectorparameter(i) (input), ith k vector parameter (a,b,g)
*   tvector(i) (input), ith component (x,y,z) of translation vector.
*     The values of tvector(i) are integers with a common denominator in
*     tvector(4).
*   irmatrix_out(k,j) (output): IR matrix
*   nr (input): number of rows in the array, irmatrix_out

      use pir_data_module

      implicit none
      integer i,j,k,n,ik,nr,irnumber,tvector(4),nd,nk,nb,kp,kvec(4,4)
      double precision kvectorparameter(3),irmatrix_out(nr,nr),kvec2(3),
     $     tvec2(3),dot,cosval,sinval
      double precision, parameter :: pi = 3.1415926535897932384d0

* check for valid selection
      if(irnumber.lt.1.or.irnumber.gt.ircount)then
        write(6,'(a)')'Error in pir_data_get_irmatrix_translation: '
     $       //'invalid value of irnumber'
        stop
      endif
      if(tvector(4).eq.0)then
        write(6,'(a)')'Error in pir_data_get_irmatrix_translation: '
     $       //'invalid value of tvector(4)'
        stop
      endif
      if(nr.lt.irdimension(irnumber))then
        write(6,'(a)')'Error in pir_data_get_irmatrix_translation: '
     $       //'invalid value of nr'
        stop
      endif

      irmatrix_out=0
* dimension of IR
      nd=irdimension(irnumber)
* number of k vectors in star of +-k
      nk=pmkcount(irnumber)
* dimension of block matrix in IR
      nb=nd/nk
* translation vector
      do j=1,3
        tvec2(j)=1d0*tvector(j)/tvector(4)
      enddo
* pointer to first k vector
      ik=kvectorpointer(irnumber)
* do each block
      do i=1,nk
* get k vector
        do j=1,4
        do k=1,4
          kvec(k,j)=kvector(ik)
          ik=ik+1
        enddo
        enddo
* evaluate k vector
        do j=1,3
          kvec2(j)=1d0*kvec(j,1)/kvec(4,1)
          if(.not.kspecial(irnumber))then
            do k=1,3
              if(kvec(4,k+1).ne.0)then
                kvec2(j)=kvec2(j)+kvectorparameter(k)*kvec(j,k+1)
     $               /kvec(4,k+1)
              endif
            enddo
          endif
        enddo
* cos and sin values
        dot=0
        do j=1,3
          dot=dot+kvec2(j)*tvec2(j)
        enddo
        cosval=dcos(2*pi*dot)
        sinval=dsin(2*pi*dot)
* construct block
        n=(i-1)*nb
* cos value on diagonal
        do j=1,nb
          irmatrix_out(n+j,n+j)=cosval
        enddo
* sin value off diagonal
        if(dabs(sinval).gt.1d-6)then
          do j=1,nb/2
            irmatrix_out(n+j,n+nb/2+j)=sinval
            irmatrix_out(n+nb/2+j,n+j)=-sinval
          enddo
        endif
      enddo
      end
*****************************************************************************
      subroutine pir_data_get_irmatrix(irnumber,
     $     kvectorparameter,operatormatrix_in,irmatrix_out,nr)
* get IR matrix for a selected IR and operator
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   kvectorparameter(i) (input), ith k vector parameter (a,b,g)
*   operatormatrix_in(j,i) (input), matrix of selected operator
*   irmatrix_out(k,j) (output): IR matrix
*   nr (input): number of rows in the array, irmatrix_out

      use pir_data_module

      implicit none
      integer nrop,nrir
      integer i,j,k,nd,nr,
     $     irnumber,operatormatrix_in(4,4),centeringmatrix(3,3,7),
     $     ctype,operatormatrix_stored(4,4),tvector(4),
     $     tvectorprimitive(4),opnumber
      double precision kvectorparameter(3),irmatrix_out(nr,nr),
     $     irmatrix_pointoperation(48,48),irmatrix_translation(48,48)
      character centeringtype(7)*1

      data centeringtype/'P','A','B','C','F','I','R'/
      data centeringmatrix/
     $     1,0,0,0,1,0,0,0,1,
     $     1,0,0,0,1,-1,0,1,1,
     $     1,0,-1,0,1,0,1,0,1,
     $     1,-1,0,1,1,0,0,0,1,
     $     -1,1,1,1,-1,1,1,1,-1,
     $     0,1,1,1,0,1,1,1,0,
     $     1,-1,0,0,1,-1,1,1,1/

* type of centering
      do ctype=1,7
        if(spacegroupsymbol(irnumber)(1:1).eq.centeringtype(ctype))exit
        if(ctype.eq.7)then
          write(6,'(a)')'Error in pir_data_get_irmatrix: centering '
     $         //'type not found'
          stop
        endif
      enddo
* identify operator
* try every operator in data base
      iloop: do i=1,operatorcount(irnumber)
        call pir_data_unpack_operatormatrix(irnumber,i,
     $     operatormatrix_stored)
* compare point operator parts
        do j=1,3
        do k=1,3
* not the same
          if(operatormatrix_in(k,j)/operatormatrix_in(4,4)
     $         .ne.operatormatrix_stored(k,j)
     $         /operatormatrix_stored(4,4))then
* last operator: operator not found
            if(i.eq.operatorcount(irnumber))then
              write(6,'(a)')'Error in pir_data_get_irmatrix: '
     $             //'operator not found'
              stop
* try next operator in data base
            else
              cycle iloop
            endif
          endif
        enddo
        enddo
* found it
        opnumber=i
        exit
      enddo iloop
* find difference between translation parts
      call pir_data_vsub(operatormatrix_in(1,4),
     $     operatormatrix_stored(1,4),tvector)
* primitive setting
      call pir_data_vmlt(centeringmatrix(1,1,ctype),tvector,
     $     tvectorprimitive)
* if difference is not a lattice vector (integers), invalid operator
      if(tvectorprimitive(4).ne.1)then
        write(6,'(a)')'Error in pir_data_get_irmatrix: translational '
     $       //'part of operator is invalid'
        stop
      endif
* add in translation part of IR
      if(.not.kspecial(irnumber))then
        call pir_data_vadd(irtranslation(1,i,irnumber),tvector,tvector)
      endif
* get translation part of IR
      call pir_data_get_irmatrix_translation(irnumber,kvectorparameter,
     $     tvector,irmatrix_translation,48)
* get point operation part of IR
      call pir_data_unpack_irmatrix_pointoperator(irnumber,
     $     opnumber,irmatrix_pointoperation,48)
* multiply them to obtain final IR matrix
      irmatrix_out=0
      nd=irdimension(irnumber)
      do i=1,nd
      do j=1,nd
      do k=1,nd
        irmatrix_out(i,j)=irmatrix_out(i,j)+irmatrix_translation(i,k)
     $       *irmatrix_pointoperation(k,j)
      enddo
      enddo
      enddo
      end
*****************************************************************************
      subroutine pir_data_test_multiplication_table(irnumber)
* get the IR matrices and check that they have the same multiplication
* table as the operators
      use pir_data_module
      implicit none
      integer i,j,k,m,n,i1,i2,j1,j2,
     $     irnumber,opmat(4,4,51),ndim,opmat2(4,4),ctype,
     $     centeringmatrix(3,3,7),centeringmatrix_denom(7)
      double precision irmat(48,48,51),irmat2(48,48),irmat3(48,48),
     $     kvectorparameter(3)
      character centeringtype(7)*1
      data centeringtype/'P','A','B','C','F','I','R'/
      data centeringmatrix/
     $     1,0,0,0,1,0,0,0,1,
     $     2,0,0,0,1,-1,0,1,1,
     $     1,0,-1,0,2,0,1,0,1,
     $     1,-1,0,1,1,0,0,0,2,
     $     0,1,1,1,0,1,1,1,0,
     $     -1,1,1,1,-1,1,1,1,-1,
     $     2,1,1,-1,1,1,-1,-2,1/
      data centeringmatrix_denom/1,2,2,2,2,2,3/
* for non-special k points, use arbitrary values for the parameters
      data kvectorparameter/0.11,0.12,0.13/
* get operator matrix for each operator
      opmat=0
      do i=1,operatorcount(irnumber)
        call pir_data_unpack_operatormatrix(irnumber,i,opmat(1,1,i))
      enddo
      n=operatorcount(irnumber)
* type of centering
      do ctype=1,7
        if(spacegroupsymbol(irnumber)(1:1).eq.centeringtype(ctype))exit
        if(ctype.eq.7)then
          write(6,'(a)')
     $         'Error in pir_data_get_test_multiplication_table:'
     $         //'centering type not found'
          stop
        endif
      enddo
* include lattice translations
      do i=1,3
        opmat(1:3,4,n+i)=centeringmatrix(1:3,i,ctype)
        opmat(4,4,n+i)=centeringmatrix_denom(ctype)
        do j=1,3
          opmat(j,j,n+i)=opmat(4,4,n+i)
        enddo
        call pir_data_factor(16,opmat(1,1,n+i))
      enddo
      n=n+3
* get IR matrices
      do i=1,n
        call pir_data_get_irmatrix(irnumber,kvectorparameter,
     $       opmat(1,1,i),irmat(1,1,i),48)
      enddo
* dimension of IR
      ndim=irdimension(irnumber)
* do each pair of operators
      do i1=1,n
      do i2=1,n
* multiply operators
        call pir_data_opmlt(opmat(1,1,i1),
     $       opmat(1,1,i2),opmat2)
* get IR for operator
        call pir_data_get_irmatrix(irnumber,kvectorparameter,
     $       opmat2,irmat2,48)
* multiply IR matrices
        call pir_data_dmatmlt(irmat(1,1,i1),irmat(1,1,i2),
     $       irmat3,ndim,ndim,ndim,48,48,48)
* compare IR matrices
        do j1=1,ndim
        do j2=1,ndim
          if(dabs(irmat3(j2,j1)-irmat2(j2,j1)).gt.1d-6)then
            write(6,*)'Error in pir_data_test_multiplication_table: '
     $           //'IR matrices do not have the same multiplication '
     $           //'table as the operators'
            stop
          endif
        enddo
        enddo
      enddo
      enddo
      end
***************************************************************************
      subroutine pir_data_get_transform2complex(irnumber,smat,smatinv,
     $     nr)
* construct matrix that transforms PIR to complex block-diagonal form
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   smat(i,j) (output): transformation matrix
*   smatinv(i,j) (output): inverse of smat
*   nr (input): number of rows in the arrays, smat and smatinv

      use pir_data_module

      implicit none
      integer i,j,k,m,n,nr,irnumber,nmod,nd,nd2,n1,n2,
     $     ipermute(4,4),ndk
      logical minuskvec
      double precision sqrt2
      complex*16 smat(nr,nr),zmat(nr,nr),c2real(2,2),smatinv(nr,nr)
      data ipermute/1,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0/
      data c2real/1,(0,1),1,(0,-1)/

      sqrt2=dsqrt(2d0)
* determine if -k is in the star of k
      if(pmkcount(irnumber).eq.kcount(irnumber))then
        minuskvec=.false.
      else
        minuskvec=.true.
      endif
* dimension of irrep matrix
      nd=irdimension(irnumber)
* number of blocks
      nmod=pmkcount(irnumber)
* size of block
      nd2=nd/nmod
* matrix to undo permutation if -k is in the star of k
      if(minuskvec.and.irtype(irnumber).ne.1)then
        smat=0
        do i=1,nmod
          do j=1,4
          do k=1,4
            n1=(j-1)*nd2/4+(i-1)*nd2
            n2=(k-1)*nd2/4+(i-1)*nd2
            do m=1,nd2/4
              smat(n1+m,n2+m)=ipermute(j,k)
            enddo
          enddo
          enddo
        enddo
      else
        smat=0
        do i=1,nd
          smat(i,i)=1
        enddo
      endif
* matrix to undo transformation to real form
      if(irtype(irnumber).ne.1.or.minuskvec)then
      zmat=0
      do i=1,nmod
        do j=1,2
        do k=1,2
          n1=(j-1)*nd2/2+(i-1)*nd2
          n2=(k-1)*nd2/2+(i-1)*nd2
          do m=1,nd2/2
            zmat(n1+m,n2+m)=dconjg(c2real(k,j))/sqrt2
          enddo
        enddo
        enddo
      enddo
      call pir_data_zmatmlt(zmat,smat,smat,nd,nd,nd,nr,nr,nr)
      endif
* matrix to permute blocks
      if(irtype(irnumber).ne.1)then
      zmat=0
      j=1
      do i=1,nmod*2
        n1=(i-1)*nd2/2
        n2=(j-1)*nd2/2
        do k=1,nd2/2
          zmat(n1+k,n2+k)=1
        enddo
        j=j+2
        if(j.gt.nmod*2)j=2
      enddo
      call pir_data_zmatmlt(zmat,smat,smat,nd,nd,nd,nr,nr,nr)
      endif
* inverse
      do i=1,nd
      do j=1,nd
        smatinv(j,i)=dconjg(smat(i,j))
      enddo
      enddo
      end
***************************************************************************
      subroutine pir_data_vadd(nv1,nv2,nv3)
* add two vectors of rational numbers: nv3=nv1+nv2
      implicit none
      integer nv1(4),nv2(4),nv3(4),j
      do j=1,3
        nv3(j)=nv1(j)*nv2(4)+nv2(j)*nv1(4)
      enddo
      nv3(4)=nv1(4)*nv2(4)
      call pir_data_factor(4,nv3)
      return
      end
****************************************************************************
      subroutine pir_data_vsub(nv1,nv2,nv3)
* subtract two vectors of rational numbers: nv3=nv1-nv2
      implicit none
      integer nv1(4),nv2(4),nv3(4),nv(4),j
      do j=1,3
        nv(j)=-nv2(j)
      enddo
      nv(4)=nv2(4)
      call pir_data_vadd(nv1,nv,nv3)
      end
*****************************************************************************
      subroutine pir_data_vmlt(mat,nv1,nv2)
* multiply a 3x3 matrix times a vector: nv2=mat*nv1
      implicit none
      integer mat(3,3),nv1(4),nv2(4),nv3(4),j,k
      nv3=0
      do j=1,3
      do k=1,3
        nv3(j)=nv3(j)+mat(j,k)*nv1(k)
      enddo
      enddo
      nv3(4)=nv1(4)
      call pir_data_factor(4,nv3)
      nv2=nv3
      end
*****************************************************************************
      subroutine pir_data_opmlt(op1,op2,op3)
* multiply two operator matrices op3=op1*op2
      implicit none
      integer i,j,k,m,n,op1(4,4),op2(4,4),op3(4,4),op(4,4)
      op=0
      do i=1,4
      do j=1,4
      do k=1,4
        op(i,j)=op(i,j)+op1(i,k)*op2(k,j)
      enddo
      enddo
      enddo
* remove common factor
      call pir_data_factor(16,op)
* copy into output
      op3=op
      end
*****************************************************************************
      subroutine pir_data_dmatmlt(x1,x2,x3,nrow1,ncol1,ncol2,nr1,nr2,
     $     nr3)
* multiply two matrices of double precision numbers, x3=x1*x2
* arguments:
*     x1,x2 (input), first and second matrix
*     x3 (output), product x1*x2
*     nrow1 (input), number of rows in x1, also the number of rows in x3
*     ncol1 (input), number of columns in x1, also the number of
*          rows in x2
*     ncol2 (input), number of columns in x2, also the number of
*          columns in x3
*     nr1 (input), number of rows in the physical array x1
*     nr2 (input), number of rows in the physical array x2
*     nr3 (input), number of rows in the physical array x3
      implicit none
      integer i,j,k,nrow1,ncol1,ncol2,nr1,nr2,nr3
      double precision x1(nr1,ncol1),x2(nr2,ncol2),x3(nr3,ncol2),
     $     x(nrow1,ncol2)
      x=0
      do i=1,ncol2
      do j=1,nrow1
      do k=1,ncol1
        x(j,i)=x(j,i)+x1(j,k)*x2(k,i)
      enddo
      enddo
      enddo
      x3(1:nrow1,1:ncol2)=x(1:nrow1,1:ncol2)
      end
*****************************************************************************
      subroutine pir_data_zmatmlt(x1,x2,x3,nrow1,ncol1,ncol2,nr1,nr2,
     $     nr3)
* multiply two matrices of complex*16 numbers, x3=x1*x2
* arguments:
*     x1,x2 (input), first and second matrix
*     x3 (output), product x1*x2
*     nrow1 (input), number of rows in x1, also the number of rows in x3
*     ncol1 (input), number of columns in x1, also the number of
*          rows in x2
*     ncol2 (input), number of columns in x2, also the number of
*          columns in x3
*     nr1 (input), number of rows in the physical array x1
*     nr2 (input), number of rows in the physical array x2
*     nr3 (input), number of rows in the physical array x3
      implicit none
      integer i,j,k,nrow1,ncol1,ncol2,nr1,nr2,nr3
      complex*16 x1(nr1,ncol1),x2(nr2,ncol2),x3(nr3,ncol2),
     $     x(nrow1,ncol2)
      x=0
      do i=1,ncol2
      do j=1,nrow1
      do k=1,ncol1
        x(j,i)=x(j,i)+x1(j,k)*x2(k,i)
      enddo
      enddo
      enddo
      x3(1:nrow1,1:ncol2)=x(1:nrow1,1:ncol2)
      end
*****************************************************************************
      subroutine pir_data_factor(n,numbers)
* remove the greatest common factor contained in n integers in numbers
      implicit none
      integer n,numbers(100000),min,i,j,factor
      factor=1
* find a nonzero integer
      do i=1,n
        if(numbers(i).ne.0)goto 2
      enddo
* all zeros
      return
* find the minimum absolute nonzero value among the integers
2     min=iabs(numbers(i))
      do i=2,n
        if(numbers(i).ne.0.and.iabs(numbers(i)).lt.min)
     +      min=iabs(numbers(i))
      enddo
* try each number from the minimum on down to 2
      do i=min,2,-1
* is i a common factor?
        do j=1,n
          if(mod(numbers(j),i).ne.0)goto 1
        enddo
* yes, divide it out
        do j=1,n
          numbers(j)=numbers(j)/i
        enddo
* save it too
        factor=factor*i
* done
        return
* try next number
1       continue
      enddo
      end





