!///////////////////////////////////////////////////////////////////////
!
!      TREE CONSTRUCTION AND TREE-ACCESS SUBROUTINES:
!
!///////////////////////////////////////////////////////////////////////
!
      subroutine SetBox
      INCLUDE 'SPH_common.h'
      real aux,vertix
#ifdef _TUBE
      real *8 scale,s1bit,ar,c,a,xperr
#else
      real scale,s1bit,ar,c,a
#endif
      common/aux/ aux(nn,dim)
      equivalence (aux,vertix)
      dimension vertix(1,dim)
!
!
!      * take the lower vertix from the trihedron enclosing the system:
!          vertix(j)=min(xp(i,j))
!      * set the major cube arest:
!          arest = max(xp(i,j))-min(xp(i,j))
!
      do j=1,dim
          vertix(1,j)=xp(1,j)
      end do
      c=xp(1,1)
      a=c
      do j=1,dim
          do i=2,n
            if(xp(i,j).gt.c)then
                c=xp(i,j)
            else
                if(xp(i,j).lt.vertix(1,j))vertix(1,j)=xp(i,j)
            endif
          end do
          if(vertix(1,j).lt.a)a=vertix(1,j)
      end do
      ar=c-a
!
!      * determina o cubo circunscrito cujas arestas sao potencias
!        inteiras de 2:
!
      s1bit=1.
      if(ar.gt.1.)then
          do while(s1bit.lt.ar)
            s1bit=s1bit+s1bit
          end do
      else
          ar=2*ar
          do while(s1bit.gt.ar)
            s1bit=0.5*s1bit
          end do
      endif
      if(s1bit.ne.0)then
          scale=1d0/s1bit
          arest(1)=0.5*s1bit
      else
          print *,
     &      '%err: particles were distributed into a singularity!!'
          stop
      endif
!
!      renormaliza as coordenadas das particulas para um cubo
!      unitario:
!
      do j=1,dim
          do i=1,n
#ifdef _TUBE
            xperr = 1e-15*(1-2*ran(iiseed))
            box(i,j)=scale*(xp(i,j)-vertix(1,j))*(1d0+xperr)
#else
            box(i,j)=scale*(xp(i,j)-vertix(1,j))
#endif
          end do
      end do
      end
!
!
      subroutine MakeTree
!      this procedure generates a tree data-structure for a given 
!      set of n mass-points in 3-dim space. the method consists on 
!      updating the particles in all corresponding nodes in each level
!      at once, which implies in a one descenting construction, 
!      discarding recursive procedures. the philosophy of tree-
!      descendent construction has the principal advantage on its easy
!      adaptation for vector machines.
!
!    data configuration:
!      nn=n
!            nnode = n - 1
!            maxoctants = eight * nnode
!
!
!
!
      integer i,i1,nsave
      real volume,dens_mean
      INCLUDE 'SPH_common.h'
      common /indices/ i,i1
!
!
!      inicializacao
!
#ifdef _VERBOSE
      print*,'maketree: start'
#endif/*_VERBOSE*/

      call SetBox
      nsave=n
      do i=1,n
          t(i,1)=0
          t(i,four)=0
          t(i,2)=i
      end do
#ifndef _N_BODY
#  ifdef _OLD_INITIAL_LENGTHS_METHOD
      if(first)then
        do i=1,n
          h(i)=0
        end do
      endif
#  endif
#endif
      node = 1
      next = node + 1
!
!      while any particle is non-leaf yet, continue building tree:
      do while(n.gt.1)
          call FindOctants(1,n)
          i=1
          do while(i.le.n)
            !make node by sharing cubic-cells
            !with respective particles:
            call MakeNode
#ifndef _N_BODY
#  ifdef _OLD_INITIAL_LENGTHS_METHOD
            !calculate sph lengths in the first run:
            if(first)call InitialLengths
#  endif
#endif
!            link node to the tree:
            call LinkTree
          end do
          call NoLeaves
      end do
!
!      finalization of tree-construction:
      maxnode=node-1
      n=nsave
#ifdef _VERBOSE
      print*,'Tree completed! n=',n
#endif/*_VERBOSE*/
!
#ifndef _N_BODY
      if(first)then
#  ifdef _OLD_INITIAL_LENGTHS_METHOD
        volume = 0
        do i=1,n
	    !** print*,h(i),i
            volume = volume + 1/h(i)
        end do
	if(n_fix.eq.0)n_fix=48
        if(tol_nfix.eq.0)tol_nfix=1
#ifdef _VERBOSE
	print*,'volume=',volume,', n_fix=',n_fix
#endif
        dens_mean = n/volume
        do i=1,n
#    ifdef _2D
            h(i)=fine/sqrt(dens_mean+h(i)/n_fix)
#    else
            h(i)=fine/(  dens_mean  +  h(i)/n_fix  )**.3333333
#    endif
        end do
#ifdef _VERBOSE
	print*,'dens=',dens_mean
#endif
#  else
        do i=1,n
          h(i)=0.3*s(i)
        end do
#  endif
          first = .false.
      endif
#endif
#if !(defined(_TUBE)||defined(_NON_SELFGRAVITATING)||defined(_2D))
      call Quadrupole
#endif
#ifdef _VERBOSE
      print*,'maketree: done'
#endif/*_VERBOSE*/
      end /*MakeTree*/



      subroutine MakeNode
      logical changing
      integer cell,ip,i1,i,j
      INCLUDE 'SPH_common.h'
      common /indices/ i,i1
      common /non_singularity/ isingularity
      isingularity=0
      call CLR
      i1=i
      do while(.not.nextnode)
          ip=t(i,2)
          cell=t(i,3)+1
          count(cell)=count(cell)+1
          cellmass(node,cell)=cellmass(node,cell)+mass(ip)
          lst(count(cell),cell)=ip
!         calculate the first-order position-momentum:
          do j=1,dim
             x_cell(node,cell,j)=x_cell(node,cell,j)+xp(ip,j)*mass(ip)
          end do
          if(count(cell).eq.1)ifirst(cell)=i
          i=i+1
          if(i.le.n)changing=t(i,1).gt.t(i-1,1)
          nextnode=changing.or.(i.gt.n)
          if(nextnode)call NonDegeneracy
      end do
      end /*MakeNode*/




      subroutine LinkTree
      integer cell,j,l,ip
      INCLUDE 'SPH_common.h'
#ifndef _TUBE
      do cell=1,vert
#else
      do cell=1,2
#endif
        np(node,cell)=count(cell)
!        make centroid:
        if(count(cell).ge.1)then
!          resume cell centroide:
          do j=1,dim
           x_cell(node,cell,j)=x_cell(node,cell,j)/cellmass(node,cell)
          end do
        endif
!
        if(count(cell).gt.1)then
            down(node,cell)=next
            arest(next)=0.5*arest(node)
            next=next+1
        else
            down(node,cell)=0
            if(count(cell).eq.1)then
                l=ifirst(cell)
                ip=t(l,2)
                s(ip)=arest(node)
                t(l,four) = 1
                label(node,cell)=ip
            endif
        endif
      end do
      node=node+1
      if(node.gt.nnode)then
        print*,'LinkTree: ERROR: tree overflow'
        print*,'node=',node,'  nnode=',nnode
        stop
      endif
      end/*LinkTree*/




      subroutine NonDegeneracy
      integer cell,i,i1,idegen
      INCLUDE 'SPH_common.h'
      common /non_singularity/ isingularity
      common /indices/ i,i1
      idegen=0
      do cell=1,vert
          if(count(cell).eq.0)idegen=idegen+1
      end do
      if(idegen.eq.(vert-1))then
          call CLR
          isingularity=isingularity+1
          call FindOctants(i1,i-1)
          arest(node)=0.5*arest(node)
          i=i1
      else
          isingularity=0
      endif
      end /*NonDegeneracy*/






      subroutine FindOctants(ia,ib)
      integer cell,i,ia,ib,j,ip
#ifdef _TUBE
      real *8 xperr
#endif
      INCLUDE 'SPH_common.h'
      common /non_singularity/ isingularity
      do i=ia,ib
          ip=t(i,2)
          cell=0
          do j=1,dim
                cell=2*cell
            box(ip,j)=2*box(ip,j)

#ifdef _CHECK_TREE_SINGULARITY
            if(isingularity.gt.23)then
              isingularity=0
#ifdef _VERBOSE
              print*,'Tree-singularity occurring: assuming dynamic round-off'
              print*,'old box: ',box(ip,j)
#endif
#ifdef _TUBE
              xperr = -1d-15*(1+ran(iiseed))
#else
              xperr = -1e-6*(1+ran(iiseed))
#endif
              box(ip,j) = box(ip,j)*(1+xperr)
#ifdef _VERBOSE
              print*,'new box: ',box(ip,j)
              print*,'x-error: ',xperr
#endif
            endif
#endif /*_CHECK_TREE_SINGULARITY*/

            if(box(ip,j).gt.1)then
                box(ip,j)=box(ip,j)-1
                cell=cell+1
            endif

            if(box(ip,j).gt.1)then
              print*,'FindOctants: error'
              print*,'box=',box(ip,j)
              stop
            endif

          end do
          t(i,3)=cell
      end do
      end /*FindOctants(ia,ib)*/




      subroutine CLR
      integer cell,j
      INCLUDE 'SPH_common.h'
#ifndef _TUBE
      do cell=1,vert
#else
      do cell=1,2
#endif
          count(cell) = 0
          cellmass(node,cell) = 0
          do j=1,dim
            x_cell(node,cell,j) = 0
          end do
      end do
      nextnode = .false.
      end
!
!
      subroutine NoLeaves
      integer k,i,j
      logical noleaf
!
      INCLUDE 'SPH_common.h'
!
!      copia somente as linhas nao-folha
!      e dispensa a coluna dos "leaf-flags"(4):
!
      k=0
      do i=1,n
          noleaf=(t(i,four).ne.1)
          if(noleaf)then
            k=k+1
            do j=1,four
                t(k,j)=t(i,j)
            end do
          endif
      end do
      n=k
      if(n.eq.0)return
!
!      calcula o novo parametro de peano promovendo os octantes para o 
!      status de no':
!
      do i=1,n
#ifndef _TUBE
          t(i,1)=vert*t(i,1)+t(i,3)
#else
          t(i,1)=2*t(i,1)+t(i,3)
#endif
      end do
!
!      ordena o array t por heap-sort dispensando a coluna dos 
!      octantes (col. 3):
!
      call hsort
!
!      reduz o parametro de peano a fim de evitar estouro de inteiro:
!
      k=0
      do i=2,n
          if(t(i,1).ne.t(i-1,1))k=k+1
          t(i-1,3)=k
      end do
      do i=2,n
          t(i,1)=t(i-1,3)
      end do
      t(1,1)=0
!
      end /*NoLeaves*/
!
!
      subroutine hsort
!
!            method used: heapsort (num. rec. 1986)
!
      integer l,ir,k,i,j
      INCLUDE 'SPH_common.h'
!
!
      l=n/2+1
      ir=n
   10      continue
          if(l.gt.1)then
            l=l-1
            do k=1,2
                row(k)=t(l,k)
                t(l,k)=t(1,k)
            end do
          else
            do k=1,2
                row(k)=t(ir,k)
                t(ir,k)=t(1,k)
            end do
            ir=ir-1
            if(ir.eq.1)then
                do k=1,2
                  t(1,k)=row(k)
                end do
                return
            endif
          endif
          i=l
          j=l+l
   20          if(j.le.ir)then
            if(j.lt.ir)then
                if(t(j,1).lt.t(j+1,1))j=j+1
            endif
            if(row(1).lt.t(j,1))then
                do k=1,2
                  t(i,k)=t(j,k)
                end do
                i=j
                j=j+j
            else
                j=ir+1
            endif
          go to 20
          endif
          do k=1,2
            t(i,k)=row(k)
          end do
      go to 10
      end
!
!
      subroutine forbid
!
      integer cell
!
      INCLUDE 'SPH_common.h'
!
#ifdef _TUBE
      do cell=1,2
#else
      do cell=1,vert
#endif
          next=down(node,cell)
          if(next.gt.0)permit(next)=.false.
      end do
      end
!
!
#if !(defined( _TUBE)||defined(_2D))
      subroutine quadrupole
      integer node1,m1,i,j,node2,m2
      real aux,D_x
      INCLUDE 'SPH_common.h'
      common /auxiliar/ aux(nn,dim)
      dimension D_x(nn,dim)
      equivalence (aux, D_x)
!
      do node1=maxnode,1,-1
        do m1=1,eight
!
          do i=1,dim
            do j=i,dim
            p_inertia(node1,m1,i,j)=0
            p_inertia(node1,m1,j,i)=0
            end do
          end do
!
          if(np(node1,m1).gt.1)then
!
            node2  =  down(node1,m1)
!
            do m2=1,eight
!
            do i=1,dim
              D_x(1,i)  =  x_cell(node1,m1,i) - x_cell(node2,m2,i)
            end do
!
!            parallel axes theorem:
            do j=1,dim
                    do i=j,dim
!
                p_inertia(node1,m1,i,j) = p_inertia(node1,m1,i,j) +
     &                                          p_inertia(node2,m2,i,j) +
     &                         D_x(1,i)*cellmass(node2,m2)*D_x(1,j)
!
                p_inertia(node1,m1,j,i) = p_inertia(node1,m1,i,j)
!
                    end do
            end do
!
            end do
!
          endif
!
        end do
      end do
      end
#endif
#ifndef _N_BODY
!
!
#ifdef _OLD_INITIAL_LENGTHS_METHOD
      subroutine InitialLengths
      integer cell,k,j
      INCLUDE 'SPH_common.h'
!     for each cell of current node do:
#ifndef _TUBE
      do cell=1,vert
#else
      do cell=1,2
#endif
        if(count(cell).gt.0)then
          !travel through the particles-in-cell list:
          do k=1,count(cell)
            !take the k-th particle of the list:
            j=lst(k,cell)
            h(j) = h(j)  +  arest(node)
          end do
          h(j)=h(j)/count(cell)
        endif
      end do
      end
#endif _OLD_INITIAL_LENGTHS_METHOD

#endif

