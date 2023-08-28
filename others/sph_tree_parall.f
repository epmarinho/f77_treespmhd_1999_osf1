! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
# 1 "SPH_tree.F"
! 
!
!      TREE CONSTRUCTION AND TREE-ACCESS SUBROUTINES:
!
! 
!
      subroutine SetBox
      INCLUDE 'SPH_common.h'
      real aux,vertix
 
 
 
      real scale,s1bit,ar,c,a
 
      common/aux/ aux(nn,dim)
      equivalence (aux,vertix)
      dimension vertix(1,dim)
!
!
!      * toma o vertice inferior do triedro circunscrito no sistema:
!          vertix(j)=min(xp(i,j))
!      * determina a maior aresta do cubo inscrito no triedro:
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
 
 
 
 
            box(i,j)=scale*(xp(i,j)-vertix(1,j))
 
          end do
      end do
      end
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
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
 
 
 
 
      call SetBox
      nsave=n
      do i=1,n
          t(i,1)=0
          t(i,four)=0
          t(i,2)=i
      end do
 
 
 
 
 
 
 
 
 
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
 
 
 
 
 
 
!            link node to the tree:
            call LinkTree
          end do
          call NoLeaves
      end do
!
!      finalization of tree-construction:
      maxnode=node-1
      n=nsave
 
 
 
!
 
      if(first)then
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
        do i=1,n
          h(i)=0.3*s(i)
        end do
 
          first = .false.
      endif
 
 
 
 
 
 
 
      end 
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
 
 
 
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
      end 
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
 
 
 
 
      subroutine LinkTree
      integer cell,j,l,ip
      INCLUDE 'SPH_common.h'
 
      do cell=1,vert
 
 
 
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
      end
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
 
 
 
 
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
      end 
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
 
 
 
 
 
 
      subroutine FindOctants(ia,ib)
      integer cell,i,ia,ib,j,ip
 
 
 
      INCLUDE 'SPH_common.h'
      common /non_singularity/ isingularity
      do i=ia,ib
          ip=t(i,2)
          cell=0
          do j=1,dim
                cell=2*cell
            box(ip,j)=2*box(ip,j)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
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
      end 
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
 
 
 
 
      subroutine CLR
      integer cell,j
      INCLUDE 'SPH_common.h'
 
      do cell=1,vert
 
 
 
          count(cell) = 0
          cellmass(node,cell) = 0
          do j=1,dim
            x_cell(node,cell,j) = 0
          end do
      end do
      nextnode = .false.
      end
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
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
 
          t(i,1)=vert*t(i,1)+t(i,3)
 
 
 
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
      end 
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
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
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 12:56:19
!### KAP/Digital_UA_F detected syntax errors in this program.
!### Check the KAP/Digital_UA_F listing for details.
!
!
      subroutine forbid
!
      integer cell
!
      INCLUDE 'SPH_common.h'
!
 
      do cell=1,vert
 
 
 
          next=down(node,cell)
          if(next.gt.0)permit(next)=.false.
      end do
      end
!
!
# 536 "SPH_tree.F"
 
 
!
!
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
