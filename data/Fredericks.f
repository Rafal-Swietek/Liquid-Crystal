       program Fredericks_transition
        integer L, k, i
        real Ef, npp1, npp2, npp3, a, b, c, h
        !mpp - 2nd differentiation of n1
        real m1, m2, m3
        real, allocatable:: n1(:), n2(:), n3(:), Etab(:)
        k = 100
        allocate (n1(k), n2(k), n3(k), Etab(k))
        open(40,file='n(E)L=50-neww-k.txt')
        open(80,file='n(E)L=50-n.txt')
        open(120,file='n(E)L=50-newwk.txt')
        open(77,file='Ef_2.txt')

        do i=1,k
         read(40,*) Etab(i), n1(i)
         read(80,*) E, n2(i)
         read(120,*) E, n3(i)
         if(E.gt.1.5) EXIT
         write(*,*) Etab(i), n1(i), n2(i), n3(i)
        enddo
        h = 0.01**2
        npp1 = abs(n1(3)+n1(1)-2*n1(2))/h
        npp2 = abs(n2(3)+n2(1)-2*n2(2))/h
        npp3 = abs(n3(3)+n3(1)-2*n3(2))/h
        write(*,*) '---------------------------------------------------'
        
        !three-point scheme for evaluating the second derivative
        
        do i=3, k-1
          a = abs(n1(i+1)+n1(i-1)-2*n1(i))/h
          b = abs(n2(i+1)+n2(i-1)-2*n2(i))/h
          c = abs(n3(i+1)+n3(i-1)-2*n3(i))/h
          if(a.lt.npp1) then
            npp1 = a
            m1 = Etab(i)
          endif
          if(b.lt.npp2) then
             npp2 = b
             m2 = Etab(i)
          endif
          if(c.lt.npp3) then
             npp3 = c
             m3 = Etab(i)
          endif
          if(Etab(i).gt.1.5) EXIT
        enddo
        Ef = (m1 + m2 + m3)/3
        write(*,*) 'E_fredericks = ' , Ef
        !write(77,*) Ef
        pause
        end
