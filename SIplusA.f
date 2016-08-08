      PROGRAM SIplusA.1.3.3




c	This program simulates the genetic evolution of a population with
c	mixed asexual and sexual reproduction and a self-incompatibility system




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DEFINITION OF VARIABLES/PARAMETERS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	INTEGER N ! population size of N diploid individuals
	INTEGER neutral_loci ! population size of N diploid individuals
	INTEGER k ! number of alleles at the k-allele model
	DOUBLE PRECISION u_neutral ! mutation rate for neutral loci
	DOUBLE PRECISION u_S ! mutation rate for S locus
	DOUBLE PRECISION u_S0 ! mut rate towards SelfCompatibility, IRREVERSIBLE
	DOUBLE PRECISION u_1 ! mutation rate for viability locus A to a
	DOUBLE PRECISION u_2 ! mutation rate for viability locus a to A
	REAL s ! selection coefficient
	REAL h ! dominance coefficient
	REAL c ! rate of asexual reproduction
	INTEGER generations ! number of generations to be run
	INTEGER burn_in ! number of generations run before calculating indices
	INTEGER sampling_interval !interval of generations for calculating indices
	INTEGER replicates ! number of simulation replicates to be run

	INTEGER saves_GenePop
	INTEGER start_as_selfcompatible
	INTEGER any_selfcompatible
	INTEGER total_selfcompatible_pop

	INTEGER parent_generation(1000,13,4) !genotypes at parental generation
	INTEGER offspring_generation(1000,12,2) !genotypes at offspring generation

	REAL freq_S(0:4000)
c	REAL freq_S0
	REAL freq_V(2)
	REAL freq_N(2,100)

	INTEGER j
	REAL T_xy,T_xz, denom
	REAL delta_xy,delta_xz
	REAL R_xy,Rxz
	REAL C_x,C_y,C_z
	REAL Hobs_x,Hobs_y,Hobs_z
	REAL F_xy,F_xz
	INTEGER	Nc_xy,Nc_xz
	REAL sum_F_xy,sum_F_xz
	REAL moy_F_xy,moy_F_xz
	REAL R_GGD_SV,R_GGD_SN1,R_GGD_N1V,R_GGD_N1N2



	INTEGER x,xx,y,yy,z,zz,w,ww

	INTEGER ovule(12,2)
	INTEGER pollen(12,2)
	INTEGER samples
	INTEGER calculate_LD

	INTEGER na,contador
	DOUBLE PRECISION ne_s

	INTEGER i,ii,iii,iiii,t,rep,using_input_file, locus

	INTEGER number_of_clones, clon, mother, father, allele
	INTEGER number_of_mutants_S,S_alleles,mutant,number_of_mutants_N

	REAL treshold_selection

	REAL is_it_time

	REAL random, poisson, binomial, random1,random2,random3
	REAL lambda

	INTEGER F00,F10,F11,F1,F0,seed


	DOUBLE PRECISION Shannon_E_n(10)
	INTEGER na_n(10)

	DOUBLE PRECISION Shannon_E, Shannon_H, sum_E

	DOUBLE PRECISION delta,Ws,Wo,L,He,Ho,Fis,Ne,pa,sum
	DOUBLE PRECISION He_n(10),Ho_n(10),Fis_n(10),ne_n(10)
	DOUBLE PRECISION Neff(12),He_s,Fis_s,Ho_s



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT OF PARAMETERS VALUES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	PRINT *,''
	PRINT *,'=====>  Welcome to SIplusA simulation   <====='
	PRINT *,'        by M. Navascues & S. Mariette'
	PRINT *,''
	PRINT *,'        this is SIplusA v. 1.3.3 beta'
	PRINT *,''


	OPEN (1,file='log.txt')
	OPEN (3,file='output_V.txt')
	OPEN (4,file='output_S.txt')
	OPEN (5,file='output_N1.txt')
	OPEN (55,file='output_N2.txt')
	OPEN (7,file='output_ef.txt')
	OPEN (10,file='S_allele_freq.txt')
c	OPEN (30,file='genetix.txt')


	WRITE (1,*),''
	WRITE (1,*),'=====>  Welcome to SIplusA simulation   <====='
	WRITE (1,*),'        by M. Navascues & S. Mariette'
	WRITE (1,*),''
	WRITE (1,*),'        this is SIplusA v. 1.3.3 beta'
	WRITE (1,*),''

5	PRINT *,'Do you want to use input file (input.txt)? Yes=1 No=0'
	READ (*,*) using_input_file
	WRITE (1,*)'Do you want to use input file (input.txt)? Yes=1 No=0'
	WRITE (1,*) using_input_file
	PRINT *,''

	if ((using_input_file.ne.1).and.(using_input_file.ne.0)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 5
	endif

c.....reading parameters from imput file

	IF (using_input_file.eq.1) then


	OPEN (2,file='input.txt')


	PRINT *,'Reading input file...'

	READ (2,*)seed
	WRITE (1,*)'Input seed (integer) for random number generator'
	WRITE (1,*)seed

	READ (2,*)N
	WRITE (1,*)'Input population size (number of diploid individuals)'
	WRITE (1,*)'(2-1000)'
	WRITE (1,*)N

	if ((N.lt.2).or.(N.gt.1000)) then
		Print *,'Invalid input file, wrong population size (2-1000)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)neutral_loci
	WRITE (1,*)'Input number of neutral loci (0-10)'
	WRITE (1,*)neutral_loci

	if ((neutral_loci.lt.0).or.(neutral_loci.gt.10)) then
		Print *,'Invalid input file, wrong number of loci (0-10)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)k
	WRITE (1,*)'Input number of alleles for neutral loci (2-100)'
	WRITE (1,*)k

	if ((k.lt.2).or.(k.gt.100)) then
		Print *,'Invalid input file, wrong number of alleles  (2-100)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)u_neutral
	WRITE (1,*)'Input mutation rate for neutral loci (0-1)'
	WRITE (1,*)u_neutral

	if ((u_neutral.lt.0).or.(u_neutral.gt.1)) then
		Print *,'Invalid input file, wrong mutation rate (0-1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)u_S
	WRITE (1,*)'Input mutation rate for S locus (0-1)'
	WRITE (1,*)u_S

	if ((u_S.lt.0).or.(u_S.gt.1)) then
		Print *,'Invalid input file, wrong mutation rate (0-1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)u_S0
	WRITE (1,*)'Input mutation rate to self-compatible allele (0-1)'
	WRITE (1,*)u_S0

	if ((u_S0.lt.0).or.(u_S0.gt.1)) then
		Print *,'Invalid input file, wrong mutation rate (0-1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif


	READ (2,*)u_1
	WRITE (1,*)'Input mutation rate for viability locus: A->a (0-1)'
	WRITE (1,*)u_1

	if ((u_1.lt.0).or.(u_1.gt.1)) then
		Print *,'Invalid input file, wrong mutation rate (0-1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)u_2
	WRITE (1,*)'Input mutation rate for viability locus: a->A (0-1)'
	WRITE (1,*)u_2

	if ((u_2.lt.0).or.(u_2.gt.1)) then
		Print *,'Invalid input file, wrong mutation rate (0-1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)s
	WRITE (1,*)'Input selection coefficient, s (0-1)'
	WRITE (1,*)s

	if ((s.lt.0).or.(s.gt.1)) then
		Print *,'Invalid input file, wrong selection coef. (0-1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)h
	WRITE (1,*)'Input dominance coefficient, h (0-1)'
	WRITE (1,*)h

	if ((h.lt.0).or.(h.gt.1)) then
		Print *,'Invalid input file, wrong dominance coef. (0-1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif
	if (h.lt.0.5) then
		Write(1,*)'deletereous allele is recesive'		
	else
		Write(1,*)'deletereous allele is dominant'		
	endif

	READ (2,*)c
	WRITE (1,*)'Input rate of asexual reproduction, c (0-1)'
	WRITE (1,*)c

	if ((c.lt.0).or.(c.gt.1)) then
		Print *,'Invalid input file, wrong asexual rep. rate(0-1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)generations
	WRITE (1,*)'Input number of generations to be run (>=1)'
	WRITE (1,*)generations

	if (generations.lt.1) then
		Print *,'Invalid input file, wrong # of generations (>=1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)burn_in
	WRITE (1,*)'Input lenght of burn-in interval (>=1)'
	WRITE (1,*)burn_in

	if (burn_in.lt.1) then
		Print *,'Invalid input file, wrong burn-in interval (>=1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)sampling_interval
	WRITE (1,*)'Input lenght of sampling interval (>=1)'
	WRITE (1,*)sampling_interval

	if (sampling_interval.lt.1) then
		Print *,'Invalid input file, wrong sampling interval (>=1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)replicates
	WRITE (1,*)'Input number of replicates to be run (>=1)'
	WRITE (1,*)replicates

	if (replicates.lt.1) then
		Print *,'Invalid input file, wrong number of replicates (>=1)'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)calculate_LD
	if (calculate_LD.eq.1) then
		WRITE (1,*)'SIplusA will calculate LD'
	else
		WRITE (1,*)'SIplusA will NOT calculate LD'
	endif

	if ((calculate_LD.ne.1).and.(calculate_LD.ne.0)) then
		Print *,'Invalid input file'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)saves_GenePop
	if (saves_GenePop.eq.1) then
		WRITE (1,*)'SIplusA will save population data'
	else
		WRITE (1,*)'SIplusA will NOT save population dat'
	endif

	if ((saves_GenePop.ne.1).and.(saves_GenePop.ne.0)) then
		Print *,'Invalid input file'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif

	READ (2,*)start_as_selfcompatible
	if (start_as_selfcompatible.eq.1) then
	  WRITE (1,*)'SIplusA will simulate a self-compatible population'
	else
	 WRITE (1,*)'SIplusA will simulate a self-incompatible population'
	endif

	if ((start_as_selfcompatible.ne.1)
     >.and.(start_as_selfcompatible.ne.0)) then
		Print *,'Invalid input file'
		Print *,'Enter parameters manually'
		Write(1,*)'Invalid input file'
		Write(1,*)'Entering parameters manually'
		using_input_file=0
		goto 20
	endif







	ENDIF

c.....introduction of parameters manually
	CLOSE (2)

20	IF (using_input_file.eq.0) then

	PRINT *,'Input seed (integer) for random number generator'
	READ (*,*)seed
	WRITE (1,*)'Input seed (integer) for random number generator'
	WRITE (1,*)seed
	PRINT *,''

22	PRINT *,'Input population size (number of diploid individuals)'
	PRINT *,'(2-1000)'
	READ (*,*)N
	WRITE (1,*)'Input population size (number of diploid individuals)'
	WRITE (1,*)'(2-1000)'
	WRITE (1,*)N
	PRINT *,''

	if ((N.lt.2).or.(N.gt.1000)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 22
	endif

24	PRINT *,'Input number of neutral loci (0-10)'
	READ (*,*)neutral_loci
	WRITE (1,*)'Input number of neutral loci (0-10)'
	WRITE (1,*)neutral_loci
	PRINT *,''

	if ((neutral_loci.lt.0).or.(neutral_loci.gt.10)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 24
	endif

26	PRINT *,'Input number of alleles for neutral loci (2-100)'
	Print *,'(uses k-allele mutation model)'
	READ (*,*)k
	WRITE (1,*)'Input number of alleles for neutral loci (2-100)'
	WRITE (1,*)k
	PRINT *,''

	if ((k.lt.2).or.(k.gt.100)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 26
	endif

28	PRINT *,'Input mutation rate for neutral loci (0-1)'
	READ (*,*)u_neutral
	WRITE (1,*)'Input mutation rate for neutral loci (0-1)'
	WRITE (1,*)u_neutral
	PRINT *,''

	if ((u_neutral.lt.0).or.(u_neutral.gt.1)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 28
	endif

30	PRINT *,'Input mutation rate for S locus (0-1)'
	READ (*,*)u_S
	WRITE (1,*)'Input mutation rate for S locus (0-1)'
	WRITE (1,*)u_S
	PRINT *,''

	if ((u_S.lt.0).or.(u_S.gt.1)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 30
	endif

31	PRINT *,'Input mutation rate to self-compatibility (0-1)'
	READ (*,*)u_S0
	WRITE (1,*)'Input mutation rate to self-compatibility (0-1)'
	WRITE (1,*)u_S0
	PRINT *,''

	if ((u_S0.lt.0).or.(u_S0.gt.1)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 31
	endif




	PRINT *,'Recomended: rate for A->a >> rate for a->A'
	PRINT *,'------------------------------------------'

32	PRINT *,'Input mutation rate for viability locus: A->a (0-1)'
	READ (*,*)u_1
	WRITE (1,*)'Input mutation rate for viability locus: A->a (0-1)'
	WRITE (1,*)u_1

	if ((u_1.lt.0).or.(u_1.gt.1)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 32
	endif

34	PRINT *,'Input mutation rate for viability locus: a->A (0-1)'
	READ (*,*)u_2
	WRITE (1,*)'Input mutation rate for viability locus: a->A (0-1)'
	WRITE (1,*)u_2
	PRINT *,''

	if ((u_2.lt.0).or.(u_2.gt.1)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 34
	endif

36	PRINT *,'Input selection coefficient, s (0-1)'
	READ (*,*)s
	WRITE (1,*)'Input selection coefficient, s (0-1)'
	WRITE (1,*)s
	PRINT *,''

	if ((s.lt.0).or.(s.gt.1)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 36
	endif

38	PRINT *,'Input dominance coefficient, h'
	Print *,'0-0.5 for recesive deletereous allele'
	Print *,'0.5-1 for dominant deletereous allele'
	READ (*,*)h
	WRITE (1,*)'Input dominance coefficient, h'
	WRITE (1,*)'0-0.5 for recesive deletereous allele'
	WRITE (1,*)'0.5-1 for dominant deletereous allele'
	WRITE (1,*)h
	PRINT *,''

	if ((h.lt.0).or.(h.gt.1)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 38
	endif
	if (h.lt.0.5) then
		Write(1,*)'deletereous allele recesive'		
	else
		Write(1,*)'deletereous allele dominant'		
	endif

40	PRINT *,'Input rate of asexual reproduction, c (0-1)'
	READ (*,*)c
	WRITE (1,*)'Input rate of asexual reproduction, c (0-1)'
	WRITE (1,*)c
	PRINT *,''

	if ((c.lt.0).or.(c.gt.1)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 40
	endif

42	PRINT *,'Input number of generations to be run (>=1)'
	READ (*,*)generations
	WRITE (1,*)'Input number of generations to be run (>=1)'
	WRITE (1,*)generations
	PRINT *,''

	if (generations.lt.1) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 42
	endif

44	PRINT *,'Input lenght of burn-in interval (>=1)'
	Print *,'Recommended >>10 if calculating linkage disequilibrium'
	READ (*,*)burn_in
	WRITE (1,*)'Input lenght of burn-in interval (>=1)'
	WRITE (1,*)burn_in
	PRINT *,''

	if (burn_in.lt.1) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 44
	endif

45	PRINT *,'Input lenght of sampling interval (>=1)'
	READ (*,*)sampling_interval
	WRITE (1,*)'Input lenght of sampling interval (>=1)'
	WRITE (1,*)sampling_interval
	PRINT *,''

	if (sampling_interval.lt.1) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 45
	endif

46	PRINT *,'Input number of replicates to be run (>=1)'
	READ (*,*)replicates
	WRITE (1,*)'Input number of replicates to be run (>=1)'
	WRITE (1,*)replicates
	PRINT *,''

	if (replicates.lt.1) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 46
	endif

47	PRINT *,'Calculate linkage disequilibrium? Yes=1, No=0'
c	Print *,'currently not working, choose "0"'
	Print *,'Recommended only with burn-in >> 10'
	READ (*,*)calculate_LD
	WRITE (1,*)'Calculate linkage disequilibrium? Yes=1, No=0'

	if ((calculate_LD.ne.1).and.(calculate_LD.ne.0)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 47
	endif


48	PRINT *,'Save population data (GenePop file)? Yes=1, No=0'
	Print *,'Not Recommended'
	READ (*,*)saves_GenePop
	WRITE (1,*)'Save population data (GenePop file)? Yes=1, No=0'

	if ((saves_GenePop.ne.1).and.(saves_GenePop.ne.0)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 48
	endif

49	PRINT *,'Start population as self-compatible? Yes=1, No=0'
	READ (*,*)start_as_selfcompatible
	WRITE (1,*)'Start population as self-compatible? Yes=1, No=0'

	if ((start_as_selfcompatible.ne.1)
     >.and.(start_as_selfcompatible.ne.0)) then
		Print *,'Invalid input'
		Write(1,*)'Invalid input'
		goto 49
	endif




	if (saves_GenePop.eq.1) then
		WRITE (1,*)'SIplusA will save population data to GenePop file'
	else
		WRITE (1,*)'SIplusA will NOT save',
     >'population data to GenePop file'
	endif




	if (calculate_LD.eq.1) then
		WRITE (1,*)'SIplusA will calculate LD'
	else
		WRITE (1,*)'SIplusA will NOT calculate LD'
	endif

	ENDIF

	idum=-seed

	samples=(generations-burn_in+sampling_interval)
     >/sampling_interval

	WRITE(3,*)replicates,'    replicates'
	WRITE(3,*)samples,'    samples'
	WRITE(3,*)'     r       t     pa    He    Ho    Fis     ',
     >'Ws     Wo      delta          L'

	WRITE(4,*)replicates,'    replicates'
	WRITE(4,*)samples,'    samples'
	WRITE(4,*)'     r       t      na           ne     He    Fis'
     >,'  ShannonE'

	WRITE(5,*)replicates,'    replicates'
	WRITE(5,*)samples,'    samples'
	WRITE(5,*)'     r      t   na1      ne1     He1    Ho1   Fis1',
     >'  ShannonE1'

	WRITE(55,*)replicates,'    replicates'
	WRITE(55,*)samples,'    samples'
	WRITE(55,*)'     r      t   na2      ne2     He2    Ho2   Fis2',
     >'  ShannonE2'



	WRITE(7,*)replicates,'    replicates'
	WRITE(7,*)samples,'    samples'
	WRITE(7,*)'     r      t     Ne_b     Ne_S',
     >'     Ne_v    Ne_n1    Ne_n2'

	WRITE(10,*)replicates,'    replicates'
	WRITE(10,*)samples,'    samples'
	WRITE(10,*)'  r      t  na     Allele frequencies'


	if (calculate_LD.eq.1) then

	OPEN (8,file='output_LD.txt')

	WRITE(8,*)replicates,'    replicates'
	WRITE(8,*)samples,'    samples'
	WRITE(8,*)'     r      t     R_GGD(S-V)',
     >'   R_GGD(S-N1)   R_GGD(N1-V)   R_GGD(N1-N2)'

	endif


	if (saves_GenePop.eq.1) then

	OPEN (20,file='GenePop.txt')

	WRITE(20,*)'SIplusA v. 1.3.1 beta:'
     >,'genetic data of simulated population'
	WRITE(20,*)'S-locus'
	WRITE(20,*)'viability_locus'
	WRITE(20,*)'neutral_locus_1'
	WRITE(20,*)'neutral_locus_2'
	endif






ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     START: LOOP REPLICATES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	DO rep=1, replicates


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     GENERATION OF INITIAL GENOTYPES (AT TIME=0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	DO i=1,N

		if (start_as_selfcompatible.eq.0) then
			parent_generation(i,1,1)=i
			parent_generation(i,1,2)=N+i
		else
			parent_generation(i,1,1)=0
			parent_generation(i,1,2)=0
		endif

		binomial=bnldev(5.0E-1,2,idum)

		if (binomial.le.1) then
			parent_generation(i,2,1)=1
		else
			parent_generation(i,2,1)=0
		endif

		if (binomial.ge.1) then
			parent_generation(i,2,2)=0
		else
			parent_generation(i,2,2)=1
		endif

		do ii=1, neutral_loci
			do iii=1,2
				random=ran2(idum)
				k_allele=1+INT(random*k)
				parent_generation(i,2+ii,iii)=k_allele
			enddo
		enddo

	ENDDO


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     START: LOOP GENERATIONS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	sampled=0.0

	if (start_as_selfcompatible.eq.0) then
		S_alleles=2*N
		any_selfcompatible=0
	else
		S_alleles=0
		any_selfcompatible=1
	endif

	DO t=1, generations

	
	Print *,'Replicate',rep,'      Generation',t



c.....sets reproductive success of parents to zero
	do i=1,N
		parent_generation(i,13,1)=0
	enddo


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NUMBER OF ASEXUAL REPRODUCTIONS DRAW FROM POISSON DISTRIBUTION
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	if ((c.eq.1).or.
     >			((S_alleles.eq.2).and.(any_selfcompatible.eq.0))) then

		number_of_clones=N

	else if (c.eq.0) then

		number_of_clones=0		

	else
		lambda=c*real(N)
		poisson=poidev(lambda,idum)
		number_of_clones=poisson

	endif

	if(number_of_clones.gt.N) number_of_clones=N

c	Print*,poisson,number_of_clones,N


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     GENERATION OF OFFSPRING GENOTYPES FROM ASEXUAL REPRODUCTION
c     BY DRAWING RANDOM PARENTAL GENOTYPES (PLUS SELECTION)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	DO offspring=1, number_of_clones




c.........gets parental clon
51		random=ran2(idum)
		clon=1+INT(random*N)

c.........fitness calculation
		if ((parent_generation(clon,2,1).eq.0)
     >       .or.(parent_generation(clon,2,2).eq.0)) then

			if ((parent_generation(clon,2,1).eq.0)
     >	       .and.(parent_generation(clon,2,2).eq.0)) then
				treshold_selection=1.0-s
			else
				treshold_selection=1.0-s*h
			endif

c.........comparison fitness & random number to decide on selection 
			random=ran2(idum)
			if (random.gt.treshold_selection) goto 51
		endif

c.........adds "gametic contribution" to parent for calculation of effective pop size
		parent_generation(clon,13,1)=parent_generation(clon,13,1)+2
	
		do i=1,2+neutral_loci
			do ii=3,4
				parent_generation(clon,i,ii)=
     >                 parent_generation(clon,i,ii)+1			
			enddo
		enddo

c.........copies genotype from parent to the cloned offspring 
		do i=1, neutral_loci+2
			do ii=1, 2
				offspring_generation(offspring,i,ii)
     >			=parent_generation(clon,i,ii)
			enddo
		enddo

	ENDDO

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     GENERATION OF OFFSPRING GENOTYPES FROM SEXUAL REPRODUCTION
c     BY DRAWING 2 RANDOM PARENTS (RECOMBINATION, SI & SELECTION)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	DO offspring=number_of_clones+1,N

		do i=1,2+neutral_loci
			do ii=1,2
				ovule(i,ii)=0
				pollen(i,ii)=0
			enddo
		enddo



c.........gets mother
52		random=ran2(idum)
		mother=1+INT(random*N)

c.........builds ovule gamete
		do i=1,neutral_loci+2
			random=ran2(idum)
			allele=1+INT(random*2)
			ovule(i,1)=parent_generation(mother,i,allele)
			ovule(i,2)=allele+2
		enddo

c.........gets father
53		random=ran2(idum)
		father=1+INT(random*N)

c.........builds pollen gamete & checks SI
		random=ran2(idum)
		allele=1+INT(random*2)
		if (parent_generation(father,1,allele).ne.0) then
			if ( (parent_generation(father,1,allele)
     >			.eq.parent_generation(mother,1,1) )
     >			.or.( parent_generation(father,1,allele)
     >			.eq.parent_generation(mother,1,2) ) )
     >			goto 53 ! SI,get other father
		endif

		pollen(1,1)=parent_generation(father,1,allele)

		do i=2,neutral_loci+2
			random=ran2(idum)
			allele=1+INT(random*2)
			pollen(i,1)=parent_generation(father,i,allele)
			pollen(i,2)=allele+2
		enddo

c.........fitness calculation
		if ((ovule(2,1).eq.0).or.(pollen(2,1).eq.0)) then

			if ((ovule(2,1).eq.0).and.(pollen(2,1).eq.0)) then
				treshold_selection=1.0-s
			else
				treshold_selection=1.0-s*h
			endif

c.........comparison fitness & random number to decide on selection 
			random=ran2(idum)
			if (random.gt.treshold_selection) goto 52
		endif

c.........builds offspring genotype from pollen and ovule haplotypes 
		do i=1, neutral_loci+2
		  offspring_generation(offspring,i,1)=ovule(i,1)
		  offspring_generation(offspring,i,2)=pollen(i,1)
		enddo

c.........adds "gametic contribution" to parents for calculation of effective pop size
	    parent_generation(mother,13,1)=
     >                        parent_generation(mother,13,1)+1
	    parent_generation(father,13,1)=
     >                        parent_generation(father,13,1)+1

		do i=1, neutral_loci+2
		   parent_generation(mother,i,ovule(i,2))=
     >                        parent_generation(mother,i,ovule(i,2))+1
	       parent_generation(father,i,pollen(i,2))=
     >                        parent_generation(father,i,pollen(i,2))+1
		enddo


		
	ENDDO

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     PACKING ALLELES AT S LOCUS TO SMALLER NUMBERS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

70	DO i=1,S_alleles
		contador=0
		do ii=1,N
			do iii=1,2
				if (offspring_generation(ii,1,iii).eq.i) then
					contador=contador+1
				endif
			enddo
		enddo
		freq_S(i)=real(contador)/real(2*N)
	ENDDO

c	freq_S0=0
c	DO i=1,N
c		do ii=1,2
c			if (offspring_generation(i,1,ii).eq.0) then
c				freq_S0=freq_S0+1
c			endif
c		enddo
c		freq_S0 = freq_S0 / real(2 * N)
c	ENDDO


	DO i=1,S_alleles
		if (freq_S(i).eq.0) then
			do ii=S_alleles,i,-1

				S_alleles=S_alleles-1

				if (freq_S(ii).ne.0) then
					do j=1,N
						do jj=1,2

						if (offspring_generation(j,1,jj).eq.ii) then
							offspring_generation(j,1,jj)=i
						endif

						enddo
					enddo
					
					goto 70

				endif

			enddo
		endif
	ENDDO

	contador=0
	do i=1,N
		do ii=1,2
			if (offspring_generation(i,1,ii).eq.0) then
				contador=contador+1
			endif
		enddo
	enddo
	freq_S(0)=real(contador)/real(2*N)


	if (freq_S(0).eq.0) then
		any_selfcompatible=0
	else
		any_selfcompatible=1
	endif




c	do i=1,S_alleles
c
c	Write(20,'(F8.4)') freq_S(i)
c
c	enddo


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CALCULATION OF INDICES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	IF (t.ge.burn_in) then
		is_it_time=(real(t)-real(burn_in))/(real(sampling_interval))

	if (is_it_time.eq.sampled) then
		sampled=sampled+1.0


c.....SAVES POPULATION DATA ON GENEPOP FILE

	if (saves_GenePop.eq.1) then

	Print *,'writing GenePop file'

	Write(20,*)'POP'
	Do i=1,N
	Write(20,'("g",I10.10,"i",I4.4,",",6x,2I3.3
     >,6x,2I3.3,6x,2I3.3,6x,2I3.3)')
     >t,i,(offspring_generation(i,1,ii),ii=1,2),
     >(offspring_generation(i,2,ii)+1,ii=1,2),
     >(offspring_generation(i,3,ii),ii=1,2),
     >(offspring_generation(i,4,ii),ii=1,2)
	enddo

	endif

c.....end of saving population data........



	Print*,'Calculating genetic indices'

	if (freq_S(0).eq.0) then
		any_selfcompatible=0
	else
		any_selfcompatible=1
	endif

c	Print*,any_selfcompatible

c.....Writing S alleles frequencies to a file

	WRITE(10,'(I4,1x,I6,1x,I3,4001F10.3)')
     >rep,t,S_alleles,(freq_S(i),i=0,S_alleles)

c.....Calculating indices for viability locus

		F00=0
		F10=0
		F11=0
		F1=0
		F0=0
		delta=0
		Ws=0
		Wo=0
		L=0
		He=0
		Ho=0
		Fis=0
		Ne=0
	    pa=0

		do i=1,N
			if ((offspring_generation(i,2,1).eq.0)
     >		   .and.(offspring_generation(i,2,2).eq.0)) then

				F00=F00+1
				F0=F0+2

			else if ((offspring_generation(i,2,1).eq.1)
     >		   .and.(offspring_generation(i,2,2).eq.1)) then

				F11=F11+1
				F1=F1+2

			else

				F10=F10+1
				F0=F0+1
				F1=F1+1

			endif
		enddo

		pa=dble(F0)/dble(2*N)

		freq_V(1)=pa
		freq_V(2)=1.0-pa
		
		He=2.0*pa*(1.0-pa)

		Ho=dble(F10)/dble(N)

		Fis=(He-Ho)/He

		if (He.eq.0) Fis=-9

		Ws= dble(F11)/dble(N)
     >        +Ho*(0.25+0.5*(1.0-dble(s*h))+0.25*(1.0-dble(s)))
     >        +(dble(F00)/dble(N))*(1.0-dble(s))

		Wo= dble(F11)/dble(N)
     >        +Ho*(1.0-dble(s*h))
     >        +(dble(F00)/dble(N))*(1.0-dble(s))

		delta=1.0-Ws/Wo

		L=1.0-Wo

		WRITE(3,'(I7,1x,I7,1x,3F6.3,1x,F6.3,2F7.4,2E11.4)')
     >rep,t,pa,He,Ho,Fis,Ws,Wo,delta,L

c.....Calculating indices for S locus

		ne_s=0
		na=0
		sum = freq_S(0)**2.0
		if (freq_S(0).eq.0) sum_E = 0.0
		else sum_E = freq_S(0)*log(freq_S(0))

c	Print*,'It gets here A'


		DO i=1,S_alleles

			sum=sum+(freq_S(i))**2.0
			if (freq_S(i).ne.0) sum_E=sum_E + freq_S(i)*log(freq_S(i))

c	Print *,'SUM',sum

		ENDDO

c	Print*,'It gets here B'

		
		na=S_alleles+any_selfcompatible
c		if (any_selfcompatible.eq.1) na=na+1
		He_s=1.0-sum
		ne_s=sum**(-1.0)
		if (any_selfcompatible.eq.0) then
			Fis_s=(He_s-1.0)/He_s
c			also it could be calculated as: Fis_s=1.0/(ne_s-1.0)
		else
			contador=0
			do i=1,N
				if( offspring_generation(i,1,1)
     >.eq.offspring_generation(i,1,2) ) then 
					contador=contador+1					
				endif
			enddo
			Ho_s=1-dble(contador)/dble(N)
			Fis_s=(He_s-Ho_s)/He_s
		endif

c	Print*,'It gets here C'


		Shannon_H=((-1.0)*sum_E)

c	Print*,s_alleles,na

		if (na.eq.1) then
			Shannon_E=-9
		else
			Shannon_E=((-1.0)*sum_E)/log(real(na))
		endif 

c	Print*,'It gets here D'


		if (He_s.eq.0) Fis_s=-9

		if (freq_S(0).eq.1) total_selfcompatible_pop=1
		else total_selfcompatible_pop=0



c	write(20,*)s_alleles,na


c	Print*,ne_s

	WRITE(4,'(I7,1x,I7,1x,I7,1x,F12.3,1x,F6.3,1x,F6.3,1x,F6.3)')
     >rep,t,na,ne_s,He_s,Fis_s,Shannon_E

c.....Calculating indices for neutral loci

		DO i=1,neutral_loci
			na_n(i)=0.0
			He_n(i)=0.0
			Ho_n(i)=0.0
			sum=0.0
			sum_E=0.0
			do ii=1,k
				contador=0
				do iii=1,N
					do iiii=1,2
					  if(offspring_generation(iii,i+2,iiii).eq.ii)then 
						contador=contador+1					
					  endif
					enddo
				enddo
				freq_N(i,ii)=dble(contador)/dble(2*N)

				sum=sum+(freq_N(i,ii))**2.0

				if (freq_N(i,ii).ne.0) then
					na_n(i)=na_n(i)+1
					sum_E=sum_E + freq_N(i,ii)*log(freq_N(i,ii))
				endif 

			enddo

			He_n(i)=1.0-sum
			ne_n(i)=sum**(-1.0)

			Shannon_E_n(i)=((-1.0)*sum_E)/log(real(na_n(i)))
			
			if (na_n(i).eq.1) Shannon_E_n(i)=-9

			contador=0
			do ii=1,N
				if (offspring_generation(ii,i+2,1)
     >				.ne.offspring_generation(ii,i+2,2)) then
					contador=contador+1
				endif
			enddo

			Ho_n(i)=dble(contador)/dble(N)

			Fis_n(i)=(He_n(i)-Ho_n(i))/He_n(i)

			if (He_n(i).eq.0) Fis_n(i)=-9

		ENDDO

	Write(5,'(2I7,I6,1x,F8.3,1x,3F7.3,1x,F6.3)')
     >rep,t,na_n(1),ne_n(1),He_n(1),Ho_n(1),Fis_n(1),Shannon_E_n(1)
 
	Write(55,'(2I7,I6,1x,F8.3,1x,3F7.3,1x,F6.3)')
     >rep,t,na_n(2),ne_n(2),He_n(2),Ho_n(2),Fis_n(2),Shannon_E_n(2)


c.....Calculating effective population size
		sum=0.0

		do i=1, N
		sum=sum+(dble(parent_generation(i,13,1))/dble(2*N))**2.0
		enddo

		Ne=sum**(-1.0)

		do i=1, 2+neutral_loci
		
			sum=0.0
	
			do ii=1,N
				do iii=3,4

	sum=sum+(dble(parent_generation(ii,i,iii))/dble(2*N))**2.0
					
				enddo
			enddo
		
			Neff(i)=sum**(-1.0)

		enddo


	WRITE(7,'(2I7,1x,F8.3,1x,F8.3,1x,F8.3,1x,F8.3,1x,F8.3)')
     >rep,t,Ne,Neff(1),Neff(2),Neff(3),Neff(4)













c.....CALCULATES LINKAGE DISEQUILIBRIUM

		if (calculate_LD.eq.1) then

	Print*,'Calculating linkage disequilibrium, please wait'


c.....Calculates linkage disequilibrium for S locus and a) viability locus b) 1st neutral locus

			Nc_xy=0
			Nc_xz=0
			sum_F_xy=0
			sum_F_xz=0
			moy_F_xy=0
			moy_F_xz=0
			R_GGD_SN1=0
			R_GGD_SV=0


			DO x=0,S_alleles
				if (freq_S(x).ne.0) then

					if (Fis.ne.-9) then
					do y=1,2
						if (freq_V(y).ne.0) then

							T_xy=0
							Hobs_x=0
							Hobs_y=0

							do i=1,N

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+2

	if ((offspring_generation(i,1,1).ne.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+1

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).ne.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+1

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,2,1).ne.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+1

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).ne.y-1)) T_xy=T_xy+1

	if ((offspring_generation(i,1,1).ne.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,2,1).ne.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+0.5

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).ne.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).ne.y-1)) T_xy=T_xy+0.5

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).ne.x).and.
     >(offspring_generation(i,2,1).ne.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+0.5

	if ((offspring_generation(i,1,1).ne.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).ne.y-1)) T_xy=T_xy+0.5


	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).eq.x)) Hobs_x=Hobs_x+1
	
	if ((offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) Hobs_y=Hobs_y+1

							enddo							

	Hobs_x=Hobs_x/real(N)
	Hobs_y=Hobs_y/real(N)

	C_x=Hobs_x-freq_S(x)**2
	C_y=Hobs_y-freq_V(y)**2

	delta_xy=(real(T_xy)/real(N))-(2*freq_S(x)*freq_V(y)) 


c	delta_xy=(real(N)/real(N-1))
c     >*((real(T_xy)/real(N))-(2*freq_S(x)*freq_V(y))) !TO BE DELETED



	denom=((freq_S(x)*(1-freq_S(x))+C_x)
     >*(freq_V(y)*(1-freq_V(y))+C_y))


	R_xy=delta_xy/(denom**0.5)


	if(R_xy.ge.1) then
c		print*,R_xy
		R_xy=0.9999999
c		print*,R_xy
c		pause
	endif


	if(R_xy.le.-1) then
c		print*,R_xy
		R_xy=-0.9999999
c		print*,R_xy
c		pause
	endif









c	if (R_xy.gt.0) then
		F_xy=0.5*(log((1+R_xy)/(1-R_xy)))
		sum_F_xy=sum_F_xy+abs(F_xy)
		Nc_xy=Nc_xy+1
c	endif

c	Print*,''
c	Print '(2I3,2x,F4.1,2x,F8.5,2x,F8.5)',x,y,T_xy,delta_xy,R_xy
c	Print*,''
c	Print*,1+R_xy
c	Print*,1-R_xy
c	Print*,(1+R_xy)/(1-R_xy)
c	Print*,log((1+R_xy)/(1-R_xy))
c	Print*,F_xy


						endif
					enddo
					endif

					if (Fis_n(1).ne.-9) then
					do z=1,k
						if (freq_N(1,z).ne.0) then


							T_xz=0
							Hobs_x=0
							Hobs_z=0

							do i=1,N


	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,3,1).eq.z).and.
     >(offspring_generation(i,3,2).eq.z)) T_xz=T_xz+2

	if ((offspring_generation(i,1,1).ne.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,3,1).eq.z).and.
     >(offspring_generation(i,3,2).eq.z)) T_xz=T_xz+1

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).ne.x).and.
     >(offspring_generation(i,3,1).eq.z).and.
     >(offspring_generation(i,3,2).eq.z)) T_xz=T_xz+1

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,3,1).ne.z).and.
     >(offspring_generation(i,3,2).eq.z)) T_xz=T_xz+1

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,3,1).eq.z).and.
     >(offspring_generation(i,3,2).ne.z)) T_xz=T_xz+1

	if ((offspring_generation(i,1,1).ne.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,3,1).ne.z).and.
     >(offspring_generation(i,3,2).eq.z)) T_xz=T_xz+0.5

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).ne.x).and.
     >(offspring_generation(i,3,1).eq.z).and.
     >(offspring_generation(i,3,2).ne.z)) T_xz=T_xz+0.5

	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).ne.x).and.
     >(offspring_generation(i,3,1).ne.z).and.
     >(offspring_generation(i,3,2).eq.z)) T_xz=T_xz+0.5

	if ((offspring_generation(i,1,1).ne.x).and.
     >(offspring_generation(i,1,2).eq.x).and.
     >(offspring_generation(i,3,1).eq.z).and.
     >(offspring_generation(i,3,2).ne.z)) T_xz=T_xz+0.5










	if ((offspring_generation(i,1,1).eq.x).and.
     >(offspring_generation(i,1,2).eq.x)) Hobs_x=Hobs_x+1

	if ((offspring_generation(i,3,1).eq.z).and.
     >(offspring_generation(i,3,2).eq.z)) Hobs_z=Hobs_z+1

							enddo							

	Hobs_x=Hobs_x/real(N)
	Hobs_z=Hobs_z/real(N)

	C_x=Hobs_x-freq_S(x)**2
	C_z=Hobs_z-freq_N(1,z)**2

	delta_xz=(real(T_xz)/real(N))-(2*freq_S(x)*freq_N(1,z))

c	delta_xz=(real(N)/real(N-1))*
c     >((real(T_xz)/real(N))-(2*freq_S(x)*freq_N(1,z))) !TO BE DELETED


	denom=((freq_S(x)*(1-freq_S(x))+C_x)
     >*(freq_N(1,z)*(1-freq_N(1,z))+C_z))

	R_xz=delta_xz/(denom**0.5)



	if(R_xz.ge.1) then
c		print*,R_xz
		R_xz=0.9999999
c		print*,R_xz
c		pause
	endif


	if(R_xz.le.-1) then
c		print*,R_xz
		R_xz=-0.9999999
c		print*,R_xz
c		pause
	endif



c	if (R_xz.gt.0) then 
		F_xz=0.5*log((1+R_xz)/(1-R_xz))
		sum_F_xz=sum_F_xz+abs(F_xz)
		Nc_xz=Nc_xz+1
c	endif

c	Print*,''
c	Print*,x,z,T_xz,R_xz


						endif
					enddo
					endif
				endif
			ENDDO

	if (Fis.ne.-9) then
		moy_F_xy=sum_F_xy/real(Nc_xy)
		a=((exp(2*moy_F_xy))-1)
		b=((exp(2*moy_F_xy))+1)
		R_GGD_SV=a/b
c	Print*,Nc_xy,moy_F_xy,a,b,R_GGD_SV
	else
		R_GGD_SV=-9
	endif

	if (Fis_n(1).ne.-9) then
		moy_F_xz=sum_F_xz/real(Nc_xz)
		R_GGD_SN1=(exp(2*moy_F_xz)-1)/(exp(2*moy_F_xz)+1)
	else
		R_GGD_SN1=-9
	endif

c.....Calculates linkage disequilibrium for 1st neutral locus and a) viability locus b) 2nd neutral locus


	IF (Fis_n(1).ne.-9) THEN


			Nc_xy=0
			Nc_xz=0
			sum_F_xy=0
			sum_F_xz=0
			moy_F_xy=0
			moy_F_xz=0
			R_GGD_N1V=0
			R_GGD_N1N2=0

			DO x=1,k
				if (freq_N(1,x).ne.0) then

					if (Fis.ne.-9) then
					do y=1,2
						if (freq_V(y).ne.0) then


							T_xy=0
							Hobs_x=0
							Hobs_y=0

							do i=1,N

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+2

	if ((offspring_generation(i,3,1).ne.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+1

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).ne.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+1

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,2,1).ne.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+1

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).ne.y-1)) T_xy=T_xy+1

	if ((offspring_generation(i,3,1).ne.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,2,1).ne.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+0.5

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).ne.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).ne.y-1)) T_xy=T_xy+0.5

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).ne.x).and.
     >(offspring_generation(i,2,1).ne.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) T_xy=T_xy+0.5

	if ((offspring_generation(i,3,1).ne.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).ne.y-1)) T_xy=T_xy+0.5



	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).eq.x)) Hobs_x=Hobs_x+1
	
	if ((offspring_generation(i,2,1).eq.y-1).and.
     >(offspring_generation(i,2,2).eq.y-1)) Hobs_y=Hobs_y+1

							enddo							



	Hobs_x=Hobs_x/real(N)
	Hobs_y=Hobs_y/real(N)



	C_x=Hobs_x-freq_N(1,x)**2
	C_y=Hobs_y-freq_V(y)**2



	delta_xy=(real(T_xy)/real(N))-(2*freq_N(1,x)*freq_V(y))

c	delta_xy=(real(N)/real(N-1))*
c     >((real(T_xy)/real(N))-(2*freq_N(1,x)*freq_V(y)))



	denom=((freq_N(1,x)*(1-freq_N(1,x))+C_x)
     >*(freq_V(y)*(1-freq_V(y))+C_y))
	


	R_xy=delta_xy/(denom**0.5)


c	Print*,x,y,T_xy,delta_xy,R_xy
c	Print*,denom,denom**0.5


	if(R_xy.ge.1) then
c		print*,R_xy
		R_xy=0.9999999
c		print*,R_xy
c		pause
	endif


	if(R_xy.le.-1) then
c		print*,R_xy
		R_xy=-0.9999999
c		print*,R_xy
c		pause
	endif


c	if (R_xy.gt.0) then
		F_xy=0.5*log((1+R_xy)/(1-R_xy))
		sum_F_xy=sum_F_xy+abs(F_xy)
		Nc_xy=Nc_xy+1
c	endif


c	Print*,''
c	Print*,x,y,T_xy,R_xy



						endif
					enddo
					endif

					if (Fis_n(2).ne.-9) then
					do z=1,k
						if (freq_N(2,z).ne.0) then


							T_xz=0
							Hobs_x=0
							Hobs_z=0

							do i=1,N



	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,4,1).eq.z).and.
     >(offspring_generation(i,4,2).eq.z)) T_xz=T_xz+2

	if ((offspring_generation(i,3,1).ne.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,4,1).eq.z).and.
     >(offspring_generation(i,4,2).eq.z)) T_xz=T_xz+1

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).ne.x).and.
     >(offspring_generation(i,4,1).eq.z).and.
     >(offspring_generation(i,4,2).eq.z)) T_xz=T_xz+1

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,4,1).ne.z).and.
     >(offspring_generation(i,4,2).eq.z)) T_xz=T_xz+1

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,4,1).eq.z).and.
     >(offspring_generation(i,4,2).ne.z)) T_xz=T_xz+1

	if ((offspring_generation(i,3,1).ne.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,4,1).ne.z).and.
     >(offspring_generation(i,4,2).eq.z)) T_xz=T_xz+0.5

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).ne.x).and.
     >(offspring_generation(i,4,1).eq.z).and.
     >(offspring_generation(i,4,2).ne.z)) T_xz=T_xz+0.5

	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).ne.x).and.
     >(offspring_generation(i,4,1).ne.z).and.
     >(offspring_generation(i,4,2).eq.z)) T_xz=T_xz+0.5

	if ((offspring_generation(i,3,1).ne.x).and.
     >(offspring_generation(i,3,2).eq.x).and.
     >(offspring_generation(i,4,1).eq.z).and.
     >(offspring_generation(i,4,2).ne.z)) T_xz=T_xz+0.5



	if ((offspring_generation(i,3,1).eq.x).and.
     >(offspring_generation(i,3,2).eq.x)) Hobs_x=Hobs_x+1

	if ((offspring_generation(i,4,1).eq.z).and.
     >(offspring_generation(i,4,2).eq.z)) Hobs_z=Hobs_z+1

							enddo							

	Hobs_x=Hobs_x/real(N)
	Hobs_z=Hobs_z/real(N)

	C_x=Hobs_x-freq_N(1,x)**2
	C_z=Hobs_z-freq_N(2,z)**2

	delta_xz=(real(T_xz)/real(N))-(2*freq_N(1,x)*freq_N(2,z))

c	delta_xz=(real(N)/real(N-1))*
c     >((real(T_xz)/real(N))-(2*freq_N(1,x)*freq_N(2,z)))




	denom=((freq_N(1,x)*(1-freq_N(1,x))+C_x)
     >*(freq_N(2,z)*(1-freq_N(2,z))+C_z))



	R_xz=delta_xz/(denom**0.5)


	if(R_xz.ge.1) then
c		print*,R_xz
		R_xz=0.9999999
c		print*,R_xz
c		pause
	endif


	if(R_xz.le.-1) then
c		print*,R_xz
		R_xz=-0.9999999
c		print*,R_xz
c		pause
	endif


c	if (R_xz.gt.0) then
		F_xz=0.5*log((1+R_xz)/(1-R_xz))
		sum_F_xz=sum_F_xz+abs(F_xz)
		Nc_xz=Nc_xz+1
c	endif

c	Print*,''
c	Print*,x,z,T_xz,R_xz



						endif
					enddo
					endif
				endif
			ENDDO

	if (Fis.ne.-9) then
		moy_F_xy=sum_F_xy/real(Nc_xy)
		R_GGD_N1V=(exp(2*moy_F_xy)-1)/(exp(2*moy_F_xy)+1)
	else
		R_GGD_N1V=-9
	endif

	if (Fis_n(2).ne.-9) then
		moy_F_xz=sum_F_xz/real(Nc_xz)
		R_GGD_N1N2=(exp(2*moy_F_xz)-1)/(exp(2*moy_F_xz)+1)
	else
		R_GGD_N1N2=-9
	endif

	ELSE
		R_GGD_N1V=-9
		R_GGD_N1N2=-9
	ENDIF


	WRITE(8,'(2I7,4x,E11.5,3x,E11.5,3x,E11.5,4x,E11.5)')
     >rep,t,R_GGD_SV,R_GGD_SN1,R_GGD_N1V,R_GGD_N1N2

c.....disequilibrilium

c			DO x=1,S_alleles
c			 if (freq_S(x).gt.0) then
c			  do xx=x,S_alleles
c				if (freq_S(xx).gt.0) then
c				 DO y=1,2
c				  if (freq_V(y).gt.0) then
c				   do yy=y,2
c					if (freq_V(yy).gt.0) then
c					 DO z=1,k
c					  if (freq_N(1,z).gt.0) then
c					   do zz=z,k
c						if (freq_N(1,zz).gt.0) then
c						 DO w=1,k
c						  if (freq_N(2,w).gt.0) then
c						   do ww=w,k
c							if (freq_N(2,ww).gt.0) then
c
c
c	Print*,x,xx,y,yy,z,zz,w,ww
c
c	Print*,freq_S(x),freq_S(xx),freq_V(y),freq_V(yy),freq_N(1,z)
c     >,freq_N(1,zz),freq_N(2,w),freq_N(2,ww)
c
c	Print*,freq_S(x)*freq_S(xx)*freq_V(y)*freq_V(yy)*freq_N(1,z)
c     >*freq_N(1,zz)*freq_N(2,w)*freq_N(2,ww)
c
c
c							endif
c						   enddo
c						  endif
c						 ENDDO
c						endif
c					   enddo
c					  endif
c					 ENDDO
c					endif
c				   enddo
c				  endif
c				 ENDDO
c				endif
c			  enddo
c			 endif
c			ENDDO

		endif

	Print*,'Linkage Disequilibium calculations finished'

c	pause !TO BE DELETED


ccccccccccccccccccccccccccccc
c	end for LD calculations
ccccccccccccccccccccccccccccc






	endif
	ENDIF





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MUTATION PROCESS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c.....Number of mutants at S locus toward self-INcompatible alleles
	DO i=1,N
	do ii=1,2 
		random=ran2(idum)
		if (offspring_generation(i,1,ii).ne.0) then
			if (random.lt.u_S) then
				S_alleles=S_alleles+1
				offspring_generation(i,1,ii)=S_alleles
			endif
		endif
	enddo
	ENDDO
c.....Number of mutants at S locus toward self-compatible alleles
	DO i=1,N
	do ii=1,2 
		random=ran2(idum)
		if (offspring_generation(i,1,ii).ne.0) then
			if (random.lt.u_S0) then
				offspring_generation(i,1,ii)=0
			endif
		endif
	enddo
	ENDDO




c	lambda=u_S*real(2*N)
c	poisson=poidev(lambda,idum)
c	number_of_mutants_S=poisson
c.....place S locus mutations at random individuals
c	DO i=1,number_of_mutants_S
c		mutant=0
c		allele=0
c		random=ran2(idum)
c		allele=1+INT(random*2)
c		if (allele.eq.1) then
c			mutant=1+INT(random*2*N)
c		else
c			mutant=1+INT(random*2*N)-N
c		endif
c		S_alleles=S_alleles+1
c		offspring_generation(mutant,1,allele)=S_alleles
c	ENDDO

c.....mutations at viability locus
	DO i=1,N
	do ii=1,2 
		random=ran2(idum)
		if (offspring_generation(i,2,ii).eq.0) then
			if (random.lt.u_2) offspring_generation(i,2,ii)=1
		else
			if (random.lt.u_1) offspring_generation(i,2,ii)=0
		endif
	enddo
	ENDDO

c.....Number of mutants at neutral loci
	lambda=u_neutral*real(2*N*neutral_loci)
	poisson=poidev(lambda,idum)
	number_of_mutants_N=poisson

	DO i=1,number_of_mutants_N
		mutant=0
		allele=0
		locus=0

		random1=ran2(idum)
		locus=3+INT(random1*neutral_loci)

		random2=ran2(idum)
		allele=1+INT(random2*2)
		if (allele.eq.1) then
			mutant=1+INT(random2*2*N)
		else
			mutant=1+INT(random2*2*N)-N
		endif

100		random3=ran2(idum)

		k_allele=1+INT(random3*k)

		if(offspring_generation(mutant,locus,allele).eq.k_allele)
     >goto 100

		offspring_generation(mutant,locus,allele)=k_allele
	ENDDO

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Copies offspring 2 parents and makes offsprint 2 zeros
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	DO i=1,N
		do ii=1,2+neutral_loci
			do iii=1,2
				parent_generation(i,ii,iii)=
     >				offspring_generation(i,ii,iii)
				offspring_generation(i,ii,iii)=0
				
			enddo
			do iii=3,4
				parent_generation(i,ii,iii)=0
			enddo
		enddo
		parent_generation(i,13,1)=0		
	ENDDO




	ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     END: LOOP GENERATIONS (above this)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     END: LOOP REPLICATES (above this)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	Write(1,*)'SIplusA finished normally'
	Print*,'SIplusA finished normally'
c	Print*,'Press RETURN to exit'
c	pause

	STOP
	END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     END OF THE PROGRAM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc








ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     FUNCTIONS FOR RANDOM NUMBER:
c     Uniform distribution
c     Poisson distribution
c     Binomial distribution
c
c	from "Numerical recipes in F77" by Press et al
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....Function poidev from "N. R. in F77" Chapter 7, pag 285
c.....USES gammln & ran2 (ran1 in the book)

	FUNCTION bnldev(pp,n,idum)
	INTEGER idum,n
	REAL bnldev,pp,PI
	PARAMETER (PI=3.141592654)
	INTEGER j,nold
	REAL am,em,en,g,oldg,p,pc,pclog,plog,pold,sq,t,y,gammln,ran2
	SAVE nold,pold,pc,plog,pclog,en,oldg
	DATA nold /-1/, pold /-1./
	if(pp.le.0.5)then
		p=pp
	else
		p=1.-pp
	endif
	am=n*p
	if (n.lt.25)then
		bnldev=0.
		do j=1,n
			if(ran2(idum).lt.p)bnldev=bnldev+1.
		enddo 
	else if (am.lt.1.) then
		g=exp(-am)
		t=1.
		do j=0,n
			t=t*ran2(idum)
			if (t.lt.g) goto 1000
		enddo
		j=n
1000		bnldev=j
	else
		if (n.ne.nold) then
			en=n
			oldg=gammln(en+1.)
			nold=n
		endif
		if (p.ne.pold) then		
			pc=1.-p
			plog=log(p)
			pclog=log(pc)
			pold=p
		endif
		sq=sqrt(2.*am*pc)
1500		y=tan(PI*ran2(idum))
		em=sq*y+am
		if (em.lt.0..or.em.ge.en+1.) goto 1500
		em=int(em)
		t=1.2*sq*(1.+y**2)*exp(oldg-gammln(em+1.)
     >	 -gammln(en-em+1.)+em*plog+(en-em)*pclog)
		if (ran2(idum).gt.t) goto 1500
		bnldev=em
	endif
	if (p.ne.pp) bnldev=n-bnldev
	return
	END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....Function poidev from "N. R. in F77" Chapter 7, pag 284
c.....USES gammln & ran2 (ran1 in the book)

	FUNCTION poidev(xm,idum)
	INTEGER idum
	REAL poidev,xm,PI
	PARAMETER (PI=3.141592654)
	REAL alxm,em,g,oldm,sq,t,y,gammln,ran2
	SAVE alxm,g,oldm,sq
	DATA oldm /-1./
	if (xm.lt.12.)then
		if (xm.ne.oldm) then
			oldm=xm
			g=exp(-xm)
		endif
		em=-1
		t=1.
2000		em=em+1.
		t=t*ran2(idum)
		if (t.gt.g) goto 2000
	else
		if (xm.ne.oldm) then
			oldm=xm
			sq=sqrt(2.*xm)
			alxm=log(xm)
			g=xm*alxm-gammln(xm+1.)
		endif
3000		y=tan(PI*ran2(idum))
		em=sq*y+xm
		if (em.lt.0.) goto 3000
		em=int(em)
		t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
		if (ran2(idum).gt.t) goto 3000
	endif
	poidev=em
	return
	END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....Function ran2 from "N R in F77" Chapter 7, pag 272

	FUNCTION ran2(idum)
	INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	REAL ran2,AM,EPS,RNMX
	PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     > IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     > IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
	INTEGER idum2,j,k,iv(NTAB),iy
	SAVE iv,iy,idum2
	DATA idum2/123456789/, iv/NTAB*0/, iy/0/
	if (idum.le.0) then
		idum=max(-idum,1)
		idum2=idum
		do j=NTAB+8,1,-1
			k=idum/IQ1
			idum=IA1*(idum-k*IQ1)-k*IR1
			if (idum.lt.0) idum=idum+IM1
			if (j.le.NTAB) iv(j)=idum
		enddo
		iy=iv(1)
	endif
	k=idum/IQ1
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1
	k=idum2/IQ2
	idum2=IA2*(idum2-k*IQ2)-k*IR2
	if (idum2.lt.0) idum2=idum2+IM2
	j=1+iy/NDIV
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+IMM1
	ran2=min(AM*iy,RNMX)
	return
	END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....Function gammln from "N. R. in F77" Chapter 6, pag 207

	FUNCTION gammln(xx)
	REAL gammln,xx
	INTEGER j
	DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
	SAVE cof,stp
	DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     > 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     > -.5395239384953d-5,2.5066282746310005d0/
	x=xx
	y=x
	tmp=x+5.5d0
	tmp=(x+0.5d0)*log(tmp)-tmp
	ser=1.000000000190015d0
	do j=1,6
		y=y+1.d0
		ser=ser+cof(j)/y
	enddo
	gammln=tmp+log(stp*ser/x)
	return
	END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
