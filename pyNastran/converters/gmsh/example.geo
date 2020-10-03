Lx= 0.500000 ;Ly= 0.070000 ;Ep= 0.030000 ;EpBS= 0.020000 ;EpBI= 0.040000 ;EcartB= 0.035000 ;REpau= 0.005250 ;RPion= 0.003000 ;//
RLocal = 0.5*(REpau+EcartB);
xR=Lx/2;
pi=3.14159;
pi=3.14159;

dyBs=Ep+5e-3;
dyBi=-5e-3;

// characteristic length
l1 = Ep*10;
l2 = REpau*10;
//
If (Ep<=5e-3)
    nh=2;
 EndIf
If ((Ep>5e-3)&&(Ep<=10e-3))
   nh=2;
EndIf
If ((Ep>10e-3)&&(Ep<=20e-3))
   nh=3;
EndIf
If (Ep>20e-3)
   nh=4;
EndIf 
//
If (EpBS<=20e-3)
   nhBS=4;
EndIf
If ((EpBS>20e-3)&&(EpBS<=50e-3))
   nhBS=5;
EndIf   
If ((EpBS>50e-3)&&(EpBS<=75e-3))
   nhBS=6;
EndIf 
If (EpBS>75e-3)
   nhBS=7;
EndIf 
//
If (EpBI<=20e-3)
   nhBI=4;
EndIf
If ((EpBI>20e-3)&&(EpBI<=50e-3))
   nhBI=5;
EndIf   
If ((EpBI>50e-3)&&(EpBI<=75e-3))
   nhBI=6;
EndIf 
If (EpBI>75e-3)
   nhBI=7;
EndIf 
n1 = 7;
n1b=5;
n1c=4;
//
If (REpau<=5e-3)
   n1d=6;
EndIf
If ((REpau>5e-3)&&(REpau<=10e-3))
   n1d=7;
EndIf     
If ((REpau>10e-3)&&(REpau<=15e-3))
   n1d=9;
EndIf
If (REpau>15e-3)
   n1d=10;
EndIf
//
n1e=4;
n2=13;
n3=11;
//
If (REpau<=5e-3)
   n4=17;
EndIf
If ((REpau>5e-3)&&(REpau<=10e-3))
   n4=16;
EndIf     
If ((REpau>10e-3)&&(REpau<=15e-3))
   n4=15;
EndIf
If (REpau>15e-3)
   n4=14;
EndIf
//
Px=1.25;
Pxa=1.35;
//
If (REpau<=5e-3)
   Py=1.45;
EndIf
If ((REpau>5e-3)&&(REpau<=10e-3))
   Py=1.38;
EndIf     
If ((REpau>10e-3)&&(REpau<=15e-3))
   Py=1.31;
EndIf
If (REpau>15e-3)
   Py=1.25;
EndIf
//

// Quadrillage
Point(1) = {0,0,0,l1};
Point(2) = {-Ly,0,0,l1};
Point(3) = {-Ly,0,-Lx,l1};
Point(4) = {0,0,-Lx,l1};
Point(5) = {0,0,-(xR-EcartB),l2};
Point(6) = {-EcartB,0,-(xR-EcartB),l2};
Point(7) = {-Ly,0,-(xR-EcartB),l2};
Point(8) = {0,0,-(xR+EcartB),l2};
Point(9) = {-EcartB,0,-(xR+EcartB),l2};
Point(10) = {-Ly,0,-(xR+EcartB),l2};
Point(11) = {-EcartB,0,0,l2};
Point(12) = {-EcartB,0,-Lx,l2};
// Points de la zone source
Point(20)={0,0,-xR,l2};
Point(21)={0,0,-(xR+REpau),l2};
Point(22)={-REpau*Sin(0.7*pi/3),0,-(xR+REpau*Cos(0.7*pi/3)),l2};
Point(23)={-REpau*Sin(0.7*pi/3),0,-(xR-REpau*Cos(0.7*pi/3)),l2};
Point(24)={0,0,-(xR-REpau),l2};
Point(25)={0,0,-(xR+RLocal),l2};
Point(26)={-RLocal*Sin(0.7*pi/3),0,-(xR+RLocal*Cos(0.7*pi/3)),l2};
Point(27)={-RLocal*Sin(0.7*pi/3),0,-(xR-RLocal*Cos(0.7*pi/3)),l2};
Point(28)={0,0,-(xR-RLocal),l2};
//
Point(31)={0,0,-(xR+RPion),l2};
Point(32)={-RPion*Sin(0.7*pi/3),0,-(xR+RPion*Cos(0.7*pi/3)),l2};
Point(33)={-RPion*Sin(0.7*pi/3),0,-(xR-RPion*Cos(0.7*pi/3)),l2};
Point(34)={0,0,-(xR-RPion),l2};
//
Point(40)={-RPion/2*Sin(0.7*pi/3),0,-(xR+RPion/2*Cos(0.7*pi/3)),l2};
Point(41)={-RPion/2*Sin(0.7*pi/3),0,-(xR-RPion/2*Cos(0.7*pi/3)),l2};
Point(42)={0,0,-(xR+RPion/2*Cos(0.7*pi/3)),l2};
Point(43)={0,0,-(xR-RPion/2*Cos(0.7*pi/3)),l2};



Line(1)={1,5};
Line(2)={5,6};
Line(3)={6,11};
Line(4)={11,1};
Line(5)={6,7};
Line(6)={7,2};
Line(7)={2,11};
Line(8)={6,9};
Line(9)={9,10};
Line(10)={10,7};
Line(12)={9,12};
Line(13)={12,3};
Line(14)={3,10};
Line(16)={8,4};
Line(17)={4,12};
Line(18)={9,8};
Line(20)={5,28};
Line(22)={27,6};
Line(24)={26,9};
Line(26)={25,8};
Line(30)={28,24};
Line(32)={23,27};
Line(34)={22,26};
Line(36)={21,25};
Line(40)={21,31};
Line(41)={42,43};
Line(42)={34,24};
Line(43)={43,41};
Line(44)={41,40};
Line(45)={40,42};
Line(46)={41,33};
Line(47)={40,32};
Line(50)={31,42};
Line(51)={31,42};
Line(52)={43,34};
Line(53)={32,22};
Line(54)={33,23};
Circle(21) = {28,20,27};
Circle(23) = {27,20,26};
Circle(25) = {26,20,25};
Circle(31) = {24,20,23};
Circle(33) = {23,20,22};
Circle(35) = {22,20,21};
Circle(60)={32,20,31};
Circle(61)={33,20,32};
Circle(62)={34,20,33};
// Construction des surfaces
Line Loop(1) = {1,2,3,4};
Line Loop(2) = {-3,5,6,7};
Line Loop(3) = {8,9,10,-5};
Line Loop(4) = {12,13,14,-9};
Line Loop(5) = {16,17,-12,18};
Line Loop(6) = {20,21,22,-2};
Line Loop(7) = {23,24,-8,-22};
Line Loop(8) = {26,-18,-24,25};
Line Loop(9) = {30,31,32,-21};
Line Loop(10) = {33,34,-23,-32};
Line Loop(11) = {36,-25,-34,35};
Line Loop(20) = {-52,43,46,-62};
Line Loop(21) = {44,47,-61,-46};
Line Loop(22) = {45,-50,-60,-47};
Line Loop(23) = {-41,-45,-44,-43};
//
Line Loop(30) = {-42,62,54,-31};
Line Loop(31) = {61,53,-33,-54};
Line Loop(32) = {60,-40,-35,-53};
//
Plane Surface(1)={1};
Plane Surface(2)={2};
Plane Surface(3)={3};
Plane Surface(4)={4};
Plane Surface(5)={5};
Plane Surface(6)={6};
Plane Surface(7)={7};
Plane Surface(8)={8};
Plane Surface(9)={9};
Plane Surface(10)={10};
Plane Surface(11)={11};
Plane Surface(20)={20};
Plane Surface(21)={21};
Plane Surface(22)={22};
Plane Surface(23)={23};
//
Plane Surface(30)={30};
Plane Surface(31)={31};
Plane Surface(32)={32};


// Definition du maillage regle
// Zone source
Transfinite Line{2,18,21,25,31,35,60,61,62} = n1;
Transfinite Line{20,22,24,26} = n1b;
Transfinite Line{30,32,34,36} = n1e;
Transfinite Line{43,45} = n1;
Transfinite Line{40,42,46,47,50,52,53,54} = n1c;
Transfinite Line{1} = n2 Using Progression 1/Pxa;
Transfinite Line{3,6} = n2 Using Progression Pxa;
Transfinite Line{4,17} = n1;

Transfinite Line{8,10,23,33,41,44,61} = n1d;

Transfinite Line{7} = n3 Using Progression 1/Py;
Transfinite Line{5,9,13} = n3 Using Progression Py;
Transfinite Line{14} = n4 Using Progression 1./Px;
Transfinite Line{12,16} = n4 Using Progression Px;




// Definition des barres de pompage sup
Point(102) = {-Ly,dyBs,0,l1};
Point(103) = {-Ly,dyBs,-Lx,l1};
Point(106) = {-EcartB,dyBs,-(xR-EcartB),l2};
Point(107) = {-Ly,dyBs,-(xR-EcartB),l2};
Point(109) = {-EcartB,dyBs,-(xR+EcartB),l2};
Point(110) = {-Ly,dyBs,-(xR+EcartB),l2};
Point(111) = {-EcartB,dyBs,0,l2};
Point(112) = {-EcartB,dyBs,-Lx,l2};
Line(103)={106,111};
Line(105)={106,107};
Line(106)={107,102};
Line(107)={102,111};
Line(108)={106,109};
Line(109)={109,110};
Line(110)={110,107};
Line(112)={109,112};
Line(113)={112,103};
Line(114)={103,110};  
Line Loop(102) = {-103,105,106,107};
Line Loop(103) = {108,109,110,-105};
Line Loop(104) = {112,113,114,-109};
Plane Surface(102)={102};
Plane Surface(103)={103};
Plane Surface(104)={104}; 
Transfinite Line{103,106} = n2 Using Progression Pxa;
Transfinite Line{108,110} = n1d;
Transfinite Line{107} = n3 Using Progression 1/Py;
Transfinite Line{105,109,113} = n3 Using Progression Py;
Transfinite Line{114} = n4 Using Progression 1./Px;
Transfinite Line{112} = n4 Using Progression Px; 


// Definition des barres de pompage inf
Point(202) = {-Ly,dyBi,0,l1};
Point(203) = {-Ly,dyBi,-Lx,l1};
Point(206) = {-EcartB,dyBi,-(xR-EcartB),l2};
Point(207) = {-Ly,dyBi,-(xR-EcartB),l2};
Point(209) = {-EcartB,dyBi,-(xR+EcartB),l2};
Point(210) = {-Ly,dyBi,-(xR+EcartB),l2};
Point(211) = {-EcartB,dyBi,0,l2};
Point(212) = {-EcartB,dyBi,-Lx,l2};
Line(203)={206,211};
Line(205)={206,207};
Line(206)={207,202};
Line(207)={202,211};
Line(208)={206,209};
Line(209)={209,210};
Line(210)={210,207};
Line(212)={209,212};
Line(213)={212,203};
Line(214)={203,210};  
Line Loop(202) = {-203,205,206,207};
Line Loop(203) = {208,209,210,-205};
Line Loop(204) = {212,213,214,-209};
Plane Surface(202)={202};
Plane Surface(203)={203};
Plane Surface(204)={204}; 
Transfinite Line{203,206} = n2 Using Progression Pxa;
Transfinite Line{208,210} = n1d;
Transfinite Line{207} = n3 Using Progression 1/Py;
Transfinite Line{205,209,213} = n3 Using Progression Py;
Transfinite Line{214} = n4 Using Progression 1./Px;
Transfinite Line{212} = n4 Using Progression Px; 

// Extrusion
//Layers{ {8,2}, {0.5,1} };  Layers{{Ceil(nh/2),Floor(nh/2)},{0.33,1}}

Transfinite Surface{9} = {24,23,27,28};
Recombine Surface{9};
cc1[]=Extrude {0,Ep,0} {
   Surface{9};  Layers{nh,1}; Recombine;
 };
Transfinite Surface{10} = {23,22,26,27};
Recombine Surface{10};
cc2[]=Extrude {0,Ep,0} {
   Surface{10}; Layers{nh,1}; Recombine;
 }; 
Transfinite Surface{11} = {21,25,26,22};
Recombine Surface{11};
cc3[]=Extrude {0,Ep,0} {
   Surface{11}; Layers{nh,1}; Recombine;
 };
 
Transfinite Surface{20} = {34,43,41,33};
Recombine Surface{20};
s1[]=Extrude {0,Ep,0} {
   Surface{20}; Layers{nh,1}; Recombine;
 } ;
Transfinite Surface{21} = {41,40,32,33};
Recombine Surface{21};
s2[]=Extrude {0,Ep,0} {
   Surface{21}; Layers{nh,1}; Recombine;
 } ;
Transfinite Surface{22} = {42,31,32,40};
Recombine Surface{22};
s3[]=Extrude {0,Ep,0} {
   Surface{22}; Layers{nh,1}; Recombine;
 } ;
Transfinite Surface{23} = {43,42,40,41};
Recombine Surface{23};
s4[]=Extrude {0,Ep,0} {
   Surface{23}; Layers{nh,1}; Recombine;
 } ;  
//
Transfinite Surface{30} = {34,33,23,24};
Recombine Surface{30};
s5[]=Extrude {0,Ep,0} {
   Surface{30}; Layers{nh,1}; Recombine;
 } ;  
Transfinite Surface{31} = {33,32,22,23};
Recombine Surface{31};
s6[]=Extrude {0,Ep,0} {
   Surface{31}; Layers{nh,1}; Recombine;
 } ;  
Transfinite Surface{32} = {31,21,22,32};
Recombine Surface{32};
s7[]=Extrude {0,Ep,0} {
   Surface{32}; Layers{nh,1}; Recombine;
 } ; 
// Autres regions  
 Transfinite Surface{6} = {5,28,27,6};
 Recombine Surface{6};
 a1[]=Extrude {0,Ep,0} {
    Surface{6}; Layers{nh,1}; Recombine;
  } ;
 Transfinite Surface{7} = {27,26,9,6};
 Recombine Surface{7};
 a2[]=Extrude {0,Ep,0} {
    Surface{7}; Layers{nh,1}; Recombine;
  };
 Transfinite Surface{8} = {25,8,9,26};
 Recombine Surface{8};
 a3[]=Extrude {0,Ep,0} {
    Surface{8}; Layers{nh,1}; Recombine;
  } ;  
// Autres regions 

Transfinite Surface{1} = {1,5,6,11};
Recombine Surface{1};
a4[]=Extrude {0,Ep,0} {
   Surface{1}; Layers{nh,1}; Recombine;
 }; 
Transfinite Surface{2} = {11,6,7,2};
Recombine Surface{2};
a5[]=Extrude {0,Ep,0} {
   Surface{2}; Layers{nh,1}; Recombine;
 };   
Transfinite Surface{3} = {6,9,10,7};
Recombine Surface{3};
a6[]=Extrude {0,Ep,0} {
   Surface{3}; Layers{nh,1}; Recombine;
 };  
Transfinite Surface{4} = {9,12,3,10};
Recombine Surface{4};
a7[]=Extrude {0,Ep,0} {
   Surface{4}; Layers{nh,1}; Recombine;
 }; 
Transfinite Surface{5} = {8,4,12,9};
Recombine Surface{5};
a8[]=Extrude {0,Ep,0} {
   Surface{5}; Layers{nh,1}; Recombine;
 } ;    
 
// Barre de pompage sup 
Transfinite Surface{102} = {111,106,107,102};
Recombine Surface{102};
bs1[]=Extrude {0,EpBS,0} {
   Surface{102}; Layers{{Ceil(nhBS/2),Floor(nhBS/2)},{0.33,1}}; Recombine;
};   
Transfinite Surface{103} = {106,109,110,107};
Recombine Surface{103};
bs2[]=Extrude {0,EpBS,0} {
   Surface{103}; Layers{{Ceil(nhBS/2),Floor(nhBS/2)},{0.33,1}}; Recombine;
 };  
Transfinite Surface{104} = {109,112,103,110};
Recombine Surface{104};
bs3[]=Extrude {0,EpBS,0} {
   Surface{104}; Layers{{Ceil(nhBS/2),Floor(nhBS/2)},{0.33,1}}; Recombine;
};  
// Barre de pompage inf 
Transfinite Surface{202} = {211,206,207,202};
Recombine Surface{202};
bi1[]=Extrude {0,-EpBI,0} {
   Surface{202}; Layers{{Ceil(nhBI/2),Floor(nhBI/2)},{0.33,1}}; Recombine;
 };   
Transfinite Surface{203} = {206,209,210,207};
Recombine Surface{203};
bi2[]=Extrude {0,-EpBI,0} {
   Surface{203}; Layers{{Ceil(nhBI/2),Floor(nhBI/2)},{0.33,1}}; Recombine;
 };  
Transfinite Surface{204} = {209,212,203,210};
Recombine Surface{204};
bi3[]=Extrude {0,-EpBI,0} {
   Surface{204}; Layers{{Ceil(nhBI/2),Floor(nhBI/2)},{0.33,1}}; Recombine;
 };     

 
 


Physical Volume(100) = {s1[1],s2[1],s3[1],s4[1],s5[1],s6[1],s7[1]}; // Source
Physical Volume(101) = {a1[1],a2[1],a3[1],a4[1],a5[1],a6[1],a7[1],a8[1],cc1[1],cc2[1],cc3[1]}; // Tou_sauf_Source
Physical Volume(102) = {bs1[1],bs2[1],bs3[1]}; // BP_SUP
Physical Volume(202) = {bi1[1],bi2[1],bi3[1]}; // BP_INF
Physical Surface(300)={543,631,697};//Entree
Physical Surface(401)={544,566,588};//T_SSUP
Physical Surface(402)={102,103,104};//BP_SSUP
Physical Surface(501)={2,3,4};//T_SINF
Physical Surface(502)={202,203,204};//BP_SINF
Physical Surface(600)={271,253,235};//MOD_LOCAL
Physical Surface(700)={301,319,341};//PION
Physical Surface(800)={30,31,32,390,412,434};//EPAULEMENTS

