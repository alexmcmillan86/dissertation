%08TTD001: Final Year Project
%SPECTRAL PROPERTIES OF A VIBRATING TIMOSHENKO BEAM
%Date: 17/12/2008
%Author: Alexander McMillan (A546840)
%Supervisor: Dr. Andrew Watson
%Loughborough University, Department of Aeronautical and Automotive Engineering
%%
%-------------------------------------------------------------------------%
%1. PROGRAM INFORMATION
%-------------------------------------------------------------------------%
%NOTE: As default, code is run in double precision
clc; %Clears command window
clear all; %Clears workspace
disp('>> LOUGHBOROUGH UNIVERSITY, DEPARTMENT OF AUTOMOTIVE AND AERONAUTICAL ENGINEERING');
disp('>> PROGRAM TO COMPUTE NATURAL FREQUENCIES FOR VIBRATING TIMOSHENKO BEAMS');
%%
%-------------------------------------------------------------------------%
%2. READING INPUT DATA
%-------------------------------------------------------------------------%
[NUMERIC]=xlsread('INPUT_DATA.xlsm'); %Extracts data from excel macro enabled input file
NP=NUMERIC(1,1); %Number of problems
TYPE=input('\nPlease Select Frequency Type (0=Dimensionless,1=Radians,2=Hertz): ');
if TYPE~=0&TYPE~=1&TYPE~=2;
warndlg('Selection not recognised!'); %Displays warning message
clc; %Clears command window
clear all; %Clears workspace
run Timoshenko_Beam_Vibration; %Creates break point
end
for NP=1:NP;
fprintf('\n\nPROBLEM NO. :%5.0f\n\n',NP);
CN=1+(NP-1)*8; %Calculates column number to read data
INPUT=NUMERIC(:,CN); %Assigns input data to an array
AA0=INPUT(5,:); %Area of LH end of beam [m^2]
AI0=INPUT(6,:); %Second moment of area [m^4]
AL=INPUT(7,:); %Length of beam [m]
TR=INPUT(8,:); %Taper ratio
E=INPUT(9,:); %Young's modulus [GPa]
G=INPUT(10,:); %Shear modulus [GPa]
RHO=INPUT(11,:); %Density [kg/m^3]
SF=INPUT(12,:); %Section shape factor (shear co-efficient)
JR=floor(INPUT(16,:)); %Number of eigenvalues required
NS=floor(INPUT(17,:)); %Number of segments
CV=floor(INPUT(18,:)); %Convergence value
RI=floor(INPUT(22,:)); %Rotary inertia
SI=floor(INPUT(23,:)); %Shear deformation
P=INPUT(24,:); %Axial load [N]
AKD1=INPUT(28,:); %Lateral spring stiffness at LH end
AKT1=INPUT(29,:); %Longitudinal spring stiffness at LH end
AKD2=INPUT(30,:); %Lateral spring stiffness at RH end
AKT2=INPUT(31,:); %Longitudinal spring stiffness at RH end
%%
-------------------------------------------------------------------------%
%3. VIBRATION THEORY
%-------------------------------------------------------------------------%
if TR==0; %Sets number of segments to 1 for no taper
NS=1;
end
if RI==1&SI==1; %Rotary inertia and shear deformation considered
disp('>> TIMOSHENKO BEAM THEORY');
disp(' ');
elseif RI==1&SI==0; %Rotary inertia considered
disp('>> RAYLEIGH BEAM THEORY');
disp(' ');
else RI==0&SI==0; %No rotary inertia or shear deformation considered
disp('>> BERNOULLI-EULER BEAM THEORY');
disp(' ');
end
%%
%-------------------------------------------------------------------------%
%4. MODELLING BEAM
%-------------------------------------------------------------------------%
tic; %Starts timing
RG=sqrt(AI0/AA0); %Radius of gyration
RGBL=RG/AL; %Radius of gyration to beam length ratio
TRN=1+TR; %Truncation ratio
fprintf('Radius of Gyration to Beam Length = %1.2f\n',RGBL);
fprintf('Truncation Ratio = %1.2f\n',TRN);
fprintf('Section Shape Factor = %1.2f\n',SF);
fprintf('Convergence Accuracy = %8.0f\n',CV);
for JR=1:JR;
FU=1E10; %Upper frequency limit [rads]
FL=0; %Lower frequency limit [rads]
TF=1; %Trial frequency value [rads]
SL=AL/(floor(NS)); %Segment length [m]
Q=0; %Counting variable
while CV*(TF-FL)>=TF; %Identifies convergence of trial frequency
Q=Q+1; %Counting term
IG=2; %No. of pivotal rows in gauss elimination
JB=0; %No. of roots passed
JE=0;
JD=0;
A(1)=AKD1; %LH end boundary condition term
A(2)=0;
A(5)=AKT1; %LH end boundary condition term
X=-SL/2; %Segment midpoint
%%
%-------------------------------------------------------------------------%
%5. DIFFERIENTIAL EQUATION SOLUTION TERMS
%-------------------------------------------------------------------------%
for IS=1:NS; %Loop to consider every segment
X=X+SL; %Midpoint of current segment
FAC=(1+TR*X/AL)^2;
AK=E*AI0*FAC*FAC/(SL*SL);
AL2=AA0*FAC*SL*SL;
B2=RHO*AL2*TF*TF/AK; %Material properties term
R2=RI*AI0*FAC*FAC/AL2; %Rotary inertia term
S2=SI*AK/(SF*AA0*FAC*G); %Shear deflection term
P2=P/AK; %Axial load term
BT=1-B2*R2*S2; %Pivotal term
%%
%-------------------------------------------------------------------------%
%6. STIFFNESS MATRIX CASE I (STANDARD SOLUTION)
%-------------------------------------------------------------------------%
if BT>=0;
SP=1-S2*P2;
DA=P2/B2+R2*SP+S2;
T1=sqrt(B2/(2*SP));
T2=sqrt(DA*DA+4*SP*BT/B2);
APA=T1*sqrt(-DA+T2); %Alpha term
BTA=T1*sqrt(DA+T2); %Beta term
Z=(SP*BTA*BTA-B2*S2)/(BTA*SL);
H=(SP*APA*APA+B2*S2)/(APA*SL);
ETA=Z/H;
S=sin(BTA);
C=cos(BTA);
SH=sinh(APA);
CH=cosh(APA);
GMA=2*ETA*(1-CH*C)+(1-ETA*ETA)*SH*S;
T1=B2/(APA*BTA);
T2=(APA+ETA*BTA)/GMA;
A(8)=T1*T2*(SH*C+ETA*CH*S)*AK/SL;
A(1)=A(1)+A(8);
A(3)=-T1*T2*(SH+ETA*S)*AK/SL;
A(10)=T2*(CH*S-ETA*SH*C)*AK*SL;
A(5)=A(5)+A(10);
A(7)=T2*(ETA*SH-S)*AK*SL;
T1=(APA-ETA*BTA)*(CH*C-1);
A(9)=-Z*SL*(T1+(BTA+ETA*APA)*SH*S)*AK/GMA;
A(2)=A(2)-A(9);
A(4)=ETA*SP*(APA*APA+BTA*BTA)*(CH-C)*AK/GMA;
A(6)=-A(4);
JB=JB+(abs(floor(BTA/pi)));
%%
%-------------------------------------------------------------------------%
%7. STIFFNESS MATRIX CASE II (ALTERNATIVE SOLUTION)
%-------------------------------------------------------------------------%
else
SP=1-S2*P2;
DA=P2/B2+R2*SP+S2;
T1=sqrt(B2/(2*SP));
T2=sqrt(DA*DA+4*SP*BT/B2);
APAP=T1*sqrt(DA-T2); %Alpha term
BTA=T1*sqrt(DA+T2); %Beta term
Z=(SP*BTA*BTA-B2*S2)/(BTA*SL);
HP=((SP*-1*APAP^2)+B2*S2)/(APAP*SL);
ETAP=Z/HP;
S=sin(BTA);
C=cos(BTA);
SAP=sin(APAP);
CAP=cos(APAP);
GMAP=2*ETAP*(1-CAP*C)+(1+ETAP*ETAP)*SAP*S;
T1=B2/(APAP*BTA);
T2=(APAP+ETAP*BTA)/GMAP;
A(8)=T1*T2*(SAP*C+ETAP*S*CAP)*AK/SL;
A(1)=A(1)+A(8);
A(3)=-1*T1*T2*(SAP+ETAP*S)*AK/SL;
A(10)=T2*(S*CAP+ETAP*C*SAP)*AK*SL;
A(5)=A(5)+A(10);
A(7)=-1*T2*(ETAP*SAP+S)*AK*SL;
T1=(APAP-ETAP*BTA)*(CAP*C-1);
A(9)=Z*SL*(T1+(BTA-ETAP*APAP)*SAP*S)*AK/GMAP;
A(2)=A(2)-A(9);
A(4)=ETAP*SP*(APAP*APAP-BTA*BTA)*(CAP-C)*AK/GMAP;
A(6)=-A(4);
JD=abs(floor(BTA/pi));
JE=abs(floor(APAP/pi))+1;
JB=JB+JE+JD;
end
%%
%-------------------------------------------------------------------------%
%8. WITTRICK-WILLIAMS ALGORITHM
%-------------------------------------------------------------------------%
if A(10)<0;
JB=JB-1;
end
if (A(10)-A(7)*A(7)/A(10))<0;
JB=JB-1;
end
if IS>=NS;
IG=3;
A(8)=A(8)+AKD2;
A(10)=A(10)+AKT2;
end
%%
%-------------------------------------------------------------------------%
%9. GAUSS ELIMINATION
%-------------------------------------------------------------------------%
for I=1:IG;
IPT=10-(4-I)*(7-I)/2;
PT=1/A(IPT);
for J=I+1:4;
IPT=IPT+1;
PIV=A(IPT)*PT;
L=IPT-1;
J1=10-(4-J)*(7-J)/2;
J2=J1+4-J;
for K=J1:J2;
L=L+1;
A(K)=A(K)-PIV*A(L);
end
end
end
if IS==NS;
IG=4;
end
for I=1:IG;
II=10-(4-I)*(7-I)/2;
if A(II)<0;
JB=JB+1;
end
end
A(1)=A(8);
A(2)=A(9);
A(5)=A(10);
end %Ends loop [FOR IS=1:NS]
%%
%-------------------------------------------------------------------------%
%10. TRIAL FREQUENCY CONVERGENCE
%-------------------------------------------------------------------------%
if JB<JR;
FL=TF; %Setting new lower limit
if FU<1E9;
TF=0.5*(FL+FU); %Resetting trial frequency
else
TF=2*TF; %Doubling of trial frequency
end
else
FU=TF; %Setting new upper limit
TF=0.5*(FL+FU); %Resetting trial frequency
end
end %Ends loop [while CV*(TF-FL >= TF)]
TF1=TF/(2*pi); %Frequency in [Hz]
LAM=TF*(sqrt((RHO*AA0*AL^4)/(E*AI0))); %Dimensionless frequency
RES(JR,:)=[JR,TF,TF1,LAM]; %Array of frequency values
end %Ends loop [for JR=1:JR]
TIME=toc; %Stops timing
fprintf('Run Time = %6.4f [seconds]\n\n',TIME);
%%
%-------------------------------------------------------------------------%
%11. WRITING OUTPUT DATA
%-------------------------------------------------------------------------%
LAM=RES(:,4); %Vector of dimensionless frequencies
TF=RES(:,2); %Vector of circular frequencies [radians]
TF1=RES(:,3); %Vector of frequencies [hertz]
FREQV=RES(:,1);
OUTPUT=INPUT; %Copying input data to output
OUTPUT([2;3;4;13;14;15;19;20;21;25;26;27])=[]; %Removing spaces in vector
OUTPUT(1,:)=NP;
format LONG G;
if TYPE==0; %Dimensionless frequencies
fprintf('>> RESULTS:\n\n');
fprintf('FREQUENCY No.\tDIMENSIONLESS FREQUENCIES\n\n');
LAM_S=[FREQV,LAM];
fprintf('%6.0f\t\t\t\t%10f\n',LAM_S');
disp(' ');
OUTPUT=[OUTPUT;NaN;NaN;LAM];
FREQ(:,NP)=[OUTPUT];
elseif TYPE==1; %Circular frequencies [radians]
fprintf('>> RESULTS:\n\n');
fprintf('FREQUENCY No.\tCIRCULAR FREQUENCY [radians]\n\n');
TF_S=[FREQV,TF];
fprintf('%6.0f\t\t\t\t%10f\n',TF_S');
disp(' ');
OUTPUT=[OUTPUT;NaN;NaN;TF];
FREQ(:,NP)=[OUTPUT];
else TYPE==2; %Frequencies [hertz]
fprintf('>> RESULTS:\n\n');
fprintf('FREQUENCY No.\tFREQUENCY [Hz]\n\n');
TF1_S=[FREQV,TF1];
fprintf('%6.0f\t\t\t\t%10f\n',TF1_S');
disp(' ');
OUTPUT=[OUTPUT;NaN;NaN;TF1];
FREQ(:,NP)=[OUTPUT];
end
end %Ends loop for number of problems
FNAME=(['OUTPUT_DATA','.xlsx']); %Output excel filename
if TYPE==0; %Dimensionless frequency
xlswrite(FNAME,FREQ,'Sheet1','A1'); %Writing output data
elseif TYPE==1; %Frequency [radians]
xlswrite(FNAME,FREQ,'Sheet1','A1'); %Writing output data
else TYPE==2; %Frequency [Hz]
xlswrite(FNAME,FREQ,'Sheet1','A1'); %Writing output data
end
disp('>> END OF PROGRAM');
disp(' ');
%--------------------------------------------------END-OF-PROGRAM--------------------------------------------------%