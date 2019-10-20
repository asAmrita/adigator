% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function p = myprogram_ADiGatorJac(x,d,m)
global ADiGator_myprogram_ADiGatorJac
if isempty(ADiGator_myprogram_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_myprogram_ADiGatorJac.myprogram_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % FIT -- Given x and d, fit() returns p
%User Line: % such that norm(V*p-d) = min, where
%User Line: % V = [1, x, x.^2, ... x.^(m-1)].
dim_x.f = size(x.f,1);
%User Line: dim_x = size(x, 1);
cadaconditional1 = lt(dim_x.f,m);
%User Line: cadaconditional1 = dim_x < m;
V.f = ones(dim_x.f,1);
%User Line: V = ones(dim_x, 1);
cada1f1 = m - 1;
cadaforvar1.f = 1:cada1f1;
%User Line: cadaforvar1 = 1 : (m-1);
V.dx = zeros(200,1);
V.f(100,3) = 0;
for cadaforcount1 = 1:2
    count.f = cadaforvar1.f(:,cadaforcount1);
    %User Line: count = cadaforvar1(:,cadaforcount1);
    cada1f1dx = count.f.*x.f(:).^(count.f-1).*x.dx;
    cada1f1dx((x.f(:) == 0 & x.dx == 0) | count.f == 0) = 0;
    cada1f1 = x.f.^count.f;
    cada1f3dx = 2.*cada1f1dx;
    cada1f3 = 2*cada1f1;
    cada1f4dx = cada1f3dx;
    cada1f4 = 1 + cada1f3;
    V.dx = V.dx(Gator1Data.Index4,1);
    V.f = V.f(:,1:2);
    cada1td1 = zeros(200,1);
    cada1td1(logical(Gator1Data.Index1(:,cadaforcount1))) = V.dx(nonzeros(Gator1Data.Index1(:,cadaforcount1)));
    cada1td1(logical(Gator1Data.Index2(:,cadaforcount1))) = cada1f4dx(nonzeros(Gator1Data.Index2(:,cadaforcount1)));
    V.dx = cada1td1;
    cada1tempf1 = [V.f(:,1:Gator1Data.Index3(cadaforcount1)),cada1f4];
    V.f = zeros(100,3);
    V.f(:,1:size(cada1tempf1,2)) = cada1tempf1;
    %User Line: V = [V, 1+2*x.^count];
end
cada1tf3 = V.f\d;
cada1td1 = sparse(Gator1Data.Index5,Gator1Data.Index6,V.dx,3,100);
cada1td1 = cada1tf3.'*cada1td1;
cada1td1 = cada1td1(:);
cada1td3 = full(cada1td1(Gator1Data.Index7));
cada1tf4 = V.f.';
cada1td1 = sparse(Gator1Data.Index8,Gator1Data.Index9,cada1td3,100,100);
cada1td1 = cada1tf4*cada1td1;
cada1td1 = cada1td1(:);
cada1td4 = full(cada1td1(Gator1Data.Index10));
cada1tf4 = (V.f*cada1tf3 - d).';
cada1td1 = sparse(Gator1Data.Index11,Gator1Data.Index12,V.dx,100,200);
cada1td1 = cada1tf4*cada1td1;
cada1td1 = cada1td1(:);
cada1td5 = full(cada1td1(Gator1Data.Index13));
cada1td3 = cada1td4;
cada1td3(Gator1Data.Index14) = cada1td3(Gator1Data.Index14) + cada1td5;
cada1tf4 = -(V.f.'*V.f);
cada1td1 = zeros(3,100);
cada1td1(Gator1Data.Index15) = cada1td3;
cada1td1 = cada1tf4\cada1td1;
cada1td1 = cada1td1(:);
p.dx = cada1td1(Gator1Data.Index16);
p.f = cada1tf3;
%User Line: p = V\d
p.dx_size = [3,100];
p.dx_location = Gator1Data.Index17;
end


function ADiGator_LoadData()
global ADiGator_myprogram_ADiGatorJac
ADiGator_myprogram_ADiGatorJac = load('myprogram_ADiGatorJac.mat');
return
end