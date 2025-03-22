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

function [ResNorm_intsum,Res_intsum] = minresCost_Y(X,U,p,T,data)
global ADiGator_minresCost_Y
if isempty(ADiGator_minresCost_Y); ADiGator_LoadData(); end
Gator1Data = ADiGator_minresCost_Y.minresCost_Y.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %costResidualMin_ModeMinRes_Adigator - cost computation for integrated residual minimization (alternating method: residual minimization) with Adigator
%User Line: %
%User Line: % Syntax:   [ ResNorm_intsum, Res_intsum ] = costResidualMin_ModeMinRes_Adigator( X,U,p,T,data)
%User Line: %
%User Line: % Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
%User Line: % The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
%User Line: % This code is published under the MIT License.
%User Line: % Department of Aeronautics and Department of Electrical and Electronic Engineering,
%User Line: % Imperial College London London  England, UK
%User Line: % ICLOCS (Imperial College London Optimal Control) Version 2.5
%User Line: % 1 Aug 2019
%User Line: % iclocs@imperial.ac.uk
dataNLP=data.dataNLP;
%User Line: dataNLP=data.dataNLP;
f=dataNLP.data.InternalDynamics;
%User Line: f=dataNLP.data.InternalDynamics;
dyn_data=data.dataNLP.data;
%User Line: dyn_data=data.dataNLP.data;
cadaconditional1 = strcmp(dataNLP.options.discretization,'globalLGR') | strcmp(dataNLP.options.discretization,'hpLGR');
%User Line: cadaconditional1 = strcmp(dataNLP.options.discretization,'globalLGR') | strcmp(dataNLP.options.discretization,'hpLGR');
    %User Line: % p/hp Transcription Method
    %User Line: n=dataNLP.sizes{3};
    %User Line: ng_eq=dataNLP.sizes{18};
    %User Line: t_0=T(1);
    %User Line: t_f=T(end);
    %User Line: delta_t=t_f-t_0;
    %User Line: U=[U;U(end,:)];
    %User Line: X_quad=data.sumInterpHMat*(data.InterpH*X);
    %User Line: U_quad=data.sumInterpHMat*(data.InterpH*U);
    %User Line: X_quad(data.interp_fixi,:)=X(data.interp_fixj,:);
    %User Line: U_quad(data.interp_fixi,:)=U(data.interp_fixj,:);
    %User Line: dX_quad=data.DT_seg_mat_d2*(data.D_mat*X_quad)/delta_t;
    %User Line: P_quad=repmat(p,data.M_quad,1);
    %User Line: T_quad=data.tau_quad*delta_t+t_0;
    %User Line: Fp=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
    %User Line: Res=(dX_quad-Fp).^2;
    %User Line: Res_int=delta_t*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
    %User Line: % h Transcription Method
    n.f = dataNLP.sizes{3};
    %User Line: n=dataNLP.sizes{3};
    nt.f = dataNLP.sizes{1};
    %User Line: nt=dataNLP.sizes{1};
    M.f = dataNLP.sizes{7};
    %User Line: M=dataNLP.sizes{7};
    ng_eq.f = dataNLP.sizes{15};
    %User Line: ng_eq=dataNLP.sizes{15};
    cadaconditional1 = nt.f;
    %User Line: cadaconditional1 = nt;
        t_0.dY = T.dY(1);
        t_0.f = T.f(1);
        %User Line: t_0=T(1);
        cada1f1 = length(T.f);
        t_f.dY = T.dY(2);
        t_f.f = T.f(cada1f1);
        %User Line: t_f=T(end);
        cada1td1 = zeros(2,1);
        cada1td1(2) = t_f.dY;
        cada1td1(1) = cada1td1(1) + -t_0.dY;
        delta_t.dY = cada1td1;
        delta_t.f = t_f.f - t_0.f;
        %User Line: delta_t=t_f-t_0;
        P.f = repmat(p,101,1);
        %User Line: P=repmat(p,M,1);
        X_col.dY = X.dY;         X_col.f = X.f;
        %User Line: X_col=X;
        U_col.dY = U.dY;         U_col.f = U.f;
        %User Line: U_col=U;
        cada1tempdY = delta_t.dY(Gator1Data.Index1);
        cada1tf1 = data.tau(Gator1Data.Index3);
        cada1f1dY = cada1tf1(:).*cada1tempdY(Gator1Data.Index2);
        cada1f1 = data.tau*delta_t.f;
        cada1tempdY = t_0.dY(Gator1Data.Index4);
        cada1td1 = zeros(201,1);
        cada1td1(Gator1Data.Index5) = cada1f1dY;
        cada1td1(Gator1Data.Index6) = cada1td1(Gator1Data.Index6) + cada1tempdY;
        T_col.dY = cada1td1;
        T_col.f = cada1f1 + t_0.f;
        %User Line: T_col=data.tau*delta_t+t_0;
        %User Line: % F_k=f(X_col,U_col,P,T_col,dyn_data);
        %User Line: % F_kph=data.DxHS_hf*X/delta_t-F_k(1:end-1,:)/2;
        %User Line: % F_kp1=data.DxHS_p1*X/delta_t+F_k(1:end-1,:);
        %User Line: % F=[F_k(1:end-1,:) F_kph F_kp1]';
        %User Line: % F=reshape(F(:),n,3*data.nps)';
        cada1f1dY = X_col.dY(Gator1Data.Index7);
        cada1f1 = X_col.f(:,1);
        cada1f2dY = X_col.dY(Gator1Data.Index8);
        cada1f2 = X_col.f(:,2);
        cada1f3dY = 2.*cada1f2dY;
        cada1f3 = 2*cada1f2;
        cada1temp1 = Gator1Data.Data1;
        cada1f4dY = cada1f3dY;
        cada1f4 = cada1temp1;
        cada1f4(:,1) = cada1f3;
        cada1f5dY = 10.*cada1f2dY;
        cada1f5 = 10*cada1f2;
        cada1f6dY = 1.44.*cada1f1dY;
        cada1f6 = 1.44*cada1f1;
        cada1f7dY = -cada1f6dY;
        cada1f7 = 0.21 - cada1f6;
        cada1td1 = zeros(202,1);
        cada1td1(Gator1Data.Index9) = cada1f7(:).*cada1f5dY;
        cada1td1(Gator1Data.Index10) = cada1td1(Gator1Data.Index10) + cada1f5(:).*cada1f7dY;
        cada1f8dY = cada1td1;
        cada1f8 = cada1f5.*cada1f7;
        cada1f9dY = 0.8.*cada1f1dY;
        cada1f9 = 0.8*cada1f1;
        cada1td1 = cada1f8dY;
        cada1td1(Gator1Data.Index11) = cada1td1(Gator1Data.Index11) + -cada1f9dY;
        cada1f10dY = cada1td1;
        cada1f10 = cada1f8 - cada1f9;
        cada1td1 = zeros(303,1);
        cada1td1(Gator1Data.Index12) = cada1f10dY;
        cada1td1(Gator1Data.Index13) = cada1td1(Gator1Data.Index13) + U_col.dY;
        cada1f11dY = cada1td1;
        cada1f11 = cada1f10 + U_col.f;
        cada1td1 = zeros(404,1);
        cada1td1(Gator1Data.Index14) = cada1f11dY;
        cada1td1(Gator1Data.Index15) = cada1f4dY(Gator1Data.Index16);
        F.dY = cada1td1;
        F.f = cada1f4;
        F.f(:,2) = cada1f11;
        %User Line: F=f(X_col,U_col,P,T_col,dyn_data);
        cada1f1 = size(X.f,1);
        cada1f2 = cada1f1 - 1;
        cada1f3 = 1:2:cada1f2;
        cada1f4dY = X.dY(Gator1Data.Index17);
        cada1f4 = X.f(cada1f3,:);
        cada1td1 = sparse(Gator1Data.Index18,Gator1Data.Index19,cada1f4dY,50,100);
        cada1td1 = data.repXend_mat*cada1td1;
        cada1td1 = cada1td1(:);
        cada1f5dY = full(cada1td1(Gator1Data.Index20));
        cada1f5 = data.repXend_mat*cada1f4;
        cada1tempdY = delta_t.dY(Gator1Data.Index21);
        cada1tf1 = data.AxHS(Gator1Data.Index23);
        cada1f6dY = cada1tf1(:).*cada1tempdY(Gator1Data.Index22);
        cada1f6 = delta_t.f*data.AxHS;
        cada1td2 = sparse(Gator1Data.Index24,Gator1Data.Index25,cada1f6dY,101,220);
        cada1td2 = F.f.'*cada1td2;
        cada1td1 = zeros(1760,1);
        cada1td1(Gator1Data.Index27) = cada1td2(Gator1Data.Index26);
        cada1td2 = sparse(Gator1Data.Index28,Gator1Data.Index29,F.dY,101,404);
        cada1td2 = cada1f6*cada1td2;
        cada1td2 = cada1td2(:);
        cada1td1(Gator1Data.Index31) = cada1td1(Gator1Data.Index31) + cada1td2(Gator1Data.Index30);
        cada1f7dY = cada1td1;
        cada1f7 = cada1f6*F.f;
        cada1td1 = zeros(1950,1);
        cada1td1(Gator1Data.Index32) = cada1f5dY;
        cada1td1(Gator1Data.Index33) = cada1td1(Gator1Data.Index33) + cada1f7dY;
        X_quad.dY = cada1td1;
        X_quad.f = cada1f5 + cada1f7;
        %User Line: X_quad=data.repXend_mat*X(1:2:end-1,:)+delta_t*data.AxHS*F;
        cada1td1 = sparse(Gator1Data.Index34,Gator1Data.Index35,U.dY,101,101);
        cada1td1 = data.AuHS*cada1td1;
        cada1td1 = cada1td1(:);
        U_quad.dY = full(cada1td1(Gator1Data.Index36));
        U_quad.f = data.AuHS*U.f;
        %User Line: U_quad=data.AuHS*U;
        cada1td1 = sparse(Gator1Data.Index37,Gator1Data.Index38,F.dY,101,404);
        cada1td1 = data.AfHS*cada1td1;
        cada1td1 = cada1td1(:);
        dX_quad.dY = full(cada1td1(Gator1Data.Index39));
        dX_quad.f = data.AfHS*F.f;
        %User Line: dX_quad=data.AfHS*F;
        P_quad.f = repmat(p,150,1);
        %User Line: P_quad=repmat(p,data.M_quad,1);
        cada1tempdY = delta_t.dY(Gator1Data.Index40);
        cada1tf1 = data.tau_quad(Gator1Data.Index42);
        cada1f1dY = cada1tf1(:).*cada1tempdY(Gator1Data.Index41);
        cada1f1 = data.tau_quad*delta_t.f;
        cada1tempdY = t_0.dY(Gator1Data.Index43);
        cada1td1 = zeros(299,1);
        cada1td1(Gator1Data.Index44) = cada1f1dY;
        cada1td1(Gator1Data.Index45) = cada1td1(Gator1Data.Index45) + cada1tempdY;
        T_quad.dY = cada1td1;
        T_quad.f = cada1f1 + t_0.f;
        %User Line: T_quad=data.tau_quad*delta_t+t_0;
        cada1f1dY = X_quad.dY(Gator1Data.Index46);
        cada1f1 = X_quad.f(:,1);
        cada1f2dY = X_quad.dY(Gator1Data.Index47);
        cada1f2 = X_quad.f(:,2);
        cada1f3dY = 2.*cada1f2dY;
        cada1f3 = 2*cada1f2;
        cada1temp1 = Gator1Data.Data2;
        cada1f4dY = cada1f3dY;
        cada1f4 = cada1temp1;
        cada1f4(:,1) = cada1f3;
        cada1f5dY = 10.*cada1f2dY;
        cada1f5 = 10*cada1f2;
        cada1f6dY = 1.44.*cada1f1dY;
        cada1f6 = 1.44*cada1f1;
        cada1f7dY = -cada1f6dY;
        cada1f7 = 0.21 - cada1f6;
        cada1tf1 = cada1f7(Gator1Data.Index48);
        cada1td1 = zeros(1290,1);
        cada1td1(Gator1Data.Index49) = cada1tf1(:).*cada1f5dY;
        cada1tf1 = cada1f5(Gator1Data.Index50);
        cada1td1(Gator1Data.Index51) = cada1td1(Gator1Data.Index51) + cada1tf1(:).*cada1f7dY;
        cada1f8dY = cada1td1;
        cada1f8 = cada1f5.*cada1f7;
        cada1f9dY = 0.8.*cada1f1dY;
        cada1f9 = 0.8*cada1f1;
        cada1td1 = cada1f8dY;
        cada1td1(Gator1Data.Index52) = cada1td1(Gator1Data.Index52) + -cada1f9dY;
        cada1f10dY = cada1td1;
        cada1f10 = cada1f8 - cada1f9;
        cada1td1 = zeros(1330,1);
        cada1td1(Gator1Data.Index53) = cada1f10dY;
        cada1td1(Gator1Data.Index54) = cada1td1(Gator1Data.Index54) + U_quad.dY;
        cada1f11dY = cada1td1;
        cada1f11 = cada1f10 + U_quad.f;
        cada1td1 = zeros(2580,1);
        cada1td1(Gator1Data.Index55) = cada1f11dY;
        cada1td1(Gator1Data.Index56) = cada1f4dY(Gator1Data.Index57);
        Fp.dY = cada1td1;
        Fp.f = cada1f4;
        Fp.f(:,2) = cada1f11;
        %User Line: Fp=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
        cada1td1 = zeros(2580,1);
        cada1td1(Gator1Data.Index58) = dX_quad.dY;
        cada1td1 = cada1td1 + -Fp.dY;
        cada1f1dY = cada1td1;
        cada1f1 = dX_quad.f - Fp.f;
        cada1tf2 = cada1f1(Gator1Data.Index59);
        Res.dY = 2.*cada1tf2(:).^(2-1).*cada1f1dY;
        Res.f = cada1f1.^2;
        %User Line: Res=(dX_quad-Fp).^2;
        cada1tempdY = delta_t.dY(Gator1Data.Index60);
        cada1tf1 = data.DT_seg_node_mat(Gator1Data.Index62);
        cada1f1dY = cada1tf1(:).*cada1tempdY(Gator1Data.Index61);
        cada1f1 = delta_t.f*data.DT_seg_node_mat;
        cada1f2dY = cada1f1dY./2;
        cada1f2 = cada1f1/2;
        cada1td1 = sparse(Gator1Data.Index63,Gator1Data.Index64,cada1f2dY,50,100);
        cada1td1 = data.sum_nps_quad.'*cada1td1;
        cada1td1 = cada1td1(:);
        cada1f3dY = full(cada1td1(Gator1Data.Index65));
        cada1f3 = cada1f2*data.sum_nps_quad;
        cada1td2 = sparse(Gator1Data.Index66,Gator1Data.Index67,cada1f3dY,150,100);
        cada1td2 = Res.f.'*cada1td2;
        cada1td1 = zeros(1100,1);
        cada1td1(Gator1Data.Index69) = cada1td2(Gator1Data.Index68);
        cada1td2 = sparse(Gator1Data.Index70,Gator1Data.Index71,Res.dY,150,610);
        cada1td2 = cada1f3*cada1td2;
        cada1td2 = cada1td2(:);
        cada1td1 = cada1td1 + full(cada1td2(Gator1Data.Index72));
        Res_int.dY = cada1td1;
        Res_int.f = cada1f3*Res.f;
        %User Line: Res_int=delta_t*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
        %User Line: t0=dataNLP.t0;
        %User Line: tf=dataNLP.tf;
        %User Line: deltat=tf-t0;
        %User Line: P=repmat(p,M,1);
        %User Line: X_col=X;
        %User Line: U_col=U;
        %User Line: T_col=data.tau*deltat+t0;
        %User Line: % F_k=f(X_col,U_col,P,T_col,dyn_data);
        %User Line: % F_kph=data.DxHS_hf*X/deltat-F_k(1:end-1,:)/2;
        %User Line: % F_kp1=data.DxHS_p1*X/deltat+F_k(1:end-1,:);
        %User Line: % F=[F_k(1:end-1,:) F_kph F_kp1]';
        %User Line: % F=reshape(F(:),n,3*data.nps)';
        %User Line: F=f(X_col,U_col,P,T_col,dyn_data);
        %User Line: X_quad=data.repXend_mat*X(1:2:end-1,:)+deltat*data.AxHS*F;
        %User Line: U_quad=data.AuHS*U;
        %User Line: dX_quad=data.AfHS*F;
        %User Line: P_quad=repmat(p,data.M_quad,1);
        %User Line: T_quad=data.tau_quad*deltat+t0;
        %User Line: Fp=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
        %User Line: Res=(dX_quad-Fp).^2;
        %User Line: Res_int=deltat*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
%User Line: % Compuation of integrated residual for each dynamics equation
cada1f1dY = Res_int.dY;
cada1f1 = Res_int.f(:);
cada1td1 = sparse(Gator1Data.Index73,Gator1Data.Index74,cada1f1dY,100,305);
cada1td1 = data.ResConstScaleMat*cada1td1;
cada1td1 = cada1td1(:);
Res_int_Const.dY = full(cada1td1(Gator1Data.Index75));
Res_int_Const.f = data.ResConstScaleMat*cada1f1;
%User Line: Res_int_Const=data.ResConstScaleMat*Res_int(:);
cada1f1 = n.f + ng_eq.f;
Res_int_Const.dY = Res_int_Const.dY;
Res_int_Const.f = reshape(Res_int_Const.f,data.nps,cada1f1);
%User Line: Res_int_Const=reshape(Res_int_Const,data.nps,n+ng_eq);
cada1td1 = sum(sparse(Gator1Data.Index76,Gator1Data.Index77,Res_int_Const.dY,50,610),1);
cada1f1dY = full(cada1td1(:));
cada1f1 = sum(Res_int_Const.f);
Res_intsum.dY = cada1f1dY;
Res_intsum.f = cada1f1.';
%User Line: Res_intsum=sum(Res_int_Const)';
%User Line: % Compuation of integrated residual norm
cada1f1dY = Res_int.dY;
cada1f1 = Res_int.f(:);
cada1td1 = sparse(Gator1Data.Index78,Gator1Data.Index79,cada1f1dY,100,305);
cada1td1 = data.ResNormScaleMat*cada1td1;
cada1td1 = cada1td1(:);
Res_int_Norm.dY = full(cada1td1(Gator1Data.Index80));
Res_int_Norm.f = data.ResNormScaleMat*cada1f1;
%User Line: Res_int_Norm=data.ResNormScaleMat*Res_int(:);
cada1f1 = n.f + ng_eq.f;
Res_int_Norm.dY = Res_int_Norm.dY;
Res_int_Norm.f = reshape(Res_int_Norm.f,data.nps,cada1f1);
%User Line: Res_int_Norm=reshape(Res_int_Norm,data.nps,n+ng_eq);
cada1td1 = sum(sparse(Gator1Data.Index81,Gator1Data.Index82,Res_int_Norm.dY,50,610),1);
cada1f1dY = full(cada1td1(:));
cada1f1 = sum(Res_int_Norm.f);
cada1td1 = zeros(2,305);
cada1td1(Gator1Data.Index83) = cada1f1dY;
cada1td1 = sum(cada1td1,1);
ResNorm_intsum.dY = cada1td1(:);
ResNorm_intsum.f = sum(cada1f1);
%User Line: ResNorm_intsum=sum(sum(Res_int_Norm));
ResNorm_intsum.dY_size = 305;
ResNorm_intsum.dY_location = Gator1Data.Index84;
Res_intsum.dY_size = [2,305];
Res_intsum.dY_location = Gator1Data.Index85;
end


function ADiGator_LoadData()
global ADiGator_minresCost_Y
ADiGator_minresCost_Y = load('minresCost_Y.mat');
return
end