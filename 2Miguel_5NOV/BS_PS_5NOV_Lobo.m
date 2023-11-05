clear 
close all
clc


%% DEFINIÇÃO DOS PARAMETROS COMUNS DO PROCESSAMENTO DOS DADOS DE TREINO E DE TESTE

%Introdução dos parametros comuns (treino e teste) a processar
disp('Introdução dos dados comuns de configuração para extração dos parâmetros das ARCs do modelo e das observações');
disp(' ');
%         angle_min=input('Ângulo mínimo: ');
        angle_min=25;
%         angle_max=input('Ângulo máximo: ');
        angle_max=55;
%         angRes=input('Resolução angular: ');
        angRes=1;
%         IntAng=input('Ângulo de "Intercept": ');
        IntAng=40;
               
        disp(' ');
        
        
        disp(['CONFIRMAÇÃO PARÂMETROS DE PROCESSAMENTO']);
        disp(['Angulo mínimo: ',num2str(angle_min),'°']);
        disp(['Angulo máximo: ',num2str(angle_max),'°']);
        disp(['Resolução angular: ',num2str(angRes),'°']);
        disp(['Angulo de intercepto: ',num2str(IntAng),'°']);
        disp(' ');
        
W=input('"Enter" para continuar');

%% PREPARAÇÃO E PROCESSAMENTO DOS DADOS BACKSCATTER A TESTAR (CLASSIFICAR)

file=input('Insira  o caminho e o nome (*.txt) do ficheiro a classificar, FMGT Processed Backscatter: ','s');
disp(' ');

fid10=fopen(file,'r');
fid20=fopen([file(1:end-4) '_ARC_VECTORS.txt'],'w'); % Ficheiro com os vetores extraídos das BS ARCs

% Se o ficheiro for criado pela versão recente do FMGT inserir mais uma coluna %f e 18 na
% Headerlines 
% file_version=input('Versao antiga(0) ou nova(1)?')
file_version=0
if file_version==0
    data_cell=textscan(fid10,'%f %f %f %f %f %f %f %f %f %f %f %f','Delimiter','\t','Headerlines',17);
else
    data_cell=textscan(fid10,'%f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter','\t','Headerlines',18);
end

% ping_date_time=[data_cell{1}];
ping_n=[data_cell{2}];
beam_n=[data_cell{3}];
back_corr=[data_cell{10}];
true_angle=[data_cell{11}];

disp(['Ficheiro:' file]);
% disp(['Date and time: ' num2str(ping_date_time(1)) ' to ' num2str(ping_date_time(end))]);
disp(['Pings: ' num2str(min(ping_n)) ' to ' num2str(max(ping_n))]);
disp(['Beams: ' num2str(min(beam_n)) ' to ' num2str(max(beam_n))]);
disp(' ');

disp('DEFINIÇÃO DIMENSÃO DO PATCH');
% patch=input('Número de pings consecutivos: ');
patch=30;
disp(' ');

%Escrita do cabeçalho no ficheiro (*_ARC_VECTORS.txt)
fprintf(fid20, '%s\n', ['Ficheiro FMGT   : ','SSMF line: ' file]);
fprintf(fid20, '%s\n', ['Pings           : ' num2str(min(ping_n)) ' to ' num2str(max(ping_n))]);
fprintf(fid20, '%s\n', ['Beams           : ' num2str(min(beam_n)) ' to ' num2str(max(beam_n))]);
fprintf(fid20, '%s\n', ['DimePatch       : ' num2str(patch) ' pings']);
fprintf(fid20, '%s\n', ['IntervAngularARC: ' num2str(angle_min),'°' ' to ' num2str(angle_max),'°']);
fprintf(fid20, '%s\n', ['ResolAngular    : ' num2str(angRes),'°']);
fprintf(fid20, '%s\n', ['AngIntercept   : ' num2str(IntAng),'°']);
fprintf(fid20, '%s\n', ('#----------------------------------------------#'));
fprintf(fid20, '%s\n', ('BORDO,NUM_PATCH,ANG_MEDIO,BS_MEDIO,INTERC,SLOPE'));



max_k=patch*(floor(max(ping_n)/patch));
start_ping=data_cell{2}(1,1);
for k=start_ping:patch:max_k
    
    ping_number_i=k;
    if k==max_k
        ping_number_f=max(ping_n);
    else
        ping_number_f=k+patch;    
    end
    
    % Numeração dos patches
    numPatch=k/patch;
    
    indices=find(ping_n>=ping_number_i & ping_n<=ping_number_f & abs(true_angle)>=angle_min & abs(true_angle)<=angle_max & back_corr<=10 & back_corr>=-70);

    back_corr_indices=back_corr(indices);
    true_angle_indices=true_angle(indices);

    back_corr_median_ps=zeros(angle_max-angle_min+1,2);

    true_angle_mean_ps=zeros(angle_max-angle_min+1,1);
%     true_angle_mean_sb=zeros(angle_max-angle_min+1,1);

    count=1;
    for i=angle_min:angRes:angle_max
%      back_corr_mean_ps(count,1)=mean(back_corr_indices(find(true_angle_indices>=-i-angRes/2 & true_angle_indices<=-i+angRes/2)));% média direta dos dBs (FMGT)
%      back_corr_mean_ps(count,1)=20*log10(mean(10.^(((back_corr_indices(find(true_angle_indices>=-i-angRes/2 & true_angle_indices<=-i+angRes/2)))/20))));% média valores lineares dos dBs
       back_corr_mean_ps(count,1)=20*log10(median(10.^((back_corr_indices(find(true_angle_indices>=-i-angRes/2 & true_angle_indices<=-i+angRes/2))/20))))+1; % mediana dos valores lineares +1dB offset   
       back_corr_mean_ps(count,2)=length(find(true_angle_indices>=-i-angRes/2 & true_angle_indices<=-i+angRes/2));
       true_angle_mean_ps(count)=mean(true_angle_indices(find(true_angle_indices>=-i-angRes/2 & true_angle_indices<=-i+angRes/2)));
       count=count+1;
    end
    
%     count=1;
%     for i=angle_min:angRes:angle_max
%        back_corr_mean_sb(count,1)=20*log10(mean(10.^(((back_corr_indices(find(true_angle_indices>=i-angRes/2 & true_angle_indices<=i+angRes/2)))/20))));
%        back_corr_mean_sb(count,2)=length(find(true_angle_indices>=i-angRes/2 & true_angle_indices<=i+angRes/2));
%        true_angle_mean_sb(count)=mean(true_angle_indices(find(true_angle_indices>=i-angRes/2 & true_angle_indices<=i+angRes/2)));
%        count=count+1;
%     end
    
    %Limpezas dos dados sem informação
    back_corr_mean_ps=back_corr_mean_ps(isnan(back_corr_mean_ps(:,1))==0);
    true_angle_mean_ps=true_angle_mean_ps(isnan(true_angle_mean_ps(:,1))==0);
%     back_corr_mean_sb=back_corr_mean_sb(isnan(back_corr_mean_sb(:,1))==0);
%     true_angle_mean_sb=true_angle_mean_sb(isnan(true_angle_mean_sb(:,1))==0);

    A2=fitlm(true_angle_mean_ps,back_corr_mean_ps(:,1),'Linear');
%     B=fitlm(true_angle_mean_sb,back_corr_mean_sb(:,1),'Linear');
    
    A2_Coef=table2array(A2.Coefficients);
%     B_Coef=table2array(B.Coefficients);
    
    A2_AngInt=((A2_Coef(2,1)*IntAng*(-1))+(A2_Coef(1,1)));
%     B_AngInt=((B_Coef(2,1)*IntAng*(-1)+(B_Coef(1,1))));
    
    disp(' ')
    
    disp('PORTSIDE');
    disp(['NÚMERO DO PATCH:    ' num2str(numPatch)]);
    disp(['ANGULO MÉDIO:       ' num2str(mean(true_angle_mean_ps(:,1)))]);
    disp(['BS MÉDIO:           ' num2str(mean(back_corr_mean_ps(:,1)))]);
    disp(['INTERCEPT A -' num2str(IntAng),'°:   ' num2str(A2_AngInt)]);
    disp(['SLOPE:              ' num2str(A2_Coef(2,1))]);

    disp(' ')

    %Escrita dos vetores no ficheiro *_ARC_VECTORS.txt
    fprintf(fid20,'%s\n',['PS',',',num2str(numPatch),',',num2str(mean(true_angle_mean_ps(:,1))),',',num2str(mean(back_corr_mean_ps(:,1))),',',num2str(A2_AngInt),',',num2str(A2_Coef(2,1))]);
%   fprintf(fid20,'%s\n',['PS',',',num2str(numPatch),',',num2str(mean(true_angle_mean_ps(:,1)),',',num2str(A_AngInt),',',num2str(A_Coef(2,1)))]);
 

%%    VISUALIZAÇÃO DE DADOS  

%         figure;
%         plot(true_angle_mean_ps,back_corr_mean_ps(:,1));
%         title(['NumPatch: ' num2str(numPatch)]);
%         xlabel('Incident angle interval average (º)');
%         ylabel('Backscatter (db)')
%         ylim([-35, -10]);
%         hold on; 
% 
%         figure;
%     
%         plot (A2);
%         title(['NumPatch: ' num2str(numPatch),'    Avrg: ' num2str(mean(back_corr_mean_ps(:,1))),'    Intercept at -' num2str(IntAng),':','  ' num2str(A2_AngInt),'     Slope: ' num2str(A2_Coef(2,1))]);
%         xlabel('Incident angle interval average (º)');
%         ylabel('Backscatter (db)')

%    2D dos ARCs referentes aos patches em análise
%     figure;
%     plot(true_angle_mean_ps,back_corr_mean_ps(:,1));
%     title(['patch: ' num2str(k/patch)]);
%     hold on;
%     
%    3D dos ARCs referentes aos patches em análise  
%     y=zeros(length(true_angle_mean_ps));
%     y(:)=k/patch;
%     plot3(true_angle_mean_ps, y, back_corr_mean_ps(:,1));
%     hold on
% 
%    plot(A);
%     title(['AngResol: ' num2str(angRes),'    Avrg: ' num2str(mean(back_corr_mean_sb(:,1))),'    Intercept at -' num2str(IntAng),': '   num2str(A_AngInt),'     Slope: ' num2str(A_Coef(2,1))]);
%     xlabel('True angle interval average: (°)');
%     ylabel('Backscatter(db)')
% 
%    figure;
%     plot(true_angle_mean_sb,back_corr_mean_sb(:,1));
%     hold on;

%    plot(B);
%     title(['AngResol: ' num2str(angRes),'    Avrg: ' num2str(mean(back_corr_mean_ps(:,1))),'    Intercept at ' num2str(IntAng),': '  num2str(B_AngInt),'     Slope: ' num2str(B_Coef(2,1))]);
%     xlabel('True angle interval average: (°)');
%     ylabel('Backscatter(db)')
    
end

%% CLASSIFICADOR ML
% 
% clear
% clc
% close all

% Importar para o Workspace os vectores BS a testar
Proc_BS_Testar = readtable("E:\BkSctt_Research\1_MatLab_MLClassif\1_ML_BkSctt Tests\#37_EstimaAM_ARC_VECTORS.txt");

% Chamar a Function_trainClassifier.m
modeloBS=Function_SVMOpt_trainClassifier;

% Testar os dados com SVM Optimizada
PS_BS_SedClass = modeloBS.predictFcn(Proc_BS_Testar);

Proc_BS_Testar.CLASSE=(PS_BS_SedClass); % Adiciona a coluna CLASSE ao ficheiro a testar

ResulClassif=Proc_BS_Testar;

disp(' ');
disp(ResulClassif);

disp('                           "RESULTADO DA CLASSIFICACAO"')

%%
% % GEOREFERENCIAÇAO DA INFORMAÇAO

% data_cell{14}=cell(length(data_cell{1}),1);

%   for i=1:length(indices(1,:))
%       data_cell{14}(cell2mat(indices(1,i)))=Proc_BS_Testar(i,7).Class;
%   end
