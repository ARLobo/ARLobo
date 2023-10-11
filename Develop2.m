% clear
% clc
% 
% % PREPARAÇÃO E PROCESSAMENTO DOS DADOS APL PARA TREINAR O CLASSIFICATION LEARNER
%%
%Introdução dos parametros comuns (treino e teste) a processar
disp('Introdução dos dados configuração de extração parametros das ARCs');
disp(' ');
%         angle_min=input('Angulo mínimo: ')
        angle_min=25;
%         angle_max=input('Angulo máximo: ')
        angle_max=55;
%         angRes=input('Resolução angular: ')
        angRes=0.5;
%         IntAng=input('Angulo de Intercept: ')
        In25tAng=40;

%%
% Define a pasta onde os ficheiros APL iniciais estão localizados
pasta_ps = 'E:\BkSctt_Research\1_MatLab_MLClassif\2_DadosAPL\ASCII_Files\APL_6col\PS';
pasta_sb = 'E:\BkSctt_Research\1_MatLab_MLClassif\2_DadosAPL\ASCII_Files\APL_6col\SB';
disp(['Pasta com ficheiros de origem : ' pasta_ps]);
disp(['Pasta com ficheiros de origem : ' pasta_sb]);

nomesFicheiros_ps = dir(fullfile(pasta_ps, 'PS_*.txt'));
nomesFicheiros_sb = dir(fullfile(pasta_sb, 'SB_*.txt'));
    
% Nome do ficheiro de destino
    nomeficheiroDestino_ps = fullfile('E:\BkSctt_Research\1_MatLab_MLClassif\2_DadosAPL\ASCII_Files\DadosDeTreino', 'PS_APLData_2TrainClassf.txt');
    nomeficheiroDestino_sb = fullfile('E:\BkSctt_Research\1_MatLab_MLClassif\2_DadosAPL\ASCII_Files\DadosDeTreino', 'SB_APLData_2TrainClassf.txt');
    
    % Abrir os ficheiros de destino para escrita
    fid_ps = fopen(nomeficheiroDestino_ps, 'w'); 
    fid_sb = fopen(nomeficheiroDestino_sb, 'w'); 
    
    if fid_ps == -1
                error(['Não foi possível criar o arquivo de destino ' nomeficheiroDestino_ps]);
    end
                
    if fid_sb == -1
                error(['Não foi possível criar o arquivo de destino ' nomeficheiroDestino_sb]);
    end  
    
    
%Introdução dos parametros comuns (treino e teste) a processar
%         angle_min=input('Angulo mínimo: ')
        angle_min=25;
%         angle_max=input('Angulo máximo: ')
        angle_max=55;
%         angRes=input('Resolução angular: ')
        angRes=0.5;
%         IntAng=input('Angulo de Intercept: ')
        IntAng=40;
              
        %Escrita do cabeçalho no ficheiro dados PS_*.txt. Esta inf é importante que acompanhe o ficheiro.
        fprintf(fid_ps, '%s\n', ('Ficheiro Modelo : APL PS'));
        fprintf(fid_ps, '%s\n', ['Intervalo ARC   : ' num2str(angle_min),'°' ' to ' num2str(angle_max),'°']);
        fprintf(fid_ps, '%s\n', ['Resol angular   : ' num2str(angRes),'°']);
        fprintf(fid_ps, '%s\n', ['Angulo Intercept: ' num2str(IntAng),'°']);
        fprintf(fid_ps, '%s\n', ('#----------------------------------------------'));
        fprintf(fid_ps, 'BORDO,BS_MEDIO,INTERC1,INTERC2,SLOPE,ANGMEDIO,CLASSE');
        fprintf(fid_ps, '%s\n', (' ') );
        
        %Escrita do cabeçalho no ficheiro dados PS_*.txt. Esta inf é importante que acompanhe o ficheiro.
        fprintf(fid_sb, '%s\n', ('Ficheiro Modelo : APL SB'));
        fprintf(fid_sb, '%s\n', ['Intervalo ARC   : ' num2str(angle_min),'°' ' to ' num2str(angle_max),'°']);
        fprintf(fid_sb, '%s\n', ['Resol angular   : ' num2str(angRes),'°']);
        fprintf(fid_sb, '%s\n', ['Angulo Intercept: ' num2str(IntAng),'°']);
        fprintf(fid_sb, '%s\n', ('#----------------------------------------------'));
        fprintf(fid_sb, 'BORDO,BS_MEDIO,INTERC1,INTERC2,SLOPE,ANGMEDIO,CLASSE');
        fprintf(fid_sb, '%s\n', (' ') );
  
% Loop para processar os ficheiros PS (numl>numero linhas)
for i = 1:numel(nomesFicheiros_ps)
    % Nome do arquivo de origem
    file = fullfile(pasta_ps, nomesFicheiros_ps(i).name);
      
    % Abrir o ficheiro de origem para leitura
    fid1 = fopen(file, 'r');
    
    if fid1 == -1
        error(['Não foi possível abrir o arquivo de origem ' file]);
    end
    
     
     % Temos os ficheiros de leitura abertos, prontos para processamento. 
       
    data_cell=textscan(fid1,'%s %f %f %f %f %s','Delimiter',' ');

        sedname=[data_cell{1}{1}]; % Converte a celula em string para poder fazer fprintf na linha 111
        phi=[data_cell{2}];
        % gama=[data_cell{3}]
        inc_ang=[data_cell{4}];
        bsdb=[data_cell{5}];
        sedclass=[data_cell{6}{1}]; % Converte a celula em string para poder fazer fprintf na linha 111

        disp(' ');
        disp(['Ficheiro: ' sedname]);
        % disp('Tipo sedimento: ' sedname(1:1));
        % disp(['Phi: ' num2str(phi(1:1))]);
        % disp(['Angulo de incidência: ' num2str(inc_ang)]);
        % disp(['BS: ' num2str(bsdb)]);
        % disp(['Classe sedimento: ' sedclass(1:1)]);
        % disp('Inserir o angulo inicial e o final: ')
        
         
        % Preparação das matrizes para guardar resultados
        indices=find(abs(inc_ang)>=angle_min & abs(inc_ang)<=angle_max); % Cria indices referentes aos ângulos entre (angle_min e angle=max)

        bsdb_indices=bsdb(indices); % Dados BS limitados aos (angle_min e angle=max)
        inc_ang_indices=inc_ang(indices); % Dados IncAng limitados aos (angle_min e angle=max)     

        bsdb_mean_ps=zeros(angle_max-angle_min+1,2); % calcula o número de elementos necessários na matriz para cobrir todos
                                                        % os ângulos entre angle_min e angle_max, incluindo esses valores extremos                                             

        inc_ang_mean_ps=zeros(angle_max-angle_min+1,1); % calcula o número de elementos necessários na matriz para cobrir todos
        % os ângulos entre angle_min e angle_max, incluindo esses valores extremos.
        % Essa matriz é usada para armazenar médias ou valores relacionados a um conjunto de ângulos, onde cada elemento da matriz
        % corresponde a um ângulo específico dentro do intervalo angle_min e angle_max.
       
        
        count=1;
        for i=angle_min:angRes:angle_max
           bsdb_mean_ps(count,1)=mean(bsdb_indices(inc_ang_indices>=-i-angRes/2 & inc_ang_indices<=-i+angRes/2));
           bsdb_mean_ps(count,2)=length(find(inc_ang_indices>=-i-angRes/2 & inc_ang_indices<=-i+angRes/2));
           inc_ang_mean_ps(count)=mean(inc_ang_indices(inc_ang_indices>=-i-angRes/2 & inc_ang_indices<=-i+angRes/2));
           count=count+1;
        end


        %Limpezas dos dados sem informação
        bsdb_mean_ps=bsdb_mean_ps(isnan(bsdb_mean_ps(:,1))==0);
        inc_ang_mean_ps=inc_ang_mean_ps(isnan(inc_ang_mean_ps(:,1))==0);
        

        A=fitlm(inc_ang_mean_ps,bsdb_mean_ps(:,1),'Linear');

        A_Coef=table2array(A.Coefficients);

        A_AngInt=((A_Coef(2,1)*IntAng+(A_Coef(1,1)))); % B_Coef(1,1) Intercept a zero graus
                                                            % B_AngInt Intercept ao angulo inserido
                                                            
                                                                                                                      
       % Escrita dos resultados em \ASCII_Files\DadosDeTreino\PS_APLData_2TrainClassf.txt 
          fprintf(fid_ps,'%s\n',['PS',',',num2str(mean(bsdb_mean_ps(:,1))),',',num2str(A_Coef(1,1)),',',num2str(A_AngInt),',',num2str(A_Coef(2,1)),',',num2str(mean(inc_ang_mean_ps(:,1))),',',sedclass]);

        % figure;
        % plot(inc_ang_mean_ps,bsdb_mean_ps(:,1));
        % hold on; 

        % plot (A);
        % title(['ResAng: ' num2str(e),'    Avrg: ' num2str(mean(bsdb_mean_ps(:,1))),'    Intercept at -' num2str(IntAng),':','  ' num2str(A_AngInt),'     Slope: ' num2str(A_Coef(2,1))]);
        % xlabel('Incident angle interval average (º)');
        % ylabel('Backscatter(db)')


        disp('PortSide')
        disp(['Mean: ' num2str(mean(bsdb_mean_ps(:,1)))])
        % disp(num2str(length(bsdb_mean_ps(:,1))))
        disp(['Intercept: ' num2str(A_Coef(1,1))]);
        disp(['Intercept at -' num2str(IntAng),':','  ' num2str(A_AngInt)])
        disp(['Slope: ' num2str(A_Coef(2,1))]);
        % disp(num2str(length(bsdb_mean_ps(:,1))))
        disp(' ');

end    

% Loop para processar os ficheiros SB (numl>numero linhas)
for i = 1:numel(nomesFicheiros_sb)
    % Nome do arquivo de origem
    file = fullfile(pasta_sb, nomesFicheiros_sb(i).name);
      
    % Abrir o ficheiro de origem para leitura
    fid2 = fopen(file, 'r');
    
    if fid2 == -1
        error(['Não foi possível abrir o arquivo de origem ' file]);
    end
    
     
     % Temos os ficheiros de leitura abertos, prontos para processamento. 
       
    data_cell=textscan(fid2,'%s %f %f %f %f %s','Delimiter',' ');

        sedname=[data_cell{1}{1}]; % Converte a celula em string para poder fazer fprintf na linha 111
        phi=[data_cell{2}];
        % gama=[data_cell{3}]
        inc_ang=[data_cell{4}];
        bsdb=[data_cell{5}];
        sedclass=[data_cell{6}{1}]; % Converte a celula em string para poder fazer fprintf na linha 111

        disp(' ');
        disp(['Ficheiro: ' sedname]);
        % disp('Tipo sedimento: ' sedname(1:1));
        % disp(['Phi: ' num2str(phi(1:1))]);
        % disp(['Angulo de incidência: ' num2str(inc_ang)]);
        % disp(['BS: ' num2str(bsdb)]);
        % disp(['Classe sedimento: ' sedclass(1:1)]);
        % disp('Inserir o angulo inicial e o final: ')
        
         
        % Preparação das matrizes para guardar resultados
        indices=find(abs(inc_ang)>=angle_min & abs(inc_ang)<=angle_max); % Cria indices referentes aos ângulos entre (angle_min e angle=max)

        bsdb_indices=bsdb(indices); % Dados BS limitados aos (angle_min e angle=max)
        inc_ang_indices=inc_ang(indices); % Dados IncAng limitados aos (angle_min e angle=max)     

        bsdb_mean_sb=zeros(angle_max-angle_min+1,2); % calcula o número de elementos necessários na matriz para cobrir todos
                                                        % os ângulos entre angle_min e angle_max, incluindo esses valores extremos                                             

        inc_ang_mean_sb=zeros(angle_max-angle_min+1,1); % calcula o número de elementos necessários na matriz para cobrir todos
        % os ângulos entre angle_min e angle_max, incluindo esses valores extremos.
        % Essa matriz é usada para armazenar médias ou valores relacionados a um conjunto de ângulos, onde cada elemento da matriz
        % corresponde a um ângulo específico dentro do intervalo angle_min e angle_max.
      

        count=1;
        for i=angle_min:angRes:angle_max
           bsdb_mean_sb(count,1)=mean(bsdb_indices(inc_ang_indices>=i-angRes/2 & inc_ang_indices<=i+angRes/2));
           bsdb_mean_sb(count,2)=length(find(inc_ang_indices>=i-angRes/2 & inc_ang_indices<=i+angRes/2));
           inc_ang_mean_sb(count)=mean(inc_ang_indices(inc_ang_indices>=i-angRes/2 & inc_ang_indices<=i+angRes/2));
           count=count+1;
        end


        %Limpezas dos dados sem informação
%         bsdb_mean_ps=bsdb_mean_ps(isnan(bsdb_mean_ps(:,1))==0);
%         inc_ang_mean_ps=inc_ang_mean_ps(isnan(inc_ang_mean_ps(:,1))==0);
        bsdb_mean_sb=bsdb_mean_sb(isnan(bsdb_mean_sb(:,1))==0);
        inc_ang_mean_sb=inc_ang_mean_sb(isnan(inc_ang_mean_sb(:,1))==0);
        


      B=fitlm(inc_ang_mean_sb,bsdb_mean_sb(:,1),'Linear');


      B_Coef=table2array(B.Coefficients);


      B_AngInt=((B_Coef(2,1)*(IntAng*-1))+(B_Coef(1,1))); % B_Coef(1,1) Intercept a zero graus
                                                            % B_AngInt Intercept ao angulo inserido
                                                            
                                                                                                                      
       % Escrita dos resultados em \ASCII_Files\DadosDeTreino\SB_APLData_2TrainClassf.txt 
        fprintf(fid_sb,'%s\n',['SB',',',num2str(mean(bsdb_mean_sb(:,1))),',',num2str(B_Coef(1,1)),',',num2str(B_AngInt),',',num2str(B_Coef(2,1)),',',num2str(mean(inc_ang_mean_sb(:,1))),',',sedclass]);

       
%         figure;
%         plot(inc_ang_mean_sb,bsdb_mean_sb(:,1));
%         hold on; 
%         
%         plot(B);
%         title (['ResAng: ' num2str(e),'    Avrg: ' num2str(mean(bsdb_mean_sb(:,1))),'    Intercept at ' num2str(IntAng),':','  ' num2str(B_AngInt),'     Slope: ' num2str(B_Coef(2,1))]);
%         xlabel('Incident angle interval average (º)');
%         ylabel('Backscatter(db)')

        disp('Starbord')
        disp(['Mean: ' num2str(mean(bsdb_mean_sb(:,1))),'dB'])
        disp(num2str(length(bsdb_mean_sb(:,1))))
        disp(['Intercept at 0°: ' num2str(B_Coef(1,1))]);
        disp(['Intercept at ' num2str(IntAng),'°',':',' ' num2str(B_AngInt)])
        disp(['Slope: ' num2str(B_Coef(2,1))]);
        disp(' ')
end    

%     close all  %fecha todas as figuras
%     fclose(fid1);
%     fclose(fid2);
    
    disp(' ');
    disp(['Os ficheiros foram processados corretamente e gravados na pasta ' nomeficheiroDestino_ps]);
    disp(' ');


    
%
% PREPARAÇÃO E PROCESSAMENTO DOS DADOS BACKSCATTER A TESTAR (CLASSIFICAR)

file=input('Insira o caminho para o ficheiro FMGT Processed Backscatter: ','s');
disp(' ');

fid10=fopen(file,'r');
fid20=fopen([file(1:end-4) '_ARC_Parameters.txt'],'w');

% Se o ficheiro for criado pela versão recente do FMGT inserir mais uma coluna %f e 18 na
% Headerlines 
data_cell=textscan(fid10,'%f %f %f %f %f %f %f %f %f %f %f %f','Delimiter',',','Headerlines',17);

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

disp('Inserir os parâmetros para definir dimensão do patch: ');

patch=input('Nº pings consecutivos: ');
disp(' ');

%Escrita do cabeçalho no ficheiro *.resultados
fprintf(fid20, '%s\n', ['Ficheiro FMGT   : ','SSMF: ' file]);
fprintf(fid20, '%s\n', ['Pings           : ' num2str(min(ping_n)) ' to ' num2str(max(ping_n))]);
fprintf(fid20, '%s\n', ['Beams           : ' num2str(min(beam_n)) ' to ' num2str(max(beam_n))]);
fprintf(fid20, '%s\n', ['DimePatch       : ' num2str(patch) ' pings']);
fprintf(fid20, '%s\n', ['IntervAngularARC: ' num2str(angle_min),'°' ' to ' num2str(angle_max),'°']);
fprintf(fid20, '%s\n', ['ResolAngular    : ' num2str(angRes),'°']);
fprintf(fid20, '%s\n', ['AngIntercept1   : ', '0°']);
fprintf(fid20, '%s\n', ['AngIntercept2   : ' num2str(IntAng),'°']);
fprintf(fid20, '%s\n', ['#----------------------------------------------#']);
fprintf(fid20, '%s\n', ['BORDO,BS_MEDIO,INTERC1,INTERC2,SLOPE,ANGMEDIO']);

max_k=patch*(floor(max(ping_n)/patch));

for k=0:patch:max_k
    
    ping_number_i=k;
    if k==max_k
        ping_number_f=max(ping_n);
    else
        ping_number_f=k+patch;    
    end
    
    
    indices=find(ping_n>=ping_number_i & ping_n<=ping_number_f & abs(true_angle)>=angle_min & abs(true_angle)<=angle_max & back_corr<=20 & back_corr>=-100);

    back_corr_indices=back_corr(indices);
    true_angle_indices=true_angle(indices);

    back_corr_mean_ps=zeros(angle_max-angle_min+1,2);

    true_angle_mean_ps=zeros(angle_max-angle_min+1,1);
    true_angle_mean_sb=zeros(angle_max-angle_min+1,1);

    count=1;
    for i=angle_min:angRes:angle_max
       back_corr_mean_ps(count,1)=mean(back_corr_indices(find(true_angle_indices>=-i-angRes/2 & true_angle_indices<=-i+angRes/2)));
       back_corr_mean_ps(count,2)=length(find(true_angle_indices>=-i-angRes/2 & true_angle_indices<=-i+angRes/2));
       true_angle_mean_ps(count)=mean(true_angle_indices(find(true_angle_indices>=-i-angRes/2 & true_angle_indices<=-i+angRes/2)));
       count=count+1;
    end
    
    count=1;
    for i=angle_min:angRes:angle_max
       back_corr_mean_sb(count,1)=mean(back_corr_indices(find(true_angle_indices>=i-angRes/2 & true_angle_indices<=i+angRes/2)));
       back_corr_mean_sb(count,2)=length(find(true_angle_indices>=i-angRes/2 & true_angle_indices<=i+angRes/2));
       true_angle_mean_sb(count)=mean(true_angle_indices(find(true_angle_indices>=i-angRes/2 & true_angle_indices<=i+angRes/2)));
       count=count+1;
    end
    
    %Limpezas dos dados sem informação
    back_corr_mean_ps=back_corr_mean_ps(isnan(back_corr_mean_ps(:,1))==0);
    true_angle_mean_ps=true_angle_mean_ps(isnan(true_angle_mean_ps(:,1))==0);
    back_corr_mean_sb=back_corr_mean_sb(isnan(back_corr_mean_sb(:,1))==0);
    true_angle_mean_sb=true_angle_mean_sb(isnan(true_angle_mean_sb(:,1))==0);

    A=fitlm(true_angle_mean_ps,back_corr_mean_ps(:,1),'Linear');
    B=fitlm(true_angle_mean_sb,back_corr_mean_sb(:,1),'Linear');
    
    A_Coef=table2array(A.Coefficients);
    B_Coef=table2array(B.Coefficients);
    
    A_AngInt=((A_Coef(2,1)*(IntAng*-1))+(A_Coef(1,1)));
    B_AngInt=((B_Coef(2,1)*IntAng+(B_Coef(1,1))));
    
    disp('PORTSIDE');
    disp(['Número de amostras: ' num2str(length(back_corr_mean_ps(:,1)))]);
    disp(['Mean: ' num2str(mean(back_corr_mean_ps(:,1)))]);
    disp(['Intercept at 0: ' num2str(A_Coef(1,1))]);
    disp(['Intercept at -' num2str(IntAng),': ' num2str(A_AngInt)]);
    disp(['Slope: ' num2str(A_Coef(2,1))]);
    disp(['True angle mean: ' num2str(mean(true_angle_mean_ps(:,1)))]);

    disp(' ')

    disp('STARBOARD');
    disp(['Número de amostras: ' num2str(length(back_corr_mean_sb(:,1)))]);
    disp(['Mean: ' num2str(mean(back_corr_mean_sb(:,1)))]);
    disp(['Intercept at 0: ' num2str(B_Coef(1,1))]);
    %disp(['Intercept angle: ' num2str(IntAng)]);
    disp(['Intercept at ' num2str(IntAng),': ' num2str(B_AngInt)]);
    disp(['Slope: ' num2str(B_Coef(2,1))]);
    disp(['True angle mean: ' num2str(mean(true_angle_mean_sb(:,1)))]);


    %Escrita dos resultados no ficheiro *_resultados.txt
    fprintf(fid20,'%s\n',['PS',',',num2str(mean(back_corr_mean_ps(:,1))),',',num2str(A_Coef(1,1)),',',num2str(A_AngInt),',',num2str(A_Coef(2,1)),',',num2str(mean(true_angle_mean_ps(:,1)))]);
%     fprintf(fid20,'%s\n',['SB',',',num2str(mean(back_corr_mean_sb(:,1))),',',num2str(B_Coef(1,1)),',',num2str(B_AngInt),',',num2str(B_Coef(2,1)),',',num2str(mean(true_angle_mean_sb(:,1)))]);
    
%     figure;
%     plot(true_angle_mean_ps,back_corr_mean_ps(:,1));
%     hold on; 
%     plot(A);
%     title(['AngResol: ' num2str(angRes),'    Avrg: ' num2str(mean(back_corr_mean_sb(:,1))),'    Intercept at -' num2str(IntAng),': '   num2str(A_AngInt),'     Slope: ' num2str(A_Coef(2,1))]);
%     xlabel('True angle interval average: (°)');
%     ylabel('Backscatter(db)')
% 
%     figure;
%     plot(true_angle_mean_sb,back_corr_mean_sb(:,1));
%     hold on;
%     plot(B);
%     title(['AngResol: ' num2str(angRes),'    Avrg: ' num2str(mean(back_corr_mean_ps(:,1))),'    Intercept at ' num2str(IntAng),': '  num2str(B_AngInt),'     Slope: ' num2str(B_Coef(2,1))]);
%     xlabel('True angle interval average: (°)');
%     ylabel('Backscatter(db)')

    
    beam_n_indices=beam_n(indices);
    count=1;
    for i=min(beam_n_indices(beam_n_indices>0)):(max(beam_n_indices(beam_n_indices<199)))
       back_corr_mean_ps_beam(count)=mean(back_corr_indices(beam_n_indices==i));
       true_angle_mean_ps_beam(count)=mean(true_angle_indices(beam_n_indices==i));
       count=count+1;
    end

    count=1;
    for i=min(beam_n_indices(beam_n_indices>199)):(max(beam_n_indices(beam_n_indices<399)))
       back_corr_mean_sb_beam(count)=mean(back_corr_indices(beam_n_indices==i));
       true_angle_mean_sb_beam(count)=mean(true_angle_indices(beam_n_indices==i));
       count=count+1;
    end

end

% close all

%%

%Classificador

%classificationLearner;

yfit = trainedModel_KNNN.predictFcn('43ProcessedBkSctt_ARC_Parameters.txt');

